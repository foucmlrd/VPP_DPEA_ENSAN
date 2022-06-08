# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:45:42 2022

@author: Cyprien Courtois
"""
import math as math
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from algo import *
from env_app import *

import scipy 

#%%

'''Fonction de calcul des forces exercées et centres de poussée velique pour chaque sailset'''

'''
ENTREE :

PARAMETRE d'ENTREE
  rho_a : masse volumique de l'eau
  AWA_deg 
  AWA_rad
  AWS_knts
  AWS_ms
  RED : Paramètre d'optimisation gérant la réduction de la voile principale et de la voile d'avant (ORC 2018)
  Mode : Si mode = 1, la réduction de la GV se fait en fixant le pourcentage des ris (Ris 1, Ris 2 et Ris 3).
  Si le mode = 0, la réduction de GV se fait par une fonction linéaire définie par le pourcentage du Ris 3

  FB average freebord

  
CARACTERISTIQUES DES RIS 
!! Guindant résiduel en pourcentage !! 
  Ris_1    
  Ris_2
  Ris_3

    
    
    
    
SORTIE :
    Fx
    Fy
    Fz
    Zce
    Yce
    Mx
    Pour :
        MF
        MJ
        MSA
        MSS

'''
NameCoefAero = ["AWA (°)", "Cl GV", "Cl Génois", "Cl Spi Sym", "Cl Spi Asy", "Cd GV", "Cd Génois", "Cl Spi Sym", "Cl Spi Asy"]
MatCoeffAero = np.array([[0, 0, 0, 0, 0, 0.0431, 0, 0, 0],
                         [7, 0.905175, 0,	0, 0, 0.02586, 0.05, 0, 0],
                         [9, 1.094825, 0,	0, 0, 0.02328, 0, 0, 0],
                         [12,	1.206895, 0, 0,	0, 0.02328,	0, 0, 0],
                         [15,	0, 1.05, 0,	0, 0, 0.032, 0, 0],
                         [20,	0, 1.425, 0, 0,	0, 0.031, 0, 0],
                         [27,	0, 1.475, 0, 0,	0, 0.037, 0, 0],
                         [28,	1.386895, 0, -0.02484, 0.0183, 0.03259,	0, 0.19152, 0.16215],
                         [40,	0, 0, 0, 0,	0, 0, 0, 0],
                         [41,	0, 0, 0.69437, 0.735, 0, 0, 0.28152, 0.25184],
                         [50,	0, 1.45, 0.90677, 0.94666, 0, 0.25, 0.35496, 0.32502],
                         [60,	1.36832, 1.25, 1.044, 1.08342, 0.11302, 0.35, 0.4392, 0.40897],
                         [67,	0, 0, 1.08, 1.10494, 0, 0, 0.4896, 0.4592],
                         [75,	0, 0, 1.08,	1.09059, 0,	0, 0.5328, 0.50225],
                         [80,	0, 0, 0, 0,	0, 0, 0, 0],
                         [90,	1.26724, 0,	0, 0, 0.3825, 0, 0, 0],
                         [100, 0, 0.4, 0.9576, 0.95427, 0, 0.73, 0.6192, 0.59839],
                         [115, 0,	0, 0.8136, 0.81077, 0, 0, 0.6588, 0.65292],
                         [120, 0.93103, 0, 0, 0, 0.96888, 0, 0, 0],
                         [130, 0, 0, 0.612, 0.60987, 0, 0, 0.6732, 0.67086],
                         [140, 0,	0, 0, 0, 0,	0, 0, 0],
                         [150, 0.38793, 0, 0.324, 0.32287, 1.31578, 0.95, 0.6732, 0.67086],
                         [160, 0,	0, 0, 0, 0,	0, 0, 0],
                         [170, 0,	0, 0.108, 0.10762, 0, 0, 0.6732, 0.67086],
                         [180, -0.11207, -0.1, 0, 0, 1.34483, 0.9, 0.6732, 0.67086]])



'''INPUT'''
#Franc bord
FB = 1.3 #(m)average freeboard
#Caractéristiques des ris -> Guindant de GV résiduel en %
Reef1 = 0.85
Reef2 = 0.70
Reef3 = 0.55

#MAIN SAIL
P = 16.7# (m) Main sail hoist
E = 5.6# (m) Foot of main sail
BAD = 0.61# (m) Height of main boom above sheer
BD = 0 #INPUT POUR LA BOME (Viens de la fonction FARDAGE)
#On import BD de fardage  

MGL = 4.2 #(m) Voir schéma en annexe
MGM = 2.8 #(m) Voir schéma en annexe
MGU = 1.45 #(m) Voir schéma en annexe
MGT = 0.7 #(m) Voir schéma en annexe
HB = 0 #(m) Voir schéma en annexe

#GENOA/ JIB
I = 16.2 #(m) Height of jib (HLU)
J = 5.1 #(m) Base of jib (Bordure)
LP = 5.4 #(m) Perpendicular of longest jib (HLP)

#Foretriangle
IFore = I #(m) Height of foretriangle
JFore = J #(m) Base of foretriangle

#SYMMETRIC SPINNAKER
SL = 18 #(m) Guindant ou chute (SLU/ SLE)
SF = 8 #(m) Borure (SHW)
SMG = 8 #(m) Largeur à mi-hauteur (SHW)

#ASYMMETRIC SPINNAKER
ALU = 18 #(m) Guindant (SLU)
ALE = 18 #(m) Chute (SLE)
ASF= 8 #(m) Bordure (SFL)
AMG = 8 #(m) Largeur à mi-hauteur (SHW)

#CARACTERISTIQUES DU MAT POUR CALCULER AR
MH = 18 #(m) Mast height above sheer
MD = 0.12 #(m) Mast diameter

'''END INPUT'''

'''PARAMETRES DE REDUCTION'''
#Flat
if RED <= 2 :
  Flat = 0.6
else : Flat = 0.38*RED - 0.14

#Reef Main  
if Mode == 0:
  if RED > 1 :
    Reef_Main = 1
  else :
    Reef_Main = (1 - Reef3)*RED + Reef3
else : #Si Mode == 1
  if RED >= 1:
    Reef_Main = 1
  elif  RED < 1 & RED >= 0.666 :
    Reef_Main = Reef1
  elif RED < 0.666 & RED >= 0.333: 
    Reef_Main = Reef2
  elif RED < 0.333 :
    Reef_Main = Reef3

#Reef Jib
if RED < 1 : 
  Reef_Jib = 0.8
else : 
  Reef_Jib = 0.2*RED+0.6
  
'''END PARAMETRES DE REDUCTION'''

'''CALCUL DE SURFACE DE GV selon ORC 2018'''
MGMH = (P/2) + E*(MGM - E/2)/P #(m)

MGLH = (MGMH/2) + ((MGL - (E + MGM)/2)/MGMH) * (E - MGM) #(m)

MGUH = (MGMH + P)/2 + MGM*(MGU - (MGM/2))/(P-MGMH) #(m)

MGTH = (MGUH + P)/2 + (MGT - (MGU/2)/(P - MGUH))*MGU #(m)

Upper_area = (P/8)*(MGL+2*MGM + 1.5*MGU +0.5*HB + MGT) #(m) Upper3/4_Aera

ROACH = ((Upper_area / (0.375*P*MGL)) - 1)/0.844 #Si ROACH > 0, il y a une corne, sinon pas.

#CALCUL DE LA SURFACE DE LA GRANDE VOILE (m²)
if BD > 0.06*E :
  Main_Sail_Area = (2*E*(BD-0.06*E))*((MGLH*(MGL+E)/2)+(((MGL+MGM)/2)*(MGMH-MGLH))+((MGM+MGU)/2)*(MGUH-MGMH)+((MGT+MGU)/2)*(MGTH-MGUH)+((MGT+HB)/2)*(P-MGTH))
else :
    
  Main_Sail_Area = (((MGL+E)/2)*MGLH)+((MGL+MGM)/2)*(MGMH-MGLH)+((MGM+MGU)/2)*(MGUH-MGMH)+(((MGT+MGU)/2)*(MGTH-MGUH))+(((MGT+HB)/2)*(P-MGTH))

'''CALCUL DE LA SURFACE NOMINALE '''
SA_N = 0.5* (P*E + J*I)#(m²) Nominal Sail Area

'''CALCUL DE L' ASPECT RATIO '''
if (AWA_deg < 45):
  AR = math.pow(1.1*(MH+FB),2)/SA_N
else : 
  AR = math.pow(1.1*MH,2)/SA_N

'''CALCUL DES SURFACES DE VOILURE'''
  
SA_For = 0.5*IFore*JFore #(m²)Foretriangle Are

SA_Jib =0.5*math.pow(math.pow(I,2) + math.pow(J,2) , 1/2)*LP #(m²) Genoa/ Jib Area

SA_SpiA = 0.5*(ALU+ALE)*(ASF+4*AMG)/6 #(m²) Asymmetric Spinnaker Area

SA_SpiS = (SL*(SF+4*SMG)/6)*math.pow(Reef_Jib,2)#(m²) Symmetric Spinnaker Area

SA_MF = Main_Sail_Area + SA_For #(m²) Main + Foretriangle

SA_MJ = Main_Sail_Area + SA_Jib #(m²) Main + Jib

SA_MSA = Main_Sail_Area + SA_SpiA #(m²) Main + Asymmetric Spinnaker Area

SA_MSS = Main_Sail_Area + SA_SpiS #(m²) Main + Symmetric Spinnaker Area

'''CALCUL DES HAUTEURS Z DES CENTRES VELIQUES PAR RAPPORT A LA FLOTTAISON''' #LARSOON, p.162

ZCEM =(0.39*P+BAD)*Reef_Main+FB #(m) ZCE Main

ZCEF =0.39*IFore*Reef_Jib+FB #(m) ZCE Foretriangle

ZCEJ = 0.39*I*Reef_Jib+FB #(m) ZCE Jib 

ZCESPIA = 0.59*I*Reef_Jib+FB #(m) ZCE Asymmetric Spinnaker ATTENTION FORMULE A VERIFIER

ZCESPIS = 0.59*I*Reef_Jib+FB #(m) ZCE Symmetric Spinnaker ATTENTION FORMULE A VERIFIER

ZCEMF = (ZCEM*Main_Sail_Area+ZCEF*SA_For)/SA_MF #(m) ZCE Main + Foretriangle

ZCEMJ = (ZCEM*Main_Sail_Area+ZCEJ*SA_Jib)/SA_MJ #(m) ZCE Main + Jib

ZCEMSA = (ZCEM*Main_Sail_Area+ZCESPIA*SA_SpiA)/SA_MSA #(m) ZCE Main + Asymmetric Spinnaler 

ZCEMSS = (ZCEM*Main_Sail_Area+ZCESPIS*SA_SpiS)/SA_MSS #(m) ZCE Main + Symmetric Spinnaler 

'''CALCUL DU COEFFICIENT DE PORTANCE CL'''

AWA_CL_GV = [0,7,9,12,28,60,90,120,150,180]
CL_GV = [0,0.905,1.095,1.207,1.387,1.368,1.267,0.931,0.388,-0.112] #[[0,0], [7,0.905], [9,1.095],[12,1.207],[28,1.387],[60,1.368],[90,1.267],[120,0.931],[180,-0.112]] #[AWA°, Cl_GV]

tck = interpolate.splrep(AWA_CL_GV, CL_GV)
print(AWA_deg)
print(interpolate.splev(AWA_deg, tck))




'''CL_M = scipy.interpolate.CubicSpline(x, y, axis=0, bc_type='not-a-knot', extrapolate=None) 

cubspline(3,AWA_deg,Coeff.Aero!B7:B31,Coeff.Aero!C7:C31)*Flat*ReefMain^2 #Main

CL_For = #Foretriangle

CL_Jib = #Genoa/Jib

CL_SpiA = #Asymmetric Spinnaker

CL_SpiS = #Symmetric Spinnaker'''