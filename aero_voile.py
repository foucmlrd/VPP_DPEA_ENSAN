# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:45:42 2022

@author: Foucauld Malard
"""
import math as math
import numpy as np
from algo import *
from env_app import *
from input import *

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

ZCEMSA = (ZCEM*Main_Sail_Area+ZCESPIA*SA_SpiA)/SA_MSA #(m) ZCE Main + Asymmetric Spinnaker 

ZCEMSS = (ZCEM*Main_Sail_Area+ZCESPIS*SA_SpiS)/SA_MSS #(m) ZCE Main + Symmetric Spinnaker 

'''CALCUL DU COEFFICIENT DE PORTANCE CL'''

AWA = [0.0, 7, 9, 12, 15, 20, 27, 28, 40, 41, 50, 60, 67, 75, 80, 90, 100, 115, 120, 130, 140, 150, 160, 170, 180]
CL_GV = [0.0, 0.905175, 1.094825, 1.206895,None, None, None, 1.386895, None, None, None, 1.36832, None, None, None, 1.26724, None, None, 0.93103, None, None, 0.38793, None, None, -0.11207]
CL_GENOIS = [None, 0, None, None, 1.050, 1.425, 1.475, None, None, None, 1.450, 1.250, None, None, None, None, 0.4, None, None, None, None, 0, None, None, -0.1]
CL_SPI_SYM = [None, None, None, None, None, None, None, -0.02484, None, 0.69437, 0.90677, 1.044, 1.080, 1.080, None, None, 0.9576, 0.8136, None, 0.612, None, 0.324, None, 0.108, 0]
CL_SPI_ASY = [None, None, None, None, None, None, None, -0.0183, None, 0.735, 0.94666, 1.08342, 1.10494, 1.09059, None, None, 0.95427, 0.81077, None, 0.60987, None, 0.32287, None, 0.10762, 0]

CD_GV = [0.0431, 0.02586, 0.02328, 0.02328,None, None, None, 0.03259, None, None, None, 0.11302, None, None, None, 0.3825, None, None, 0.96888, None, None, 1.31578, None, None, 1.34483]
CD_GENOIS = [None, 0.05, None, None, 0.032, 0.031, 0.037, None, None, None, 0.25, 0.35, None, None, None, None, 0.73, None, None, None, None, 0.950, None, None, 0.9]
CD_SPI_SYM = [None, None, None, None, None, None, None, 0.19152, None, 0.28152, 0.35496, 0.4392, 0.4896, 0.5328, None, None, 0.6192, 0.6588, None, 0.6732, None, 0.6732, None, 0.6732, 0.6732]
CD_SPI_ASY = [None, None, None, None, None, None, None, 0.16215, None, 0.25184, 0.32502, 0.40897, 0.4592, 0.50225, None, None, 0.59839, 0.65292, None, 0.67086, None, 0.67086, None, 0.67086, 0.67086]


'''tck = interpolate.splrep(AWA_CL_GV, CL_GV)
print(AWA_deg)
print(interpolate.splev(AWA_deg, tck))'''


def Cubic_Spline(methode, xi, x_points, y_points):
  x=[]
  y=[]
  if methode == 1:
    j=0
  else :
    j=-1

  i=0
  while (i<len(x_points)):
    if(y_points[i]!= None):
      x.append(x_points[i])
      y.append(y_points[i])
    i=i+1

  if (methode ==1):
  #NR cubic spline
  #Get y2
    y2 = Spline(x, y,len(x), pow(10,30), pow(10,30))
    yi = splint (x, y, y2, len(x), xi)
  elif methode ==3 :
    #Own cubic spline
    yi = SplineX3(xi, x, y)
  return yi

  
#PAS UTILISE ICI..
def Spline(x_points, y_points, N, ypl, ypn):
  '''Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., y i = f(xi), with
  x1<x2< :::<xN , and given values yp1 and ypn for the first derivative of the interpolating function at points 1     and n, respectively, this routine returns an array y2(1:n) of length n which contains the second derivatives of the   interpolating function at the tabulated points xi. If yp1 and/or ypn are equal to 1 * 10^30 or larger, the routine   is signaled to set the corresponding boundary condition for a natural spline, with zero second derivative on that   boundary. 
  '''

  u=[]
  y2_points = [0]
  #The lower boundary condition is set either to be natural
  if (ypl > 9.9e29):
    y2_points[0]=0
    u[0]=0
  else:
  #or else to have a specicied first derivative.
    y2_points[0] = -0.5
    u.append((3 / (x_points[1] - x_points[0])) * ((y_points[1] - y_points[0])/ (x_points[1] - x_points[0]) - ypl))

#This is the decomposition loop of the tridiagonal algorithm. y2_points and u are used for temporary storage of the decomposed factors.
  for i in range(1,N-1,1):
    sig = (x_points[i]-x_points[i-1])/(x_points[i+1] - x_points[i-1])
    p = sig * y2_points[i-1] + 2
    y2_points.append((sig - 1)/p)
    u.append( (6 * ((y_points[i+1]- y_points[i])/(x_points[i+1]-x_points[i]) - (y_points[i] - y_points[i-1]) / (x_points[i] - x_points[i-1])) / (x_points[i+1] - x_points[i-1]) - sig*u[i-1]) / p)
  #The upper boundary condition is set either to be natural
  if (ypn > 9.9e29):
    qn = 0
    un = 0
  else :
    #or else to have a specified first derivative.
    qn = 0.5
    un = (3/ (x_points[N-1] - x_points[N-2]))*(ypn - (y_points[N-1] - y_points[N - 2])/ (x_points[N-1] - x_points[N-2]))
  y2_points.append((un - qn*u[N-2])/(qn * y2_points[N-2]+1))
  #This is the backsubstitution loop of the tridiagonal algorithm

  for k in range(N - 2, -1, -1):
    y2_points[k] = y2_points[k]*y2_points[k+1] + u[k]
  return y2_points

#PAS UTILISE ICI..
def splint (xa, ya, y2a, N, x):
  '''Given the arrays xa(1:n) and ya(1:n)   of length n, which tabulate a function   (with the       xai 's in order), and given   the array y2a(1:n), which is the output   from spline above, and    given a value of x, this routine returns a cubic-spline interpolated value y.'''
  klo=1
  khi = N
  while khi - klo > 1:
    k = int((khi + klo)/2)
    
    if xa[k] > x :
      khi = k
    else :
      klo = k
  #klo and khi now bracket the input value of x.
  h = xa[khi] - xa[klo]
  if(h==0):
    print("bad xa input in splint")
  #Cubic spline polynomial is now evaluated
  A = (xa[khi]-x)/h
  B = (x - xa[klo])/h
  y = A * ya[klo] + B * ya[khi] + ((pow(A,3)-A)*y2a[klo]+(pow(B,3)-B) * y2a[khi]) * pow(h,2)/6
  return y

def dxx(x1, x0):
  dxxx = x1 - x0
  if dxxx == 0:
    dxxx = pow(10,30)
  return dxxx

def SplineX3 (x, xx, yy):
  '''Function returns y value for a corresponding x value, based on cubic spline.
   Will never oscillates or overshoot. No need to solve matrix.
   Also calculate constants for cubic in case needed (for integration).

   xx(0 to No_of_lines) is x values
    * Must be unique (no two consequetive ones the same)
    * Must be in ascending order
    * No of lines = Number of points - 1
   yy(0 to No_of_lines) is y values

   Uses function dxx to prevent div by zero.

  Developer: C Kruger, Guildford, UK
  Date: December 2001'''
  #1st and 2nd derivative for left and right ends of line
  gxx = [0,0] #Initiation puis ils sont remplis
  ggxx = [0,0]
  Num = 0
  
  #Number of lines = points - 1
  Nmax = len(xx)
  
  #(1a) Find LineNumber or segment. Linear extrapolate if outside range.
  if x < xx[0] or x > xx[Nmax-1]:
    
  #X outisde range. Linear interpolate - Below min or max?
    if x < xx[0]:
      Num=1
    else:
      Num = Nmax - 1
    B = (yy[Num] - yy[Num - 1])/dxx(xx[Num],xx[Num - 1])
    A = yy[Num] - B * xx[Num]
    return A+ B*x 
    
    #1b) Find LineNumber or segment. Linear extrapolate if outside range.
  else :
    #X in range. Get line.
    for i in range(Nmax):
      if x <= xx[i]:
        Num = i
        break
        
  #(2) Calc first derivative (slope) for intermediate points   
  for j in range(2):
    if (Num == 0): 
      i = Num + j
    elif Num == Nmax - 1:
      i = Num - 2
    else :
      i = Num - 1 + j #Two points around line

    if i==0 or i == Nmax -1:
      #Set very large slope at ends
      gxx[j] = pow(10,30)
    elif yy[i+1] - yy[i] == 0 or yy[i] - yy[i-1] == 0:
      #Only check for 0 dy. dx assumed NEVER equals 0 !
      gxx[j]=0
    elif ((xx[i+1] - xx[i])/(yy[i+1] - yy[i]) + (xx[i] - xx[i-1])/(yy[i] - yy[i-1])) == 0 :
      #Pos PLUS neg slope is 0. Prevent div by zero.
      gxx[j]=0
    elif (yy[i+1]-yy[i]) * (yy[i] - yy[i-1])<0 :
      #Pos AND neg slope, assume slope = 0 to prevent overshoot
      gxx[j]=0
    else : 
      #Calculate an average slope for point based on    connecting lines
      gxx[j]= 2 / (dxx(xx[i+1], xx[i])/(yy[i+1] - yy[i]) + dxx(xx[i], xx[i-1])/(yy[i]-yy[i-1]))
      
  # (3) Reset first derivative (slope) at first and last point
  if Num == 0 :
    #First point has 0 2nd derivative
    gxx[0] = 3 / 2*(yy[Num +1] - yy[Num])/dxx(xx[Num+1],xx[Num]) - gxx[1]/2
  if Num == Nmax - 1:
    #Last point has 0 2nd derivative
    gxx[1] = 3 / 2*(yy[Num] - yy[Num - 1])/dxx(xx[Num],xx[Num-1]) - gxx[0]/2
    
  #(4) Calc second derivative at points
  ggxx[0] = -2 *(gxx[1] + 2* gxx[0])/dxx(xx[Num],xx[Num-1]) + 6 * (yy[Num]- yy[Num -1])/pow(dxx(xx[Num], xx[Num-1]),2)
  ggxx[1] = 2 *(2*gxx[1] + gxx[0])/dxx(xx[Num],xx[Num-1]) - 6 * (yy[Num]- yy[Num -1])/pow(dxx(xx[Num], xx[Num-1]),2)

  #(5) Calc constants for cubic
  D = 1/ 6 *(ggxx[1]-ggxx[0])/dxx(xx[Num], xx[Num-1])
  C = 1/2* (xx[Num]*ggxx[0] - xx[Num-1]*ggxx[1])/dxx(xx[Num], xx[Num-1])
  B = (yy[Num] - yy[Num -1] - C*(pow(xx[Num],2) - pow(xx[Num-1],2)) - D * (pow(xx[Num],3) - pow(xx[Num-1],3)))/dxx(xx[Num], xx[Num-1])
  A = yy[Num -1] - B * xx[Num - 1] - C*pow(xx[Num - 1],2) - D* pow(xx[Num-1],3)
  
  #Return function
  Result = A + B * x + C * pow(x,2) + D * pow(x,3)
  return Result
  


#y_test = Cubic_Spline(3, 20, AWA_CL_GV, CL_GV)
#print(y_test)*

#Calcul du coefficient de portance Cl
  
CL_M = Cubic_Spline(3, AWA_deg, AWA, CL_GV)*Flat*pow(Reef_Main,2) #Main

CL_For = Cubic_Spline(3, AWA_deg, AWA, CL_GENOIS)*Flat*pow(Reef_Jib,2) #Foretriangle

CL_Jib = Cubic_Spline(3, AWA_deg, AWA, CL_GENOIS)*Flat*pow(Reef_Jib,2) #Genoa/Jib

CL_SpiA = Cubic_Spline(3, AWA_deg, AWA, CL_SPI_ASY)*Flat*pow(Reef_Jib,2) #Asymmetric Spinnaker

CL_SpiS = Cubic_Spline(3, AWA_deg, AWA, CL_SPI_SYM)*Flat*pow(Reef_Jib,2) #Symmetric Spinnaker

#Calcul du coefficient de trainée parasite Cdp
Cdp_M = Cubic_Spline(3, AWA_deg, AWA, CD_GV)*pow(Reef_Main,2) #Main

Cdp_For = Cubic_Spline(3, AWA_deg, AWA, CD_GENOIS)*pow(Reef_Jib,2) #Foretriangle

Cdp_Jib = Cubic_Spline(3, AWA_deg, AWA, CD_GENOIS)*pow(Reef_Jib,2) #Genoa/Jib

Cdp_SpiA = Cubic_Spline(3, AWA_deg, AWA, CD_SPI_ASY)*pow(Reef_Jib,2) #Asymmetric Spinnaker

Cdp_SpiS = Cubic_Spline(3, AWA_deg, AWA, CD_SPI_SYM)*pow(Reef_Jib,2) #Symmetric Spinnaker

#Calcul du coefficient de trainée induite Cdi
Cdi_M =(pow(CL_M,2))/(math.pi*AR+0.005) #Main

Cdi_For = (pow(CL_For,2))/(math.pi*AR+0.005) #Foretriangle

Cdi_Jib = (pow(CL_Jib,2))/(math.pi*AR+0.005) #Genoa/Jib

Cdi_SpiA = (pow(CL_SpiA,2))/(math.pi*AR+0.005) #Asymmetric Spinnaker

Cdi_SpiS = (pow(CL_SpiS,2))/(math.pi*AR+0.005) #Symmetric Spinnaker

#Calcul de la portance
lift_M = 0.5*rho_a*pow(AWS_ms,2)*Main_Sail_Area*CL_M #(N) Main

lift_For = 0.5*rho_a*pow(AWS_ms,2)*SA_For*CL_For #(N) Foretriangle

lift_Jib = 0.5*rho_a*pow(AWS_ms,2)*SA_Jib*CL_Jib #(N) Genoa/Jib

lift_SpiA = 0.5*rho_a*pow(AWS_ms,2)*SA_SpiA*CL_SpiA #(N) Asymmetric Spinnaker

lift_SpiS = 0.5*rho_a*pow(AWS_ms,2)*SA_SpiS*CL_SpiS #(N) Symmetric Spinnaker

#Calcul de la trainée parasite
dp_M = 0.5*rho_a*pow(AWS_ms,2)*Main_Sail_Area*Cdp_M #(N) Main

dp_For = 0.5*rho_a*pow(AWS_ms,2)*SA_For*Cdp_For #(N) Foretriangle

dp_Jib = 0.5*rho_a*pow(AWS_ms,2)*SA_Jib*Cdp_Jib #(N) Genoa/Jib

dp_SpiA = 0.5*rho_a*pow(AWS_ms,2)*SA_SpiA*Cdp_SpiA #(N) Asymmetric Spinnaker

dp_SpiS = 0.5*rho_a*pow(AWS_ms,2)*SA_SpiS*Cdp_SpiS #(N) Symmetric Spinnaker

#Calcul de la trainée induite
di_M = 0.5*rho_a*pow(AWS_ms,2)*Main_Sail_Area*Cdi_M #(N) Main

di_For = 0.5*rho_a*pow(AWS_ms,2)*SA_For*Cdi_For #(N) Foretriangle

di_Jib = 0.5*rho_a*pow(AWS_ms,2)*SA_Jib*Cdi_Jib #(N) Genoa/Jib

di_SpiA = 0.5*rho_a*pow(AWS_ms,2)*SA_SpiA*Cdi_SpiA #(N) Asymmetric Spinnaker

di_SpiS = 0.5*rho_a*pow(AWS_ms,2)*SA_SpiS*Cdi_SpiS #(N) Symmetric Spinnaker

#Calcul de la trainée totale
drag_M = dp_M + di_M #(N) Main

drag_For = dp_For + di_For #(N) Foretriangle

drag_Jib = dp_Jib + di_Jib #(N) Genoa/Jib

drag_SpiA = dp_SpiA + di_SpiA #(N) Asymmetric Spinnaker

drag_SpiS = dp_SpiS + di_SpiS #(N) Symmetric Spinnaker

#Calcul de la portance des différents sailsets
lift_MF = lift_M + lift_For #(N) Main + Foretriangle

lift_MJ = lift_M + lift_Jib #(N) Main + Jib

lift_MSA = lift_M + lift_SpiA #(N) Main + Asymmetric Spinnaker 

lift_MSS = lift_M + lift_SpiS #(N) Main + Symmetric Spinnaker 

#Calcul de la trainée des différents sailsets
drag_MF = drag_M + drag_For #(N) Main + Foretriangle

drag_MJ = drag_M + drag_Jib #(N) Main + Jib

drag_MSA = drag_M + drag_SpiA #(N) Main + Asymmetric Spinnaker 

drag_MSS = drag_M + drag_SpiS #(N) Main + Symmetric Spinnaker 

#Driving force     #Repère : Réf : X - axe, Y - babord, Z - verticale
driving_force_MF = (lift_MF * math.sin(AWA_rad) - drag_MF * math.cos(AWA_rad))*math.cos(Derive_rad) #(N) Main + Foretriangle

driving_force_MJ = (lift_MJ * math.sin(AWA_rad) - drag_MJ * math.cos(AWA_rad))*math.cos(Derive_rad) #(N) Main + Jib

driving_force_MSA = (lift_MSA * math.sin(AWA_rad) - drag_MSS * math.cos(AWA_rad))*math.cos(Derive_rad) #(N) Main + Asymmetric Spinnaker 

driving_force_MSS = (lift_MSS * math.sin(AWA_rad) - drag_MSS * math.cos(AWA_rad))*math.cos(Derive_rad) #(N) Main + Symmetric Spinnaker 

#Side force     #Repère : Réf : X - axe, Y - babord, Z - verticale
side_force_MF = (lift_MF * math.cos(AWA_rad) + drag_MF * math.sin(AWA_rad))*math.cos(Heeling_rad)*math.cos(Derive_rad) #(N) Main + Foretriangle

side_force_MJ = (lift_MJ * math.cos(AWA_rad) + drag_MJ * math.sin(AWA_rad))*math.cos(Heeling_rad)*math.cos(Derive_rad) #(N) Main + Jib

side_force_MSA = (lift_MSA * math.cos(AWA_rad) + drag_MSA * math.sin(AWA_rad))*math.cos(Heeling_rad)*math.cos(Derive_rad) #(N) Main + Asymmetric Spinnaker 

side_force_MSS = (lift_MSS * math.cos(AWA_rad) + drag_MSS * math.sin(AWA_rad))*math.cos(Heeling_rad)*math.cos(Derive_rad) #(N) Main + Symmetric Spinnaker

#Vertical force FZ     #Repère : Réf : X - axe, Y - babord, Z - verticale
FZ_MF = (lift_MF * math.cos(AWA_rad) + drag_MF * math.sin(AWA_rad))*math.sin(Derive_rad) #(N) Main + Foretriangle

FZ_MJ = (lift_MJ * math.cos(AWA_rad) + drag_MJ * math.sin(AWA_rad))*math.sin(Derive_rad) #(N) Main + Jib

FZ_MSA = (lift_MSA * math.cos(AWA_rad) + drag_MSS * math.sin(AWA_rad))*math.sin(Derive_rad) #(N) Main + Asymmetric Spinnaker 

FZ_MSS = (lift_MSS * math.cos(AWA_rad) + drag_MSS * math.sin(AWA_rad))*math.sin(Derive_rad) #(N) Main + Symmetric Spinnaker 

#Valeurs de sortie
#FX
FX_MF = max(0, driving_force_MF) #(N) Main + Foretriangle
FX_MJ = max(0, driving_force_MJ) #(N) Main + Jib
FX_MSA = max(0, driving_force_MSA) #(N) Main + Asymmetric Spinnaker 
FX_MSS = max(0, driving_force_MSS) #(N) Main + Symmetric Spinnaker 

#FY
if FX_MF == 0 : FY_MF = 0 
else : FY_MF = side_force_MF #(N) Main + Foretriangle
  
if FX_MJ == 0 : FY_MJ = 0 
else : FY_MJ = side_force_MJ #(N) Main + Jib
  
if FX_MSA == 0 : FY_MSA = 0 
else : FY_MSA = side_force_MSA #(N) Main + Asymmetric Spinnaker 
  
if FX_MSS == 0 : FY_MSS = 0 
else : FY_MSS = side_force_MSS #(N) Main + Symmetric Spinnaker 

#Mx (N.m)
MX_MF = FY_MF * ZCEMF + FZ_MF *(ZCEMF * math.sin(Heeling_rad)) #(N) Main + Foretriangle ZCEMF * math.sin(Heeling_rad) c'est Yce
MX_MJ = FY_MJ * ZCEMJ + FZ_MJ *(ZCEMJ * math.sin(Heeling_rad)) #(N) Main + Jib
MX_MSA = FY_MSA * ZCEMSA + FZ_MSA *(ZCEMSA * math.sin(Heeling_rad)) #(N) Main + Asymmetric Spinnaker 
MX_MSS = FY_MSS * ZCEMSS + FZ_MSS *(ZCEMSS * math.sin(Heeling_rad)) #(N) Main + Symmetric Spinnaker 
