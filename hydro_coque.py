# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:37:43 2022

@author: Cyprien Courtois
"""
import math as math
import numpy as np
import matplotlib.pyplot as plt


#%%

''' Fonction de calcul de la résistance hydrodynamique de la coque '''

'''
ENTREE : 
    Paramètres extérieurs : 
        rho, 
        Vb (vitesse boat), 
        Heeling (angle de gîte)
    Paramètres géométriques : 
        Vc (déplacement), 
        Lwl (waterline length), 
        Bwl (waterline beam), 
        B (beam), 
        Fb (average freeboard) #INPUT DANS AERO_VOILE, 
        T (total dragt including keel),
        Tc (draft carène)
        Lcb (longitudinal center of buyancy from forward perpendicular)
        Lcf (longitudianl center of flotation from froward perpendicular)
        Cp (Prismatic coefficient)
        Cm (Midship section coefficient)
        AW (Waterplane area)
        Sw_h (Wetted surface)
        Kg (center of gravity above baseline)
        
SORTIE : 
    Repère sailing : X-heading, Y-ship, Z-ship
    R_coque 
    Zcc CLR
    Repère réf : X-axe, Y-babord, Z-vertiale
    Fx
    Fy
    Zce
    Mx 
    Tracé Rrh fct Fr
'''

# Type : SansGite, AvecGite, Visqueux

def methode_larsson(Vc,rho,g,Lwl,Lcb,Cp,Aw,Bwl,Lcf,Tc,Cm, Sw_h, Type = "None"):
    if Type == "None":
        return("Entrer un type de calcul")
    
    if Type == "SansGite" :
        LarssonMat = np.array([[-0.0005, 0.0023, -0.0086, -0.0015, 0.0061, 0.001, 0.0001, 0.0052],
                      [-0.0003, 0.0059, -0.0064, 0.007, 0.0014, 0.0013, 0.0005, -0.002],
                      [-0.0002, -0.0156, 0.0031, -0.0021, -0.007, 0.0148, 0.001, -0.0043],
                      [-0.0009, 0.0016,0.0337,-0.0285,-0.0367,0.0218,0.0015,-0.0172],
                      [-0.0026, -0.0567,0.0446,-0.1091,-0.0707,0.0914,0.0021,-0.0078],
                      [-0.0064, -0.4034,-0.125,0.0273,-0.1341,0.3578,0.0045,0.1115],
                      [-0.0218, -0.5261,-0.2945,0.2485,-0.2428,0.6293,0.0081,0.2086],
                      [-0.0388, -0.5968,-0.3038,0.6033,-0.043,0.8332,0.0106,0.1336],
                      [-0.0347, -0.4764,-0.2361,0.8726,0.4219,0.899,0.0096,-0.2272],
                      [-0.0361, 0.0037,-0.296,0.9661,0.6123,0.7534,0.01,-0.3352],
                      [0.0008, 0.3728,-0.3667,1.3957,1.0343,0.323,0.0072,-0.4632],
                      [0.0108, -0.1238,-0.2026,1.1282,1.1836,0.4973,0.0038,-0.4477],
                      [0.1023, 0.7726, 0.504, 1.7867, 2.1934, -1.5479, -0.0115, -0.0977]])
                      
        K = np.array([Lwl/(Vc**(1/3)), Lcb/Lwl, Cp, Vc**(2/3)/Aw, Bwl/Lwl, Lcb/Lcf, Bwl/Tc, Cm])
        
        k1 = Vc*rho*g
        k2 = Vc**(1/3)/Lwl   
        
        X = np.linspace(0.15,0.75,13)
        Y = []
        for i in range (0,13):
            Y += [k1*k2*np.dot(K,LarssonMat[i])]
        plt.xlabel("Fr")
        plt.ylabel("Rrh")
        plt.plot(X,Y)
        plt.show()
        
        return()
        
        
    if Type == "AvecGite" :
        
        LarssonMat0 = np.array([[-0.0268, -0.0014, -0.0057, 0.0016, -0.007, -0.0017],
                               [0.6628,	-0.0632, -0.0699, 0.0069, 0.0459, -0.0004],
                               [1.6433,	-0.2144, -0.164, 0.0199, -0.054, -0.0268],
                               [-0.8659, -0.0354, 0.2226, 0.0188, -0.58, -0.1133],
                               [-3.2715, 0.1372, 0.5547, 0.0268, -1.0064, -0.2026],
                               [-0.1976, -0.148, -0.6593, 0.1862, -0.7489, -0.1648],
                               [1.5873, -0.3749, -0.7105, 0.2146, -0.4818, -0.1174]])
        LarssonMat = LarssonMat0/1000
        
        k1 = Vc*rho*g
        
        K = np.array([1,Lwl/Bwl, Bwl/Tc, (Bwl/Tc)**2, (Lcb-Lwl/2)/Lwl, ((Lcb-Lwl/2)/Lwl)**2])
        X = np.linspace(0.25,0.55,7)
        Y = []
        for i in range (0,7):
            Y += [k1*np.dot(K,LarssonMat[i])]
        plt.xlabel("Fr")
        plt.ylabel("Delta_Rh_20deg")
        plt.plot(X,Y, label= "Avec Gite")
        plt.show()
        return()
    
    if Type == "Visqueux" :

        LarssonMat = np.array([[-4.112, 0.054, -0.027, 6.329],
                               [-4.522, -0.132, -0.077, 8.738],
                               [-3.291, -0.389, -0.118, 8.949],
                               [1.85, -1.2, -0.109, 5.364],
                               [6.51, -2.305, -0.066, 3.443],
                               [12.334, -3.911, 0.024, 1.767],
                               [14.648, -5.182, 0.102, 3.497]])        
        
        
        k1 = 0.01
        
        K = np.array([1, Bwl/Tc, (Bwl/Tc)**2, Cm])
        X = np.linspace(0,35,8)
        print(X)
        Y = [Sw_h]
        for i in range (0,7):
            Y += [Sw_h*(1+k1*np.dot(K,LarssonMat[i]))]
        print(Y)
        plt.xlabel("Heel angle")
        plt.ylabel("Sw_hH")
        plt.plot(X,Y, label= "frottement visqueux")
        plt.show()
        return()
    else:
        
        return("choisir avec ou sans gite")
    
methode_larsson(6.3,1025,9.81,11.9,6.47,0.56,26.34,3.18,6.81,0.4,0.72,27.34, Type="SansGite")
#Test = np.array([6.443,0.544,0.560,0.13,0.267,0.95,7.95,0.715])
methode_larsson(6.3,1025,9.81,11.9,6.47,0.56,26.34,3.18,6.81,0.4,0.72,27.34, Type="AvecGite")
methode_larsson(6.3,1025,9.81,11.9,6.47,0.56,26.34,3.18,6.81,0.4,0.72,27.34, Type="Visqueux")

