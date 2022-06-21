import math as math
VB = 0 #Vitesse bateau en knot

RED = 3

Mode = 0 #INPUT PAR L'UTILISATEUR

Heeling = 0 #(deg) Angle de gite du navire

Derive = 0 #(deg) Angle de dérive du navire

#Transformation d'angle de degres au radians


def degtorad(a):
    return(a*math.pi/180)

'''
Transformation d'angle de radians en degres
'''

def radtodeg(a):
  if a < 0:
    return(a*180/math.pi+180)
  else: 
    return(a*180/math.pi)

def knottoms(a):
  return 1852*a/3600

def mstoknot (a):
  return 3600*a/1852

#Passage en unité SI pour les calculs

Vb_ms = knottoms(VB) #vitesse bateau en m/s

Heeling_rad = radtodeg(Heeling)
Derive_rad = radtodeg(Derive)

