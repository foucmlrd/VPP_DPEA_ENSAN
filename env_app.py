#Objectif : Calculer la vitesse et l'angle du vent apparant

import math as math
from input import *
from algo import *

TWSx = TWS_ms*math.cos(TWA_rad) #(m/s) True wind speed en x 

TWSy = TWS_ms*math.sin(TWA_rad) #(m/s) True wind speed en x 


AWSx = Vb_ms + TWSx #(m/s) Apparent wind speed en x

AWSy = TWSy*math.cos(Heeling_rad) #(m/s) Apparent wind speed en y

AWS = pow(pow(AWSx,2) + pow(AWSy,2),1/2) #(m/s) Apparent wind speed

AWA = math.acos(AWSx/AWS) #(rad) Apparent wind angle

'''VALEURS DE SORTIE
PASSAGE DES DONNEES EN DEGRE ET EN NOEUDS
'''

AWS_knts = mstoknot(AWS) #Apparent wind speed en noeuds

AWA_deg = radtodeg(AWA) #Apparent wind angle en degré
