import math as math
from algo import *

rho_eau = 1025 #kg/m3 Masse volumique de l'eau
g = 9.81 #m/s² constante de pesanteur
rho_air = 1.225 #kg/m3 Masse volumique de l'air
viscosite = 0.000001 #m²/s Viscosité de l'eau de mer

TWS = 8 #Faire une boucle de calcul avec toutes les wind speed et les angles #[4, 8, 10] #True Wind Speed en noeuds

TWS_ms = knottoms(TWS)

TWA = 60 #[40, 60, 80, 100, 120, 140] #True Wind Angle en degré

TWA_rad = degtorad(TWA)



#REMPLIR LE MODE DANS LA FONCTION ALGO

#Remplir les inputs des voiles dans la fonction aero_voile


