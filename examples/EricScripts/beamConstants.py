import numpy as np
um=1e-6
mm=1e-3

#Twiss from BPM93 calculation
betaX = 1006.80 #[m]
alphX = -60.44
nemtX = 0.11315 #[mm-mrad]
betaY = 129.72 #[m]
alphY = -7.72
nemtY = 0.12155 #[mm-mrad]

partA = 938.27209 #[MeV/c2]
partZ = 1
energy = 570 #[MeV]
gamma_rel = 1 + energy/partA #from PrintTwissParameters
beta_rel = np.sqrt(gamma_rel*gamma_rel -1 )/gamma_rel
gemtX = nemtX*um / (beta_rel * gamma_rel) #[m]
gemtY = nemtY*um / (beta_rel * gamma_rel) #[m]