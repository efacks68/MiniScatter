import numpy as np


betax = 941.25 #[m]
alphx = -58.81
nemtx = 0.11315 #[um-mrad]
betay = 120.03 #[m]
alphy = -7.46
nemty = 0.12155 #[um-mrad]

partA = 938.27209 #[MeV/c2]
partZ = 1
um=1e-6
mm=1e-3
energy = 570
gamma_rel = 1 + energy/partA #from PrintTwissParameters
beta_rel = np.sqrt(gamma_rel*gamma_rel -1 )/gamma_rel
gemtx = nemtx*um / (beta_rel * gamma_rel)
gemty = nemty*um / (beta_rel * gamma_rel)