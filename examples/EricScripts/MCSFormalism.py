#MCSFormalism.py
from plotFit import calcEq8#,calcEq16
import numpy as np

radLenAl = 88.97 #[mm] Al
radLenH2O = 360.8 #[mm] liquid Water
partZ = 1 #[C]?
partA = 938.27209 #[MeV/c2]
um = 1e-6 #[m]

Inemx = 0.11315 #[um-mrad]
Ibetax = 341.25 #[m]
Ialphx = -18.81
Inemy = 0.12155 #[um-mrad]
Ibetay = 320.03 #[m]
Ialphy = -17.46

energy = 570
#energy = float(input("What is the beam Energy? "))
gamma_rel = 1 + energy/partA #from PrintTwissParameters
beta_rel = np.sqrt(gamma_rel*gamma_rel -1 )/gamma_rel

#Get rid of Normalized Emittance!
Igemx = Inemx/(beta_rel*gamma_rel)
Igemy = Inemy/(beta_rel*gamma_rel)
#Make Twiss array
TwissIx   = [Ibetax,Ialphx,Igemx*um,((1+Ialphx*Ialphx)/Ibetax)] #Initial Twiss
TwissIy   = [Ibetay,Ialphy,Igemy*um,((1+Ialphy*Ialphy)/Ibetay)]


##Highland Equation Radiation Length Calculation
p = np.sqrt((energy+partA)**2 - (partA)**2) #[MeV/c] #derived with Kyrre 15.6.22
betap = beta_rel*p #Eq 5

#m1Len = 1.0
m2Len = 2.0
m3Len = 1.25
m1Len = float(input("What is the first Al plate thickness? "))
#m2Len = float(input("What is the water thickness? "))
#m3Len = float(input("What is the second Al plate thickness? "))

#Al Front contribution
thetasqAl1 = 13.6 * partZ / betap * np.sqrt(m1Len/radLenAl) * (1 + 0.038 * np.log(m1Len/radLenAl))
Twisse8xAl1 = calcEq8(thetasqAl1, TwissIx,m1Len,beta_rel,gamma_rel)
Twisse8yAl1 = calcEq8(thetasqAl1, TwissIy,m1Len,beta_rel,gamma_rel)

#H2O contribution
thetasqH2O = 13.6 * partZ / betap * np.sqrt(m2Len/radLenH2O) * (1 + 0.038 * np.log(m2Len/radLenH2O))
Twisse8xH2O = calcEq8(thetasqH2O, Twisse8xAl1,m2Len,beta_rel,gamma_rel)
Twisse8yH2O = calcEq8(thetasqH2O, Twisse8yAl1,m2Len,beta_rel,gamma_rel)

#Al Back contribution
thetasqAl2 = 13.6 * partZ / betap * np.sqrt(m3Len/radLenAl) * (1 + 0.038 * np.log(m3Len/radLenAl))
Twisse8x = calcEq8(thetasqAl2, Twisse8xH2O,m3Len,beta_rel,gamma_rel)
Twisse8y = calcEq8(thetasqAl2, Twisse8yH2O,m3Len,beta_rel,gamma_rel)

print(Twisse8x)

from plotFit import toTarget,getMoments
e8TargxReal = toTarget(Twisse8x,"e8XReal")
e8TargyReal = toTarget(Twisse8y,"e8YReal")

#print(e8TargxReal)
Mex = getMoments(e8TargxReal)
Mey = getMoments(e8TargyReal)
print("The beam size is x {:.2f}mm, y {:.2f}mm".format(Mex[0],Mey[0]))