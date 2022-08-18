#MCSFormalism.py
import numpy as np
import sys

def calcEq8(thetasq,Twiss,thick,beta_rel,gamma_rel):
  m=1
  #Twiss=[beta,alpha,gemt,gamma]
  #Calculations from Eq 7 and 8 from Twiss MCS Formalism Calculations from https://cds.cern.ch/record/499590
  e8dgem = 0.5 * Twiss[0]*m * thetasq * thetasq #[m*rad^2]
  e8alph = Twiss[2] * Twiss[1] / (Twiss[2] + e8dgem)
  e8beta = Twiss[2] * Twiss[0]*m / (Twiss[2] + e8dgem) #[m]
  e8gamma = (Twiss[2] * Twiss[3] + thetasq * thetasq ) / (Twiss[2] + e8dgem) #[m^-1]
  e8gemt = Twiss[2] + e8dgem
  #print("e8",thetasq,e8beta,e8alph,e8gemt,e8gamma)
  return [e8beta,e8alph,e8gemt,e8gamma]

def getMoments(Twiss):
  mm=1e-3
  #[beta,alpha,gemt,gamma]
  sigma = np.sqrt(Twiss[2]*Twiss[0]) /mm #sigma x rms = beam size = sqrt(beta*epsG)
  sigmapx = np.sqrt(Twiss[2]*Twiss[3]) #mrad actually
  #print(sigma,sigmapx,sigpx)
  return [sigma,sigmapx]

def toTarget(Twiss,label):
  #Extend the distributions to the Target Location
  #Twiss=[beta,alpha,gemt,gamma]
  PBWexitBetaMx = np.array([[Twiss[0],-Twiss[1]],[-Twiss[1],Twiss[3]]])

  d_PBW_Targ = 5
  drift_PBW_Targ = np.array([ [1, d_PBW_Targ],[ 0, 1]])
  Calc_PBW_Targ = np.linalg.multi_dot([drift_PBW_Targ,PBWexitBetaMx,np.transpose(drift_PBW_Targ)])
  #print(label,"PBWexit:",PBWexitBetaMx,"\n",drift_PBW_Targ,"\n",Calc_PBW_Targ)
 
  return [Calc_PBW_Targ[0][0],-Calc_PBW_Targ[1][0],Twiss[2],Calc_PBW_Targ[1][1]]

radLenAl = 88.97 #[mm] Al
radLenH2O = 360.8 #[mm] liquid Water
partZ = 1 #[C]
partA = 938.27209 #[MeV/c2]
um = 1e-6 #[m]

Inemtx = 0.11315 #[um-mrad]
Ibetax = 141.25 #[m]
Ialphx = -8.881
Inemty = 0.12155 #[um-mrad]
Ibetay = 920.03 #[m]
Ialphy = -57.46

print("Enter Twiss as arguments, beta, alpha, normalized emittance")
if len(sys.argv) > 1:
  Ibetax = float(sys.argv[1])
  Ialphx = float(sys.argv[2])
  Inemtx = float(sys.argv[3])
  Ibetay = float(sys.argv[4])
  Ialphy = float(sys.argv[5])
  Inemty = float(sys.argv[6])

energy = 570
#energy = float(input("What is the beam Energy? "))
gamma_rel = 1 + energy/partA #from PrintTwissParameters
beta_rel = np.sqrt(gamma_rel*gamma_rel -1 )/gamma_rel

#Get rid of Normalized Emittance!
Igemtx = Inemtx/(beta_rel*gamma_rel)
Igemty = Inemty/(beta_rel*gamma_rel)

#Make Twiss array = [beta,alpha,gemt,gamma]
TwissIx   = [Ibetax,Ialphx,Igemtx*um,((1+Ialphx*Ialphx)/Ibetax)] #Initial Twiss
TwissIy   = [Ibetay,Ialphy,Igemty*um,((1+Ialphy*Ialphy)/Ibetay)]
print(TwissIx,"\n",TwissIy)

##Highland Equation Radiation Length Calculations
p = np.sqrt((energy+partA)**2 - (partA)**2) #[MeV/c] #derived with Kyrre 15.6.22
betap = beta_rel*p #Eq 5

m1Len = 1.0
m2Len = 2.0
m3Len = 1.25
#m1Len = float(input("What is the first Al plate thickness? "))
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
#print(Twisse8x)

e8TargxReal = toTarget(Twisse8x,"e8XReal")
e8TargyReal = toTarget(Twisse8y,"e8YReal")
print("Twiss array = [beta [m],alpha,gemt [m],gamma [?]]\n",e8TargxReal)

Mex = getMoments(e8TargxReal)
Mey = getMoments(e8TargyReal)
print("The beam size is x {:.2f}mm, y {:.2f}mm".format(Mex[0],Mey[0]))