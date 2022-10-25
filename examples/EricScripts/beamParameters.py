import numpy as np

energy = float(input("What is the energy of the proton beam in MeV? "))
#partA = input("What is the mass of the particle? ")

partA = 938.27209 #[MeV/c2]
partZ = 1
um=1e-6
mm=1e-3
gamma_rel = 1 + energy/partA #from PrintTwissParameters
beta_rel = np.sqrt(gamma_rel*gamma_rel -1 )/gamma_rel

print("The gamma and beta are: ",gamma_rel,beta_rel)

nemt = float(input("What is the normalized emittance in um? "))

gemt = nemt*um / (beta_rel * gamma_rel)

print("The Geometric Emittance is ",gemt/um," [um]")

beta = float(input("What is the beta in m? "))

sigma  = np.sqrt(gemt * beta)

print("The beam size is ", sigma/mm,"mm")

p = np.sqrt((energy+partA)**2 - (partA)**2) #[MeV/c] #derived with Kyrre 15.6.22
betap = beta_rel*p #Eq 5

radLen = float(input("What is the radiation Length in mm? "))
thick = float(input("what is the thickness of the plate? "))

thetasq = 13.6 * partZ / betap * np.sqrt(thick/radLen) * (1 + 0.038 * np.log(thick/radLen))

print("Your average sqrt(angle) is ", thetasq)