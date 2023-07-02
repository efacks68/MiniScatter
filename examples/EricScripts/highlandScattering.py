###highlandScattering.py
##Simple model for getting distribution of particles at Target with Highland Scattering only
import ROOT
import numpy as np
from numpy.random import default_rng
from os.path import isfile
scratchPath = "/scratch2/ericdf/PBWScatter/"
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"

##Set particle characteristic values
nParticles = 1e7
energy = 570 #[MeV]
partA = 938.27209 #[MeV/c2]
partZ = 1
gamma_rel = (energy+partA)/partA #from PrintTwissParameters
beta_rel = np.sqrt(gamma_rel*gamma_rel - 1)/gamma_rel

##Get Highland Equation Radiation Length Calculation variables
p = np.sqrt((energy+partA)**2 - (partA)**2) #[MeV/c] #derived with Kyrre 15.6.22
#pAppleby = beta_rel * gamma_rel * partA / c #Appleby2022 Eqn 2.16
betap = beta_rel*p #Eq 5
radLenAl = 88.97 #[mm]
radLenH2O = 360.8 #[mm] liquid Water

##Need to account for Y bending!
thickAl1 = 1.0
thickH2O = 2.0
thickAl2 = 1.25
L=3565 #thickAl1+thickH2O+thickAl2

##Calculate sigmas for each layer
hSigAl1 = 13.6 * partZ / betap * np.sqrt(thickAl1/radLenAl) * (1 + 0.038 * np.log(thickAl1/radLenAl))
hSigH2O = 13.6 * partZ / betap * np.sqrt(thickH2O/radLenAl) * (1 + 0.038 * np.log(thickH2O/radLenAl))
hSigAl2 = 13.6 * partZ / betap * np.sqrt(thickAl2/radLenAl) * (1 + 0.038 * np.log(thickAl2/radLenAl))

print("Highland Angles: Al1: {:.3e}, H2O: {:.3e}, Al2: {:.3e}, Total: {:.3e}".format(hSigAl1,hSigH2O,hSigAl2,(hSigAl1+hSigH2O+hSigAl2)))

###If file present, load it and print 
filename = "HighlandScatteredDistribution"
rootFile="/scratch2/ericdf/PBWScatter/"+filename+"_{:.0e}".format(nParticles)+".root"
if isfile(rootFile):
    print("Loading data from rootFile!")
    myFile = ROOT.TFile.Open(rootFile,"r")
    hTarget = myFile.Get("Highland Scattered Proton Distribution").Clone("Highland Scattered Proton Distribution")

else:
###If file not present, make it
##Make ROOT histogram to hold
    print("Making New Data")
    width=200
    hTarget = ROOT.TH1D("haTarget","Proton Distribution at Target",4000,-width,width)

    ##Make pencil beam of particles
    ##Send through sheets
    random = default_rng()
    for i in range(int(nParticles)):
        xp1 = random.normal(scale=hSigAl1)
        xp2 = random.normal(scale=hSigH2O)
        xp3 = random.normal(scale=hSigAl2)
        xTarget = L*(xp1+xp2+xp3)
        ##Fill in ROOT histogram
        hTarget.Fill(xTarget)
                    
    integral = hTarget.Integral(hTarget.GetXaxis().FindBin(-width),hTarget.GetXaxis().FindBin(width))
    hTarget.Scale(1/integral)
    sumPri = hTarget.Integral()

    ##Add writing to file
    myFile = ROOT.TFile.Open("/scratch2/ericdf/PBWScatter/"+filename+"_{:.0e}".format(nParticles)+".root","RECREATE")
    print("/scratch2/ericdf/PBWScatter/"+filename+"_{:.0e}".format(nParticles)+".root")
    myFile.WriteObject(hTarget,"Highland Scattered Proton Distribution")

##Fit distribution
maxim = 50
f1 = ROOT.TF1('Highland','gaus',-maxim,maxim)
f1.SetNpx(5000)
f1_res = hTarget.Fit(f1, 'RSQ')

##Show distribution
c1 = ROOT.TCanvas("Highland Scattered Proton Distribution","Highland Scattered Proton Distribution",0,0,500*8,400*8)
ROOT.gStyle.SetOptStat("iRMe")
ROOT.gStyle.SetOptFit(101)
c1.SetLogy()
f1.Draw()
hTarget.Draw()
hTarget.GetXaxis().SetRangeUser(-maxim,maxim)
c1.Print(picpath+"HighlandScattering.png")

