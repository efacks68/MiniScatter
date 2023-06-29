###convolve.py
###convolve the pencil fit with the erf function of the raster
##load in pencil fit function

##Gaussian function for 
def gaussian(x,amplitude,mu,sigma):
    from numpy import exp as npExp
    return amplitude * npExp( - (x - mu) ** 2 / (2 * sigma ** 2))

##Lorentzian function, deviates from ROOT as "lambda" seems to be reserved in Python
def lorentz(x,alpha,gamma): 
    return alpha / ( x * x + gamma * gamma)

#from sys import path as sysPath
#from os import chdir
#MiniScatter_path="../../MiniScatter/build/."
#sysPath.append(MiniScatter_path)
#chdir("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/")
#import miniScatterDriver
#print("else GaussFit")
paths = {"scratchPath":"/scratch2/ericdf/PBWScatter/"}
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
#Nparts = 1e5
#TRYLOAD = True 
#G4_FBERTZ = {'PHYS': 'FTF_BIC_EMZ', 'ZOFFSET': '*-2', 'WORLDSIZE': 1000.0, 'QUICKMODE': False, 'MINIROOT': True, 
#                    'ANASCATTER': True, 'EDEP_DZ': 1.0, 'CUTOFF_RADIUS': 1000.0, 'CUTOFF_ENERGYFRACTION': 0.99, 'POSLIM': 1000.0, 
#                    'DIST': [3565], 'BEAM': 'proton', 'ENERGY': 570, 'COVAR': (0.0001, 0.15, 0, 0.0001, 0.15, 0), 'N': Nparts, 
#                    'THICK': 0.0, 'MAGNET': [{'type': 'PBW', 'length': 0, 'gradient': 0.0, 'keyval': {'material': 'G4_Al', 'radius': 88.0, 
#                                                        'al1Thick': 1.0, 'waterThick': 2.0, 'al2Thick': 1.25}, 'pos': 24.125}], 
#                    'OUTFOLDER': '/scratch2/ericdf/PBWScatter/ESS/'}
#outname = "PBW_{:.0f}MeV_eX{:.2f},eY{:.2f}um_bX{:.2f},bY{:.2f}m_aX{:.2f},aY{:.2f}_N{:.0e}".format(G4_FBERTZ['ENERGY'],G4_FBERTZ['COVAR'][0],
#                G4_FBERTZ['COVAR'][3],G4_FBERTZ['COVAR'][1],G4_FBERTZ['COVAR'][4],G4_FBERTZ['COVAR'][2],G4_FBERTZ['COVAR'][5],G4_FBERTZ["N"])+"_FBZ"
#G4_FBERTZ['OUTNAME'] = outname
#(twiss_G4_FBERTZ, numPart_G4_FBERTZ, objects_G4_FBERTZ) = miniScatterDriver.getData_tryLoad(G4_FBERTZ, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212"])
import ROOT
import scipy.special
import numpy as np
import matplotlib.pyplot as plt
from plotFit import gaussianFit,converter,fitGaussians
EMZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08EMonly_EMZ.root"
fEMZ = ROOT.TFile(EMZ_file,"r")
G4_EMZ_TH2 = fEMZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
nGaus=4
Lorentz=True
_,_,coeffsEMy,_,_,_ = gaussianFit(G4_EMZ_TH2,"y",100,500,picpath+"convolve_G4_EMZ",2,25,True,True,nGaus,Lorentz)
print(coeffsEMy)

##Open file with ROOT
QBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
fQBERTZ = ROOT.TFile(QBERTZ_file,"r")
##Get TH2D and clone in
G4_QBERTZ_TH2 = fQBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
#Send to fit function, returns coefficients
_,_,coeffsQBERTZy,_,_,_ = gaussianFit(G4_QBERTZ_TH2,"y",5,500,picpath+"physListComp_G4_QBERTZ",2,25,True,True,nGaus,True)


##Make Numpy map out of TH2D
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--Nb",        type=int,   default=100,help="Number of macroparticles per beamlet. Default=10")
parser.add_argument("--nP",        type=float, default=1e3,help="Numper of beamlets in pulse. Default=1e3")
parser.add_argument("--reBin",     type=int,   default=1,  help="Number of bins to make into 1 in 2D histogram for smoothing")
parser.add_argument("--Nbeamlet",  type=float, default=5e8)
parser.add_argument("--sim",       type=str,   default="beamlet")
args = parser.parse_args()
(ImgGeant4, xax, yax) = converter(G4_EMZ_TH2,False,picpath+"concolveTry",paths,args,False) #convert from TH2D to numpy map
(ImgGeant4_QBERTZ, xax, yax) = converter(G4_QBERTZ_TH2,False,picpath+"concolveTry",paths,args,False) #convert from TH2D to numpy map
Error=np.ones(np.shape(ImgGeant4))
coeffsNP = fitGaussians(ImgGeant4,Error,picpath,100,"y",4)

##Make projection of numpy map
xlim=300
a=1
lenG = np.shape(ImgGeant4)[0]
projG = ImgGeant4[round(lenG/2)-a:round(lenG/2)+a,:]
projG = (projG[0,:] + projG[1,:]) #/ (a*2+1)
g = np.linspace(-xlim,xlim,len(projG))

#projG = projG[round(lenG/2)-round(len(projM)/2):round(lenG/2)+round(len(projM)/2)]/0.04
lenQ = np.shape(ImgGeant4)[0]
projQ = ImgGeant4_QBERTZ[round(lenQ/2)-a:round(lenQ/2)+a,:]
projQ = (projQ[0,:] + projQ[1,:]) #/ (a*2+1)
q = np.linspace(-xlim,xlim,len(projQ))

lenx= 10000
x = np.linspace(-xlim,xlim,lenx)
a = np.linspace(-xlim,xlim,2*lenx)
#b = np.linspace(-300,300,3*lenx-1)##used for the convolution if used mode='full'
ax = 48.7867
offsetx = 0
ay = 16.3991
offsety = 0
sigmax = np.sqrt(0.118980737408497*1085.63306926394)
sigmay = np.sqrt(0.123632934174567*136.062409365455)
#y = gaussian(x,coeffsEMy[2],0,coeffsEMy[3])+gaussian(x,coeffsEMy[4],0,coeffsEMy[5])+gaussian(x,coeffsEMy[6],0,coeffsEMy[7])+lorentz(x,coeffsEMy[0],coeffsEMy[1])
y = gaussian(x,coeffsNP[0],0,coeffsNP[1])+gaussian(x,coeffsNP[2],0,coeffsNP[3])+gaussian(x,coeffsNP[4],0,coeffsNP[5])
z = 0.01 * (scipy.special.erf((a+ax-offsetx) /sigmax/np.sqrt(2)) -   scipy.special.erf((a-ax-offsetx) /sigmax/np.sqrt(2))) 
##convolve
c = 10 * np.convolve(y,z,mode='same')
###compare the convolved with the actual raster
##plot on top of each other
plt.plot(x,y,label="Fit Pencil")
plt.plot(a,z,label="Reference Raster")
plt.plot(a,c,label="Convolution")
#plt.plot(g,projG,label="Pencil Proj")
#plt.plot(q,projQ,label="Pencil Proj QBERTZ")

plt.yscale("log")
plt.ylim(1e-10,2e-1)
#plt.xlim(-50,50)
plt.grid(True)
plt.legend(loc="upper left")
plt.title("Scattered Beams")
plt.xlabel("Horizontal [mm]")
plt.ylabel(r"Density [mm$^{-2}$]")
plt.savefig(picpath+"npConvolve_G4_EMZ.png")
print(picpath+"npConvolve_G4_EMZ.png")



"""
##plot the fit
maxim=200
if nGaus == 1:
    f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1]))',-maxim,maxim)
elif nGaus ==2:
    if Lorentz:
        f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] / (x**2 + [3]**2) + [4]',-maxim,maxim) #[2] / (x**2 + [3]**2) #[2] / x ** 3 + [3]
    else:
        f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3]))',-maxim,maxim)
elif nGaus == 3:
    if Lorentz:
        f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] / (x**2 + [3]**2) + [4] * exp(-x*x/(2*[5]*[5])) ',-maxim,maxim)
    else:
        f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3])) + [4] * exp(-x*x/(2*[5]*[5])) ',-maxim,maxim)
elif nGaus == 4:
    if Lorentz:
        f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] / (x**2 + [3]**2) + [4] * exp(-x*x/(2*[5]*[5])) + [6] * exp(-x*x/(2*[7]*[7]))',-maxim,maxim)
        f2.SetParameter(0,coeffsEMy[2]); f2.SetParameter(1,coeffsEMy[3]); f2.SetParameter(2,coeffsEMy[0]); f2.SetParameter(3,coeffsEMy[1]); f2.SetParameter(4,coeffsEMy[4])
        f2.SetParameter(5,coeffsEMy[5]); f2.SetParameter(6,coeffsEMy[6]); f2.SetParameter(7,coeffsEMy[7])
    else:
        f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3])) + [4] * exp(-x*x/(2*[5]*[5])) + [6] * exp(-x*x/(2*[7]*[7]))',-maxim,maxim)
elif nGaus == 5:
    if Lorentz:
        f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] / (x**2 + [3]**2) + [4] * exp(-x*x/(2*[5]*[5])) + [6] * exp(-x*x/(2*[7]*[7])) + [8] * exp(-x*x/(2*[9]*[9]))',-maxim,maxim)
    else:
        f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3])) + [4] * exp(-x*x/(2*[5]*[5])) + [6] * exp(-x*x/(2*[7]*[7])) + [8] * exp(-x*x/(2*[9]*[9]))',-maxim,maxim)

f2.SetNpx(10000)

fRaster = ROOT.TF1("fRaster"," 0.5 * (erf( ( x+[0]-0) /[1]/sqrt(2)) - erf( ( x-[0]-0) /[1]/sqrt(2))) * 0.5 * (erf( ( y+[2]-0) /[3]/sqrt(2)) - erf( ( Y-[2]-0) /[3]/sqrt(2)))")
fRaster.SetParameter(0,48.7867); fRaster.SetParameter(1,14); fRaster.SetParameter(2,16.3991); fRaster.SetParameter(3,5.3)
fRaster.SetNpx(10000)

c1 = ROOT.TCanvas("Scattered Beams","Scattered Beams",0,0,300*8,200*8)
c1.SetLogy()
#f2.Draw()
fRaster.Draw()
c1.Print(picpath+"convolveFunc_G4_EMZ.png")
"""

