###convolve.py
##improvement 16 Aug
###########
##Plan:
## 1. Load in Pencil, beamlet and Rastered beams
## 2. Fit Pencil with 1,3,5 Gaussians
## 3. Plot all together
## 4. Convolve Pencil fit with Raster function
###########

##Function for fitting projections to n Gaussians
def fit_projection(proj,nGauss,p0,p2):
    r=0.01
    ##Set function
    if nGauss == 2:
        func = ROOT.TF1('func','([0] / sqrt(2*pi*[1]) ) * exp(-x*x/(2*[1]*[1])) + ([2] / sqrt(2*pi*[3]) ) * exp(-x*x/(2*[3]*[3]))',-xlim,xlim)
    elif nGauss == 3:
        func = ROOT.TF1('func','([0] / sqrt(2*pi*[1]) ) * exp(-x*x/(2*[1]*[1])) + ([2] / sqrt(2*pi*[3]) ) * exp(-x*x/(2*[3]*[3])) + ([4] / sqrt(2*pi*[5]) ) * exp(-x*x/(2*[5]*[5]))',-xlim,xlim)
    elif nGauss == 4:
        func = ROOT.TF1('func','([0] / sqrt(2*pi*[1]) ) * exp(-x*x/(2*[1]*[1])) + ([2] / sqrt(2*pi*[3]) ) * exp(-x*x/(2*[3]*[3])) + ([4] / sqrt(2*pi*[5]) ) * exp(-x*x/(2*[5]*[5])) + ([6] / sqrt(2*pi*[7]) ) * exp(-x*x/(2*[7]*[7]))',-xlim,xlim)
    elif nGauss == 5:
        func = ROOT.TF1('func','([0] / sqrt(2*pi*[1]) ) * exp(-x*x/(2*[1]*[1])) + ([2] / sqrt(2*pi*[3]) ) * exp(-x*x/(2*[3]*[3])) + ([4] / sqrt(2*pi*[5]) ) * exp(-x*x/(2*[5]*[5])) + ([6] / sqrt(2*pi*[7]) ) * exp(-x*x/(2*[7]*[7])) + ([8] / sqrt(2*pi*[9]) ) * exp(-x*x/(2*[9]*[9]))',-xlim,xlim)

    ##Set parameters. Setting more than needed for a function does nothing
    func.SetParameter(0,p0)                 #A1
    func.SetParLimits(0,1e-4,p0*100)        #A1
    func.SetParameter(1,p2)                 #sigma1
    func.SetParLimits(1,p2,p2*(1+r))        #sigma1

    func.SetParameter(2,p0*(r**2))          #A2
    func.SetParLimits(2,0, p0*5)            #A2
    func.SetParameter(3,p2*2)               #sigma2
    func.SetParLimits(3,p2*(1+2*r), p2*2.6) #sigma2

    func.SetParameter(4,p0*r**2)            #A3
    func.SetParLimits(4,0, p0*5)            #A3
    func.SetParameter(5,p2*7)               #sigma3
    func.SetParLimits(5,p2*2.5, p2*12)      #sigma3

    func.SetParameter(6,p0*r**4)            #4
    func.SetParLimits(6,0, p0)              #A4
    func.SetParameter(7,p2*10)              #sigma4
    func.SetParLimits(7,p2*8, p2*23)        #sigma4

    func.SetParameter(8,p0*r**5)            #A5
    func.SetParLimits(8,0, p0)              #A5
    func.SetParameter(9,p2*25)              #sigma5
    func.SetParLimits(9,p2*15, p2*30)       #sigma5

    ##Number of fit points and Fit. Return function and fit results
    func.SetNpx(10000)
    func_res = proj.Fit(func,'RSQ')
    return func, func_res

import ROOT
import numpy as np
#useful variables
paths = {"scratchPath":"/scratch2/ericdf/PBWScatter/"}
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
filename="convolutions"
##plot settings
nGauss=5
xlim = 60
width=10
ylim = [1e-9,1.5e2]
leg = [0.13,0.6,0.4,0.9]
Pencil=True
Raster=True
ReferenceRaster=True
Beamlet=True
Convolution=True

##Create Canvas and Legend
canvas = ROOT.TCanvas("Scattered Beams","Scattered Beams",0,0,400*8,250*8)
leg = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])
ROOT.gPad.SetLogy()
ROOT.gStyle.SetLegendTextSize(0.04)
leg.SetMargin(0.12)

##Add Pencil Slice
if Pencil:
    ##Open file with ROOT
    QBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
    fQBERTZ = ROOT.TFile(QBERTZ_file,"r")
    ##Get TH2D and clone in
    G4_QBERTZ_TH2 = fQBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("QBZ")
    ##Send to fit function, returns coefficients
    #_,_,coeffsQBERTZy,_,_,_ = gaussianFit(G4_QBERTZ_TH2,"y",5,500,picpath+"physListComp_G4_QBERTZ",2,25,True,True,nGauss,True)
    sliceY_G4_QBERTZ = G4_QBERTZ_TH2.ProjectionY("Y",G4_QBERTZ_TH2.GetXaxis().FindBin(-width),G4_QBERTZ_TH2.GetXaxis().FindBin(width),"e")
    sliceY_G4_QBERTZ.SetName("QBERTZ Y Slice")
    integral = sliceY_G4_QBERTZ.Integral(sliceY_G4_QBERTZ.GetXaxis().FindBin(-xlim),sliceY_G4_QBERTZ.GetXaxis().FindBin(xlim))
    sliceY_G4_QBERTZ.Scale(1/integral)
    sum = sliceY_G4_QBERTZ.Integral()
    #print("Sum:",sum)

    ##Draw slice
    sliceY_G4_QBERTZ.Draw()
    sliceY_G4_QBERTZ.GetXaxis().SetRangeUser(-xlim,xlim)
    sliceY_G4_QBERTZ.GetYaxis().SetRangeUser(ylim[0],ylim[1])
    sliceY_G4_QBERTZ.SetMarkerStyle(34)
    sliceY_G4_QBERTZ.SetMarkerColor(1)
    sliceY_G4_QBERTZ.SetMarkerSize(2)
    sliceY_G4_QBERTZ.SetStats(False)
    sliceY_G4_QBERTZ.SetTitle("Scattered Pencil Beam at Target")
    leg.AddEntry(sliceY_G4_QBERTZ,"Y Slice")

    ##Fit N Gaussians
    a1y = ROOT.TF1('a1','gaus',-xlim,xlim)
    a1y.SetNpx(10000)
    a1y_res = sliceY_G4_QBERTZ.Fit(a1y, 'RSQ')
    p0 = a1y.GetParameter(0) #A
    p2 = a1y.GetParameter(2) #sigma
    a1y.Draw("SAME")
    leg.AddEntry(a1y,"1 Gaussian Fit")
    a3y, a3y_res = fit_projection(sliceY_G4_QBERTZ,3,p0,p2)
    a3y.Draw("SAME")
    a3y.SetLineColor(4)
    a3y.SetLineWidth(2)
    leg.AddEntry(a3y,"3 Gaussian Fit")
    a5y, a5y_res = fit_projection(sliceY_G4_QBERTZ,5,p0,p2)
    a5y.Draw("SAME")
    a5y.SetLineColor(3)
    a5y.SetLineWidth(4)
    leg.AddEntry(a5y,"5 Gaussian Fit")

##Add Beamlet Twiss
if Beamlet:
    Twiss = [0.118980737408497,1085.63306926394,-65.1638312921654,0.123632934174567,136.062409365455,-8.12599512314246]
    twissFile = "/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.12,eY0.12um_bX1085.63,bY136.06m_aX-65.16,aY-8.13_N1e+08_miniR_QBZ.root"
    fTwiss = ROOT.TFile(twissFile,"r")
    G4_Twiss = fTwiss.Get("tracker_cutoff_xy_PDG2212").Clone("Twiss")
    slice_Twiss = G4_Twiss.ProjectionY("Y",G4_Twiss.GetXaxis().FindBin(-width),G4_Twiss.GetXaxis().FindBin(width),"e")
    slice_Twiss.SetName("QBERTZ Y Slice")
    integral = slice_Twiss.Integral(slice_Twiss.GetXaxis().FindBin(-xlim),slice_Twiss.GetXaxis().FindBin(xlim))
    slice_Twiss.Scale(1/integral)
    sum = slice_Twiss.Integral()
    #print("Twiss Sum:",sum)
    slice_Twiss.Draw("SAME")
    slice_Twiss.SetMarkerColor(4)
    slice_Twiss.SetMarkerStyle(41)
    slice_Twiss.SetMarkerSize(1)
    leg.AddEntry(slice_Twiss,"Twiss Scattered")

##Add Convolution of Beamlet and Pencil? No, should be pencil and convolution with covariance? Or something like that...

##Add Reference Raster Slice
if ReferenceRaster:
    #ax = 48.7867;   ay = 16.3991
    #sigmax = 11.365272681993753;   sigmay = 4.101438150297074
    refRaster = ROOT.TF1("RefRaster", "0.5 * (TMath::Erf((x+16.3991)/4.101/sqrt(2)) - TMath::Erf((x-16.3991)/4.101/sqrt(2)))", -xlim, xlim)
    refRaster.Draw("SAME")
    refRaster.SetLineColor(6)
    leg.AddEntry(refRaster,"Reference Raster")

##Add Raster Slice
if Raster:
    rasterFile = "/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ.root"
    fRaster = ROOT.TFile(rasterFile,"r")
    G4_Raster = fRaster.Get("tracker_cutoff_xy_PDG2212").Clone("Raster")
    slice_Raster = G4_Raster.ProjectionY("Y",G4_Raster.GetXaxis().FindBin(-width),G4_Raster.GetXaxis().FindBin(width),"e")
    slice_Raster.SetName("QBERTZ Y Slice")
    integral = slice_Raster.Integral(slice_Raster.GetXaxis().FindBin(-xlim),slice_Raster.GetXaxis().FindBin(xlim))
    slice_Raster.Scale(1/integral)
    sum = slice_Raster.Integral()
    #print("Raster Sum:",sum)
    slice_Raster.Draw("SAME HIST P") #removes error bars, as low bar was obscuring everything below it
    slice_Raster.SetMarkerColor(6)
    slice_Raster.SetMarkerStyle(41)
    slice_Raster.SetMarkerSize(2)
    leg.AddEntry(slice_Raster,"Raster Scattered")

##Add Convolution of Pencil and Ref Raster
convG1Raster = ROOT.TF1Convolution(a1y,refRaster,-100,100,True) #not showing up...
convG1Raster.SetNofPointsFFT(100000)
convG1Raster.SetRange(-xlim,xlim)
convG1Raster.Draw("SAME")
#convG1Raster.SetMarkerColor(7)
leg.AddEntry(convG1Raster,"1 Gaussian Raster Convolution")

#ROOT.gPad.RedrawAxis("g") #supposed to make a grid, but it doesn't seem to...
leg.Draw("SAME")
canvas.Update()

canvas.Draw()
canvas.Print(picpath+filename+".png")
#canvas.Print(picpath+filename".pdf")




#EMZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08EMonly_EMZ.root"
#fEMZ = ROOT.TFile(EMZ_file,"r")
#G4_EMZ_TH2 = fEMZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
#_,_,coeffsEMy,_,_,_ = gaussianFit(G4_EMZ_TH2,"y",100,500,picpath+"convolve_G4_EMZ",2,25,True,True,nGauss,Lorentz)
#print(coeffsEMy)
"""
import scipy.special
import numpy as np
import matplotlib.pyplot as plt
from plotFit import gaussianFit,converter,fitGaussians
Lorentz=False

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
coeffsNP = fitGaussians(ImgGeant4,Error,picpath,100,"y",nGauss)

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

