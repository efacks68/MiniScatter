####physListCompare.py
###Load in MCNPX
import ROOT
import numpy as np
import csv

#Read in particle density map
def decode(path,file,lenx,leny,ind):
    import csv
    from numpy import zeros
    Img = zeros((leny,lenx))
    Error = zeros((leny,lenx))
    j=0
    i=0
    if file[-1] != "v":
        print(file)
        file+=".csv"
        print(path+file)
    with open(path+file,mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        next(csv_reader)
        for row in csv_reader:
            #print(j,i,ind,row)
            Img[j,i] = float(row[ind])
            Error[j,i] = float(row[ind+1])
            #print(j,i,Img[j,i])
            i+=1
            if i == lenx:
                i=0
                j+=1
            if j==leny+1:
                print(j,i,Img[j,i])
                print("ending",j)
                break
    csv_file.close()
    print(Img.shape)
    return Img,Error

#path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/examples/EricScripts/"
path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
file = "MCNP_ProtonDensity_570MeV"
filename = "MCNP_ProtonDensityMap_570MeV"
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
ind = 2
lenx=77
lowx = -76
highx = 76
lowy = -200
highy = 200
leny = 201
Img,Error = decode(path,file,lenx,leny,ind)
print(np.shape(Img))

#Convert to ROOT
MC_BERT = ROOT.TH2D("protonDensity","Scattered Pencil Beam at Target;Horizontal [mm];Vertical [mm]",lenx,lowx,highx,leny,lowy,highy)
#Fill Image map with 2D histogram values
#Help from Kyrre's converter function and https://root-forum.cern.ch/t/create-and-fill-a-th2f-using-arrays-for-binning/27161/2
for xi in range(lenx):
    for yi in range(leny):
        #bx = hIn.GetBin(xi+1,yi+1)
        #MC_BERT.GetBin(yi,xi) = Img[yi,xi]
        bin = MC_BERT.GetBin(xi+1,yi+1)
        MC_BERT.SetBinContent(bin,Img[yi,xi])
        MC_BERT.SetBinError(bin,Error[yi,xi])
        #if xi == 38 and yi == 100:
            #print(Img[yi,xi],MC_BERT.GetBinContent(MC_BERT.GetBin(xi+1,yi+1)))

#c1 = ROOT.TCanvas()
#MC_BERT.SetOption("COLZ")
#c1.SetLogz()
#MC_BERT.GetXaxis().SetRangeUser(-200,200)
#MC_BERT.Draw()
#c1.Print(picpath+filename+"_ROOT2D.png")




###Load in all scattered beams of different phys lists
from sys import path as sysPath
from os import chdir
MiniScatter_path="../../MiniScatter/build/."
sysPath.append(MiniScatter_path)
chdir("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/")
import miniScatterDriver
print("else GaussFit")
paths = {"scratchPath":"/scratch2/ericdf/PBWScatter/"}
Nparts = 1e5
TRYLOAD = True 
G4_QBERTZ = {'PHYS': 'QGSP_BERT_EMZ', 'ZOFFSET': '*-2', 'WORLDSIZE': 1000.0, 'QUICKMODE': False, 'MINIROOT': True, 
                    'ANASCATTER': True, 'EDEP_DZ': 1.0, 'CUTOFF_RADIUS': 1000.0, 'CUTOFF_ENERGYFRACTION': 0.99, 'POSLIM': 1000.0, 
                    'DIST': [3565], 'BEAM': 'proton', 'ENERGY': 570, 'COVAR': (0.0001, 0.15, 0, 0.0001, 0.15, 0), 'N': Nparts, 
                    'THICK': 0.0, 'MAGNET': [{'type': 'PBW', 'length': 0, 'gradient': 0.0, 'keyval': {'material': 'G4_Al', 'radius': 88.0, 
                                                        'al1Thick': 1.0, 'waterThick': 2.0, 'al2Thick': 1.25}, 'pos': 24.125}], 
                    'OUTFOLDER': '/scratch2/ericdf/PBWScatter/ESS/'}
outname = "PBW_{:.0f}MeV_eX{:.2f},eY{:.2f}um_bX{:.2f},bY{:.2f}m_aX{:.2f},aY{:.2f}_N{:.0e}".format(G4_QBERTZ['ENERGY'],G4_QBERTZ['COVAR'][0],
                G4_QBERTZ['COVAR'][3],G4_QBERTZ['COVAR'][1],G4_QBERTZ['COVAR'][4],G4_QBERTZ['COVAR'][2],G4_QBERTZ['COVAR'][5],G4_QBERTZ["N"])+"_QBZ"
G4_QBERTZ['OUTNAME'] = outname
(twiss_QBERTZ, numPart_QBERTZ, objects_G4_QBERTZ) = miniScatterDriver.getData_tryLoad(G4_QBERTZ, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212"])

G4_FBERTZ = {'PHYS': 'FTF_BIC_EMZ', 'ZOFFSET': '*-2', 'WORLDSIZE': 1000.0, 'QUICKMODE': False, 'MINIROOT': True, 
                    'ANASCATTER': True, 'EDEP_DZ': 1.0, 'CUTOFF_RADIUS': 1000.0, 'CUTOFF_ENERGYFRACTION': 0.99, 'POSLIM': 1000.0, 
                    'DIST': [3565], 'BEAM': 'proton', 'ENERGY': 570, 'COVAR': (0.0001, 0.15, 0, 0.0001, 0.15, 0), 'N': Nparts, 
                    'THICK': 0.0, 'MAGNET': [{'type': 'PBW', 'length': 0, 'gradient': 0.0, 'keyval': {'material': 'G4_Al', 'radius': 88.0, 
                                                        'al1Thick': 1.0, 'waterThick': 2.0, 'al2Thick': 1.25}, 'pos': 24.125}], 
                    'OUTFOLDER': '/scratch2/ericdf/PBWScatter/ESS/'}
outname = "PBW_{:.0f}MeV_eX{:.2f},eY{:.2f}um_bX{:.2f},bY{:.2f}m_aX{:.2f},aY{:.2f}_N{:.0e}".format(G4_FBERTZ['ENERGY'],G4_FBERTZ['COVAR'][0],
                G4_FBERTZ['COVAR'][3],G4_FBERTZ['COVAR'][1],G4_FBERTZ['COVAR'][4],G4_FBERTZ['COVAR'][2],G4_FBERTZ['COVAR'][5],G4_FBERTZ["N"])+"_FBZ"
G4_FBERTZ['OUTNAME'] = outname
(twiss_G4_FBERTZ, numPart_G4_FBERTZ, objects_G4_FBERTZ) = miniScatterDriver.getData_tryLoad(G4_FBERTZ, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212"])

G4_EMZ = {'PHYS': 'EMonly_EMZ', 'ZOFFSET': '*-2', 'WORLDSIZE': 1000.0, 'QUICKMODE': False, 'MINIROOT': True, 
                    'ANASCATTER': True, 'EDEP_DZ': 1.0, 'CUTOFF_RADIUS': 1000.0, 'CUTOFF_ENERGYFRACTION': 0.99, 'POSLIM': 1000.0, 
                    'DIST': [3565], 'BEAM': 'proton', 'ENERGY': 570, 'COVAR': (0.0001, 0.15, 0, 0.0001, 0.15, 0), 'N': Nparts, 
                    'THICK': 0.0, 'MAGNET': [{'type': 'PBW', 'length': 0, 'gradient': 0.0, 'keyval': {'material': 'G4_Al', 'radius': 88.0, 
                                                        'al1Thick': 1.0, 'waterThick': 2.0, 'al2Thick': 1.25}, 'pos': 24.125}], 
                    'OUTFOLDER': '/scratch2/ericdf/PBWScatter/ESS/'}
outname = "PBW_{:.0f}MeV_eX{:.2f},eY{:.2f}um_bX{:.2f},bY{:.2f}m_aX{:.2f},aY{:.2f}_N{:.0e}".format(G4_EMZ['ENERGY'],G4_EMZ['COVAR'][0],
                G4_EMZ['COVAR'][3],G4_EMZ['COVAR'][1],G4_EMZ['COVAR'][4],G4_EMZ['COVAR'][2],G4_EMZ['COVAR'][5],G4_EMZ["N"])+G4_EMZ['PHYS']
G4_EMZ['OUTNAME'] = outname
(twiss_G4_EMZ, numPart_G4_EMZ, objects_G4_EMZ) = miniScatterDriver.getData_tryLoad(G4_EMZ, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212"])

###Plot them all together with different colors
#make projections, normalized.
width=100
integral = MC_BERT.Integral(MC_BERT.GetXaxis().FindBin(-100),MC_BERT.GetXaxis().FindBin(100),MC_BERT.GetYaxis().FindBin(-100),MC_BERT.GetYaxis().FindBin(100))
MC_BERT.Scale(1/integral)
sum = MC_BERT.Integral()
proj_MC_BERT = MC_BERT.ProjectionY("y",MC_BERT.GetXaxis().FindBin(-width),MC_BERT.GetXaxis().FindBin(width),"e")

hist = objects_G4_QBERTZ["tracker_cutoff_xy_PDG2212"]
integral = hist.Integral(hist.GetXaxis().FindBin(-100),hist.GetXaxis().FindBin(100),hist.GetYaxis().FindBin(-100),hist.GetYaxis().FindBin(100))
hist.Scale(1/integral)
sum = hist.Integral()
proj_G4_QBERTZ = hist.ProjectionY("y",hist.GetXaxis().FindBin(-width),hist.GetXaxis().FindBin(width),"e")

hist = objects_G4_FBERTZ["tracker_cutoff_xy_PDG2212"]
integral = hist.Integral(hist.GetXaxis().FindBin(-100),hist.GetXaxis().FindBin(100),hist.GetYaxis().FindBin(-100),hist.GetYaxis().FindBin(100))
hist.Scale(1/integral)
sum = hist.Integral()
proj_G4_FBERTZ = hist.ProjectionY("y",hist.GetXaxis().FindBin(-width),hist.GetXaxis().FindBin(width),"e")

hist = objects_G4_EMZ["tracker_cutoff_xy_PDG2212"]
integral = hist.Integral(hist.GetXaxis().FindBin(-100),hist.GetXaxis().FindBin(100),hist.GetYaxis().FindBin(-100),hist.GetYaxis().FindBin(100))
hist.Scale(1/integral)
sum = hist.Integral()
proj_G4_EMZ = hist.ProjectionY("y",hist.GetXaxis().FindBin(-width),hist.GetXaxis().FindBin(width),"e")

#now try to draw them together
xlim = 50
#hs = ROOT.THStack("hs","Unstacked Histograms")
c2 = ROOT.TCanvas("Scattered Beams","Scattered Beams",0,0,500*8,300*8)

proj_MC_BERT.Draw()
#proj_MC_BERT.SetFillColor(ROOT.kBlack)
#proj_MC_BERT.SetMarkerStyle(21)
proj_MC_BERT.SetLineColor(ROOT.kBlack)
proj_MC_BERT.GetXaxis().SetRangeUser(-xlim,xlim)
#hs.Add(proj_MC_BERT)

proj_G4_QBERTZ.Draw("SAME")
#proj_G4_QBERTZ.SetFillColor(ROOT.kRed)
#proj_G4_QBERTZ.SetMarkerStyle(21)
proj_G4_QBERTZ.SetLineColor(ROOT.kRed)
proj_G4_QBERTZ.GetXaxis().SetRangeUser(-xlim,xlim)
#hs.Add(proj_G4_QBERTZ)

proj_G4_FBERTZ.Draw("SAME")
#proj_G4_FBERTZ.SetFillColor(ROOT.kGreen)
#proj_G4_FBERTZ.SetMarkerStyle(21)
proj_G4_FBERTZ.SetLineColor(ROOT.kGreen)
proj_G4_FBERTZ.GetXaxis().SetRangeUser(-xlim,xlim)
#hs.Add(proj_G4_FBERTZ)

proj_G4_EMZ.Draw("SAME")
#proj_G4_EMZ.SetFillColor(ROOT.kBlue)
#proj_G4_EMZ.SetMarkerStyle(21)
proj_G4_EMZ.SetLineColor(ROOT.kBlue)
proj_G4_EMZ.GetXaxis().SetRangeUser(-xlim,xlim)
#hs.Add(proj_G4_EMZ)

#c2 = ROOT.TCanvas("Scattered Beams","Scattered Beams",0,0,500*8,300*8)
#hs.Draw("nostack")
#hs.GetXaxis().SetRangeUser(-xlim,xlim)
c2.Print(picpath+filename+"BeamOverlap.png") #only one is being drawn.

###Renormalize them to be 1
###Export that to a histogram