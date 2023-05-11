#MCNPXCompare.py

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
        print(file)
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
    #print(Img.shape)
    return Img,Error

#path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/examples/EricScripts/"
path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
file = "ProtonDensityMaps_570MeV"
filename = "ProtonDensityMap2_570MeV"
ind = 2
lenx=77
lowx = -76
highx = 76
lowy = -200
highy = 200
leny = 201
Img,Error = decode(path,file,lenx,leny,ind)
print(np.shape(Img))
#optional: Visualize PDM or write in better format to CSV
def plot2DMPL(Img,ImgName,picpath):
    from matplotlib.pyplot import subplots,pcolormesh,close,tight_layout,savefig,setp,title,xlabel,ylabel
    from matplotlib.axes import Axes
    from matplotlib.patches import Rectangle
    from matplotlib.colors import LogNorm
    xax = np.linspace(lowx,highx,lenx)
    yax = np.linspace(lowy,highy,leny)
    #print(xax[-1])
    X, Y = np.meshgrid(xax,yax) #Set mesh of plot from the axes given from converter function
    close() #make sure no prior plotting messes up

    fig,ax = subplots(dpi=300)
    #print(datetime.now().strftime("%H-%M-%S"))

    #Set maximum value depending on the maximum current density
    from math import log10,ceil
    maxim = Img.max() + 1e-2
    minim = 10**ceil(log10(1e-2))
    #cbarVals  = [minim,minim*10,minim*100,maxim] #make array for color bar values
    #cbarLabels = ["{:.2f}".format(cbarVals[0]),"{:.1f}".format(cbarVals[1]),"{:.1f}".format(cbarVals[2]),"{:.0f}".format(cbarVals[3])]
    cbarLabel = "Proton Density"
    #print("Max current Density: ",Img.max(),"/",maxim,datetime.now().strftime("%H-%M-%S"))

    #Use pcolor to show density map, use log scale
    c = ax.pcolormesh(X,Y,Img,shading='auto',norm=LogNorm(vmin=minim, vmax=maxim), cmap='viridis',rasterized=True) #viridis or magma are perceptually uniform

    lw=1
    fs=15

    #Set Plot Properties
    xlim=50
    #setp(ax,xlim=([-xlim,xlim]),ylim=([-xlim,xlim]))
    ax.set_title("Scattered Pencil Beam at Target",fontsize=fs+2)
    ax.set_ylabel("Vertical [mm]",fontsize=fs)
    ax.set_xlabel("Horizontal [mm]",fontsize=fs)
    cbar = fig.colorbar(c, ax=ax,pad=0.01)
    cbar.set_label(cbarLabel,fontsize=fs-1)
    #cbar.set_ticks(cbarVals)
    #cbar.set_ticklabels(cbarLabels)
    tight_layout()
    savefig(picpath+ImgName+".png")#
    close(fig)
    close()
    print(picpath+ImgName+".png")    

picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
filename = file+"_2DMCNPXData"
#plotMPL(Img,filename,picpath)

with open(path+filename,mode='w') as csv_file:
    csv_writer = csv.writer(csv_file,delimiter = ',')
    csv_writer.writerows(Img)
    csv_file.close()


#print(Img[100,38])
#Convert to ROOT
hist = ROOT.TH2D("protonDensity","Scattered Pencil Beam at Target;Horizontal [mm];Vertical [mm]",lenx,lowx,highx,leny,lowy,highy)
#Fill Image map with 2D histogram values
#Help from Kyrre's converter function and https://root-forum.cern.ch/t/create-and-fill-a-th2f-using-arrays-for-binning/27161/2
for xi in range(lenx):
    for yi in range(leny):
        #bx = hIn.GetBin(xi+1,yi+1)
        #hist.GetBin(yi,xi) = Img[yi,xi]
        bin = hist.GetBin(xi+1,yi+1)
        hist.SetBinContent(bin,Img[yi,xi])
        hist.SetBinError(bin,Error[yi,xi])
        #if xi == 38 and yi == 100:
            #print(Img[yi,xi],hist.GetBinContent(hist.GetBin(xi+1,yi+1)))

#c1 = ROOT.TCanvas()
#hist.Draw()
#hist.SetOption("COLZ")
#c1.SetLogz()
#c1.Print(picpath+filename+"_ROOT2D.png")

#Fit
yBinSize = 200
width=yBinSize
xBinSize = 200
maximx = 760
maximy = 2000
savename = picpath+filename+"GaussiansFit"
from plotFit import gaussianFit,fitGaussians
def plot2DROOT(proj,ext,name):
        ca = ROOT.TCanvas()
        proj.Draw()
        #proj.GetXaxis().SetRangeUser(-20,20)
        #ca.SetLogy()
        ca.Print(name+"_ROOT2D_"+ext+".png")
        print(name+"_ROOT2D_"+ext+".png")

def plot1DMPL(x,Img,ImgName,picpath,lab):
    import matplotlib.pyplot as plt
    plt.scatter(x,Img,label=lab)
    plt.xlim(-40,40)
    plt.legend()
    plt.title(ImgName)
    plt.savefig(picpath+ImgName+".png")
    print(picpath+ImgName+".png")
    if lab == "GEANT4": plt.close()

def compare(ImgMCNP,ImgName,picpath,lim,axis):
    from plotFit import converter
    from sys import path as sysPath
    from os import chdir
    MiniScatter_path="../../MiniScatter/build/."
    sysPath.append(MiniScatter_path)
    chdir("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/")
    import miniScatterDriver
    print("else GaussFit")
    paths = {"scratchPath":"/scratch2/ericdf/PBWScatter/"}
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--Nb",        type=int,   default=100,   help="Number of macroparticles per beamlet. Default=10")
    parser.add_argument("--nP",        type=float, default=1e3,   help="Numper of beamlets in pulse. Default=1e3")
    parser.add_argument("--reBin",     type=int, default=1,  help="Number of bins to make into 1 in 2D histogram for smoothing")
    parser.add_argument("--Nbeamlet",type=float, default=1e7)
    parser.add_argument("--sim",       type=str,   default="beamlet")
    args = parser.parse_args()
    
    simSetup_simple1 = {'PHYS': 'QGSP_BERT_EMZ', 'ZOFFSET': '*-1', 'WORLDSIZE': 1000.0, 'QUICKMODE': False, 'MINIROOT': True, 
                        'ANASCATTER': True, 'EDEP_DZ': 1.0, 'CUTOFF_RADIUS': 1000.0, 'CUTOFF_ENERGYFRACTION': 0.9, 'POSLIM': 1000.0, 
                        'DIST': [3565], 'BEAM': 'proton', 'ENERGY': 570, 'COVAR': (0.0001, 0.15, 0, 0.0001, 0.15, 0), 'N': 10000000.0, 
                        'THICK': 0.0, 'MAGNET': [{'type': 'PBW', 'length': 0, 'gradient': 0.0, 'keyval': {'material': 'G4_Al', 'radius': 88.0, 
                                                            'al1Thick': 1.0, 'waterThick': 2.0, 'al2Thick': 1.25}, 'pos': 24.125}], 
                        'OUTFOLDER': '/scratch2/ericdf/PBWScatter/ESS/', 'OUTNAME': 'PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N1e+07_QBZ'}
    TRYLOAD = True 
    (twiss_PBW, numPart_PBW, objects_PBW) = miniScatterDriver.getData_tryLoad(simSetup_simple1, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212"])
    (ImgGeant4, xax, yax) = converter(objects_PBW["tracker_cutoff_xy_PDG2212"],True,"PencilBeamHist",paths,args,True) #convert from TH2D to numpy map
    from numpy import linspace
    a=1
    if axis in {"y","Y"}:
        lenM = round(np.shape(ImgMCNP)[1])
        projM = ImgMCNP[:,round(lenM/2)-a:round(lenM/2)+a]
        projM = (projM[:,0] + projM[:,1]) / (a*2+1)
    
        lenG = np.shape(ImgGeant4)[1]
        projG = ImgGeant4[:,round(lenG/2)-a:round(lenG/2)+a]
        projG = (projG[:,0] + projG[:,1]) / (a*2+1)
        projG = projG[round(lenG/2)-round(len(projM)/2):round(lenG/2)+round(len(projM)/2)]
    elif axis in {"x","X"}:
        lenM = round(np.shape(ImgMCNP)[0])
        projM = ImgMCNP[round(lenM/2)-a:round(lenM/2)+a,:]
        projM = (projM[0,:] + projM[1,:]) / (a*2+1)
    
        lenG = np.shape(ImgGeant4)[0]
        projG = ImgGeant4[round(lenG/2)-a:round(lenG/2)+a,:]
        projG = (projG[0,:] + projG[1,:]) / (a*2+1)
        projG = projG[round(lenG/2)-round(len(projM)/2):round(lenG/2)+round(len(projM)/2)]
    xM = linspace(-lim,lim,len(projM))
    xG = linspace(-lim,lim,len(projG))
    print("M",lenM/2-a,len(xM),len(projM))#,"\n",projM)
    print("G",round(lenG/2)-a,len(xG),len(projG))#,"\n",projG)
    plot1DMPL(xM,projM,ImgName+axis,picpath,"MCNP")
    plot1DMPL(xG,projG,ImgName+axis,picpath,"GEANT4")
    


axis="x"
#if   axis == "y" or axis == "Y": 
#    proj = hist.ProjectionY(axis,hist.GetYaxis().FindBin(-width),hist.GetYaxis().FindBin(width))
#elif axis == "x" or axis == "X":
#    proj = hist.ProjectionX(axis,hist.GetXaxis().FindBin(-width),hist.GetXaxis().FindBin(width))
#plot2DROOT(proj,"proj"+axis,picpath+filename)
#fitGaussians(Img,Error,picpath,100,"y",1)#y
#fitGaussians(Img,Error,picpath,highx,"x",3)#y
compare(Img,"ComboImages",picpath,highx,"X")
compare(Img,"ComboImages",picpath,highy,"Y")
#diffNy,diffPy,coeffsy, differenceNLy,differencePLy,coeffsLy = gaussianFit(hist,"y",yBinSize,maximx,savename,2,25,True,True)
#diffNx,diffPx,coeffsx, differenceNLx,differencePLx,coeffsLx = gaussianFit(hist,"x",xBinSize,maximy,savename,3,10,True,True)


