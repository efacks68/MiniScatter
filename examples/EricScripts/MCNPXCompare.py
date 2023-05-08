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
            Img[j,i] = float(row[ind])*1e5
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

path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/examples/EricScripts/"
file = "ProtonDensityMap2"
filename = "ProtonDensityMap2_570MeV"
ind = 2
lenx=77
lowx = -76
highx = 76
lowy = -200
highy = 200
leny = 201
Img,Error = decode(path,file,lenx,leny,ind)
#print(Img)
#optional: Visualize PDM or write in better format to CSV
def plot(Img,refImgName,picpath):
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
    savefig(picpath+refImgName+".png")#
    close(fig)
    close()
    print(picpath+refImgName+".png")    

picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
filename = file+"_2DMCNPXData"
plot(Img,filename,picpath)

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

c1 = ROOT.TCanvas()
hist.Draw()
hist.SetOption("COLZ")
c1.SetLogz()
c1.Print(picpath+filename+"_ROOT2D.png")

#Fit
yBinSize = 200
width=yBinSize
xBinSize = 200
maximx = 380
maximy = 1000
savename = picpath+filename+"GaussiansFit"
from plotFit import gaussianFit
def plot(proj,ext,name):
        ca = ROOT.TCanvas()
        proj.Draw()
        #proj.GetXaxis().SetRangeUser(-20,20)
        #ca.SetLogy()
        ca.Print(name+"_ROOT2D_"+ext+".png")
        print(name+"_ROOT2D_"+ext+".png")
axis="x"
#if   axis == "y" or axis == "Y": 
#    proj = hist.ProjectionY(axis,hist.GetYaxis().FindBin(-width),hist.GetYaxis().FindBin(width))
#elif axis == "x" or axis == "X":
#    proj = hist.ProjectionX(axis,hist.GetXaxis().FindBin(-width),hist.GetXaxis().FindBin(width))
#plot(proj,"proj"+axis,picpath+filename)
#diffNy,diffPy,coeffsy, differenceNLy,differencePLy,coeffsLy = gaussianFit(hist,"y",yBinSize,maximx,savename,2,25,True,True)
diffNx,diffPx,coeffsx, differenceNLx,differencePLx,coeffsLx = gaussianFit(hist,"x",xBinSize,maximy,savename,3,10,True,True)

