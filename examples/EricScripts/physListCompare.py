####physListCompare.py
import ROOT
##which axis to project?
axis = "X" #"X","Y" 

###Load in MCNPX
##Function to read in particle densities as a 2D map:
def decode2DMap(path,file,lenx,leny,column):
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
            Img[j,i] = float(row[column])
            Error[j,i] = float(row[column+1])
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
#Settings to read in
path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
file = "MCNP_ProtonDensity_570MeV"
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
column = 2 #column with p>564MeV
lenx  = 77 #indices
lowx  = -76  #[mm]
highx = 76   #[mm]
lowy  = -200 #[mm]
highy = 200  #[mm]
leny  = 201 #indices
##Call decoder to 2D Map function
Img,Error = decode2DMap(path,file,lenx,leny,column)

###Convert MCNPX 2D map to ROOT
##Make TH2D histogram to read the 2D Map into, same size as 2D Map
MCNP_BERT = ROOT.TH2D("MCBERT_TargetProtons","Scattered Pencil Beam at Target;Horizontal [mm];Vertical [mm]",lenx,lowx,highx,leny,lowy,highy)
#Help from Kyrre's converter function and https://root-forum.cern.ch/t/create-and-fill-a-th2f-using-arrays-for-binning/27161/2
##Get each bin and set bin content and error
for xi in range(lenx):
    for yi in range(leny):
        #bx = hIn.GetBin(xi+1,yi+1)
        #MCNP_BERT.GetBin(yi,xi) = Img[yi,xi]
        bin = MCNP_BERT.GetBin(xi+1,yi+1)
        MCNP_BERT.SetBinContent(bin,Img[yi,xi])
        MCNP_BERT.SetBinError(bin,Error[yi,xi])
        #if xi == 38 and yi == 100:
            #print(Img[yi,xi],MCNP_BERT.GetBinContent(MCNP_BERT.GetBin(xi+1,yi+1)))


###Load in all scattered beams of different GEANT4 phys lists and fit
g=5 #number of Gaussians to fit to data
from plotFit import gaussianFit

##Open file with ROOT
QBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
fQBERTZ = ROOT.TFile(QBERTZ_file,"r")
##Get TH2D and clone in
G4_QBERTZ_TH2 = fQBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
#Send to fit function, returns coefficients
_,_,coeffsQBERTZy,_,_,_ = gaussianFit(G4_QBERTZ_TH2,"y",100,500,picpath+"physListComp_G4_QBERTZ",2,25,True,True,g,True)

EMZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08EMonly_EMZ.root"
fEMZ = ROOT.TFile(EMZ_file,"r")
G4_EMZ_TH2 = fEMZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
_,_,coeffsEMy,_,_,_ = gaussianFit(G4_EMZ_TH2,"y",100,500,picpath+"physListComp_G4_EMZ",2,25,True,True,g,True)

FBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
fFBERTZ = ROOT.TFile(FBERTZ_file,"r")
G4_FBERTZ_TH2 = fFBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
_,_,coeffsQBERTZy,_,_,_ = gaussianFit(G4_FBERTZ_TH2,"y",100,500,picpath+"physListComp_G4_FBERTZ",2,25,True,True,g,True)


##Normalize TH2Ds and make Projections in the selcted axis
width=500 #width of map to sum over
##MCNP_BERT
integral = MCNP_BERT.Integral(MCNP_BERT.GetXaxis().FindBin(-width),MCNP_BERT.GetXaxis().FindBin(width),MCNP_BERT.GetYaxis().FindBin(-width),MCNP_BERT.GetYaxis().FindBin(width))
MCNP_BERT.Scale(1/integral)
sum = MCNP_BERT.Integral()
##Make Projection
if axis in {"y","Y"}:
    proj_MCNP_BERT = MCNP_BERT.ProjectionY(axis,MCNP_BERT.GetXaxis().FindBin(-width),MCNP_BERT.GetXaxis().FindBin(width),"e")
else:
    proj_MCNP_BERT = MCNP_BERT.ProjectionX(axis,MCNP_BERT.GetYaxis().FindBin(-width),MCNP_BERT.GetYaxis().FindBin(width),"e")
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
## The names are NOT inherited from their TH2D, like I had expected.
proj_MCNP_BERT.SetName("proj_MCBERT")

##G4_QBERTZ (QGSP_BERT_EMZ)
integral = G4_QBERTZ_TH2.Integral(G4_QBERTZ_TH2.GetXaxis().FindBin(-width),G4_QBERTZ_TH2.GetXaxis().FindBin(width),G4_QBERTZ_TH2.GetYaxis().FindBin(-width),G4_QBERTZ_TH2.GetYaxis().FindBin(width))
G4_QBERTZ_TH2.Scale(1/integral)
sum = G4_QBERTZ_TH2.Integral()
if axis in {"y","Y"}:
    proj_G4_QBERTZ = G4_QBERTZ_TH2.ProjectionY(axis,G4_QBERTZ_TH2.GetXaxis().FindBin(-width),G4_QBERTZ_TH2.GetXaxis().FindBin(width),"e")
else:
    proj_G4_QBERTZ = G4_QBERTZ_TH2.ProjectionX(axis,G4_QBERTZ_TH2.GetYaxis().FindBin(-width),G4_QBERTZ_TH2.GetYaxis().FindBin(width),"e")
proj_G4_QBERTZ.SetName("proj_G4QBERTZ")

##G4_FBERTZ (FTFP_BERT_EMZ)
integral = G4_FBERTZ_TH2.Integral(G4_FBERTZ_TH2.GetXaxis().FindBin(-width),G4_FBERTZ_TH2.GetXaxis().FindBin(width),G4_FBERTZ_TH2.GetYaxis().FindBin(-width),G4_FBERTZ_TH2.GetYaxis().FindBin(width))
G4_FBERTZ_TH2.Scale(1/integral)
sum = G4_FBERTZ_TH2.Integral()
if axis in {"y","Y"}:
    proj_G4_FBERTZ = G4_FBERTZ_TH2.ProjectionY(axis,G4_FBERTZ_TH2.GetXaxis().FindBin(-width),G4_FBERTZ_TH2.GetXaxis().FindBin(width),"e")
else:
    proj_G4_FBERTZ = G4_FBERTZ_TH2.ProjectionX(axis,G4_FBERTZ_TH2.GetYaxis().FindBin(-width),G4_FBERTZ_TH2.GetYaxis().FindBin(width),"e")
proj_G4_FBERTZ.SetName("proj_G4FBERTZ")

#G4_EMZ (EMonly_EMZ)
integral = G4_EMZ_TH2.Integral(G4_EMZ_TH2.GetXaxis().FindBin(-width),G4_EMZ_TH2.GetXaxis().FindBin(width),G4_EMZ_TH2.GetYaxis().FindBin(-width),G4_EMZ_TH2.GetYaxis().FindBin(width))
G4_EMZ_TH2.Scale(1/integral)
sum = G4_EMZ_TH2.Integral()
if axis in {"y","Y"}:
    proj_G4_EMZ = G4_EMZ_TH2.ProjectionY(axis,G4_EMZ_TH2.GetXaxis().FindBin(-width),G4_EMZ_TH2.GetXaxis().FindBin(width),"e")
else:
    proj_G4_EMZ = G4_EMZ_TH2.ProjectionX(axis,G4_EMZ_TH2.GetYaxis().FindBin(-width),G4_EMZ_TH2.GetYaxis().FindBin(width),"e")
proj_G4_EMZ.SetName("proj_G4EMZ")

###Plot them all together with different colors
##Make Canvas and Lengend
c2 = ROOT.TCanvas("Scattered Beams","Scattered Beams",0,0,500*8,400*8)
c2.SetLogy()
leg = ROOT.TLegend(0.65,0.75,0.9,0.9)
xlim = 125

##Draw 1st projection and configure plot settings
proj_MCNP_BERT.Draw()
proj_MCNP_BERT.SetLineColor(1)
proj_MCNP_BERT.SetMarkerStyle(34)
proj_MCNP_BERT.SetMarkerColor(1)
proj_MCNP_BERT.SetMarkerSize(3)
proj_MCNP_BERT.SetStats(False)
proj_MCNP_BERT.SetTitle("Scattered Pencil Beam Density at Target")
proj_MCNP_BERT.GetXaxis().SetRangeUser(-xlim,xlim)
proj_MCNP_BERT.GetYaxis().SetRangeUser(1e-6,1.2e-1)
##Make Legend Entry!
leg.AddEntry(proj_MCNP_BERT,"MCNP BERTINI")

##2nd Projection
proj_G4_QBERTZ.Draw("SAME")
proj_G4_QBERTZ.SetLineColor(2)
proj_G4_QBERTZ.SetMarkerStyle(21)
proj_G4_QBERTZ.SetMarkerColor(2)
proj_G4_QBERTZ.SetMarkerSize(4)
leg.AddEntry(proj_G4_QBERTZ,"GEANT4 QGSP_BERTINI_EMZ")

proj_G4_FBERTZ.Draw("SAME")
proj_G4_FBERTZ.SetLineColor(3)
proj_G4_FBERTZ.SetMarkerStyle(20)
proj_G4_FBERTZ.SetMarkerColor(3)
proj_G4_FBERTZ.SetMarkerSize(3)
leg.AddEntry(proj_G4_FBERTZ,"GEANT4 FTFP_BERTINI_EMZ")

proj_G4_EMZ.Draw("SAME")
proj_G4_EMZ.SetLineColor(4)
proj_G4_EMZ.SetMarkerStyle(47)
proj_G4_EMZ.SetMarkerColor(4)
proj_G4_EMZ.SetMarkerSize(3)
leg.AddEntry(proj_G4_EMZ,"GGEANT4 EMZ")

##Draw and Print Canvas and Legend
leg.Draw()
c2.Draw()
c2.Print(picpath+"PhysLists_BeamOverlap"+axis+".png") #only one is being drawn.

###Add together the different lists. Renormalize them to be 1
###Export that to a histogram