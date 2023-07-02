####physListCompare.py
import ROOT
from plotFit import gaussianFit
##which axis to project?
axis = "Y" #"X","Y" 
column = 2 #columns: 2= p>564MeV; 4=p<564MeV; 10=total
picFormat = ".png" #".png",".pdf", ".jpeg"
xlim = 500
g=3 #number of Gaussians to fit to data
##Settings for reading
path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
file = "MCNP_ProtonDensity_570MeV_10June"
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
scratchPath = "/scratch2/ericdf/PBWScatter/"
if axis in {"x","X"}:
    ylim = [1e-5,1.2e-1]
    mcnpMarkerSize = 5
else:
    ylim = [1e-7,1.2e-1]
    mcnpMarkerSize = 4

###Load in MCNP data for Total and Primary beams
fMCNP_BERT_Tot = ROOT.TFile(scratchPath+"MCNP_Pencil_570MeV_Total.root","r")
MCNP_BERT_Tot = fMCNP_BERT_Tot.Get("MCBERT_TargetProtons").Clone("MCNP_BERT_Total")
_,_,coeffsMCNPToty,_,_,_ = gaussianFit(MCNP_BERT_Tot,axis,100,500,picpath+"physListComp_MCNP_BERT_Tot",2,25,True,True,1,True)

##Primary Beam (E>564MeV)
fMCNP_BERT_Pri = ROOT.TFile(scratchPath+"MCNP_Pencil_570MeV_Primary.root","r")
MCNP_BERT_Pri = fMCNP_BERT_Pri.Get("MCBERT_TargetProtons").Clone("MCNP_BERT_Primary")
_,_,coeffsMCNPPriy,_,_,_ = gaussianFit(MCNP_BERT_Pri,axis,100,500,picpath+"physListComp_MCNP_BERT_Pri",2,25,True,True,1,True)


###Load in all scattered beams of different GEANT4 phys lists and fit
##Open file with ROOT
##QGSP_BERT_EMZ
QBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
fQBERTZ = ROOT.TFile(QBERTZ_file,"r")
##Get TH2D and clone in
G4_QBERTZ_TH2 = fQBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
#Send to fit function, returns coefficients
_,_,coeffsQBERTZy,_,_,_ = gaussianFit(G4_QBERTZ_TH2,axis,100,500,picpath+"physListComp_G4_QBERTZ",2,25,True,True,g,True)

##EMonly_EMZ
EMZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08EMonly_EMZ.root"
fEMZ = ROOT.TFile(EMZ_file,"r")
G4_EMZ_TH2 = fEMZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
_,_,coeffsEMy,_,_,_ = gaussianFit(G4_EMZ_TH2,axis,100,500,picpath+"physListComp_G4_EMZ",2,25,True,True,g,True)

##FTFP_BERT_EMZ
#FBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
#fFBERTZ = ROOT.TFile(FBERTZ_file,"r")
#G4_FBERTZ_TH2 = fFBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
#_,_,coeffsFBERTZy,_,_,_ = gaussianFit(G4_FBERTZ_TH2,axis,100,500,picpath+"physListComp_G4_FBERTZ",2,25,True,True,g,True)

##QGSP_FTFP_BERT__SS
QFBERTSS_file = "/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N1e+06QGSP_FTFP_BERT__SS.root"
fQFBERTSS = ROOT.TFile(QFBERTSS_file,"r")
G4_QFBERTSS_TH2 = fQFBERTSS.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
_,_,coeffsQFBERTSSy,_,_,_ = gaussianFit(G4_QFBERTSS_TH2,axis,100,500,picpath+"physListComp_G4_QFBERTSS",2,25,True,True,g,True)


##Normalize TH2Ds and make Projections in the selcted axis
width=500 #width of map to sum over
##MCNP_BERT_Tot
integral = MCNP_BERT_Tot.Integral(MCNP_BERT_Tot.GetXaxis().FindBin(-width),MCNP_BERT_Tot.GetXaxis().FindBin(width),MCNP_BERT_Tot.GetYaxis().FindBin(-width),MCNP_BERT_Tot.GetYaxis().FindBin(width))
MCNP_BERT_Tot.Scale(1/integral)
sum = MCNP_BERT_Tot.Integral()
##Make Projection
if axis in {"y","Y"}:
    proj_MCNP_BERT_Tot = MCNP_BERT_Tot.ProjectionY(axis,MCNP_BERT_Tot.GetXaxis().FindBin(-width),MCNP_BERT_Tot.GetXaxis().FindBin(width),"e")
else:
    proj_MCNP_BERT_Tot = MCNP_BERT_Tot.ProjectionX(axis,MCNP_BERT_Tot.GetYaxis().FindBin(-width),MCNP_BERT_Tot.GetYaxis().FindBin(width),"e")
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
## The names are NOT inherited from their TH2D, like I had expected.
proj_MCNP_BERT_Tot.SetName("proj_MCBERT_Tot")

##MCNP_BERT_Pri
integral = MCNP_BERT_Pri.Integral(MCNP_BERT_Pri.GetXaxis().FindBin(-width),MCNP_BERT_Pri.GetXaxis().FindBin(width),MCNP_BERT_Pri.GetYaxis().FindBin(-width),MCNP_BERT_Pri.GetYaxis().FindBin(width))
MCNP_BERT_Pri.Scale(1/integral)
sum = MCNP_BERT_Pri.Integral()
##Make Projection
if axis in {"y","Y"}:
    proj_MCNP_BERT_Pri = MCNP_BERT_Pri.ProjectionY(axis,MCNP_BERT_Pri.GetXaxis().FindBin(-width),MCNP_BERT_Pri.GetXaxis().FindBin(width),"e")
else:
    proj_MCNP_BERT_Pri = MCNP_BERT_Pri.ProjectionX(axis,MCNP_BERT_Pri.GetYaxis().FindBin(-width),MCNP_BERT_Pri.GetYaxis().FindBin(width),"e")
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
## The names are NOT inherited from their TH2D, like I had expected.
proj_MCNP_BERT_Pri.SetName("proj_MCBERT_Pri")

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
#integral = G4_FBERTZ_TH2.Integral(G4_FBERTZ_TH2.GetXaxis().FindBin(-width),G4_FBERTZ_TH2.GetXaxis().FindBin(width),G4_FBERTZ_TH2.GetYaxis().FindBin(-width),G4_FBERTZ_TH2.GetYaxis().FindBin(width))
#G4_FBERTZ_TH2.Scale(1/integral)
#sum = G4_FBERTZ_TH2.Integral()
#if axis in {"y","Y"}:
#    proj_G4_FBERTZ = G4_FBERTZ_TH2.ProjectionY(axis,G4_FBERTZ_TH2.GetXaxis().FindBin(-width),G4_FBERTZ_TH2.GetXaxis().FindBin(width),"e")
#else:
#    proj_G4_FBERTZ = G4_FBERTZ_TH2.ProjectionX(axis,G4_FBERTZ_TH2.GetYaxis().FindBin(-width),G4_FBERTZ_TH2.GetYaxis().FindBin(width),"e")
#proj_G4_FBERTZ.SetName("proj_G4FBERTZ")

#G4_EMZ (EMonly_EMZ)
integral = G4_EMZ_TH2.Integral(G4_EMZ_TH2.GetXaxis().FindBin(-width),G4_EMZ_TH2.GetXaxis().FindBin(width),G4_EMZ_TH2.GetYaxis().FindBin(-width),G4_EMZ_TH2.GetYaxis().FindBin(width))
G4_EMZ_TH2.Scale(1/integral)
sum = G4_EMZ_TH2.Integral()
if axis in {"y","Y"}:
    proj_G4_EMZ = G4_EMZ_TH2.ProjectionY(axis,G4_EMZ_TH2.GetXaxis().FindBin(-width),G4_EMZ_TH2.GetXaxis().FindBin(width),"e")
else:
    proj_G4_EMZ = G4_EMZ_TH2.ProjectionX(axis,G4_EMZ_TH2.GetYaxis().FindBin(-width),G4_EMZ_TH2.GetYaxis().FindBin(width),"e")
proj_G4_EMZ.SetName("proj_G4EMZ")

##G4_QFBERTSS (QGSP_FTFP_BERT_SS)
integral = G4_QFBERTSS_TH2.Integral(G4_QFBERTSS_TH2.GetXaxis().FindBin(-width),G4_QFBERTSS_TH2.GetXaxis().FindBin(width),G4_QFBERTSS_TH2.GetYaxis().FindBin(-width),G4_QFBERTSS_TH2.GetYaxis().FindBin(width))
G4_QFBERTSS_TH2.Scale(1/integral)
sum = G4_QFBERTSS_TH2.Integral()
if axis in {"y","Y"}:
    proj_G4_QFBERTSS = G4_QFBERTSS_TH2.ProjectionY(axis,G4_QFBERTSS_TH2.GetXaxis().FindBin(-width),G4_QFBERTSS_TH2.GetXaxis().FindBin(width),"e")
else:
    proj_G4_QFBERTSS = G4_QFBERTSS_TH2.ProjectionX(axis,G4_QFBERTSS_TH2.GetYaxis().FindBin(-width),G4_QFBERTSS_TH2.GetYaxis().FindBin(width),"e")
proj_G4_QFBERTSS.SetName("proj_G4QFBERTSS")

###Plot them all together with different colors
##Make Canvas and Lengend
c2 = ROOT.TCanvas("Scattered Beams","Scattered Beams",0,0,500*8,400*8)
c2.SetLogy()
leg = ROOT.TLegend(0.65,0.7,0.98,0.9)

##Draw 1st projection and configure plot settings
proj_MCNP_BERT_Tot.Draw()
#proj_MCNP_BERT_Tot.SetLineColor(1)
proj_MCNP_BERT_Tot.SetMarkerStyle(34)
proj_MCNP_BERT_Tot.SetMarkerColor(1)
proj_MCNP_BERT_Tot.SetMarkerSize(mcnpMarkerSize)
proj_MCNP_BERT_Tot.SetStats(False)
proj_MCNP_BERT_Tot.SetTitle(axis+" Scattered Pencil Beam Density at Target")
proj_MCNP_BERT_Tot.GetXaxis().SetRangeUser(-xlim,xlim)
proj_MCNP_BERT_Tot.GetYaxis().SetRangeUser(ylim[0],ylim[1])
##Make Legend Entry!
leg.AddEntry(proj_MCNP_BERT_Tot,"MCNP BERTINI Total Energy Beam")

##2nd Projection
proj_MCNP_BERT_Pri.Draw("SAME")
#proj_G4_QBERTZ.SetLineColor(2)
proj_MCNP_BERT_Pri.SetMarkerStyle(41)
proj_MCNP_BERT_Pri.SetMarkerColor(6)
proj_MCNP_BERT_Pri.SetMarkerSize(4)
leg.AddEntry(proj_MCNP_BERT_Pri,"MCNP BERTINI E > 564MeV Beam")

##2nd Projection
proj_G4_QBERTZ.Draw("SAME")
#proj_G4_QBERTZ.SetLineColor(2)
proj_G4_QBERTZ.SetMarkerStyle(21)
proj_G4_QBERTZ.SetMarkerColor(2)
proj_G4_QBERTZ.SetMarkerSize(4)
leg.AddEntry(proj_G4_QBERTZ,"GEANT4 QGSP_BERTINI_EMZ")

proj_G4_QFBERTSS.Draw("SAME")
#proj_G4_FBERTZ.SetLineColor(3)
proj_G4_QFBERTSS.SetMarkerStyle(20)
proj_G4_QFBERTSS.SetMarkerColor(3) #8 is mid green, 3 is light green, 7 is light blue
proj_G4_QFBERTSS.SetMarkerSize(4)
leg.AddEntry(proj_G4_QFBERTSS,"GEANT4 QF_BERTINI__SS 1e6")

proj_G4_EMZ.Draw("SAME")
#proj_G4_EMZ.SetLineColor(4)
proj_G4_EMZ.SetMarkerStyle(47)
proj_G4_EMZ.SetMarkerColor(4)
proj_G4_EMZ.SetMarkerSize(3)
leg.AddEntry(proj_G4_EMZ,"GEANT4 EMZ")

##Draw and Print Canvas and Legend
leg.Draw()
c2.Draw()
c2.Print(picpath+"PhysLists_BeamOverlap"+axis+"_TotPri"+picFormat) #only one is being drawn.

###Add together the different lists. Renormalize them to be 1
###Export that to a histogram