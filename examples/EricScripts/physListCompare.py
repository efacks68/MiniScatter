####physListCompare.py
import ROOT
from plotFit import gaussianFit
##which axis to project?
axis = "X" #"X","Y" 
column = 2 #columns: 2= p>564MeV; 4=p<564MeV; 10=total
picFormat = ".png" #".png",".pdf", ".jpeg"
xlim = 500
width=10 #width of map to sum over
g=3 #number of Gaussians to fit to data
##Settings for reading
path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
file = "MCNP_ProtonDensity_570MeV_10June"
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
scratchPath = "/scratch2/ericdf/PBWScatter/"
if axis in {"x","X"}:
    ylim = [1e-9,2e-1]
    mcnpMarkerSize = 5
else:
    ylim = [1e-9,2e-1]
    mcnpMarkerSize = 4

###Load in MCNP data for Total and Primary beams
fMCNP_BERT_Tot = ROOT.TFile(scratchPath+"MCNP_Pencil_570MeV_Total.root","r")
MCNP_BERT_Tot = fMCNP_BERT_Tot.Get("MCBERT_Target-TotalProtons").Clone("MCNP_BERT_Total")
#_,_,coeffsMCNPToty,_,_,_ = gaussianFit(MCNP_BERT_Tot,axis,100,500,picpath+"physListComp_MCNP_BERT_Tot",2,25,True,True,1,True)

##Primary Beam (E>564MeV)
fMCNP_BERT_Pri = ROOT.TFile(scratchPath+"MCNP_Pencil_570MeV_Primary.root","r")
MCNP_BERT_Pri = fMCNP_BERT_Pri.Get("MCBERT_Target-PrimaryProtons").Clone("MCNP_BERT_Primary")
#_,_,coeffsMCNPPriy,_,_,_ = gaussianFit(MCNP_BERT_Pri,axis,100,500,picpath+"physListComp_MCNP_BERT_Pri",2,25,True,True,1,True)


###Load in all scattered beams of different GEANT4 phys lists and fit
##Open file with ROOT
##QGSP_BERT_EMZ
QBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
fQBERTZ = ROOT.TFile(QBERTZ_file,"r")
##Get TH2D and clone in
G4_QBERTZ_TH2 = fQBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("QBERTZ")
#Send to fit function, returns coefficients
#_,_,coeffsQBERTZy,_,_,_ = gaussianFit(G4_QBERTZ_TH2,axis,100,500,picpath+"physListComp_G4_QBERTZ",2,25,True,True,g,True)

##EMonly_EMZ
EMZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08EMonly_EMZ.root"
fEMZ = ROOT.TFile(EMZ_file,"r")
G4_EMZ_TH2 = fEMZ.Get("tracker_cutoff_xy_PDG2212").Clone("EMZ")
#_,_,coeffsEMy,_,_,_ = gaussianFit(G4_EMZ_TH2,axis,100,500,picpath+"physListComp_G4_EMZ",2,25,True,True,g,True)

##QGSP_FTFP_BERT__SS
QFBERTSS_file = "/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N1e+06QGSP_FTFP_BERT__SS.root"
fQFBERTSS = ROOT.TFile(QFBERTSS_file,"r")
G4_QFBERTSS_TH2 = fQFBERTSS.Get("tracker_cutoff_xy_PDG2212").Clone("QFBSS")
#_,_,coeffsQFBERTSSy,_,_,_ = gaussianFit(G4_QFBERTSS_TH2,axis,100,500,picpath+"physListComp_G4_QFBERTSS",2,25,True,True,g,True)

##QGSP_FTFP_VERT_PEN
QFBPEN_file = "/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08QGSP_FTFP_BERT_PEN.root"
fQFBPEN = ROOT.TFile(QFBPEN_file,"r")
G4_QFBPEN_TH2 = fQFBPEN.Get("tracker_cutoff_xy_PDG2212").Clone("QFBPEN")

##QGSP_FTFP_VERT_LIV
QFBLIV_file = "/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08QGSP_FTFP_BERT_LIV.root"
fQFBLIV = ROOT.TFile(QFBLIV_file,"r")
G4_QFBLIV_TH2 = fQFBLIV.Get("tracker_cutoff_xy_PDG2212").Clone("QFBLIV")

###Make & Normalize Projections in the selcted axis
##MCNP_BERT_Tot
#integral = MCNP_BERT_Tot.Integral(MCNP_BERT_Tot.GetXaxis().FindBin(-xlim),MCNP_BERT_Tot.GetXaxis().FindBin(xlim),MCNP_BERT_Tot.GetYaxis().FindBin(-xlim),MCNP_BERT_Tot.GetYaxis().FindBin(xlim))
#MCNP_BERT_Tot.Scale(1/integral)
#sum = MCNP_BERT_Tot.Integral()
##Make Projection
if axis in {"y","Y"}:
    proj_MCNP_BERT_Tot = MCNP_BERT_Tot.ProjectionY(axis,MCNP_BERT_Tot.GetXaxis().FindBin(-width),MCNP_BERT_Tot.GetXaxis().FindBin(width),"e")
    integral = proj_MCNP_BERT_Tot.Integral(proj_MCNP_BERT_Tot.GetXaxis().FindBin(-xlim),proj_MCNP_BERT_Tot.GetXaxis().FindBin(xlim),option="width")
    proj_MCNP_BERT_Tot.Scale(1/integral)
    sum = proj_MCNP_BERT_Tot.Integral()
else:
    proj_MCNP_BERT_Tot = MCNP_BERT_Tot.ProjectionX(axis,MCNP_BERT_Tot.GetYaxis().FindBin(-width),MCNP_BERT_Tot.GetYaxis().FindBin(width),"e")
    integral = proj_MCNP_BERT_Tot.Integral(proj_MCNP_BERT_Tot.GetXaxis().FindBin(-xlim),proj_MCNP_BERT_Tot.GetXaxis().FindBin(xlim),option="width")
    proj_MCNP_BERT_Tot.Scale(1/integral)
    sum = proj_MCNP_BERT_Tot.Integral()
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
## The names are NOT inherited from their TH2D, like I had expected.
proj_MCNP_BERT_Tot.SetName("proj_MCBERT_Tot")

##MCNP_BERT_Pri
#integral = MCNP_BERT_Pri.Integral(MCNP_BERT_Pri.GetXaxis().FindBin(-xlim),MCNP_BERT_Pri.GetXaxis().FindBin(xlim),MCNP_BERT_Pri.GetYaxis().FindBin(-xlim),MCNP_BERT_Pri.GetYaxis().FindBin(xlim))
#MCNP_BERT_Pri.Scale(1/integral)
#sum = MCNP_BERT_Pri.Integral()
##Make Projection
if axis in {"y","Y"}:
    proj_MCNP_BERT_Pri = MCNP_BERT_Pri.ProjectionY(axis,MCNP_BERT_Pri.GetXaxis().FindBin(-width),MCNP_BERT_Pri.GetXaxis().FindBin(width),"e")
    integral = proj_MCNP_BERT_Pri.Integral(proj_MCNP_BERT_Pri.GetXaxis().FindBin(-xlim),proj_MCNP_BERT_Pri.GetXaxis().FindBin(xlim),option="width")
    proj_MCNP_BERT_Pri.Scale(1/integral)
    sum = proj_MCNP_BERT_Pri.Integral()
else:
    proj_MCNP_BERT_Pri = MCNP_BERT_Pri.ProjectionX(axis,MCNP_BERT_Pri.GetYaxis().FindBin(-width),MCNP_BERT_Pri.GetYaxis().FindBin(width),"e")
    integral = proj_MCNP_BERT_Pri.Integral(proj_MCNP_BERT_Pri.GetXaxis().FindBin(-xlim),proj_MCNP_BERT_Pri.GetXaxis().FindBin(xlim),option="width")
    proj_MCNP_BERT_Pri.Scale(1/integral)
    sum = proj_MCNP_BERT_Pri.Integral()
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
## The names are NOT inherited from their TH2D, like I had expected.
proj_MCNP_BERT_Pri.SetName("proj_MCBERT_Pri")

##G4_QBERTZ (QGSP_BERT_EMZ)
if axis in {"y","Y"}:
    proj_G4_QBERTZ = G4_QBERTZ_TH2.ProjectionY(axis,G4_QBERTZ_TH2.GetXaxis().FindBin(-width),G4_QBERTZ_TH2.GetXaxis().FindBin(width),"e")
    integral = proj_G4_QBERTZ.Integral(proj_G4_QBERTZ.GetXaxis().FindBin(-xlim),proj_G4_QBERTZ.GetXaxis().FindBin(xlim),option="width")
    proj_G4_QBERTZ.Scale(1/integral)
    sum = proj_G4_QBERTZ.Integral()
else:
    proj_G4_QBERTZ = G4_QBERTZ_TH2.ProjectionX(axis,G4_QBERTZ_TH2.GetYaxis().FindBin(-width),G4_QBERTZ_TH2.GetYaxis().FindBin(width),"e")
    integral = proj_G4_QBERTZ.Integral(proj_G4_QBERTZ.GetXaxis().FindBin(-xlim),proj_G4_QBERTZ.GetXaxis().FindBin(xlim),option="width")
    proj_G4_QBERTZ.Scale(1/integral)
    sum = proj_G4_QBERTZ.Integral()
proj_G4_QBERTZ.SetName("proj_G4QBERTZ")

#G4_EMZ (EMonly_EMZ)
if axis in {"y","Y"}:
    proj_G4_EMZ = G4_EMZ_TH2.ProjectionY(axis,G4_EMZ_TH2.GetXaxis().FindBin(-width),G4_EMZ_TH2.GetXaxis().FindBin(width),"e")
    integral = proj_G4_EMZ.Integral(proj_G4_EMZ.GetXaxis().FindBin(-xlim),proj_G4_EMZ.GetXaxis().FindBin(xlim),option="width")
    proj_G4_EMZ.Scale(1/integral)
    sum = proj_G4_EMZ.Integral()
else:
    proj_G4_EMZ = G4_EMZ_TH2.ProjectionX(axis,G4_EMZ_TH2.GetYaxis().FindBin(-width),G4_EMZ_TH2.GetYaxis().FindBin(width),"e")
    integral = proj_G4_EMZ.Integral(proj_G4_EMZ.GetXaxis().FindBin(-xlim),proj_G4_EMZ.GetXaxis().FindBin(xlim),option="width")
    proj_G4_EMZ.Scale(1/integral)
    sum = proj_G4_EMZ.Integral()
proj_G4_EMZ.SetName("proj_G4EMZ")

##G4_QFBERTSS (QGSP_FTFP_BERT_SS)
if axis in {"y","Y"}:
    proj_G4_QFBERTSS = G4_QFBERTSS_TH2.ProjectionY(axis,G4_QFBERTSS_TH2.GetXaxis().FindBin(-width),G4_QFBERTSS_TH2.GetXaxis().FindBin(width),"e")
    integral = proj_G4_QFBERTSS.Integral(proj_G4_QFBERTSS.GetXaxis().FindBin(-xlim),proj_G4_QFBERTSS.GetXaxis().FindBin(xlim),option="width")
    proj_G4_QFBERTSS.Scale(1/integral)
    sum = proj_G4_QFBERTSS.Integral()
else:
    proj_G4_QFBERTSS = G4_QFBERTSS_TH2.ProjectionX(axis,G4_QFBERTSS_TH2.GetYaxis().FindBin(-width),G4_QFBERTSS_TH2.GetYaxis().FindBin(width),"e")
    integral = proj_G4_QFBERTSS.Integral(proj_G4_QFBERTSS.GetXaxis().FindBin(-xlim),proj_G4_QFBERTSS.GetXaxis().FindBin(xlim),option="width")
    proj_G4_QFBERTSS.Scale(1/integral)
    sum = proj_G4_QFBERTSS.Integral()
proj_G4_QFBERTSS.SetName("proj_G4QFBERTSS")

##G4_QFBPEN (QGSP_FTFP_BERT_PEN)
if axis in {"y","Y"}:
    proj_G4_QFBPEN = G4_QFBPEN_TH2.ProjectionY(axis,G4_QFBPEN_TH2.GetXaxis().FindBin(-width),G4_QFBPEN_TH2.GetXaxis().FindBin(width),"e")
    integral = proj_G4_QFBPEN.Integral(proj_G4_QFBPEN.GetXaxis().FindBin(-xlim),proj_G4_QFBPEN.GetXaxis().FindBin(xlim),option="width")
    proj_G4_QFBPEN.Scale(1/integral)
    sum = proj_G4_QFBPEN.Integral()
else:
    proj_G4_QFBPEN = G4_QFBPEN_TH2.ProjectionX(axis,G4_QFBPEN_TH2.GetYaxis().FindBin(-width),G4_QFBPEN_TH2.GetYaxis().FindBin(width),"e")
    integral = proj_G4_QFBPEN.Integral(proj_G4_QFBPEN.GetXaxis().FindBin(-xlim),proj_G4_QFBPEN.GetXaxis().FindBin(xlim),option="width")
    proj_G4_QFBPEN.Scale(1/integral)
    sum = proj_G4_QFBPEN.Integral()
proj_G4_QFBPEN.SetName("proj_G4QFBERTSS")

##G4_QFBLIV (QGSP_FTFP_BERT_LIV)
if axis in {"y","Y"}:
    proj_G4_QFBLIV = G4_QFBLIV_TH2.ProjectionY(axis,G4_QFBLIV_TH2.GetXaxis().FindBin(-width),G4_QFBLIV_TH2.GetXaxis().FindBin(width),"e")
    integral = proj_G4_QFBLIV.Integral(proj_G4_QFBLIV.GetXaxis().FindBin(-xlim),proj_G4_QFBLIV.GetXaxis().FindBin(xlim),option="width")
    proj_G4_QFBLIV.Scale(1/integral)
    sum = proj_G4_QFBLIV.Integral()
else:
    proj_G4_QFBLIV = G4_QFBLIV_TH2.ProjectionX(axis,G4_QFBLIV_TH2.GetYaxis().FindBin(-width),G4_QFBLIV_TH2.GetYaxis().FindBin(width),"e")
    integral = proj_G4_QFBLIV.Integral(proj_G4_QFBLIV.GetXaxis().FindBin(-xlim),proj_G4_QFBLIV.GetXaxis().FindBin(xlim),option="width")
    proj_G4_QFBLIV.Scale(1/integral)
    sum = proj_G4_QFBLIV.Integral()
proj_G4_QFBLIV.SetName("proj_G4QFBERTSS")

###Plot them all together with different colors
##Make Canvas and Lengend
c2 = ROOT.TCanvas("Scattered Beams","Scattered Beams",0,0,500*8,400*8)
c2.SetLogy()
leg = ROOT.TLegend(0.56,0.73,0.98,0.9)
ROOT.gStyle.SetLegendTextSize(0.0265)

##Draw 1st projection and configure plot settings
proj_G4_QBERTZ.Draw()
#proj_MCNP_BERT_Tot.SetLineColor(1)
proj_G4_QBERTZ.SetMarkerStyle(21)
proj_G4_QBERTZ.SetMarkerColor(2)
proj_G4_QBERTZ.SetMarkerSize(4)
proj_G4_QBERTZ.SetStats(False)
proj_G4_QBERTZ.SetTitle(axis+" Scattered Pencil Beam Density at Target")
proj_G4_QBERTZ.GetXaxis().SetRangeUser(-xlim,xlim)
proj_G4_QBERTZ.GetYaxis().SetRangeUser(ylim[0],ylim[1])
##Make Legend Entry!
leg.AddEntry(proj_G4_QBERTZ,"GEANT4 QGSP_BERTINI_EMZ")

##2nd Projection
proj_MCNP_BERT_Tot.Draw("SAME")
#proj_G4_QBERTZ.SetLineColor(2)
proj_MCNP_BERT_Tot.SetMarkerStyle(41)
proj_MCNP_BERT_Tot.SetMarkerColor(6)
proj_MCNP_BERT_Tot.SetMarkerSize(4)
proj_MCNP_BERT_Tot.SetStats(False)
leg.AddEntry(proj_MCNP_BERT_Tot,"MCNP BERTINI Total Energy Beam")

##2nd Projection
proj_MCNP_BERT_Pri.Draw("SAME")
#proj_G4_QBERTZ.SetLineColor(2)
proj_MCNP_BERT_Pri.SetMarkerStyle(34)
proj_MCNP_BERT_Pri.SetMarkerColor(1)
proj_MCNP_BERT_Pri.SetMarkerSize(mcnpMarkerSize)
proj_MCNP_BERT_Pri.SetStats(False)
leg.AddEntry(proj_MCNP_BERT_Pri,"MCNP BERTINI E > 564MeV Beam")

proj_G4_QFBERTSS.Draw("SAME")
#proj_G4_FBERTZ.SetLineColor(3)
proj_G4_QFBERTSS.SetMarkerStyle(20)
proj_G4_QFBERTSS.SetMarkerColor(3) #8 is mid green, 3 is light green, 7 is light blue
proj_G4_QFBERTSS.SetMarkerSize(4)
proj_G4_QFBERTSS.SetStats(False)
leg.AddEntry(proj_G4_QFBERTSS,"GEANT4 QF_BERTINI__SS 1e6")

proj_G4_EMZ.Draw("SAME")
#proj_G4_EMZ.SetLineColor(4)
proj_G4_EMZ.SetMarkerStyle(47)
proj_G4_EMZ.SetMarkerColor(4)
proj_G4_EMZ.SetMarkerSize(3)
proj_G4_EMZ.SetStats(False)
leg.AddEntry(proj_G4_EMZ,"GEANT4 EMZ")

proj_G4_QFBPEN.Draw("SAME")
#proj_G4_EMZ.SetLineColor(4)
proj_G4_QFBPEN.SetMarkerStyle(47)
proj_G4_QFBPEN.SetMarkerColor(5)
proj_G4_QFBPEN.SetMarkerSize(3)
proj_G4_QFBPEN.SetStats(False)
leg.AddEntry(proj_G4_QFBPEN,"GEANT4 QF_BERT_PEN")

proj_G4_QFBLIV.Draw("SAME")
#proj_G4_EMZ.SetLineColor(4)
proj_G4_QFBLIV.SetMarkerStyle(47)
proj_G4_QFBLIV.SetMarkerColor(7)
proj_G4_QFBLIV.SetMarkerSize(3)
proj_G4_QFBLIV.SetStats(False)
leg.AddEntry(proj_G4_QFBLIV,"GEANT4 QF_BERT_LIV")

##Draw and Print Canvas and Legend
leg.Draw()
c2.Draw()
c2.Print(picpath+"PhysLists_BeamOverlap"+axis+"_Many"+picFormat) #only one is being drawn.
