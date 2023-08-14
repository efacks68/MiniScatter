####physListCompare.py
import ROOT
from plotFit import gaussianFit
##which axis to project?
column = 2 #columns: 2= p>564MeV; 4=p<564MeV; 10=total
picFormat = ".png" #".png",".pdf", ".jpeg"
xlim = 550
g=3 #number of Gaussians to fit to data
##Settings for reading
path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
file = "MCNP_ProtonDensity_570MeV_10June"
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
scratchPath = "/scratch2/ericdf/PBWScatter/"
ylim = [1e-7,1e-1]
markerSize = 3
width=500 #width of map to sum over

###Load in MCNP data
##Primary Beam (E>564MeV)
fMCNP_BERT_Pri = ROOT.TFile(scratchPath+"MCNP_Pencil_570MeV_Primary.root","r")
MCNP_BERT_Pri = fMCNP_BERT_Pri.Get("MCBERT_Target-PrimaryProtons").Clone("MCNP_BERT_Primary")


###Load in all scattered beams of different GEANT4 phys lists and fit
##Open file with ROOT
##QGSP_BERT_EMZ
QBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
fQBERTZ = ROOT.TFile(QBERTZ_file,"r")
##Get TH2D and clone in
G4_QBERTZ_TH2 = fQBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("QBERTZ")

##MCNP_BERT_Pri
##Make Projection and Normalize
projY_MCNP_BERT_Pri = MCNP_BERT_Pri.ProjectionY("Y",MCNP_BERT_Pri.GetXaxis().FindBin(-width),MCNP_BERT_Pri.GetXaxis().FindBin(width),"e")
integralY = projY_MCNP_BERT_Pri.Integral(projY_MCNP_BERT_Pri.GetXaxis().FindBin(-width),projY_MCNP_BERT_Pri.GetXaxis().FindBin(width),option="width")
projY_MCNP_BERT_Pri.Scale(1/integralY)
sum = projY_MCNP_BERT_Pri.Integral()
projX_MCNP_BERT_Pri = MCNP_BERT_Pri.ProjectionX("X",MCNP_BERT_Pri.GetYaxis().FindBin(-width),MCNP_BERT_Pri.GetYaxis().FindBin(width),"e")
integralX = projX_MCNP_BERT_Pri.Integral(projX_MCNP_BERT_Pri.GetXaxis().FindBin(-width),projX_MCNP_BERT_Pri.GetXaxis().FindBin(width),option="width")
projX_MCNP_BERT_Pri.Scale(1/integralX)
sum = projX_MCNP_BERT_Pri.Integral()
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
## The names are NOT inherited from their TH2D, like I had expected.
projY_MCNP_BERT_Pri.SetName("projY_MCBERT_Pri")
projX_MCNP_BERT_Pri.SetName("projX_MCBERT_Pri")

##G4_QBERTZ (QGSP_BERT_EMZ)
projY_G4_QBERTZ = G4_QBERTZ_TH2.ProjectionY("Y",G4_QBERTZ_TH2.GetXaxis().FindBin(-width),G4_QBERTZ_TH2.GetXaxis().FindBin(width),"e")
integralY = projY_G4_QBERTZ.Integral(projY_G4_QBERTZ.GetXaxis().FindBin(-width),projY_G4_QBERTZ.GetXaxis().FindBin(width),option="width")
projY_G4_QBERTZ.Scale(1/integralY)
sum = projY_G4_QBERTZ.Integral()
projX_G4_QBERTZ = G4_QBERTZ_TH2.ProjectionX("X",G4_QBERTZ_TH2.GetYaxis().FindBin(-width),G4_QBERTZ_TH2.GetYaxis().FindBin(width),"e")
integralX = projX_G4_QBERTZ.Integral(projX_G4_QBERTZ.GetXaxis().FindBin(-width),projX_G4_QBERTZ.GetXaxis().FindBin(width),option="width")
projX_G4_QBERTZ.Scale(1/integralX)
sum = projX_G4_QBERTZ.Integral()
projY_G4_QBERTZ.SetName("projY_G4QBERTZ")
projX_G4_QBERTZ.SetName("projX_G4QBERTZ")


###Plot them all together with different colors
##Make Canvas and Lengend
c2 = ROOT.TCanvas("Scattered Beams","Scattered Beams",0,0,500*8,400*8)
c2.SetLogy()
#leg = ROOT.TLegend(0.58,0.75,0.98,0.9)
#ROOT.gStyle.SetLegendTextSize(0.024)

##Draw 1st projection and configure plot settings
projY_G4_QBERTZ.Draw()
#proj_G4_QBERTZ.SetLineColor(2)
projY_G4_QBERTZ.SetMarkerStyle(34)#39
projY_G4_QBERTZ.SetMarkerColor(4)#8
projY_G4_QBERTZ.SetMarkerSize(markerSize+2)
projY_G4_QBERTZ.SetStats(False)
projY_G4_QBERTZ.SetTitle("Scattered Pencil Beam at Target Comparison")
projY_G4_QBERTZ.GetXaxis().SetRangeUser(-xlim,xlim)
projY_G4_QBERTZ.GetYaxis().SetRangeUser(ylim[0],ylim[1])
projY_G4_QBERTZ.GetXaxis().SetTitle("Position [mm]")
##Make Legend Entry!
#leg.AddEntry(projY_G4_QBERTZ,"GEANT4 QGSP_BERTINI_EMZ Y")

##2nd Projection
#projX_MCNP_BERT_Pri.Draw("SAME")
##proj_G4_QBERTZ.SetLineColor(2)
#projX_MCNP_BERT_Pri.SetMarkerStyle(41)
#projX_MCNP_BERT_Pri.SetMarkerColor(2)
#projX_MCNP_BERT_Pri.SetMarkerSize(markerSize)
#leg.AddEntry(projX_MCNP_BERT_Pri,"MCNP BERTINI E > 564MeV Beam X")

##3rd Projection
#projY_MCNP_BERT_Pri.Draw("SAME")
##proj_MCNP_BERT_Tot.SetLineColor(1)
#projY_MCNP_BERT_Pri.SetMarkerStyle(34)
#projY_MCNP_BERT_Pri.SetMarkerColor(1)
#projY_MCNP_BERT_Pri.SetMarkerSize(markerSize+1)
#leg.AddEntry(projY_MCNP_BERT_Pri,"MCNP BERTINI E > 564MeV Beam Y")

##4th Projection
#projX_G4_QBERTZ.Draw("SAME")
##proj_G4_QBERTZ.SetLineColor(2)
#projX_G4_QBERTZ.SetMarkerStyle(41)
#projX_G4_QBERTZ.SetMarkerColor(4)
#projX_G4_QBERTZ.SetMarkerSize(markerSize)
#leg.AddEntry(projX_G4_QBERTZ,"GEANT4 QGSP_BERTINI_EMZ X")

##Draw and Print Canvas and Legend
#leg.Draw()
c2.Draw()
c2.Print(picpath+"Figure6.pdf")#"PhysLists_BeamOverlap_MCNPPri-G4QBZ"+picFormat) #only one is being drawn.
c2.Print(picpath+"Figure6.png")
