###mcnpFits.py
###Showing that the MCNP curves are purely Gaussian and compare with Highland scattering estimate.

###Import ROOT and configure
import ROOT
axis = "Y"
xlim = 80
nofits=True
legS= [0.66,0.57,0.985,0.9]
legTextSize=0.043
ylim = [1e-7,1.2e-1]
scratchPath = "/scratch2/ericdf/PBWScatter/"
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"

###Load file and TH2
priFile = ROOT.TFile.Open(scratchPath + "MCNP_Pencil_570MeV_Primary.root","r")
totFile = ROOT.TFile.Open(scratchPath + "MCNP_Pencil_570MeV_Total.root","r")
Tot_TH2 = totFile.Get("MCBERT_Target-TotalProtons").Clone("MCNP_BERT_Total")
Pri_TH2 = priFile.Get("MCBERT_Target-PrimaryProtons").Clone("MCNP_BERT_Primary")

###Test that getting correct TH2s... yes, easy to see since area covered is different!
#c0 = ROOT.TCanvas("MCNP Test","MCNP Test",0,0,100*8,100*8)
#Pri_TH2.Draw()
#Tot_TH2.Draw()
#c0.Print(picpath+"test.png")

###Normalize the TH2s if not already.
width = 200
maxim = 200
#integralTot = Tot_TH2.Integral(Tot_TH2.GetXaxis().FindBin(-width),Tot_TH2.GetXaxis().FindBin(width),Tot_TH2.GetYaxis().FindBin(-width),Tot_TH2.GetYaxis().FindBin(width),option="width")
#Tot_TH2.Scale(1/integralTot)
#sumTot = Tot_TH2.Integral()

#integralPri = Pri_TH2.Integral(Pri_TH2.GetXaxis().FindBin(-width),Pri_TH2.GetXaxis().FindBin(width),Pri_TH2.GetYaxis().FindBin(-width),Pri_TH2.GetYaxis().FindBin(width),option="width")
#Pri_TH2.Scale(1/integralPri)
#sumPri = Pri_TH2.Integral()
#print("Pri",Pri_TH2.Integral(Pri_TH2.GetXaxis().FindBin(-width),Pri_TH2.GetXaxis().FindBin(width),Pri_TH2.GetYaxis().FindBin(-width),Pri_TH2.GetYaxis().FindBin(width),option="width"))

#GEANT4 too
QBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
fQBERTZ = ROOT.TFile(QBERTZ_file,"r")
##Get TH2D and clone in
G4_QBERTZ_TH2 = fQBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("QBERTZ")

###Open Canvas and configure, select stats output
c1 = ROOT.TCanvas("MCNP Beam Fits","MCNP Beam Fits",0,0,400*8,250*8)
ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetOptFit(000)
#ROOT.gStyle.SetStatW(0.15)
#ROOT.gStyle.SetStatH(0.1)
c1.SetLogy()
leg = ROOT.TLegend(legS[0],legS[1],legS[2],legS[3])#0.13,0.65,0.43,0.9)
ROOT.gStyle.SetLegendTextSize(legTextSize)
leg.SetMargin(0.12)


###Project onto the desired axis and Normalize
if axis in {"Y","y"}:
    proj_Pri = Pri_TH2.ProjectionY(axis,Pri_TH2.GetXaxis().FindBin(-width),Pri_TH2.GetXaxis().FindBin(width),"e")
    integralPri = proj_Pri.Integral(proj_Pri.GetXaxis().FindBin(-width),proj_Pri.GetXaxis().FindBin(width),option="width")
    proj_Pri.Scale(1/integralPri)
else:
    proj_Pri = Pri_TH2.ProjectionX(axis,Pri_TH2.GetYaxis().FindBin(-width),Pri_TH2.GetYaxis().FindBin(width),"e")
    integralPri = proj_Pri.Integral(proj_Pri.GetXaxis().FindBin(-width),proj_Pri.GetXaxis().FindBin(width),option="width")
    proj_Pri.Scale(1/integralPri)
proj_Pri.SetName("proj_PrimaryVeryDifferent")
proj_Pri.Draw()
proj_Pri.SetMarkerStyle(41)
proj_Pri.SetMarkerColor(2)
proj_Pri.SetMarkerSize(4)
proj_Pri.SetTitle(axis+" Scattered Pencil Beam Density at Target")
proj_Pri.GetXaxis().SetRangeUser(-xlim,xlim)
proj_Pri.GetYaxis().SetRangeUser(ylim[0],ylim[1])

#MCNP_Tot
##Apparently they need to be separated to not conflict, even though all variables are separate and THs get unique names...
if axis in {"y","Y"}:
    proj_Tot = Tot_TH2.ProjectionY(axis,Tot_TH2.GetXaxis().FindBin(-width),Tot_TH2.GetXaxis().FindBin(width),"e")
    integralTot = proj_Tot.Integral(proj_Tot.GetXaxis().FindBin(-width),proj_Tot.GetXaxis().FindBin(width),option="width")
    proj_Tot.Scale(1/integralTot)
else:
    proj_Tot = Tot_TH2.ProjectionX(axis,Tot_TH2.GetYaxis().FindBin(-width),Tot_TH2.GetYaxis().FindBin(width),"e")
    integralTot = proj_Tot.Integral(proj_Tot.GetXaxis().FindBin(-width),proj_Tot.GetXaxis().FindBin(width),option="width")
    proj_Tot.Scale(1/integralTot)
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
## The names are NOT inherited from their TH2D, like I had expected.
proj_Tot.SetName("proj_Tot")
#proj_Tot.Draw("SAME")
###Configure plots
proj_Tot.SetMarkerStyle(34)
proj_Tot.SetMarkerColor(1)
proj_Tot.SetMarkerSize(4)

if axis in {"Y","y"}:
    proj_QBZ = G4_QBERTZ_TH2.ProjectionY(axis,G4_QBERTZ_TH2.GetXaxis().FindBin(-width),G4_QBERTZ_TH2.GetXaxis().FindBin(width),"e")
    integralQBZ = proj_QBZ.Integral(proj_QBZ.GetXaxis().FindBin(-width),proj_QBZ.GetXaxis().FindBin(width),option="width")
    proj_QBZ.Scale(1/integralQBZ)
else:
    proj_QBZ = G4_QBERTZ_TH2.ProjectionX(axis,G4_QBERTZ_TH2.GetYaxis().FindBin(-width),G4_QBERTZ_TH2.GetYaxis().FindBin(width),"e")
    integralQBZ = proj_QBZ.Integral(proj_QBZ.GetXaxis().FindBin(-width),proj_QBZ.GetXaxis().FindBin(width),option="width")
    proj_QBZ.Scale(1/integralQBZ)
proj_QBZ.SetName("proj_QBERTZ")
proj_QBZ.Draw("SAME")
proj_QBZ.SetMarkerStyle(39)
proj_QBZ.SetMarkerColor(1)
proj_QBZ.SetMarkerSize(4)

###Fit projections to Gaussians
if nofits:
    fitOption='RSQ0'
    ext="tight_nofits"
else:
    fitOptions='RSQ'
    ext="tight_simplified"
f1 = ROOT.TF1('f1','gaus',-maxim,maxim)
f1.SetNpx(5000)
f1.SetLineColor(1)
#f1Tot_res = proj_Tot.Fit(f1, 'RSQ')

f2 = ROOT.TF1('f2','gaus',-maxim,maxim)
f2.SetLineColor(2)
f2.SetLineWidth(3)
f2.SetNpx(5000)
f2Pri_res = proj_Pri.Fit(f2, 'RSQ')

f3 = ROOT.TF1('f3','gaus',-maxim,maxim)
f3.SetLineColor(1)
f3.SetLineWidth(3)
f3.SetNpx(5000)
f3Pri_res = proj_QBZ.Fit(f3, 'RSQ')

###Import data from Highland Scattering script
nParticles=1e7
filename = "HighlandScatteredDistribution"
rootFile="/scratch2/ericdf/PBWScatter/"+filename+"_{:.0e}".format(nParticles)+".root"
myFile = ROOT.TFile.Open(rootFile,"r")
hTarget = myFile.Get("Highland Scattered Proton Distribution").Clone("Highland Scattered Proton Distribution")
integralTar = hTarget.Integral(hTarget.GetXaxis().FindBin(-width),hTarget.GetXaxis().FindBin(width),option="width")
hTarget.Scale(1/integralTar)
#print("Targ",hTarget.Integral(hTarget.GetXaxis().FindBin(-width),hTarget.GetXaxis().FindBin(width),option="width"))
hTarget.Draw("SAME")
hTarget.SetMarkerStyle(5)
hTarget.SetMarkerColor(4)
hTarget.SetMarkerSize(4)

##Fit distribution
f4 = ROOT.TF1('Highland','gaus',-maxim,maxim)
f4.SetLineColor(4)
f4.SetLineWidth(3)
f4.SetNpx(5000)
f4_res = hTarget.Fit(f4, 'RSQ')

#Set legend entries
#leg.AddEntry(proj_Tot,"MCNP_BERTINI Total Energy Beam")#; #sigma = {:.3f}".format(f1.GetParameter(2)))
leg.AddEntry(proj_QBZ,"GEANT4 QGSP_BERTINI")#; #sigma = {:.3f}".format(f3.GetParameter(2)))
leg.AddEntry(ROOT.nullptr,"  E>564MeV, #sigma={:.1f}mm".format(f3.GetParameter(2)),"")
leg.AddEntry(proj_Pri,"MCNP_BERTINI") #; #sigma = {:.3f}".format(f2.GetParameter(2)))
leg.AddEntry(ROOT.nullptr,"  E>564MeV, #sigma={:.1f}mm".format(f2.GetParameter(2)),"")
leg.AddEntry(hTarget,"Highland Scattered")#; #sigma = {:.3f}".format(f4.GetParameter(2)))
leg.AddEntry(ROOT.nullptr,"  E=570MeV, #sigma={:.1f}mm".format(f4.GetParameter(2)),"")

###Draw and Print
c1.Draw()
leg.Draw()
c1.Print(picpath+"MCNPFits"+axis+ext+".png") #nofits