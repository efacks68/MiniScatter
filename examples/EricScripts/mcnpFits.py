###mcnpFits.py
###Showing that the MCNP curves are purely Gaussian and compare with Highland scattering estimate.

###Import ROOT and configure
import ROOT
axis = "X"
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
integralTot = Tot_TH2.Integral(Tot_TH2.GetXaxis().FindBin(-width),Tot_TH2.GetXaxis().FindBin(width),Tot_TH2.GetYaxis().FindBin(-width),Tot_TH2.GetYaxis().FindBin(width))
Tot_TH2.Scale(1/integralTot)
sumTot = Tot_TH2.Integral()

integralPri = Pri_TH2.Integral(Pri_TH2.GetXaxis().FindBin(-width),Pri_TH2.GetXaxis().FindBin(width),Pri_TH2.GetYaxis().FindBin(-width),Pri_TH2.GetYaxis().FindBin(width))
Pri_TH2.Scale(1/integralPri)
sumPri = Pri_TH2.Integral()

###Open Canvas and configure, select stats output
c1 = ROOT.TCanvas("MCNP Beam Fits","MCNP Beam Fits",0,0,500*8,400*8)
ROOT.gStyle.SetOptStat(000)
ROOT.gStyle.SetOptFit(000)
#ROOT.gStyle.SetStatW(0.15)
#ROOT.gStyle.SetStatH(0.1)
c1.SetLogy()
leg = ROOT.TLegend(0.65,0.6,0.98,0.9)#0.13,0.65,0.43,0.9)
xlim = 200
ylim = [1e-7,1.2e-1]

###Project onto the desired axis
if axis in {"y","Y"}:
    proj_Tot = Tot_TH2.ProjectionY(axis,Tot_TH2.GetXaxis().FindBin(-width),Tot_TH2.GetXaxis().FindBin(width),"e")
else:
    proj_Tot = Tot_TH2.ProjectionX(axis,Tot_TH2.GetYaxis().FindBin(-width),Tot_TH2.GetYaxis().FindBin(width),"e")
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
## The names are NOT inherited from their TH2D, like I had expected.
proj_Tot.SetName("proj_Tot")
proj_Tot.Draw()
leg.AddEntry(proj_Tot,"MCNP_BERTINI Total Energy Beam")

##Apparently they need to be separated to not conflict, even though all variables are separate and THs get unique names...
if axis in {"Y","y"}:
    proj_Pri = Pri_TH2.ProjectionY(axis,Pri_TH2.GetXaxis().FindBin(-width),Pri_TH2.GetXaxis().FindBin(width),"e")
else:
    proj_Pri = Pri_TH2.ProjectionX(axis,Pri_TH2.GetYaxis().FindBin(-width),Pri_TH2.GetYaxis().FindBin(width),"e")
proj_Pri.SetName("proj_PrimaryVeryDifferent")
proj_Pri.Draw("SAME")
leg.AddEntry(proj_Pri,"MCNP_BERTINI Beam E > 564MeV")

###Fit projections to Gaussians
f1 = ROOT.TF1('f1','gaus',-maxim,maxim)
f1.SetNpx(5000)
f1.SetLineColor(ROOT.kBlue)
f1Tot_res = proj_Tot.Fit(f1, 'RSQ')
##Can add parameter value to the legend entry!
leg.AddEntry(f1,"Total Energy Beam Fit; #sigma = {:.3f}".format(f1.GetParameter(2)))

f2 = ROOT.TF1('f2','gaus',-maxim,maxim)
f2.SetLineColor(ROOT.kGreen)
f2.SetNpx(5000)
f2Pri_res = proj_Pri.Fit(f2, 'RSQ')
leg.AddEntry(f2,"Beam E > 564MeV Fit; #sigma = {:.3f}".format(f2.GetParameter(2)))

###Configure plots
proj_Tot.SetMarkerStyle(34)
proj_Tot.SetMarkerColor(1)
proj_Tot.SetMarkerSize(3)
proj_Tot.SetTitle(axis+" Scattered Pencil Beam Density at Target")
proj_Tot.GetXaxis().SetRangeUser(-xlim,xlim)
proj_Tot.GetYaxis().SetRangeUser(ylim[0],ylim[1])

proj_Pri.SetMarkerStyle(41)
proj_Pri.SetMarkerColor(2)
proj_Pri.SetMarkerSize(3)

###Import data from Highland Scattering script
nParticles=1e7
filename = "HighlandScatteredDistribution"
rootFile="/scratch2/ericdf/PBWScatter/"+filename+"_{:.0e}".format(nParticles)+".root"
myFile = ROOT.TFile.Open(rootFile,"r")
hTarget = myFile.Get("Highland Scattered Proton Distribution").Clone("Highland Scattered Proton Distribution")
hTarget.Draw("SAME")
leg.AddEntry(hTarget,"Highland Scattered Proton Distribution")
##Fit distribution
maxim = 50
f3 = ROOT.TF1('Highland','gaus',-maxim,maxim)
f3.SetNpx(5000)
f3_res = hTarget.Fit(f3, 'RSQ')
leg.AddEntry(f3,"Highland Scattered Fit; #sigma = {:.3f}".format(f3.GetParameter(2)))

###Draw and Print
c1.Draw()
leg.Draw()
c1.Print(picpath+"MCNPFits"+axis+".png")