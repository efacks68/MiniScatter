###xyFits.py
###Plot X and Y projections next to each other for paper Figure 4, with n Gaussian fits

##Function for fitting projections to n Gaussians
def fit_projection(proj,nGauss,p0,p2):
    r=0.1
    ##Set function
    if nGauss == 2:
        func = ROOT.TF1('func','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3]))',-width,width)
    elif nGauss == 3:
        func = ROOT.TF1('func','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3])) + [4] * exp(-x*x/(2*[5]*[5]))',-width,width)
    elif nGauss == 4:
        func = ROOT.TF1('func','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3])) + [4] * exp(-x*x/(2*[5]*[5])) + [6] * exp(-x*x/(2*[7]*[7]))',-width,width)
    elif nGauss == 5:
        func = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3])) + [4] * exp(-x*x/(2*[5]*[5])) + [6] * exp(-x*x/(2*[7]*[7])) + [8] * exp(-x*x/(2*[9]*[9]))',-width,width)

    ##Set parameters. Setting more than needed for a function does nothing
    func.SetParameter(0,p0)               #weight1
    func.SetParLimits(0,p0*(1-r),p0*(1+r))#weight1
    func.SetParameter(1,p2)               #sigma1
    func.SetParLimits(1,p2*(1-r),p2*(1+r))#sigma1
    func.SetParameter(2,p0*r**2)   #weight2
    func.SetParLimits(2,0, p0)     #weight2
    func.SetParameter(3,p2*(1+2*r))#sigma2
    func.SetParLimits(3,p2, p2*20)#sigma2
    func.SetParameter(4,p0*r**2)   #weight2
    func.SetParLimits(4,0, p0)     #weight2
    func.SetParameter(5,p2*(1+3*r))#sigma2
    func.SetParLimits(5,p2, p2*50)#sigma2
    func.SetParameter(6,p0*r**4)   #weight2
    func.SetParLimits(6,0, p0)     #weight2
    func.SetParameter(7,p2*(1+4*r))#sigma2
    func.SetParLimits(7,p2, p2*70)#sigma2
    func.SetParameter(8,p0*r**5)   #weight2
    func.SetParLimits(8,0, p0)     #weight2
    func.SetParameter(9,p2*(1+5*r))#sigma2
    func.SetParLimits(9,p2, p2*1e2)#sigma2

    ##Number of fit points and Fit. Return function and fit results
    func.SetNpx(10000)
    func_res = proj.Fit(func,'RSQ')
    return func, func_res


###Imports and settings
import ROOT
width = 700
xlim = 700
ylim = [1e-7,1.2e-1]
leg = [0.13,0.75,0.475,0.9]
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"

###Get TH2s
##Open file with ROOT
QBERTZ_file="/scratch2/ericdf/PBWScatter/ESS/PBW_570MeV_eX0.00,eY0.00um_bX0.15,bY0.15m_aX0.00,aY0.00_N5e+08_QBZ.root"
fQBERTZ = ROOT.TFile(QBERTZ_file,"r")
##Get TH2D and clone in
G4_QBERTZ_TH2 = fQBERTZ.Get("tracker_cutoff_xy_PDG2212").Clone("QBERTZ")

##Normalize TH2Ds and make Projections
##G4_QBERTZ (QGSP_BERT_EMZ)
integral = G4_QBERTZ_TH2.Integral(G4_QBERTZ_TH2.GetXaxis().FindBin(-width),G4_QBERTZ_TH2.GetXaxis().FindBin(width),G4_QBERTZ_TH2.GetYaxis().FindBin(-width),G4_QBERTZ_TH2.GetYaxis().FindBin(width))
G4_QBERTZ_TH2.Scale(1/integral)
sum = G4_QBERTZ_TH2.Integral()

canvas = ROOT.TCanvas("Scattered Beams","Scattered Beams",0,0,500*8,400*8)
canvas.Divide(2,1) #Divide the canvas into 2 in Horizontal direction
canvas.SetTopMargin(0.3)
#ROOT.gStyle.SetOptStat("iRMe")
#ROOT.gStyle.SetOptFit(000)
ROOT.gStyle.SetStatY(0.9)
ROOT.gStyle.SetStatX(0.9)
#ROOT.gStyle.SetStatW(0.15)
#ROOT.gStyle.SetStatH(0.1)
ROOT.gStyle.SetLegendTextSize(0.017)

##Use TPaveText because it can center the text in the box you define, other options needed manual centering
##Add an overall title for whole figure
ptTitle = ROOT.TPaveText(0,0.92,1,1,"NB")
ptTitle.SetTextAlign(22)
ptTitle.SetTextSize(0.05)
ptTitle.AddText("Scattered Pencil Beam At Target")
ptTitle.SetFillColor(ROOT.kWhite)
ptTitle.Draw()
##Set Title for X Projection
ptX = ROOT.TPaveText(0,0.9,0.5,0.92,"NB")
ptX.SetTextAlign(22)
ptX.SetTextFont(42)
ptX.SetTextSize(0.04)
ptX.AddText("X Projection")
ptX.SetFillColor(ROOT.kWhite)
ptX.Draw()
##Set Title for Y Projection
ptY = ROOT.TPaveText(0.5,0.9,1,0.92,"NB")
ptY.SetTextAlign(22)
ptY.SetTextFont(42)
ptY.SetTextSize(0.04)
ptY.AddText("Y Projection")
ptY.SetFillColor(ROOT.kWhite)
ptY.Draw()


##Left Plot
canvas.cd(1)
legend1 = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])
ROOT.gPad.SetLogy()
projX_G4_QBERTZ = G4_QBERTZ_TH2.ProjectionX("X",G4_QBERTZ_TH2.GetYaxis().FindBin(-width),G4_QBERTZ_TH2.GetYaxis().FindBin(width),"e")
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
projX_G4_QBERTZ.SetName("QBERTZ X Projection")

projX_G4_QBERTZ.Draw()
projX_G4_QBERTZ.GetXaxis().SetRangeUser(-xlim,xlim)
projX_G4_QBERTZ.GetYaxis().SetRangeUser(ylim[0],ylim[1])
projX_G4_QBERTZ.SetMarkerStyle(34)
projX_G4_QBERTZ.SetMarkerColor(1)
projX_G4_QBERTZ.SetMarkerSize(2)
legend1.AddEntry(projX_G4_QBERTZ,"X Projection")

##Fit N Gaussians
## First with gaus and extract parameters for nGauss use
a1x = ROOT.TF1('a1','gaus',-width,width)
a1x.SetNpx(10000)
a1x_res = projX_G4_QBERTZ.Fit(a1x, 'RSQ')
p0 = a1x.GetParameter(0) #weight
p2 = a1x.GetParameter(2) #sigma
a1x.Draw("SAME")
a1x.SetLineColor(1)
##Output sigma of Gaussian in Legend entry
legend1.AddEntry(a1x,"1 Gaussian Fit #sigma = {:.0f}".format(a1x.GetParameter(2)))

##Get N Gaussians
a2x, a2x_res = fit_projection(projX_G4_QBERTZ,2,p0,p2)
a2x.Draw("SAME")
a2x.SetLineColor(2)
legend1.AddEntry(a2x,"2 Gaussian Fit #sigma = {:.0f}, {:.0f}".format(a2x.GetParameter(1),a2x.GetParameter(3)))

a3x, a3x_res = fit_projection(projX_G4_QBERTZ,3,p0,p2)
a3x.Draw("SAME")
a3x.SetLineColor(6)

a4x, a4x_res = fit_projection(projX_G4_QBERTZ,4,p0,p2)
a4x.Draw("SAME")
a4x.SetLineColor(4)

a5x, a5x_res = fit_projection(projX_G4_QBERTZ,5,p0,p2)
a5x.Draw("SAME")
a5x.SetLineColor(3)
legend1.AddEntry(a3x,"3 Gaussian Fit #sigma = {:.0f}, {:.0f}, {:.0f}".format(a3x.GetParameter(1),a3x.GetParameter(5),
                                                                             a3x.GetParameter(3)))
legend1.AddEntry(a4x,"4 Gaussian Fit #sigma = {:.0f}, {:.0f}, {:.0f}, {:.0f}".format(a4x.GetParameter(1),a4x.GetParameter(5),
                                                                                     a4x.GetParameter(7),a4x.GetParameter(3)))
legend1.AddEntry(a5x,"5 Gaussian Fit #sigma = {:.0f}, {:.0f}, {:.0f}, {:.0f}, {:.0f}".format(a5x.GetParameter(1),a5x.GetParameter(5),
                                                                                             a5x.GetParameter(3),a5x.GetParameter(7),
                                                                                             a5x.GetParameter(9)))
legend1.Draw("SAME")
canvas.Update()


##Right Plot
canvas.cd(2)
legend2 = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])
ROOT.gPad.SetLogy()
projY_G4_QBERTZ = G4_QBERTZ_TH2.ProjectionY("Y",G4_QBERTZ_TH2.GetXaxis().FindBin(-width),G4_QBERTZ_TH2.GetXaxis().FindBin(width),"e")
##IMPORTANT: In order for multiple projections to be printed together, each must have a unique name!
projY_G4_QBERTZ.SetName("QBERTZ Y Projection")

projY_G4_QBERTZ.Draw()
projY_G4_QBERTZ.GetXaxis().SetRangeUser(-xlim,xlim)
projY_G4_QBERTZ.GetYaxis().SetRangeUser(ylim[0],ylim[1])
projY_G4_QBERTZ.SetMarkerStyle(34)
projY_G4_QBERTZ.SetMarkerColor(1)
projY_G4_QBERTZ.SetMarkerSize(2)
legend2.AddEntry(projY_G4_QBERTZ,"Y Projection")

##Fit N Gaussians
a1y = ROOT.TF1('a1','gaus',-width,width)
a1y.SetNpx(5000)
a1y_res = projY_G4_QBERTZ.Fit(a1y, 'RSQ')
p0 = a1y.GetParameter(0) #weight
p2 = a1y.GetParameter(2) #sigma
a1y.Draw("SAME")
a1y.SetLineColor(1)
legend2.AddEntry(a1y,"1 Gaussian Fit #sigma = {:.0f}".format(a1y.GetParameter(2)))

a2y, a2y_res = fit_projection(projY_G4_QBERTZ,2,p0,p2)
a2y.Draw("SAME")
a2y.SetLineColor(2)
legend2.AddEntry(a2y,"2 Gaussian Fit #sigma = {:.0f}, {:.0f}".format(a2y.GetParameter(1),a2y.GetParameter(3)))

a3y, a3y_res = fit_projection(projY_G4_QBERTZ,3,p0,p2)
a3y.Draw("SAME")
a3y.SetLineColor(6)

a4y, a4y_res = fit_projection(projY_G4_QBERTZ,4,p0,p2)
a4y.Draw("SAME")
a4y.SetLineColor(4)

a5y, a5y_res = fit_projection(projY_G4_QBERTZ,5,p0,p2)
a5y.Draw("SAME")
a5y.SetLineColor(3)
legend2.AddEntry(a3y,"3 Gaussian Fit #sigma = {:.0f}, {:.0f}, {:.0f}".format(a3y.GetParameter(1),a3y.GetParameter(5),
                                                                             a3y.GetParameter(3)))
legend2.AddEntry(a4y,"4 Gaussian Fit #sigma = {:.0f}, {:.0f}, {:.0f}, {:.0f}".format(a4y.GetParameter(1),a4y.GetParameter(5),
                                                                                     a4y.GetParameter(3),a4y.GetParameter(7)))
legend2.AddEntry(a5y,"5 Gaussian Fit #sigma = {:.0f}, {:.0f}, {:.0f}, {:.0f}, {:.0f}".format(a5y.GetParameter(1),a5y.GetParameter(5),
                                                                                             a5y.GetParameter(3),a5y.GetParameter(7),
                                                                                             a5y.GetParameter(9)))

legend2.Draw("SAME")
canvas.Update()

canvas.Draw()
canvas.Print(picpath+"Figure4.png")

