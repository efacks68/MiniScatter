#plotFit.py
#29 March- Present

#This holds the functions to plot figures for RasterScatter

##For simple material-variant plots
def plotFit(xs,ys,savename,xlim,ylim,material,thick): 
    import numpy as np
    import matplotlib.pyplot as plt
    #from plotFit import findFit
    #For creating the PDF of particle distribution from findFit function optimized mu and sigma.
    from scipy.stats import norm 

    fs = 14 #set the axis label font size early

    #Use Scipy.optimize.curve_fit in my findFit function to get mus and sigmas:
    mux, sigmax, amplx, xinterval = findFit(xs,[0.01,0.1,10],(0,[1e3,3e5,3e5])) #dynamically gets parameters AND histogram intervals!
    muy, sigmay, amply, yinterval = findFit(ys,[0.01,0.1,10],(0,[1e3,3e5,3e5]))

    #Find range of particles that are outside 3 sigma
    sigx=np.abs(xs)>3*sigmax# and np.abs(xs)<10*sigma)
    sigy=np.abs(ys)>3*sigmay

    #Find % of particles outside 3 sigma
    pOut3sigx = len(xs[sigx])/len(xs)*100
    pOut3sigy = len(ys[sigy])/len(ys)*100
    #print(pOut3sigx,"% outside 3 sigma")
    #print(pOut3sigy,"% outside 3 sigma")

    #Create the fig with 2 plots side by side
    plt.clf()
    fig = plt.figure(figsize=(15,5))
    plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
    s1 = fig.add_subplot(1,2,1)
    s2 = fig.add_subplot(1,2,2)
    
    #Define material short name for plot labeling
    if material == "G4_Galactic":
        mat = "Vacuum"
    elif material == "G4_Al":
        mat ="Al"
    elif material == "G4_AIR":
        mat ="Air"
    elif material == "G4_Au":
        mat = "Au"

    if thick == 0:
        thick = 4.25
        mat = "Real"

    #Make the histogram of the full energy distrubtion for X with findFit intervals
    nx, binsx, patchesx = s1.hist(xs, xinterval, density=True, facecolor="green", alpha=0.75)
    
    #Add the "best fit" line using the earlier findFit mus and sigmas for X
    y1 = norm.pdf(binsx, mux, sigmax)
    l1 = s1.plot(binsx, y1, "r--", linewidth=2) #important to keep it as l# or it won't work

    #Make the histogram of the full energy distrubtion for Y with findFit intervals
    ny, binsy, patchesy = s2.hist(ys, yinterval, density=True, facecolor="green", alpha=0.75)

    #Add the "best fit" line using the earlier findFit mus and sigmas for Y
    y2 = norm.pdf(binsy, muy, sigmay)
    l2 = s2.plot(binsy, y2, "r--", linewidth=2)

    #change limits based on received parameters
    if ylim !=0: #for halo zoom in
        s1.set_ylim([0,ylim])
        s2.set_ylim([0,ylim])
        savename=savename+"_halo"
    if xlim <= 10:
        s1.set_xlim([-xlim*sigmax,xlim*sigmax])
        s2.set_xlim([-xlim*sigmax,xlim*sigmax])
        savename = savename + "_" + str(xlim) + "sigma"
    elif xlim > 10: #allows dynamic, expected case of ~4-500
        s1.set_xlim([-xlim,xlim])
        s2.set_xlim([-xlim,xlim])
        savename = savename + "_extended"
    else: #for standard case
        s1.set_xlim([-xlim,xlim])
        s2.set_xlim([-xlim,xlim])

    #Set Plot characterists
    #s1.set_title("After "+mat+", Fitted X Distribution \n"+rf"$\mu= {{:.1e}}, \sigma = {{:.3f}}$,".format(mux,sigmax)+
    #        " {:.3f}% outisde".format(pOut3sigx)+r" 3$\sigma$",fontsize=fs)
    #s1.set_title("At {:.2f}mm ".format(thick)+mat+" PBW Exit, X Distribution \n"+rf"$\sigma = {{:.3f}}$,".format(sigmax)+
    #        " {:.3f}% outside".format(pOut3sigx)+r" 3$\sigma$",fontsize=fs)
    #s1.set_xlabel("X Position [mm]",fontsize=fs)
    #s1.set_ylabel("Probability Density",fontsize=fs)
    plt.setp(s1,title="At {:.2f}mm ".format(thick)+mat+" PBW Exit, X Distribution \n"+rf"$\sigma = {{:.3f}}$,".format(sigmax)+
            " {:.3f}% outside".format(pOut3sigx)+r" 3$\sigma$",xlabel="X Position [mm]",ylabel="Probability Density")

    #s2.set_title("After "+mat+", Fitted Y Distribution \n"+rf"$\mu= {{:.1e}}, \sigma = {{:.3f}}$,".format(muy,sigmay)+
    #        " {:.3f}% outisde".format(pOut3sigy)+r" 3$\sigma$",fontsize=fs)
    #s2.set_title("At {:.2f}mm ".format(thick)+mat+" PBW Exit, Y Distribution \n"+rf"$\sigma = {{:.3f}}$,".format(sigmay)+
    #        " {:.3f}% outside".format(pOut3sigy)+r" 3$\sigma$",fontsize=fs)
    #s2.set_xlabel("Y Position [mm]",fontsize=fs)
    #s2.set_ylabel("Probability Density",fontsize=fs)
    plt.setp(s2,title="At {:.2f}mm ".format(thick)+mat+" PBW Exit, Y Distribution \n"+rf"$\sigma = {{:.3f}}$,".format(sigmay)+
            " {:.3f}% outside".format(pOut3sigy)+r" 3$\sigma$",xlabel="Y Position [mm]",ylabel="Probability Density")

    #Can date stamp the multi plot for easier tracking of changes, if necessary
    #from datetime import datetime
    #dt = datetime.now()

    #Update Name
    #name = savename+"_"+dt.strftime("%H-%M-%S")+".png"  
    name = savename+".png" #I'd prefer to have it overwrite old files for now
    print(name) #show so it is easy to find the file
    plt.tight_layout()
    plt.savefig(name,bbox_inches='tight',dpi=500)
    #plt.show()
    plt.close() #be sure to close the plot













##Gaussian function for findFit
def gaussian(x,amplitude,mu,sigma):
    from numpy import exp as npExp
    return amplitude * npExp( - (x - mu) ** 2 / (2 * sigma ** 2))













##Uses SciPy.optimize.curve_fit to fit a Gaussian to a curve
def findFit(data,guess,lims,nBins):
    ##with help from MSeifert in stackoverflow fit a curve to a histogram in python
    ## https://stackoverflow.com/questions/35544233/fit-a-curve-to-a-histogram-in-python
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    #from plotFit import gaussian

    #print(data)
    #print("here",np.shape(data))
    #bin_heights, bin_borders, _ = plt.hist(data,bins=nBins,label="histogram")
    bin_heights, bin_borders = np.histogram(data,bins=nBins)
    pwd="/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    #plt.savefig(pwd+"hist.png"); print(pwd+"hist.png")
    #print(len(bin_borders))
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    #should weight with 1/sqrt(N) check errors should be rel or abs. curve fit will take either cov matrix of input (sigma^2) could take diagonal if no correlation (1/sigma^2(weight)).
    #  need to make array of sigmas 
    bin_sigmas=np.sqrt(bin_heights)
    #for idx, val in np.ndenumerate(bin_sigmas):
    #    if val == 0: bin_sigmas[idx] = 1e-6

    ##Exclude bins with no data, like ROOT chi2 fitting    
    bin_nonzero_idxs = bin_heights > 0
    bin_heights = bin_heights[bin_nonzero_idxs]
    bin_centers = bin_centers[bin_nonzero_idxs]
    bin_sigmas  = bin_sigmas [bin_nonzero_idxs]

    popt,_ = curve_fit(gaussian, bin_centers, bin_heights,sigma=bin_sigmas, p0=guess,bounds=lims) #p0 should be good start
    x_interval_for_fit = np.linspace(bin_borders[0],bin_borders[-1],len(bin_borders))
    #plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit,*popt),label="fit")
    #plt.legend()
    #plt.xlim([bin_borders[0],bin_borders[-1]])
    #plt.show()
    #pwd="/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    #plt.savefig(pwd+"fit.png")
    #plt.close() #be sure to close this or else these show up in the multi plots!
    #print(popt)
    #return the mean, abs(sigma), interval #sigma can sometimes be - so must abs() it.
        #   mu,     sigma,       ampl
    return(popt[1],abs(popt[2]),popt[0],x_interval_for_fit) #for Vac and Air sigma sometimes give - number...












##Take in x,y 1D arrays and calculate covariance, returning Twiss in form [beta,alpha,gemt,gamma]
def calcTwiss(labelx,labely,x,y):
    #Purely get covariance of the distributions
    #use generic x,y but usually x, x' or y,y'!
    import numpy as np
    #beta_rel = 0.948
    #gamma_rel = 3.132
    #get Twiss from np.cov
    xcov = np.cov(x,y)
    #print(labelx,labely,xcov,"\n")
    #calculated Twiss parameters
    gemt = np.sqrt(np.linalg.det(xcov)) 
    beta = xcov[0,0]/gemt
    alpha = -xcov[0,1]/gemt
    gamma = xcov[1,1]/gemt
    #nemt = gemt*beta_rel*gamma_rel
    #print(labelx,labely,":","\nGeometric Emittance: {:.3e}m, Beta: {:.3e}m, Alpha: {:.3e}\n".format(gemt,beta,alpha,gamma))
    return [beta,alpha,gemt,gamma]












##Calculate new Twiss after a thin scatterer using Mueller Eq 8
def calcEq8(thetasq,Twiss,thick,beta_rel,gamma_rel):
    #import numpy as np
    m=1
    um=1e-6
    mm=1e-3
    #print("init",Twiss)
    #Twiss=[beta,alpha,gemt,gamma]
    #Calculations from Eq 7 and 8 from Twiss MCS Formalism Calculations from https://cds.cern.ch/record/499590
    e8dgem = 0.5 * Twiss[0]*m * thetasq * thetasq #[m*rad^2]
    e8alph = Twiss[2] * Twiss[1] / (Twiss[2] + e8dgem)
    e8beta = Twiss[2] * Twiss[0]*m / (Twiss[2] + e8dgem) #[m]
    e8gamma = (Twiss[2] * Twiss[3] + thetasq * thetasq ) / (Twiss[2] + e8dgem) #[m^-1]
    e8gemt = Twiss[2] + e8dgem
    #e8nemt = e8nemt/(beta_rel*gamma_rel)
    #print("e8",thick,thetasq,e8dgem,e8beta,e8alph,e8gemt)
    #28.7-supposed to have thetasq*thetasq in gamma. Previously did NOT have it! Now numbers are great!
    return [e8beta,e8alph,e8gemt,e8gamma]











##Calculate new Twiss after a thin scatterer using Mueller Eq 16; Not working Oct22-Jun23
def calcEq16(thetasq,Twiss,thick,beta_rel,gamma_rel):
    #import numpy as np
    m=1
    um=1e-6
    mm=1e-3
    #print("init",Twiss)
    #Twiss=[beta,alpha,gemt,gamma]
    #Calculations from Eq 15 and 16 from Twiss MCS Formalism Calculations from https://cds.cern.ch/record/499590
    e16dgem = 0.5 * thetasq * thetasq * (Twiss[0]*m + thick*mm * Twiss[1] + thick*mm*thick*mm/3 * Twiss[3]) #[m*rad^2]
    e16alph = (Twiss[2] * Twiss[1] - thick*mm * 0.5 * thetasq * thetasq ) / (Twiss[2] + e16dgem)
    e16beta = (Twiss[2] * Twiss[0]*m + thick*mm * thick*mm / 3 * thetasq * thetasq ) / (Twiss[2] + e16dgem) #[m]
    e16gamma = (Twiss[2] * Twiss[3] + thetasq * thetasq ) / (Twiss[2] + e16dgem) #[m^-1]
    e16gemt = Twiss[2]*um + e16dgem
    #print("e16",thick,thetasq,e16dgem,e16beta,e16alph,e16gemt)
    #28.7-supposed to have thetasq*thetasq in gamma, alph and beta! Previously did NOT have it! Now numbers are great!
    return [e16beta,e16alph,e16gemt,e16gamma]












##Calculates RMS from the given Twiss, returns position and momentum RMSes
def getMoments(Twiss):#[beta,alpha,gemt,gamma]
    from numpy import sqrt
    mm=1e-3
    sigma = sqrt(Twiss[2]*Twiss[0]) /mm #sigma x rms = beam size = sqrt(beta*epsG)
    sigmapx = sqrt(Twiss[2]*Twiss[3]) #mrad actually
    #print(sigma,sigmapx,sigpx)
    return [sigma,sigmapx]












##Plots position and momentum projections of the xs, pxs arrays and fits Mueller Eqns to them (obsolete)
def plotTwissFit(xs,pxs,savename,mat,titledescr,axis,thick,thetasq,beta_rel,gamma_rel,TwissI):
    #16.8.22 No longer useful, as it doesn't use layer by layer #19.2.2023-?
    #For plotting particle Distribution at Target Exit with PDF from Twiss parameters.
    import numpy as np
    import matplotlib.pyplot as plt
    #from plotFit import findFit, calcTwiss, getMoments
    from scipy.stats import norm 
    mm = 1e-3 #the position numbers are in mm so must convert the PDF parameters to mm!
    urad = 1e-6
    fs = 14

    #Recalculate Twiss, 
    #calcTwiss => [beta,alpha,gemt,gamma]
    calcTwf = calcTwiss("Recheck "+axis,"Recheck "+axis+"'",xs*mm,pxs) 
    #Calculate Moments for calculated Twiss
    #getMoments => [sigma,sigmapx]
    Mcal = getMoments(calcTwf) #this is sample variance, more sensitive to outliers.
    Me8 = getMoments(calcEq8(thetasq,TwissI,thick,beta_rel,gamma_rel)) #must send ITwiss!!!!
    PDFlabele8 = "Eq 8"
    Me16 = getMoments(calcEq16(thetasq,TwissI,thick,beta_rel,gamma_rel))
    PDFlabele16 = "Eq 16"
    #print(Mcal,Me8,Me16)

    print("Recalculated Twiss:","Betax: {:.2f}m, alphax: {:.2f}, Geo Emitt x: {:.3e}m".format(calcTwf[0],calcTwf[1],calcTwf[2]))

    #Get intervals and info about distribution
    mux,sigmax,amplx,xinterval =    findFit(xs, [0.01,0.1,10],(0,[1,10,100]),"auto")
    mupx,sigmapx,amply,pxinterval = findFit(pxs,[0.01,0.1,10],(0,[1,10,100]),"auto") #this is a least square fit to Gaussian, so less sensitive to outliers.
    print("Histogram   sigmaTwx: {:.2f}mm, sigmaTwX': {:.2f}mrad".format(sigmax,sigmapx/mm))
    print("Calculated  sigmaTwx: {:.2f}mm, sigmaTwX': {:.2f}mrad".format(Mcal[0],Mcal[1]/mm))
    print("Equation 8  sigmaTwx: {:.2f}mm, sigmaTwX': {:.2f}mrad".format(Me8[0],Me8[1]/mm))
    print("Equation 16 sigmaTwx: {:.2f}mm, sigmaTwX': {:.2f}mrad".format(Me16[0],Me16[1]/mm))

    #Create the fig with 2 plots side by side
    plt.clf()
    fig = plt.figure(figsize=(16,6.0))
    plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
    s1 = fig.add_subplot(1,2,1)
    s2 = fig.add_subplot(1,2,2)

    #Make the histogram of the full energy distrubtion for X with findFit intervals
    nx, binsx, patchesx = s1.hist(xs, xinterval, density=True, facecolor="green", alpha=0.75,label="MiniScatter "+axis+" Histogram")
    npx,binspx, patchespx = s2.hist(pxs,pxinterval,density=True, facecolor="green", alpha=0.75,label="MiniScatter "+axis+"' Histogram")

    if "Pre" in titledescr:
        #only have hist and fit
        y1a = norm.pdf(binsx, mux, sigmax)
        l1a = s1.plot(binsx, y1a, "k--", linewidth=1,label=axis+" Histogram Least Square") #least Square
        y1b = norm.pdf(binsx, mux, Mcal[0])
        l1b = s1.plot(binsx, y1b, "r--", linewidth=2,label=axis+" Twiss RMS") #RMS

        y2a = norm.pdf(binspx, mupx, sigmapx)
        l2a = s2.plot(binspx, y2a, "k--", linewidth=1,label=axis+"' Histogram Least Square")
        y2b = norm.pdf(binspx, mupx, Mcal[1])
        l2b = s2.plot(binspx, y2b, "r--", linewidth=2,label="Filtered "+axis+"' Twiss RMS")

        order = [2,0,1]#for no Eq fits for initial Hist plot
    else:
        #Add the "best fit" line using the findFit mu and sigma for x and display sigma
        y1a = norm.pdf(binsx, mux, sigmax)
        l1a = s1.plot(binsx, y1a, "k--", linewidth=1,label=axis+" Histogram Least Square")
        y1b = norm.pdf(binsx, mux, Mcal[0])
        l1b = s1.plot(binsx, y1b, "r--", linewidth=2,label=axis+" Twiss RMS")
        y1d = norm.pdf(binsx, mux, Me8[0])
        l1d = s1.plot(binsx, y1d, "y--", linewidth=1,label=PDFlabele16+" "+axis+" Twiss")

        y2a = norm.pdf(binspx, mupx, sigmapx)
        l2a = s2.plot(binspx, y2a, "k--", linewidth=1,label=axis+"' Histogram Least Square")
        y2b = norm.pdf(binspx, mupx, Mcal[1])
        l2b = s2.plot(binspx, y2b, "r--", linewidth=2,label=axis+"' Twiss RMS")
        y2c = norm.pdf(binspx, mupx, Me8[1])
        l2c = s2.plot(binspx, y2c, "b.", linewidth=1,label=PDFlabele8+" "+axis+"' Twiss RMS")
        y2d = norm.pdf(binspx, mupx, Me16[1])
        l2d = s2.plot(binspx, y2d, "y--", linewidth=1,label=PDFlabele16+" "+axis+"' Twiss RMS")

        order=[4,0,1,2,3] #order of labels in legend

    titlestart = "At {:.2f}mm ".format(thick) + mat + " " + titledescr
    #Finish setting various plot characteristics
    s1.set_title(titlestart+", Twiss Fit "+axis+" Distribution",fontsize=fs) #+"\n"+rf"$\sigma_D=${{:.1f}}mm".format(sigmax)
    s1.set_xlabel(axis+" Position [mm]",fontsize=fs)
    s1.set_ylabel("Probability Density",fontsize=fs)
    handles1, labels1 = plt.gca().get_legend_handles_labels()
    s1.legend([handles1[idx] for idx in order],[labels1[idx] for idx in order],fontsize=9)
    ylim1=s1.get_ylim() #dynamically get the ylimits
    xlim1=s1.get_xlim()
    if axis=="X":
        sigmatextx = r"$\sigma_{X}$:"
        sigmatextpx = r"$\sigma_{X'}$:"
    elif axis=="Y":
        sigmatextx = r"$\sigma_{Y}$:"
        sigmatextpx = r"$\sigma_{Y'}$:"
    sigmatextx +="\nHistogram = "+"{:.2f}".format(sigmax)+"mm"
    sigmatextx +="\nTwiss RMS = "+"{:.2f}".format(Mcal[0])+"mm"
    sigmatextx +="\n"+PDFlabele8+" RMS = "+"{:.2f}".format(Me8[0])+"mm"
    sigmatextx +="\n"+PDFlabele16+" RMS = "+"{:.2f}".format(Me16[0])+"mm"
    sigmatextpx +="\nHistogram = "+"{:.2f}".format(sigmapx/mm)+"mrad"
    sigmatextpx +="\nTwiss RMS = "+"{:.2f}".format(Mcal[1]/mm)+"mrad"
    sigmatextpx +="\n"+PDFlabele8+" RMS = "+"{:.2f}".format(Me8[1]/mm)+"mrad"
    sigmatextpx +="\n"+PDFlabele16+" RMS = "+"{:.2f}".format(Me16[1]/mm)+"mrad"
    s1.text(xlim1[0]*0.97,ylim1[1]*0.7,sigmatextx,fontsize=fs-2)

    ylim=5/(10**(5+0)) #change if not N=10^5
    #s2.set_ylim([0,ylim])
    s2.set_xlim(-sigmapx*8,sigmapx*8)
    s2.set_title(titlestart+", Twiss Fit  "+axis+"' Distribution",fontsize=fs)
    s2.set_xlabel(axis+"' Angle [rad]",fontsize=fs)
    s2.set_ylabel("Probability Density",fontsize=fs)
    handles2, labels2 = plt.gca().get_legend_handles_labels()
    s2.legend([handles2[idx] for idx in order],[labels2[idx] for idx in order],fontsize=9)
    ylim2=s2.get_ylim() #dynamically get the ylimits
    xlim2=s2.get_xlim()
    s2.text(xlim2[0]*0.97,ylim2[1]*0.7,sigmatextpx,fontsize=fs-2)

    #Can date stamp the multi plot for easier tracking of changes, if necessary
    from datetime import datetime
    dt = datetime.now()

    name = savename+"_TwissFit"+axis+"'"+dt.strftime("%H-%M-%S") +".png"##
    #plt.show()
    plt.tight_layout()
    plt.savefig(name,bbox_inches='tight',dpi=500)
    plt.close()
    print(name)












##Does matrix multiplication to beam Twiss to drift the beam the driftLength
def Drift(TwPm,label,driftLength):
    #Extend the distributions to the Target Location
    import numpy as np
    #Twiss=[beta,alpha,gemt,gamma]
    PBWexitBetaMx = np.array([[TwPm[0],-TwPm[1]],[-TwPm[1],TwPm[3]]])

    drift_PBW_Targ = np.array([ [1, driftLength/1000],[ 0, 1]]) #[m]
    Calc_PBW_Targ = np.linalg.multi_dot([drift_PBW_Targ,PBWexitBetaMx,np.transpose(drift_PBW_Targ)])
    #print(label,"PBWexit:",PBWexitBetaMx,"\n",drift_PBW_Targ,"\n",Calc_PBW_Targ)

    return [Calc_PBW_Targ[0][0],-Calc_PBW_Targ[1][0],TwPm[2],Calc_PBW_Targ[1][1]]












##Compares MiniScatter Target distribution (targxTwissf) to initTarg, exitTarg, e8Targ and e16Targ PDFs
def compareTargets(targx,targy,targTwx,targTwy,fitTwx,fitTwy,fitlabel,savename,mat,PBWTwx,PBWTwy,args):
    import numpy as np
    #from plotFit import findFit, getMoments
    from scipy.stats import norm 
    import matplotlib.pyplot as plt
    mm = 1e-3
    fs=18

    ##Find fit to histogram                                     #ampl,mu,sigma
    print("compareTargs")
    mux,sigmax,amplx,xinterval = findFit(targx,[0.01,0.1,15],(0,[5e8,100,100]),"auto") #to accomodate 1e7 protons
    muy,sigmay,amply,yinterval = findFit(targy,[0.01,0.1,15],(0,[5e8,100,100]),"auto")
    print("findfit: {:.3f}, {:.3f}; np: {:.3f}, {:.3f}".format(sigmax,sigmay,np.std(targx),np.std(targy)))

    ##Find range of particles that are outside 3 sigma
    sigx=np.abs(targx)>3*sigmax# and np.abs(xs)<10*sigma)
    sigy=np.abs(targy)>3*sigmay

    ##Find % of particles outside 3 sigma
    pOut3sigx = len(targx[sigx])/len(targx)*100
    pOut3sigy = len(targy[sigy])/len(targy)*100
    print("Distribution at Target: {:0.3f}% outside 3 sigma X".format(pOut3sigx))
    print("Distribution at Target: {:0.3f}% outside 3 sigma Y".format(pOut3sigy))

    ##Get sigma for position and angle
    Mfx = getMoments(fitTwx)
    Mfy = getMoments(fitTwy)
    Mhx = getMoments(targTwx)
    Mhy = getMoments(targTwy)
    #print(fitlabel,Mhx,Mfx)

    ##Sigmas of Histogram vs fit
    print("\nDistribution at Target \t sigma X: {:.2f}mm".format(sigmax))
    print("PDF at Target \t\t sigma X: {:.2f}mm".format(Mhx[0]))
    print(fitlabel,"\t\t sigma X: {:.2f}mm".format(Mfx[0]))

    print("\nDistribution at Target \t sigma Y: {:.2f}mm".format(sigmay))
    print("PDF at Target \t\t sigma Y: {:.2f}mm".format(Mhy[0]))
    print(fitlabel,"\t\t sigma Y: {:.2f}mm".format(Mfy[0]))

    plt.clf()
    xPlot = False
    if xPlot:
        ##Create the fig with 2 plots side by side
        fig = plt.figure(figsize=(16,6.0))
        plt.subplots_adjust(wspace=0.2) #increase width space to not overlap
        s1 = fig.add_subplot(1,2,1)
        s2 = fig.add_subplot(1,2,2)
    else:
        fig = plt.figure(figsize = (8,5))
        s2 = fig.add_subplot(1,1,1)

    if xPlot:
        ##Make the histogram of the full energy distrubtion for X with findFit intervals
        nx, binsx, patchesx = s1.hist(targx, xinterval, density=True, log=True, facecolor="green", alpha=0.75,label="MiniScatter Target X Histogram")

        ##Add the "best fit" line using the findFit mu and sigma for x and display sigma
        y1a = norm.pdf(binsx, mux, sigmax)
        l1a = s1.plot(binsx, y1a, "k--", linewidth=1,label="Histogram Least Square")
        #y1b = norm.pdf(binsx, mux, Mhx[0])
        #l1b = s1.plot(binsx, y1b, "r--", linewidth=2,label="Target Twiss RMS") #sample RMS
        y1c = norm.pdf(binsx, mux, Mfx[0])
        l1c = s1.plot(binsx, y1c, "b--", linewidth=2,label=fitlabel+" Twiss RMS")

        ##Set various plot variables, often found through trial and error
        ##if don't set ylim, log goes to e-58
        plt.setp(s1,xlim=([-args.xlim,args.xlim]),ylim=([1e-6,3]))
        s1.set_title("X Distribution At Target",fontsize=fs+2) #+"\n"+rf"$\sigma_D=${{:.1f}}mm".format(sigmax)
        s1.set_xlabel("X Position [mm]",fontsize=fs)
        s1.set_ylabel("Probability Density",fontsize=fs)
        #s1.set_xlim([-xlim,xlim])
        ##s1.set_ylim([1e-6,0.1])
        #s1.set_ylim([1e-6,1])  #if don't set, then log goes to e-58

        ##Write texts to display sigmas
        import re
        if re.search(",",fitlabel): #Not sure what this is for...
            label = re.sub(".+(?<=(,))","",fitlabel)
            label = label.replace(" ","")
        else:
            label = "Müller Eqn 8"
        sigmatextx = r"$\sigma_{X}$:"
        sigmatextx +="\nHistogram = "+"{:.2f}".format(sigmax)+"mm"
        #sigmatextx +="\nTwiss RMS = "+"{:.2f}".format(Mhx[0])+"mm"
        sigmatextx +="\n"+label+" RMS = "+"{:.2f}".format(Mfx[0])+"mm"
        #sigmatextx +="\n\n{:.3f}% outside".format(pOut3sigx)+r" 3$\sigma$"
        sigmatextx += "\nHistogram/Müller={:.2f}".format((sigmax/Mfx[0]-1)*100)+"%"

        xlim1 = s1.get_xlim()
        ylim1 = s1.get_ylim()
        s1.text(xlim1[0]*0.97,5e-2,sigmatextx,fontsize=fs-5)
        s1.text(xlim1[0]*0.97,1.2e-2,"Beam Twiss at PBW:",fontsize=fs-5)
        s1.text(xlim1[0]*0.97,6.0e-3,r"$\epsilon_{Nx}$ = "+"{:.3f}".format(PBWTwx[2])+r"$_{[mm \cdot mrad]}$",fontsize=fs-5)
        s1.text(xlim1[0]*0.97,2.5e-3,r"$\beta_{x}$ = "+"{:.1f}".format(PBWTwx[0])+r"$_{[m]}$",fontsize=fs-5)
        s1.text(xlim1[0]*0.97,1.0e-3,r"$\alpha_{x}$ = "+"{:.2f}".format(PBWTwx[1]),fontsize=fs-5)

        handles1, labels1 = plt.gca().get_legend_handles_labels()
        order1=[0,1,2]
        s1.legend([handles1[idx] for idx in order1],[labels1[idx] for idx in order1],fontsize=fs-6,loc="upper right")

    ##Make the histogram of the full energy distrubtion for Y with findFit intervals
    ny, binsy, patchesy = s2.hist(targy, yinterval, density=True, log=True, facecolor="green", alpha=0.75,label="MiniScatter Target Y Histogram")

    ##Add the "best fit" line using the findFit mu and sigma for Y and display sigma
    y2a = norm.pdf(binsy, muy, sigmay)
    l2a = s2.plot(binsy, y2a, "k--", linewidth=1,label="Histogram Least Square")
    #y2b = norm.pdf(binsy, muy, Mhy[0])
    #l2b = s2.plot(binsy, y2b, "r--", linewidth=2,label="Target Twiss RMS")
    y2c = norm.pdf(binsy, muy, Mfy[0])
    l2c = s2.plot(binsy, y2c, "b--", linewidth=2,label=fitlabel+" Twiss RMS")

    ##Write texts to display sigmas
    import re
    if re.search(",",fitlabel): #Not sure what this is for...
        label = re.sub(".+(?<=(,))","",fitlabel)
        label = label.replace(" ","")
    else:
        label = "Müller Eqn 8"

    sigmatexty = r"$\sigma_{Y}$:"
    sigmatexty +="\nHistogram = "+"{:.2f}".format(sigmay)+"mm"
    #sigmatexty +="\nTwiss RMS = "+"{:.2f}".format(Mhy[0])+"mm"
    sigmatexty +="\n"+label+" RMS = "+"{:.2f}".format(Mfy[0])+"mm"
    #sigmatexty +="\n\n{:.3f}% outside".format(pOut3sigy)+r" 3$\sigma$"
    sigmatexty += "\nHistogram/Müller={:.2f}".format((sigmay/Mfy[0]-1)*100)+"%"

    ##Set various plot variables, often found through trial and error
    ##Set s2
    plt.setp(s2,xlim=([-args.xlim,args.xlim]),ylim=([1e-6,3]))
    s2.set_title("Y Distribution At Target",fontsize=fs+2)
    #if don't set ylim, log goes to e-58
    #s2.set_xlim([-xlim,xlim])
    #s2.set_ylim([1e-6,1]) #if don't set, then log goes to e-58
    s2.set_xlabel("Y Position [mm]",fontsize=fs)
    s2.set_ylabel("Probability Density",fontsize=fs)

    handles2, labels2 = plt.gca().get_legend_handles_labels()
    order2=[0,1,2]
    s2.legend([handles2[idx] for idx in order2],[labels2[idx] for idx in order2],fontsize=fs-6,loc="upper right")
    ylim2=s2.get_ylim() #dynamically get the graph limits
    xlim2=s2.get_xlim()

    s2.text(xlim2[0]*0.97,5e-2,sigmatexty,fontsize=fs-4)
    #PBWTwx = [Ibetax,Ialphx,Inemtx]
    s2.text(xlim2[0]*0.97,1.2e-2,"Beam Twiss at PBW:",fontsize=fs-5)
    s2.text(xlim2[0]*0.97,6.0e-3,r"$\epsilon_{Ny}$ = "+"{:.3f}".format(PBWTwy[2])+r"$_{[mm \cdot mrad]}$",fontsize=fs-5)
    s2.text(xlim2[0]*0.97,2.5e-3,r"$\beta_{y}$ = "+"{:.1f}".format(PBWTwy[0])+r"$_{[m]}$",fontsize=fs-5)
    s2.text(xlim2[0]*0.97,1.0e-3,r"$\alpha_{y}$ = "+"{:.2f}".format(PBWTwy[1]),fontsize=fs-5)

    ##Can date stamp the multi plot for easier tracking of changes, if necessary
    from datetime import datetime
    dt = datetime.now()
    print(dt.strftime("%H-%M-%S"))

    name = savename+"_compTargs_"+"."+args.picFormat##+dt.strftime("%H-%M-%S")
    #plt.show()
    if args.savePics:
        #plt.tight_layout()
        plt.savefig(name,bbox_inches='tight',dpi=args.dpi)
    plt.close()
    print(name)
    return sigmax,sigmay












##Prints particles in form for importParticles option of MiniScatter
def printParticles(savename,xinit,pxinit,yinit,pyinit,Einit):
    import csv
    from datetime import datetime
    dt = datetime.now()
    z = -1

    fname = savename+'.csv'
    with open(fname,mode = 'w') as part_file:
        part_writer = csv.writer(part_file,delimiter = ',')
        for i in range(len(Einit)):
            part_writer.writerow(["proton",xinit[i],pxinit[i],yinit[i],pyinit[i],z,Einit[i]])
    part_file.close()
    print(fname)













##Plots the projection of particle positions for a rastered beam (obsolete)
def plot1DRaster(targx,targy,fitlabel,savename,mat,position):
    # To Do:
    # --Add sobel filter for finding edge of Raster
    # --Add XY view with 99% box, beam edges, current density calculation, 
    import numpy as np
    #from plotFit import findFit, getMoments
    #from scipy.stats import norm 
    import matplotlib.pyplot as plt
    mm = 1e-3
    fs = 14
    #99% box widths from Beam Requirements
    xline = 160
    yline = 64

    #Find % outside Open Space Box
    boxX = np.abs(targx) > xline
    boxY = np.abs(targy) > yline
    pOutBoxX = len(targx[boxX])/len(targx)*100
    pOutBoxY = len(targy[boxY])/len(targy)*100
    print(fitlabel,pOutBoxX,"% outside 99% box X")
    print(fitlabel,pOutBoxY,"% outside 99% box Y")

    #Useful to get intervals for nice plotting
    mux,sigmax,amplx,xinterval = findFit(targx,[0.01,1,1],(0,[1e3,3e5,3e5]),"auto")
    muy,sigmay,amply,yinterval = findFit(targy,[0.01,1,1],(0,[1e3,3e5,3e5]),"auto")

    #Create the fig with 2 plots side by side
    plt.clf()
    fig = plt.figure(figsize=(15,6.0))
    plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
    s1 = fig.add_subplot(1,2,1)
    s2 = fig.add_subplot(1,2,2)

    #Make the histogram of the full energy distrubtion for X with findFit intervals
    nx, binsx, patchesx = s1.hist(targx, xinterval, density=True, log=True, facecolor="green", alpha=0.75,label="MiniScatter Target X Histogram")

    #Make the histogram of the full energy distrubtion for Y with findFit intervals
    ny, binsy, patchesy = s2.hist(targy, yinterval, density=True, log=True, facecolor="green", alpha=0.75,label="MiniScatter Target Y Histogram")

    #General range:
    xlim1 = 300
    ylim1 = [1e-6,1e-1]
    xlim2 = 300
    ylim2 = [1e-6,1e-1]
    #to zoom in on peaks:
    #xlim1 = 60
    #ylim1 = [5e-3,1.3e-2]
    #xlim2 = 20
    #ylim2 = [2.1e-2,2.7e-2]

    plt.setp(s1,title="X Raster Distribution At "+position,xlabel="X Position [mm]",
            ylabel="Probability Density",xlim=([-xlim1,xlim1]),ylim=(ylim1))
    #s1.set_title("X Raster Distribution At "+position,fontsize=fs) #+"\n"+rf"$\sigma_D=${{:.1f}}mm".format(sigmax)
    #s1.set_xlabel("X Position [mm]",fontsize=fs)
    #s1.set_ylabel("Probability Density",fontsize=fs)
    #s1.set_xlim([-xlim1,xlim1])
    #s1.set_ylim(ylim1) #if don't set, then log goes to e-58

    plt.setp(s2,title="Y Raster Distribution At "+position,xlabel="X Position [mm]",
            ylabel="Probability Density",xlim=([-xlim2,xlim2]),ylim=(ylim2))
    #s2.set_title("Y Raster Distribution At "+position,fontsize=fs)
    #s2.set_xlabel("Y Position [mm]",fontsize=fs)
    #s2.set_ylabel("Probability Density",fontsize=fs)
    #s2.set_xlim([-xlim2,xlim2])
    #s2.set_ylim(ylim2) #if don't set, then log goes to e-100

    #99% Box label. Comment out placement if zooming in!
    textX = "{:.3f}% outside 99% Box".format(pOutBoxX)
    textY = "{:.3f}% outside 99% Box".format(pOutBoxY)
    s1.text(-xlim1*0.97,3e-2,textX,fontsize=fs-2)
    s2.text(-xlim2*0.97,3e-2,textY,fontsize=fs-2)

    #Makes 99% box lines
    s1.vlines(-xline,ylim1[0],ylim1[1]/10)
    s1.vlines(xline,ylim1[0],ylim1[1]/10)
    s2.vlines(-yline,ylim2[0],ylim2[1]/10)
    s2.vlines(yline,ylim2[0],ylim2[1]/10)
    
    #To focus on small area and have fine grid
    ##https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large/44079725#44079725
    #import matplotlib.ticker as tick
    #locmin = tick.LogLocator(base=10.0, subs=(1e-3,5e-2,1e-2,1e-1 ))
    #s2.yaxis.set_minor_locator(locmin)
    #s1.grid(visible=True,which='both')
    #s2.grid(visible=True, which='both')

    #Can date stamp the multi plot for easier tracking of changes, if necessary
    from datetime import datetime
    dt = datetime.now()

    name = savename+"_"+fitlabel+"_"+dt.strftime("%H-%M-%S")##
    #plt.show()
    plt.tight_layout
    plt.savefig(name+".pdf",bbox_inches='tight',dpi=500)
    plt.close()
    print(name)












#Counts number of lines in a CSV
def numLines(filename):
    import csv
    if filename[-1] != "v":
        pass
    else:
        print(filename)
        filename-=".csv"
        print(filename)
    with open(filename+".csv") as csv_file:
        return sum(1 for line in csv_file) 
        #extreme improvement from https://stackoverflow.com/questions/16108526/how-to-obtain-the-total-numbers-of-rows-from-a-csv-file-in-python
        #csv_reader = csv.reader(csv_file, delimiter=',')
        #N = len(list(csv_reader))
        #N = 0
        #for row in csv_reader:
        #    N += 1
        #csv_file.close()
    #from datetime import datetime
    #dt = datetime.now()
    #print("There are ",N,"lines "+dt.strftime("%H-%M-%S"))
    #return N






##WIP:Attempt at a 2D fit to a scattered beamlet
def gaussian2DFit(hist,axis,width,maxim,name,y1,y2,saveFits,density):
    #from Kyrre's doubleGaussian.ipynb example
    import ROOT
    #def plot(proj,ext):
    #    ca = ROOT.TCanvas()
    #    proj.Draw()
    #    #proj.GetXaxis().SetRangeUser(-maxim,maxim)
    #    #ca.SetLogy()
    #    picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
    #    ca.Print(picpath+"MCNPXDataROOT2D_"+ext+"2D.png")
    #    print(picpath+"MCNPXDataROOT2D_"+ext+"2D.png")

    ROOT.gStyle.SetOptStat("iouRMe")
    ROOT.gStyle.SetOptFit(100)
    ROOT.gStyle.SetStatW(0.15)
    ROOT.gStyle.SetStatH(0.1)
    #plot(hist,"test")
    #Project center slice 
    integral = hist.Integral(hist.GetXaxis().FindBin(-2),hist.GetXaxis().FindBin(2),hist.GetYaxis().FindBin(-2),hist.GetYaxis().FindBin(2))
    hist.Scale(1/integral)
    sum = hist.Integral()
    #print(integral,sum)
    #print(f"axis={axis}","binY",hist.GetYaxis().FindBin(-width), hist.GetYaxis().FindBin(width),"binX",hist.GetXaxis().FindBin(-width), hist.GetXaxis().FindBin(width))

    #Define a 2D gaussian function and fit it to the projection
    f1 = ROOT.TF2('f1','bigaus',-maxim,maxim,-maxim,maxim)
    #TF2 *f2 = new TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",0,10,0,10); f2-<SetParameters(1,3,1,4,1); nGauss->Fit(f2);
    f1.SetNpx(500)
    f1_res = hist.Fit(f1, 'RSQ')
    #ca = ROOT.TCanvas("Scattered Beam Fit","Scattered Beam Fit",0,0,300*8,300*8)
    #hist.Draw("E0")
    #f1.SetLineColor(ROOT.kRed)
    #f1_res.Draw('same')
    #ca.SetLogz()
    #ca.Print(name+"_init2D.png")
    print(f1_res)
    p0 = f1.GetParameter(0) #weight x 
    p1 = f1.GetParameter(1) #mu x 
    p2 = f1.GetParameter(2) #sigma x 
    p3 = f1.GetParameter(3) #weight x 
    p4 = f1.GetParameter(4) #mu x 
    p5 = f1.GetParameter(5) #sigma x 
    print("initial 2D Gaussian fit parameters:",p0,p1,p2,p3,p4,p5)



##ROOT fitting of scattered beamlet. Fits to multiple Gaussians, nGauss, or with one Lorentz
def gaussianFit(hist,axis,width,maxim,name,y1,y2,saveFits,density,nGauss,Lorentz):
    #from Kyrre's doubleGaussian.ipynb example
    import ROOT
    #def plot(proj,ext):
    #    ca = ROOT.TCanvas()
    #    proj.Draw()
        #proj.GetXaxis().SetRangeUser(-maxim,maxim)
        #ca.SetLogy()
    #    picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
    #    ca.Print(picpath+"MCNPXDataROOT2D_"+ext+".png")
    #    print(picpath+"MCNPXDataROOT2D_"+ext+".png")

    ROOT.gStyle.SetOptStat("iRMe")
    ROOT.gStyle.SetOptFit(101)
    #plot(hist,"test")
    #Project center slice 
    integral = hist.Integral(hist.GetXaxis().FindBin(-100),hist.GetXaxis().FindBin(100),hist.GetYaxis().FindBin(-100),
                             hist.GetYaxis().FindBin(100))
    hist.Scale(1/integral)
    sum = hist.Integral()
    #print(integral,sum)
    
    if   axis == "y" or axis == "Y": 
        proj = hist.ProjectionY(axis,hist.GetXaxis().FindBin(-width),hist.GetXaxis().FindBin(width))
    elif axis == "x" or axis == "X":
        proj = hist.ProjectionX(axis,hist.GetYaxis().FindBin(-width),hist.GetYaxis().FindBin(width)) #why?
    #plot(proj,"proj"+axis)
    print(f"axis={axis}","binY",hist.GetYaxis().FindBin(-width), hist.GetYaxis().FindBin(width),"binX",hist.GetXaxis().FindBin(-width),
          hist.GetXaxis().FindBin(width))
    differenceNG = 100
    total = 100
    differencePG = 100
    coeffsG = [100]

    #Define a gaussian function and fit it to the projection
    f1 = ROOT.TF1('f1','gaus',-maxim,maxim)
    f1.SetNpx(5000)
    f1_res = proj.Fit(f1, 'RSQ')
    #ca = ROOT.TCanvas("Scattered Beam Fit","Scattered Beam Fit",0,0,300*8,300*8)
    #proj.Draw("E0")
    #f1.SetLineColor(ROOT.kRed)
    #f1_res.Draw('same')
    #ca.SetLogy()
    #ca.Print(name+"_init.png")
    #print(f1_res)
    p0 = f1.GetParameter(0) #weight
    p1 = f1.GetParameter(1) #mu
    p2 = f1.GetParameter(2) #sigma
    #print("initial Gaussian fit parameters:",p0,p1,p2)

    r=0.01
    ##Set function
    if nGauss == 2:
        f2 = ROOT.TF1('f2','([0] / sqrt(2*pi*[1]) ) * exp(-x*x/(2*[1]*[1])) + ([2] / sqrt(2*pi*[3]) ) * exp(-x*x/(2*[3]*[3]))',-maxim,maxim)
    elif nGauss == 3:
        f2 = ROOT.TF1('f2','([0] / sqrt(2*pi*[1]) ) * exp(-x*x/(2*[1]*[1])) + ([2] / sqrt(2*pi*[3]) ) * exp(-x*x/(2*[3]*[3])) + ([4] / sqrt(2*pi*[5]) ) * exp(-x*x/(2*[5]*[5]))',-maxim,maxim)
    elif nGauss == 4:
        f2 = ROOT.TF1('f2','([0] / sqrt(2*pi*[1]) ) * exp(-x*x/(2*[1]*[1])) + ([2] / sqrt(2*pi*[3]) ) * exp(-x*x/(2*[3]*[3])) + ([4] / sqrt(2*pi*[5]) ) * exp(-x*x/(2*[5]*[5])) + ([6] / sqrt(2*pi*[7]) ) * exp(-x*x/(2*[7]*[7]))',-maxim,maxim)
    elif nGauss == 5:
        f2 = ROOT.TF1('f2','([0] / sqrt(2*pi*[1]) ) * exp(-x*x/(2*[1]*[1])) + ([2] / sqrt(2*pi*[3]) ) * exp(-x*x/(2*[3]*[3])) + ([4] / sqrt(2*pi*[5]) ) * exp(-x*x/(2*[5]*[5])) + ([6] / sqrt(2*pi*[7]) ) * exp(-x*x/(2*[7]*[7])) + ([8] / sqrt(2*pi*[9]) ) * exp(-x*x/(2*[9]*[9]))',-maxim,maxim)

    ##Set parameters. Setting more than needed for a function does nothing
    f2.SetParameter(0,p0)                 #A1
    f2.SetParLimits(0,1e-4,p0*100)        #A1
    f2.SetParameter(1,p2)                 #sigma1
    f2.SetParLimits(1,p2,p2*(1+r))        #sigma1

    f2.SetParameter(2,p0*(r**2))          #A2
    f2.SetParLimits(2,0, p0*5)            #A2
    f2.SetParameter(3,p2*2)               #sigma2
    f2.SetParLimits(3,p2*(1+2*r), p2*2.6) #sigma2

    f2.SetParameter(4,p0*r**2)            #A3
    f2.SetParLimits(4,0, p0*5)            #A3
    f2.SetParameter(5,p2*7)               #sigma3
    f2.SetParLimits(5,p2*2.5, p2*12)      #sigma3

    f2.SetParameter(6,p0*r**4)            #4
    f2.SetParLimits(6,0, p0)              #A4
    f2.SetParameter(7,p2*10)              #sigma4
    f2.SetParLimits(7,p2*8, p2*23)        #sigma4

    f2.SetParameter(8,p0*r**5)            #A5
    f2.SetParLimits(8,0, p0)              #A5
    f2.SetParameter(9,p2*25)              #sigma5
    f2.SetParLimits(9,p2*15, p2*30)       #sigma5

    ##Number of fit points and Fit. Return function and fit results
    f2.SetNpx(10000)
    f2_res = proj.Fit(f2,'RSQ')
    if nGauss==1:
        print("1 Gaussian: {:.7f} {:.7f}".format(f2.GetParameter(0),f2.GetParameter(1)))
    elif nGauss==2:
        print("2 Gaussians: {:.7f} {:.7f} {:.7f} {:.7f}".format(f2.GetParameter(0),f2.GetParameter(1),
                                                             f2.GetParameter(2),f2.GetParameter(3)))
    elif nGauss==3:
        print("3 Gaussians: {:.7f} {:.7f} {:.7f} {:.7f} {:.7f} {:.7f}".format(f2.GetParameter(0),f2.GetParameter(1),
                                                                           f2.GetParameter(2),f2.GetParameter(3),
                                                                           f2.GetParameter(4),f2.GetParameter(5)))
    elif nGauss==4:
        print("4 Gaussians: {:.7f} {:.7f} {:.7f} {:.7f} {:.7e} {:.7f} {:.7e} {:.7f}".format(f2.GetParameter(0),f2.GetParameter(1),
                                                                                         f2.GetParameter(2),f2.GetParameter(3),
                                                                                         f2.GetParameter(4),f2.GetParameter(5),
                                                                                         f2.GetParameter(6),f2.GetParameter(7)))
    elif nGauss==5:
        print("5 Gaussians: {:.7f} {:.7f} {:.7f} {:.7f} {:.7e} {:.7f} {:.7e} {:.7f} {:.7e} {:.7f}".format(f2.GetParameter(0),f2.GetParameter(1),
                                                                                                       f2.GetParameter(2),f2.GetParameter(3),
                                                                                                       f2.GetParameter(4),f2.GetParameter(5),
                                                                                                       f2.GetParameter(6),f2.GetParameter(7),
                                                                                                       f2.GetParameter(8),f2.GetParameter(9)))

    #print(f2_res)

    if saveFits: #if uncommented it opens a canvas even if false...?

        c1 = ROOT.TCanvas("Scattered Beam Fit","Scattered Beam Fit",0,0,500*8,300*8)
        proj.Draw("E0")
        proj.SetMarkerStyle(34)
        proj.SetMarkerColor(1)
        proj.SetMarkerSize(3)
        f2.SetLineColor(ROOT.kRed)
        f2_res.Draw('same')
        #f1.SetLineColor(ROOT.kRed)
        #f1.Draw('same')
        proj.GetXaxis().SetRangeUser(-maxim,maxim)
        proj.GetYaxis().SetRangeUser(2e-10,2e-1)
        c1.SetLogy()

        
        if nGauss==1:
            leg = ROOT.TLegend(0.125,0.75,0.3,0.85) #https://root.cern.ch/root/htmldoc/guides/users-guide/Graphics.html
            leg.AddEntry(proj,"Particle Density")
            leg.AddEntry(f2,"A_{1} e^{#frac{-x^{2}}{2#sigma_{1}^{2}}}")
        elif nGauss==2:
            leg = ROOT.TLegend(0.125,0.672,0.32,0.85) #https://root.cern.ch/root/htmldoc/guides/users-guide/Graphics.html
            leg.AddEntry(proj,"Particle Density")
            leg.AddEntry(f2,"A_{1} e^{#frac{-x^{2}}{2#sigma_{1}^{2}}} + A_{2} e^{#frac{-x^{2}}{2#sigma_{2}^{2}}}")##frac{#lambda}{x^{2} + #gamma^{2}}")
        elif nGauss==3:
            leg = ROOT.TLegend(0.1,0.7,0.37,0.85)
            leg.AddEntry(proj,"Particle Density")
            leg.AddEntry(f2,"A_{1} e^{#frac{-x^{2}}{2#sigma_{1}^{2}}} + A_{2} e^{#frac{-x^{2}}{2#sigma_{2}^{2}}} + A_{3} e^{#frac{-x^{2}}{2#sigma_{3}^{2}}}")#frac{#lambda}{x^{2} + #gamma^{2}}")
            #leg.AddEntry(f2,"A_{1} e^{#frac{-x^{2}}{#sigma_{1}^{2}}} + A_{2} e^{#frac{-x^{2}}{#sigma_{2}^{2}}} + #frac{#lambda}{x^{3}} + #gamma")
        elif nGauss==4:
            leg = ROOT.TLegend(0.1,0.675,0.425,0.85)
            leg.AddEntry(proj,"Particle Density")
            fourGaussL = "A_{1} e^{#frac{-x^{2}}{2#sigma_{1}^{2}}} + A_{2} e^{#frac{-x^{2}}{2#sigma_{2}^{2}}} + "
            fourGaussL += "A_{3} e^{#frac{-x^{2}}{2#sigma_{3}^{2}}} + A_{4} e^{#frac{-x^{2}}{2#sigma_{4}^{2}}}"##frac{#lambda}{x^{2} + #gamma^{2}}"
            leg.AddEntry(f2,fourGaussL)
            #leg.AddEntry(f2,"A_{1} e^{#frac{-x^{2}}{#sigma_{1}^{2}}} + A_{2} e^{#frac{-x^{2}}{#sigma_{2}^{2}}} + 
            # A_{3} e^{#frac{-x^{2}}{#sigma_{3}^{2}}} + #frac{#lambda}{x^{3}} + #gamma")elif nGauss==4:
        elif nGauss==5:
            leg = ROOT.TLegend(0.1,0.675,0.45,0.85)
            leg.AddEntry(proj,"Particle Density")
            fiveGaussL = "A_{1} e^{#frac{-x^{2}}{2#sigma_{1}^{2}}} + A_{2} e^{#frac{-x^{2}}{2#sigma_{2}^{2}}} + "
            fiveGaussL += "A_{3} e^{#frac{-x^{2}}{2#sigma_{3}^{2}}} + A_{4} e^{#frac{-x^{2}}{2#sigma_{4}^{2}}} + "
            fiveGaussL += "A_{5} e^{#frac{-x^{2}}{2#sigma_{5}^{2}}}"#"#frac{#lambda}{x^{2} + #gamma^{2}}"
            leg.AddEntry(f2,fiveGaussL)
            #leg.AddEntry(f2,"A_{1} e^{#frac{-x^{2}}{#sigma_{1}^{2}}} + A_{2} e^{#frac{-x^{2}}{#sigma_{2}^{2}}} + 
            # A_{3} e^{#frac{-x^{2}}{#sigma_{3}^{2}}} + #frac{#lambda}{x^{3}} + #gamma")
        #"A e^{#frac{-x^{2}}{#sigma^{2}}} + #frac{B}{x^{2} + #gamma^{2}}") #"A e^{#frac{-x^{2}}{#sigma^{2}}} + #frac{B}{x^{3}} + C"
        #leg.SetHeader("Fits")
        leg.SetMargin(0.1)
        leg.Draw("E0")

        #ta = ROOT.char(6)
        
        #t1 = ROOT.TText(0.1,0.65,f2.GetParameter("w1"))
        #t2 = ROOT.TText(0.1,0.60,ROOT.Form("#sigma = %.3f", f2.GetParameter("s1")))
        #t3 = ROOT.TText(0.1,0.55,ROOT.Form("B = %.3f",f2.GetParameter("w2")))
        #t4 = ROOT.TText(0.1,0.50,ROOT.Form("#gamma = %.3f", f2.GetParameter("s2")))
        #t1.Draw()
        #t2.Draw()
        #t3.Draw()
        #t4.Draw()
        c1.Print(name+"_"+axis+str(nGauss)+"GaussFit"+str(maxim)+".png")

    #    import numpy as np
    #    import matplotlib.pyplot as plt
    #    plt.close()
    #    w = np.zeros(maxim)
    #    diff = np.zeros(maxim)
    #    for i in range(maxim):
    #        w[i] = i
    #        diff[i] = proj.Integral(proj.GetXaxis().FindBin(-i),proj.GetXaxis().FindBin(i),'width')  -  f2.Integral(-i,i)
    #        if i == maxim-1: print(diff[i])
    #    plt.plot(diff)
    #    plt.legend()
    #    plt.title("Y Axis Difference")
    #    plt.yscale("log")
    #    plt.savefig(name+"FitDifference.png")
    #    plt.close()

    #total=f2_res.GetNDF()
    Gchi2 = f2_res.Chi2()
    differenceNG = 0#f2_res.GetProb()
    #print(ndf,chisq,prob)

    #differenceNG = proj.Integral(proj.GetXaxis().FindBin(-maxim),proj.GetXaxis().FindBin(maxim),'width')  -  f2.Integral(-maxim,maxim)
    total = proj.Integral(proj.GetXaxis().FindBin(-maxim),proj.GetXaxis().FindBin(maxim),'width') #need better name
    #differencePG = differenceNG/total*100
    coeffsG = [f2.GetParameter(0),f2.GetParameter(1),f2.GetParameter(2),f2.GetParameter(3),f2.GetParameter(4),f2.GetParameter(5),
               f2.GetParameter(6),f2.GetParameter(7),f2.GetParameter(8),f2.GetParameter(9)]
    # = [f2.GetParameter("#lambda"),f2.GetParameter("#gamma"),f2.GetParameter("A_{1}"),f2.GetParameter("#sigma_{1}"),f2.GetParameter("A_{2}"),f2.GetParameter("#sigma_{2}"),
    #            f2.GetParameter("A_{3}"),f2.GetParameter("#sigma_{3}")]#,f2.GetParameter("A_{4}"),f2.GetParameter("#sigma_{4}")]
    #print(axis," difference: {:.0f},\t total: {:.0f},\tchi^2: {:.3f} /".format(differenceNG,total,Gchi2),f2_res.Ndf())#,"\n")#,coeffs)
    
    """
    #Now do for Lorentzian
    #Define a gaussian function and fit it to the projection
    f1 = ROOT.TF1('f1','gaus',-maxim,maxim)
    f1_res = proj.Fit(f1, 'RSQ')
    #print(f1_res)
    p0 = f1.GetParameter(0)
    p1 = f1.GetParameter(1)
    p2 = f1.GetParameter(2)
    print("initial Gaussian fit parameters:",p0,p1,p2)

    #https://root.cern/doc/master/fitConvolution_8py.html
    f3 = ROOT.TF1('f3L','[0] * [1]/(x**2 + [1]**2) + [2] * exp(-x*x/(2*[3]*[3])) ',-maxim,maxim)#+ [4] * exp(-x*x/(2*[5]*[5]))
    #f3L = ROOT.TF1('f3L','[0] * [1]/(x**2 + [1]**2)',-maxim,maxim)
    #f3G = ROOT.TF1('f3G','[2] * exp(-x*x/(2*[3]*[3]))',-maxim,maxim)
    #f3C = ROOT.TF1Convolution('f3L','f3G',-maxim,maxim,True)
    #f3C.SetRange(-maxim,maxim)
    #f3C.SetNofPointsFFT(10000)
    #f3 = ROOT.TF1('f3',f3C,-maxim,maxim,f3C.GetNpar()) #not ever finding a fit!


    #if axis in {"y","Y"}:
    f3.SetParameters(p2,p2,p2,p2)
    #elif axis in {"x","X"}:
    #    f3.SetParameters(p2,p2,p2,p2)
    f3.SetParLimits(0,0,p0*y2*y2)
    f3.SetParLimits(1,0,p0*y2*y2)
    f3.SetParLimits(2,0,p0*y2*y2)
    f3.SetParLimits(3,0,p0*y2*y2)
    f3.SetParLimits(4,0,p0*y2*y2)
    f3.SetParLimits(5,0,p0*y2*y2)
    f3_res = proj.Fit(f3, 'RSQ')
    Vchi2 = f3_res.Chi2()
    """
    differenceNL = 0#proj.Integral(proj.GetXaxis().FindBin(-maxim),proj.GetXaxis().FindBin(maxim),'width')  -  f3.Integral(-maxim,maxim)
    total = 1#proj.Integral(proj.GetXaxis().FindBin(-maxim),proj.GetXaxis().FindBin(maxim),'width') #need better name
    differencePL = differenceNL/total*100
    coeffsL = [0,0,0]#[f3.GetParameter(0),f3.GetParameter(1),f3.GetParameter(2),f3.GetParameter(3)]
    #print("Voigt",coeffsL)
    #print(axis," L difference: {:.0f},\t total: {:.0f},\tchi^2: {:.3f}".format(differenceNL,total,Vchi2))#,"\n")#,coeffs)
    """
    if saveFits: #if uncommented it opens a canvas even if false...?
        c2 = ROOT.TCanvas()
        proj.Draw('E0')
        f3.SetLineColor(ROOT.kRed)
        f3.Draw('same')
        #print("hope this works")
        #f1.SetLineColor(ROOT.kRed)
        #f1.Draw('same')
        #if axis == "y" or axis == "Y": proj.GetXaxis().SetRangeUser(-maxim,maxim)
        #elif axis =="x" or axis == "X": proj.GetXaxis().SetRangeUser(-maxim,maxim)
        proj.GetXaxis().SetRangeUser(-maxim,maxim)
        c2.SetLogy()
        c2.Print(name+"_"+axis+"LorentzFit"+str(maxim)+".pdf")
    """
    #    import numpy as np
    #    import matplotlib.pyplot as plt
    #    plt.close()
    #    w = np.zeros(maxim)
    #    diff = np.zeros(maxim)
    #    for i in range(maxim):
    #        w[i] = i
    #        diff[i] = proj.Integral(proj.GetXaxis().FindBin(-i),proj.GetXaxis().FindBin(i),'width')  -  f2.Integral(-i,i)
    #        if i == maxim-1: print(diff[i])
    #    plt.plot(diff)
    #    plt.legend()
    #    plt.title("Y Axis Difference")
    #    plt.yscale("log")
    #    plt.savefig(name+"FitDifference.png")
    #    plt.close()

    return differenceNG, differencePG, coeffsG, differenceNL,differencePL,coeffsL







def rasterImage(savename,position,histogram2D,parts,args,Twiss,options,boxes,paths,sample):
    import numpy as np
    name=savename+"_"+position+"Image_rB"+str(args.reBin)
    #um = 1e-6
    from datetime import datetime
    #from plotFit import converter

    #print(datetime.now().strftime("%H-%M-%S"))
    (Img, xax, yax) = converter(histogram2D,args.saveHist,name,paths,args,False) #convert from TH2D to numpy map
    #print(datetime.now().strftime("%H-%M-%S"))
    xBinSize = xax[round(len(xax)*0.5)+1]-xax[round(len(xax)*0.5)]
    yBinSize = yax[round(len(xax)*0.5)+1]-yax[round(len(xax)*0.5)]
    sumTot=np.sum(Img)
    ImgJ = Img #PMAS converts whole array into J, so send this one instead of making it inside PMAS.
    #print(xBinSize,yBinSize,sumTot)
    I_pulse = 62.5*1e3 #[uA]

    #99% box #for multiple boxes, use arrays
    #From Dmitriy Work on https://stackoverflow.com/questions/432112/is-there-a-numpy-function-to-return-the-first-index-of-something-in-an-array
    #  np.ndenumerate scales the best of the methods compared to 10^4 elements.
    #Find % outside the 99% Box area
    #print(Img.shape)
    #sumCharge = np.sum(ImgJ)
    widths = np.zeros(len(boxes))
    heights = np.zeros(len(boxes))
    pLargers = np.zeros(len(boxes)) #percent larger than initial box(?)
    #coreImax = np.zeros(len(boxes))
    #coreMeanI = np.zeros(len(boxes))
    widths[0] = 160
    heights[0] = 64
    boxLLxs = np.zeros(len(boxes))
    boxLLys = np.zeros(len(boxes))
    pOutsideBoxes = np.zeros(len(boxes))
    Pprotons = sumTot / parts * 100 #this is constant
    #print("Protons:",parts,sumTot,Pprotons,"\nCharge:",sumCharge,sumCharge/I_pulse*100)
    cols = ["r","cyan","gold","lime","k"]
    rdVal = args.reBin*2
    #if args.reBin > 1:
    #    print(args.reBin%2)
    #    if args.reBin % 2 == 0:
    #        reBinOffset = 0
    #    else: 
    #        reBinOffset = 0
    for i in range(len(boxes)): 
        pLargers[i] = boxes[i]
        widths[i] = round(widths[0]*(1+pLargers[i]) / 2) * 2
        heights[i] = round(heights[0]*(1+pLargers[i]) / 2) * 2 #+ reBinOffset
        boxLLxs[i] = -round((widths[i]/2) / 2) * 2 #-90
        boxLLys[i] = -round((heights[i]/2) / 2) * 2 #+ reBinOffset #-32
        #print(len(xax),boxLLxs[i],boxLLys[i],widths[i],heights[i],round(boxLLys[i]/rdVal)*rdVal)

        #11 Mar Verified the Box is the correct locations, the beam is just displaced when reBinned!
            #But Edges (tight) are correct as well
        boxBInd=460;boxTInd=540;boxLInd=455;boxRInd=545
        for idx, val in np.ndenumerate(xax):
            if int(val) == round(boxLLxs[i]/rdVal)*rdVal:
                #print(idx[0])
                boxLInd = idx[0] #455
            if int(val) == round(boxLLxs[i]/rdVal)*rdVal+round(widths[i]/rdVal)*rdVal:
                #print(idx[0])
                boxRInd = idx[0] #545
        for idx, val in np.ndenumerate(yax):
            #print(idx[0],int(val),round(boxLLys[i]/rdVal)*rdVal)
            if int(val) == round(boxLLys[i]/rdVal)*rdVal:
                boxBInd = idx[0] #460
            if int(val) == round(boxLLys[i]/rdVal)*rdVal+round(heights[i]/rdVal)*rdVal:
                #print(idx[0])
                boxTInd = idx[0] #540
        #Find % outside the 99% Box area
        core = Img[boxBInd:boxTInd,boxLInd:boxRInd]
        sumCore = np.sum(core)+1
        pOutsideBoxes[i] = (sumTot-sumCore)/sumTot*100
        #print(sumCore,pOutsideBoxes[i])
        #coreN = Img[boxBInd:boxTInd,boxLInd:boxRInd] #why a 2nd?
        #coreImax[i] = core.max() #[uA/cm^2]
        #coreMeanI[i] = np.mean(core) #[uA/cm^2]

    #from plotFit import PMAS
    sPMAS=datetime.now()#.strftime("%H-%M-%S"))
    jMax,pOutsideBox,rValue,edges,EI,beamArea,centX,centY,chi2 = PMAS(ImgJ,args,parts,xax,yax,name,paths,sample) #,dispY,dispX
    fPMAS = datetime.now()
    #print(EI)#"finished PMAS in",datetime.now()-sPMAS)#.strftime("%H-%M-%S"))

    #Warning print is outside PMAS because it takes a long time
    jMaxLim = 53
    pOutLim = 4
    centYLim = 3
    centXLim = 4
    rValLim = 0.10
    nomArea = 3200 #~32*100, nominal gives 34*104=3536
    chi2Lim = 5e3
    #print warnings if exceed limit of value:
    """if jMax >= jMaxLim:
        print("\tCurrent Density: {:0.1f}".format(jMax),"uA/cm^2 greater than",jMaxLim,"uA/cm^2")
    if pOutsideBox >= pOutLim:
        print("\t% Outside Box:",pOutsideBox,"% greater than",pOutLim,"%")
    if args.reBin == 1:
        if beamArea <= nomArea:
            print("\tBeam Area:",beamArea,"mm^2, less than",nomArea,"mm^2")
        if centY >= centYLim:
            print("\tVertical Displacement:",centY,"mm greater than",centYLim)
        if centY <= -centYLim:
            print("\tVertical Displacement:",centY,"mm greater than",centYLim)
        if centX >= centXLim:
            print("\tVertical Displacement:",centX,"mm greater than",centXLim)
        if centX <= -centXLim:
            print("\tHorizontal Displacement:",centX,"mm greater than",centXLim)
        chi2Lim = 1e4
    if rValue >= rValLim:
        print("\tR Value: {:.3e} above Limit: {:.3e}".format(rValue,rValLim))
    if chi2 >= chi2Lim:
        print("\tChi^2: {:.2e} above Limit: {:.1e}".format(chi2,chi2Lim))
    """
    #print(jMax,pOutsideBox,rValue,chi2,beamArea,dispX,dispY)
    #core = Img[EI[1]:EI[0],EI[2]:EI[3]]
    #coreJMax = core.max()
    coreJMean = np.mean(Img[EI[1]:EI[0],EI[2]:EI[3]]) #decrease by 1 index
    coreArea = (edges[0]-edges[1])*(edges[3]-edges[2])
    #print("J in core ",coreArea,"mm^2 is {:.2f} uA/cm^2".format(coreJMean),sep="")#np.sum(Img[idxMinY:idxMaxY,idxMinX:idxMaxX])))
    #print("Beam Core Area: ",edges[3]-edges[2],"mm wide x ",edges[0]-edges[1],"mm high",sep="")
    #print("Beam Area:",beamArea,"mm^2, nominal:",128*42,"mm^2")
    #print("Beam moved",dispY,"vertical and",dispX,"horizontally")

    #Normalize to full current, see 
    #print(sumTot,parts,sumCore,pOutsideBox,Img.max())
    #Img[boxBInd:boxTInd,boxLInd:boxRInd] = 0
    C = I_pulse / parts / (xBinSize * yBinSize * 1e-2) * 0.04 #[uA/cm^2]: /mm^2->/cm^2 = /1e-2, A->uA = 1e6, 4% duty cycle
    #Protmax = Img.max() #protons/mm^2
    #for i in range(len(Img)):
    #    for j in range(len(Img[i])):
    #        Img[i][j] = Img[i][j] * C #[uA/cm^2]
    #jMax = Img.max()
    jMin = 0.9 * C * parts/3e5 #background (0 hits) will be un-colored, scaled to 1e5 parts
    print(sample,"PMAS:",fPMAS-sPMAS,"s: Jmax {:.1f}, J in core {:.0f}mm^2: {:.2f}uA/cm^2, RVal {:.5f}, Chi2 {:.2f}, % Out TA {:.2f}"\
                .format(jMax,coreArea,coreJMean,rValue,chi2,pOutsideBoxes[0])) #sample here is not sample in rasterScatter...
    ##print("Beam Top: {:.1f}mm, Bottom: {:.1f}mm, Left: {:.1f}mm, Right: {:.1f}mm".format(edges[0],edges[1],edges[2],edges[3]),\
    ##            "Beam Center at (X,Y): ({:.0f},{:.0f})".format(centX,centY))

    if args.gaussFit:
        Lorentz=True
                                           #(hist,     axis, width,maxim,options,name,   y1,y2,saveFits)
        diffNy,diffPy,coeffsy, differenceNLy,differencePLy,coeffsLy = gaussianFit(histogram2D,"y",yBinSize,500,savename,
                                                                                  2,30,args.saveFits,False,args.gauss2Fit,Lorentz)
        diffNx,diffPx,coeffsx, differenceNLx,differencePLx,coeffsLx = gaussianFit(histogram2D,"x",xBinSize,500,savename,
                                                                                  3,20,args.saveFits,False,args.gauss2Fit,Lorentz)
        #add minimize function for these?

    if args.savePics:
        from matplotlib.pyplot import subplots,pcolormesh,close,tight_layout,savefig,setp,title,xlabel,ylabel
        from matplotlib.axes import Axes
        from matplotlib.patches import Rectangle
        from matplotlib.colors import LogNorm
        X, Y = np.meshgrid(xax,yax) #Set mesh of plot from the axes given from converter function
        close() #make sure no prior plotting messes up

        fig,ax = subplots(dpi=args.dpi)
        #print(datetime.now().strftime("%H-%M-%S"))
    
        #Set maximum value depending on the maximum current density
        from math import log10,ceil
        maxim = ceil(jMax / 10) * 10
        if args.maxim != 0: maxim = args.maxim #user provided maximum
        minim = 10**ceil(log10(jMin))
        cbarVals  = [minim,minim*10,minim*100,maxim] #make array for color bar values
        cbarLabels = ["{:.2f}".format(cbarVals[0]),"{:.1f}".format(cbarVals[1]),"{:.1f}".format(cbarVals[2]),"{:.0f}".format(cbarVals[3])]
        cbarLabel = r"$\mu A / cm^2$"
        #print("Max current Density: ",Img.max(),"/",maxim,datetime.now().strftime("%H-%M-%S"))

        #Use pcolor to show density map, use log scale
        c = ax.pcolormesh(X,Y,Img,shading='auto',norm=LogNorm(vmin=minim, vmax=maxim), cmap='viridis',rasterized=True) #viridis or magma are perceptually uniform
        lw=1
        fs=15

        #Set Plot Properties
        setp(ax,xlim=([-args.xlim,args.xlim]),ylim=([-args.ylim-10,args.ylim]))
        ax.set_title("Macro-Particle Beam Distribution at "+position,fontsize=fs+2)
        ax.set_ylabel("Vertical [mm]",fontsize=fs)
        ax.set_xlabel("Horizontal [mm]",fontsize=fs)
        cbar = fig.colorbar(c, ax=ax,pad=0.01)
        cbar.set_label(cbarLabel,fontsize=fs-1)
        cbar.set_ticks(cbarVals)
        cbar.set_ticklabels(cbarLabels)

        if position =="Target":
            #name = name+"_2212only"
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

        #Display beam characteristics
        bgdbox=dict(pad=1,fc='w',ec='none')
        propsR = dict(horizontalalignment="right",verticalalignment="bottom", backgroundcolor = 'w',bbox=bgdbox,fontsize=fs-4.5, transform=ax.transAxes)
        propsL = dict(fontsize=fs-4.5, backgroundcolor = 'w',bbox=bgdbox,transform=ax.transAxes)

        #Show 99% box
        if not args.noBox: #for user clarity, call is noBox, include % outside text for this
            ax.add_patch(Rectangle((boxLLxs[0],boxLLys[0]),widths[0],heights[0],linewidth=lw,edgecolor=cols[0],fill=False))
            ax.text(0.01, 0.925, "{:.2f}".format(pOutsideBox)+r"% Outside 160x64mm$^2$ Target", color=cols[0], 
                    fontsize = fs-2, fontweight='bold', backgroundcolor = 'w',bbox=bgdbox,transform=ax.transAxes)
        else: #yes no Boxes
            name+="_noBox"

        if not args.noText: #for user clarity, call is noText
            ax.set_title("Distribution at "+position+"\n{:.2f}% of {:.2e} Total Macro-Particles".format(Pprotons,parts),fontsize=fs)
            ax.text(0.01, 0.175, "Beam Twiss at PBW:", fontsize=fs-3, backgroundcolor = 'w',bbox=dict(pad=1.5,fc='w',ec='none'),transform=ax.transAxes)
            ax.text(0.01, 0.13, r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[mm \cdot mrad]}$", propsL)
            ax.text(0.01, 0.08, r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$", propsL)
            ax.text(0.01, 0.03, r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]), propsL)
            #ax.text(0.99, 0.21, " r = {:.5f}".format(rValue), propsR)
            ax.text(0.99, 0.14, r"Core $\langle\bf{J}\rangle$: "+"{:.1f} ".format(coreJMean)+r"$\mu$A/cm$^2$", propsR)
            ax.text(0.99, 0.09, r"Peak $\langle\bf{J}\rangle$: "+"{:.1f} ".format(jMax)+r"$\mu$A/cm$^2$", propsR)
            ax.text(0.99, 0.05, "RM Amplitude: {:.1f}, {:.1f}".format(args.aX,args.aY)+r"$_{[mm]}$",propsR,fontsize=fs-6)
            ax.text(0.99, 0.01, options['physList'], propsR,fontsize=fs-6)

            if chi2 > 1e3:
                ax.text(0.99, 0.195, r"$\chi^2$"+"={:.2e}".format(chi2), propsR)
            else:
                ax.text(0.99, 0.195, r"$\chi^2$"+"={:.0f}".format(chi2), propsR)

            for i in range(1,len(boxes)): #make multiple boxes
                ax.add_patch(Rectangle((boxLLxs[i],boxLLys[i]),widths[i],heights[i],linewidth=lw,edgecolor=cols[i],fill=False))
                ax.text(xlim[0]*0.90, ylim[1]*(0.85-i*0.1), "{:.2f}".format(pOutsideBoxes[i])+"% Outside {:.0f}% Larger Box".format(pLargers[i]*100), 
                                  color=cols[i], fontweight='bold',fontsize = fs-2, backgroundcolor = 'w',bbox=dict(pad=1))
                                    #,path_effects=[path_effects.withStroke(linewidth=1, foreground='k')])
        else: #yes no text
            name+="_noText"


        if args.saveEdges:
            edgeCol = 'k'
            ax.hlines(edges[0],edges[2],edges[3],colors=edgeCol,linewidths=lw)
            ax.hlines(edges[1],edges[2],edges[3],colors=edgeCol,linewidths=lw)
            ax.vlines(edges[2],edges[1],edges[0],colors=edgeCol,linewidths=lw)
            ax.vlines(edges[3],edges[1],edges[0],colors=edgeCol,linewidths=lw)
            #print("Edges printed!!")
            name+="_Edges"
            ax.text(0.01, 0.88,"Beam Center:({:.1f},{:.1f}) [mm]".format(centX,centY),propsL)

        dt = datetime.now()
        #from os.path import isfile
        #if isfile(name+"*.png"):
        #  print("already present")
        tight_layout()
        savefig(name+"."+args.picFormat,dpi=args.dpi)#
        close(fig)
        close()
        print(name+"."+args.picFormat)#"_"+dt.strftime("%H-%M-%S")+
    #dt = datetime.now()
    #print(dt-start)

    #passing lists to decrease length of argument calls
    return [jMax,pOutsideBox,beamArea,coreJMean,centX,centY,rValue,chi2]












##Converts ROOT TH2D to 2D Numpy Image, includes to save histogram as a CSV 2D map
def converter(hIn,saveHist,name,paths,args,density):
    import ROOT
    from numpy import zeros,linspace

    if args.sim == "beamlet":
        nParts= args.Nbeamlet
    elif args.sim == "raster":
        t_pulse = round(0.04 * 1/14 * 1e6) # mus
        pulse_start = 10 #why have a delay?
        t_end = (2* pulse_start + t_pulse) /1000# - 2 #2.857 ms
        time_length = round(t_end * args.nP) #number of pulses, nominal = 2.86e3
        t = linspace(0,t_end,time_length) #array of steps of length time_length
        N_t = len(t)
        nParts = N_t*10*args.Nb

    if args.reBin > 1:
        #print("Rebinning")
        oldName=hIn.GetName()
        hIn = hIn.Rebin2D(args.reBin,args.reBin,oldName+"_"+str(args.reBin))
        #print("Shape of old:",hIn.GetXaxis().GetNbins(),"\nShape of New:",hIn.GetXaxis().GetNbins())

    ##Get X axis from ROOT
    xax = zeros(hIn.GetXaxis().GetNbins()+1)
    for i in range(len(xax)-1):
        xax[i] = hIn.GetXaxis().GetBinLowEdge(i+1) #Set elements as bin edges
    xax[-1] = hIn.GetXaxis().GetBinUpEdge(hIn.GetXaxis().GetNbins())

    ##Get Y axis from ROOT
    yax = zeros(hIn.GetYaxis().GetNbins()+1)
    for i in range(len(yax)-1):
        yax[i] = hIn.GetYaxis().GetBinLowEdge(i+1) #Set elements as bin edges
    yax[-1] = hIn.GetYaxis().GetBinUpEdge(hIn.GetYaxis().GetNbins())
    
    hOut = zeros( (len(yax), len(xax)) ) #2D Image map

    ##Scales to a density rather than particle count
    if density:
        integral = hIn.Integral()
        hIn.Scale(1/integral)
        print(integral,1/integral)

    ##Fill Image map with 2D histogram values
    for xi in range(hIn.GetXaxis().GetNbins()):
        for yi in range(hIn.GetYaxis().GetNbins()):
            bx = hIn.GetBin(xi+1,yi+1)
            if density:
                hOut[yi,xi] = hIn.GetBinContent(bx)#/nParts#*reBin #rebin now is summed, kind of... this messes up the counts if
                #           reBin!=1 is used, so this doesn't work.
            else:
                hOut[yi,xi] = hIn.GetBinContent(bx)
    ##Must add Overflow options!!!

    if saveHist: #for rValue calculations
        import os,csv,re
        if os.uname()[1] in {"tensor.uio.no","heplab01.uio.no", "heplab04.uio.no","heplab03.uio.no"}:
            ##Remove upper directories that may have come with name for appending outname to scratch folder
            #print(name)
            if re.search("/PBW_",name):
                #print("\n",outname,"\n")
                name = re.sub(".+(?=(PBW_))","",name) #substitutes "" for all preceeding PBW_
                #print("Histogram-removed",name)
            elif re.search("/HEBT-A2T_",name):
                #print("\n",outname,"\n")
                name = re.sub(".+(?=(HEBT-A2T_))","",name) #substitutes "" for all preceeding PBW_
                #print("Histogram-removed",name)
            elif re.search("/Vac_",name):
                #print("\n",outname,"\n")
                name = re.sub(".+(?=(Vac_))","",name) #substitutes "" for all preceeding Vac_
              #print("Histogram-removed",name)
        else:
            if re.search("/PBW_",name):
                #print("\n",outname,"\n")
                name = re.sub(".+(?=(PBW_))","",name) #substitutes "" for all preceeding PBW_
                #print("Histogram-removed",name)
            if not os.path.isdir(paths['scratchPath']+"targetCSVs/"):
                os.mkdir(paths['scratchPath']+"targetCSVs/")
        if not os.path.isfile(paths['scratchPath']+"targetCSVs/"+name+".csv"):
            if not os.path.isdir(paths['scratchPath']+"targetCSVs/"):
                os.path.mkdir(paths['scratchPath']+"targetCSVs/")
            with open(paths['scratchPath']+"targetCSVs/"+name+".csv",mode = 'w',newline=None) as hist_file:
                hist_writer = csv.writer(hist_file,delimiter = ',')
                hist_writer.writerows(hOut)
            hist_file.close()
            print("saved Distribution at Target",paths['scratchPath']+"targetCSVs/"+name+".csv")
        else: print("Distribution already present at:",paths['scratchPath']+"targetCSVs/"+name+".csv")

    return (hOut,xax,yax)















##WIP: Add function to check if ROOT file is complete, else delete and make new
def findRoot(savename,paths):
    from os.path import isfile
    locations = ["/scratch2/ericdf/PBWScatter/ESS/","/scratch2/ericdf/PBWScatter/pencil/","/scratch2/ericdf/PBWScatter/2Dmaps/"]
    pwd="/scratch2/ericdf/PBWScatter/"
    for i in range(len(locations)):
        if isfile(locations[i]+savename+".root"):
            pwd = locations[i]
        else: 
            print("Root file not found in",locations[i])
    print(pwd,savename,".root",sep="")

    ##need to add completeness check
    return pwd














##Image Comparison Measurements: returns Rvalue and chi^2 
def imgCompare(Im,paths,args):
    import numpy as np
    import os

    lenx = np.shape(Im)[0]-1
    leny = np.shape(Im)[1]-1

    ##Find reference file
    if os.uname()[1] in {"tensor.uio.no","heplab01.uio.no","heplab04.uio.no","heplab03.uio.no"}:
        if args.Nb == 100: #do fill in from args
            #print("imgCompare Nb 100")
            if args.reBin == 5:
                refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ_TargetImage_rB5.csv"
            elif args.reBin == 4:
                refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ_TargetImage_rB4.csv"
            elif args.reBin == 1:
                if args.beamClass == "jitter":
                    import re
                    if re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(Jitter)))",args.twissFile)[1] == "-04":
                        refImgName = "HEBT-A2T_100pctField_1.0e-04Jitter_250x_570MeV_OrigbX1085.63,bY136.06m_N2.9e+06_NpB100_runW_QBZ_rB1_refImg.csv"
                    elif re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(Jitter)))",args.twissFile)[1] == "-03":
                        refImgName = "HEBT-A2T_100pctField_1.0e-03Jitter_200x_570MeV_OrigbX1085.63,bY136.06m_N2.9e+06_NpB100_runW_QBZ_rB1_refImg.csv"
                elif args.beamClass == "ESS":
                    refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ_rB1_refImg.csv"
                else:
                    refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ_rB1_refImg.csv"
                
                #refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ_rB1_refImg.csv"
        elif args.Nb == 200:
            print("imgCompare Nb 200")
            refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N5.8e+06_NpB200_runW_QBZ_TargetImage_rB4.csv"
            if args.reBin == 1:#don't care about this case, just getting pics from it, need axis to be the same
                refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ_TargetImage_rB1.csv"
        elif args.Nb == 500:
            refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N5.8e+06_NpB200_runW_QBZ_TargetImage_rB4.csv"
        elif args.Nb == 50:
            refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N1.4e+06_NpB50_runW_QBZ_TargetImage_rB4.csv"
            if args.reBin == 1:#don't care about this case, just getting pics from it, need axis to be the same
                refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ_TargetImage_rB1.csv"
        else: #Nb=10
            refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+05_NpB10_runW_QBZ_TargetImage_rB4.csv"
            if args.reBin == 1:#don't care about this case, just getting pics from it, need axis to be the same
                refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ_TargetImage_rB1.csv"
            #Iref = np.genfromtxt(open(paths['scratchPath']+"Vac_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03_run_QBZ_TargetImage.csv"),delimiter=",")
            #Iref = np.genfromtxt(open("/scratch2/ericdf/PBWScatter/PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03_runW_QBZ_TargetImage.csv"),delimiter=",")
        
        if os.path.isfile(paths['scratchPath']+"targetCSVs/"+refImgName):
            Iref = np.genfromtxt(open(paths['scratchPath']+"targetCSVs/"+refImgName),delimiter=",")
        else:
            Iref = np.zeros((leny,lenx))
    else:
        Iref = np.zeros((leny,lenx))

    divSum = np.zeros((leny,lenx))
    ndof = 0 ## of pixels

    ##Compute Rvalue
    for i in range(lenx):
        for j in range(leny):
            #Iref[j,i] += offset
            if Iref[j,i] == 0: continue #required to not crash 15.1.23
            if Iref[j,i] <= args.threshold: continue #doesn't count noise
            divSum[j,i] = ( ( ( (Im[j,i]) / (Iref[j,i]) - 1) ** 2 ) / (leny * lenx) )
            #if div[j,i] != 0:
            #    print(j,i, div[j,i])
    rDiv = np.sqrt(np.sum(divSum))

    ##Compute chi^2 and pull values
    ##print("Reference image is:",refImgName)
    chi2 = 0
    pull_array = []
    for i in range(lenx):
        for j in range(leny):
            if Iref[j,i] == 0: continue
            if Iref[j,i] < args.threshold: continue
            ndof += 1
            pull = ( ( (Im[j,i] - Iref[j,i])  ) / np.sqrt(Iref[j,i]) )
            pull_array.append(pull)
            chi2 += pull * pull
            #print(pull,chi2)

    return rDiv,chi2,pull_array,ndof

















##Finds the edges of the beam using a Sobel filter
def findEdges(Img,jMax,graph,savename,xax,yax,args):
    from numpy import array
    from scipy.signal import convolve2d as scipySigConv2d

    gradXX = array([[1,0,-1],[2,0,-2],[1,0,-1]])
    gradYY = array([[1,2,1],[0,0,0],[-1,-2,-1]])

    gradX = scipySigConv2d(Img,gradXX,'same')
    gradY = scipySigConv2d(Img,gradYY,'same')

    ##Prominence of gradient to determine if actual edge or not
    promX = 3
    promY = 4

    EI,edges = localExtrema(gradX,promX,gradY,promY,xax,yax)
    #print(edges)

    ##Only do check for Top and Left?
    #if EI[0] == 0:
    #    promY -= 1
    #    print("Not strongly defined beam on Top",promY)
    #    EI,edges = localExtrema(gradX,promX,gradY,promY,xax,yax)
    #    if EI[1] == 0:
    #        promY -= 1
    #        print("Not strongly defined beam on Top",promY)
    #        EI,edges = localExtrema(gradX,promX,gradY,promY,xax,yax)
    #        if EI[1] == 0:
    #            promY -= 1
    #            print("Not strongly defined beam on Top",promY)
    #            EI,edges = localExtrema(gradX,promX,gradY,promY,xax,yax)
    if EI[2] == 0:
        promX -= 1
        print("Not strongly defined beam on Left",promX)
        EI,edges = localExtrema(gradX,promX,gradY,promY,xax,yax)
        if EI[1] == 0:
            print("Major Error Left Edge")

    ##Graph the gradients and sum of gradients per column/row for x,y, respectively
    if graph:
        import matplotlib.pyplot as plt
        from numpy import meshgrid, sum as npSum
        plt.close()
        X,Y = meshgrid(xax,yax)
        fig = plt.figure(figsize=(12,7))
        plt.subplots_adjust(wspace=0.35,hspace=0.35) #increase width space to not overlap
        ax1 = fig.add_subplot(2,2,1)
        ax2 = fig.add_subplot(2,2,3)
        ax3 = fig.add_subplot(2,2,2)
        ax4 = fig.add_subplot(2,2,4)
        a = ax1.pcolor(X,Y,gradX,shading='auto')
        b = ax2.pcolor(X,Y,gradY,shading='auto')
        c = ax3.scatter(xax,npSum(gradX,axis=0))#,shading='auto')
        d = ax4.scatter(yax,npSum(gradY,axis=1))#,shading='auto')
        #[topInd,botInd,lefInd,rigInd]
        ax3.vlines(edges[2],ax3.get_ylim()[0],ax3.get_ylim()[1],color='m')
        ax3.vlines(edges[3],ax3.get_ylim()[0],ax3.get_ylim()[1],color='m')
        ax4.vlines(edges[0],ax4.get_ylim()[0],ax4.get_ylim()[1],color='m')
        ax4.vlines(edges[1],ax4.get_ylim()[0],ax4.get_ylim()[1],color='m')
        plt.setp(ax1,xlim=(-args.xlim,args.xlim),ylim=(-args.ylim,args.ylim),title="(a) Gradient X",xlabel="X [mm]",ylabel="Y [mm]")
        plt.setp(ax2,xlim=(-args.xlim,args.xlim),ylim=(-args.ylim,args.ylim),title="(c) Gradient Y",xlabel="X [mm]",ylabel="Y [mm]")
        plt.setp(ax3,xlim=(-args.xlim,args.xlim),title="(b) Sum of Gradient X",xlabel="X [mm]",ylabel="Gradient Sum")#,ylim=(-500,500))
        plt.setp(ax4,xlim=(-args.xlim,args.xlim),title="(d) Sum of Gradient Y",xlabel="Y [mm]",ylabel="Gradient Sum")#,ylim=(-500,500))
        ac = fig.colorbar(a, ax=ax1,pad=0.01)
        ac.set_label(r"d$\bf{J}$/dx",labelpad=3,fontsize=12)
        bc = fig.colorbar(b, ax=ax2,pad=0.01)
        bc.set_label(r"d$\bf{J}$/dy",labelpad=3,fontsize=12)
        #plt.setp(ax3,xlim=(-100,100),ylim=(-100,100))
        plt.tight_layout()
        plt.savefig(savename+str(promY)+"_GradPics."+args.picFormat,bbox_inches='tight',dpi=args.dpi)
        print(savename,promY,"_GradPics.",args.picFormat,sep="")
        plt.close()

    return EI,edges















##Finds Edge Indices from Gradient maps
def localExtrema(gradX,nx,gradY,n,xax,yax):
    from numpy import sum as npSum,argpartition,max as npMax,min as npMin,ndenumerate#,ndenumerate,argmax,argmin,amax,amin,
    from math import ceil

    sumGradX = npSum(gradX,axis=0)
    sumGradY = npSum(gradY,axis=1)
    #print(len(sumGradY),npMin(sumGradY))
    a=1 #decrease core size by a

    for idx, val in ndenumerate(sumGradY):
        if val == npMin(sumGradY):
            topInd = idx[0]-a
            #print("top",val,idx[0])
        if val == npMax(sumGradY):
            botInd = idx[0]+a
            #print("bot",val,idx[0])
    for idx,val in ndenumerate(sumGradX):
        if val == npMin(sumGradX):
            rigInd = idx[0]-a
            #print("rig",val,idx[0])
        if val == npMax(sumGradX):
            lefInd = idx[0]+a
            #print("lef",val,idx[0])

    lim=round(len(xax-1)*0.15)
    #print(lim+round(len(xax)*0.5),-lim+round(len(xax)*0.5))
    #print(topInd,botInd,lefInd,rigInd)
    if topInd > lim+round(len(xax)*0.5):
        topInd = 0
        #print(avgMinSumGradY,avgMinSumGradY/n,avgMaxSumGradY)
        print("\n\nError in Gradient calculation of top Index!\n\n")
    elif botInd < -lim+round(len(xax)*0.5):
        botInd = 0
        #print(maxSumGradYInds,avgMaxSumGradY,avgMaxSumGradY/n)
        print("\n\nError in Gradient calculation of bottom Index!\n\n")
    elif lefInd < -lim+round(len(xax)*0.5):
        lefInd = 0
        print("\n\nError in Gradient calculation of left Index!\n\n")
    elif rigInd > lim+round(len(xax)*0.5): 
        rigInd = 0
        print("\n\nError in Gradient calculation of right Index!\n\n")

    ##Add a check in localExtrema that the gradient of the row/column on either side of the max is at least half(?) of the max. 
        ##Could also check if the index is high (<100,>900) and tell it to recalculate.
    ##ideally this would then delete the row for this reconing if False is returned and find the next highest...

    lMaxX = xax[lefInd]
    lMinX = xax[rigInd] 
    lMaxY = yax[botInd]
    lMinY = yax[topInd]

    return [topInd,botInd,lefInd,rigInd], [lMinY,lMaxY,lMaxX,lMinX] 













##Detection Algorithm!!!
def PMAS(Img,args,parts,xax,yax,name,paths,sample):
    import numpy as np
    from numpy import ndenumerate,sum as npSum,max,abs

    xBinSize = xax[round(len(xax)*0.5)+1]-xax[round(len(xax)*0.5)]
    yBinSize = yax[round(len(xax)*0.5)+1]-yax[round(len(xax)*0.5)]

    ##Find Target Area box (99% within)
    width = 160
    height = 64
    boxLLx = round(-width/2 / 2) * 2 #-90
    boxLLy = round(-height/2 / 2) * 2 #-80
    #print(len(xax),boxLLx,boxLLy,width,height)
    rdVal=args.reBin*2
    boxBInd=460;boxTInd=540;boxLInd=455;boxRInd=545
    for idx, val in ndenumerate(xax):
        if int(val) == round(boxLLx/rdVal)*rdVal:
            boxLInd = idx[0] #bc it is rounding
        if int(val) == round(boxLLx/rdVal)*rdVal+round(width/rdVal)*rdVal:
            boxRInd = idx[0] #545
    for idx, val in ndenumerate(yax):
        if int(val) == round(boxLLy/rdVal)*rdVal:
            boxBInd = idx[0] #bc it is rounding
        if int(val) == round(boxLLy/rdVal)*rdVal+round(height/rdVal)*rdVal:
            boxTInd = idx[0] #540
    ##Find % outside the 99% Box area
    core = Img[boxBInd:boxTInd,boxLInd:boxRInd]
    sumTot = npSum(Img)+1
    sumCore = npSum(core)+1
    pOut = (sumTot-sumCore)/sumTot*100
    #print(pOut)

    ##R value for algorithm. Works when use Current Density, not Nprotons
    rValue=0;chi2=0;ndof = 0;pull=0
    #addRefImg(Img,paths,args,sample)
    rValue,chi2,pull,ndof = imgCompare(Img,paths,args)
    if args.savePull:
        fs=16
        X, Y = np.meshgrid(xax,yax)
        import matplotlib.pyplot as plt
        #from scipy.stats import chi2 as chiSq
        mu,sigma,ampl,bins = findFit(pull,[0.1,1,1],(0,[1e3,3e5,3e5]),30)
        plt.close()
        fig,ax = plt.subplots(dpi=args.dpi)
        ax.hist(pull,bins)
        ax.plot(bins, gaussian(bins,ampl,mu,sigma), "r--", linewidth=2)

        bgdbox=dict(pad=1,fc='w',ec='none')
        propsR = dict(horizontalalignment="right",verticalalignment="top", backgroundcolor = 'w',bbox=bgdbox,fontsize=fs-4, transform=ax.transAxes)
        ax.text(0.98, 0.95,r"$\mu$="+"{:.3f}".format(mu)+"\n"+r"$\sigma$="+"{:.3f}".format(sigma), propsR)

        if max(pull) < 5: xlim = 5
        elif max(pull) <10: xlim = 10
        elif max(pull) <15: xlim = 15
        elif max(pull) <20: xlim = 20
        else: xlim=max(pull)
        plt.setp(ax,xlim=([-xlim,xlim]))
        ax.set_title("Pull Values for Beam\nThreshold: "+str(args.threshold)+" particles; "+str(len(pull))+" Pixels Counted",fontsize=fs+2)
        ax.set_ylabel("Counts",fontsize=fs)
        ax.set_xlabel("Pull Values",fontsize=fs)
        plt.savefig(name+"_pull"+str(args.threshold)+".png")
        print("pull",name+"_pull"+str(args.threshold)+".png")

    ##Normalize to full current
    I_pulse = 62.5*1e3 #[uA]
    C = I_pulse / parts / (xBinSize * yBinSize * 1e-2) * 0.04 #[uA/cm^2]: /mm^2->/cm^2 = /1e-2, A->uA = 1e6, 4% duty cycle
    for i in range(len(Img)):
        for j in range(len(Img[i])):
            Img[i][j] = Img[i][j] * C #[uA/cm^2] #redefinition was messing with proton sum, but not necessary for final algorithm
    jMax = max(Img)

    ##Find Edges
    #from plotFit import findEdges
    EI,edges = findEdges(Img,jMax,args.saveGrads,name,xax,yax,args) #[topInd,botInd,lefInd,rigInd]
    beamArea = (abs(edges[0])+abs(edges[1]))*(abs(edges[3])+abs(edges[2]))
    #Need to display beam Center index in X,Y, not displacement 
    centX = (edges[2] + edges[3]) / 2 #+ round(xBinSize*0.5)
    centY = (edges[0] + edges[1]) / 2 #+ round(yBinSize*0.5)
    #print("Beam Center at ",centX,",",centY)

    return jMax,pOut,rValue,edges,EI,beamArea,centX,centY,chi2 #make this a list!











##Save run values in stats CSV
def saveStats(statsPWD,Twiss,rasterBeamFile,jMax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue,chi2,fileName,reBin,args):
    import csv
    from os import path as osPath
    from re import sub
    rasterBeamFile = sub(".+(?<=(/CSVs/))","",rasterBeamFile) #all before CSVs gets removed
    if fileName != "":  
            statsName = fileName+".csv"
    else:
        statsName = statsPWD+args.statsFile+".csv" #don't forget to include statsPWD!!

    #print("Is file present?",osPath.isfile(statsPWD+statsName))
    writeMore = False
    if not osPath.isfile(statsName): #does this work?
        with open(statsName,mode = 'w') as csv_file:
            csv_writer = csv.writer(csv_file,delimiter = ',')
            csv_writer.writerow(["BeamFile","Emittance X [mm-mrad]","Beta X [m]","Alpha X","Emittance Y [mm-mrad]","Beta Y [m]",
                            "Alpha Y","Peak Current Desnity [uA/cm2]","% Outside Boxes","Core Area","Core J Mean","Center X",
                            "Center Y","R Value","Chi^2"])
            csv_file.close()
        writeMore = True
        print("new file made",statsName)
    elif osPath.isfile(statsName):
        #print("saving values")
        found = False
        foundRow = 0
        i = 0
        with open(statsName,mode='r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                i += 1
                if row[0] == args.beamFile: #No duplicates
                    found = True
                    foundRow = i 
                    break
            csv_file.close()
        if found:
            print("Values already in row",foundRow)
        else:
            writeMore = True
    else:
        Exception("Something weird with CSV?")
    if writeMore: #if found, won't enter
        with open(statsName,mode = 'a') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter = ',')
            csv_writer.writerow([args.beamFile,Twiss[0],Twiss[1],Twiss[2],Twiss[3],Twiss[4],Twiss[5],jMax,pOutsideBoxes,beamArea,
                                 coreJMean,centX,centY,rValue,chi2])
        csv_file.close()
        print("saved",args.beamFile,"{:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.0f}, {:.0f}, {:.3f}, {:.1f}".format(jMax,pOutsideBoxes,
                                                                                                            beamArea,coreJMean,
                                                                                                            centX,centY,rValue,chi2))#,chi2)
        #,statsName,"\n"









##Finds the matching data for a set of runs in the statsCSV and plots it
def plotSpread(args,Twiss,statsPWD,paramName,ind,unit,paramLabel,pFitLims,paramBins,pHistLims,origBX,origBY,beamFile):
    import csv,re
    from matplotlib.pyplot import subplots,hist,savefig,close,plot,xlim,ylim,text,title,xlabel,ylabel,tight_layout,xticks
    from numpy import greater,zeros,mean,std
    #from plotFit import findFit, gaussian
    from math import floor,log10

    read = zeros(args.samples)
    read.fill(-750) #allow better filter than nonzero?
    #print(Twiss,beamFile)
    i=0
    lenbreak=False

    ##How to select which entries? Make keys to search for
    pbwKey = "PBW_{:.0f}MeV".format(args.energy)
    nBkey = "{:.0f}_NpB{:.0f}".format(floor(log10(2.88e4*args.Nb)),args.Nb)
    betaKey = "beta{:.2f},{:.2f}m".format(Twiss[1],Twiss[4])
    origKey = "OrigbX{:.2f},bY{:.2f}".format(origBX,origBY)
    failKey = "failure{:.0f}-{:.0f}f".format(args.failure,args.magFails)
    qpKey = "QP"+args.qpNum
    ##Depending on the options specified, read the lines that match options. This can get tricky...
    with open(statsPWD+args.statsFile+".csv",mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        if args.beamClass == "ESS" and args.failure != 0: #RM fails, though only "PBW..." files at this time
            for row in csv_reader:
                #print(args.beamFile,row[0])
                if re.search(nBkey,row[0]) and re.search(failKey,row[0]) and re.search(pbwKey,row[0]) and re.search(betaKey,row[0]) :
                    #print(args.beamFile,row[0])
                    read[i] = float(row[ind]) #[mm-mrad]
                    if i == args.samples-1:
                        lenbreak = True
                        break
                    i+=1
        elif args.beamClass == "jitter": #for jitter studies from OXAL
            pctKey = re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(pct)))",args.twissFile)[1]+"pctField"
            jitterKey = re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(Jitter)))",args.twissFile)[1]+"Jitter"
            #print(pctKey,nBkey,jitterKey,origKey)
            if args.physList == "FTFP_BERT_EMZ":
                print("in FTFP")
                for row in csv_reader:
                    #if i in {0,1,10,50,100,150,199}: #for finding indices out of 7sigma
                    #    print("index",i,row[0],args.beamFile,re.search(args.beamFile,row[0]))
                    #print("nominal",re.search(failKey,row[0]))#re.search(nBkey,row[0]) , re.search(pbwKey,row[0]) , (re.search(betaKey,row[0])))
                    if re.search(nBkey,row[0]) and re.search(args.twissFile,row[0]) and re.search(pctKey,row[0]) and re.search("FBZ",row[0]):
                        read[i] = float(row[ind]) #[mm-mrad]
                        if i == args.samples-1:
                            lenbreak = True
                            break
                        i+=1
            else:
                for row in csv_reader:
                    #print(args.beamFile,row[0])
                    #if i in {1,10,50,100,150,199}: #for finding indices out of 7sigma
                    #    print("index",i,"file:",row[0],row[ind])
                    #print(re.search(nBkey,row[0]) , re.search(jitterKey,row[0]) , (re.search(origKey,row[0])))
                    if re.search(nBkey,row[0]) and re.search(args.twissFile,row[0]) and re.search(pctKey,row[0]):
                        read[i] = float(row[ind])
                        if i == args.samples-1:
                            lenbreak = True
                            break
                        i+=1
        elif args.beamClass == "qpFail": #for individual QP spread given in qpNum from OXAL twissFile 
            pctKey = re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(pct)))",args.twissFile)[1]+"pctField"
            jitterKey = re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(Jitter)))",args.twissFile)[1]+"JitterJL"
            print("qpFail",pctKey,nBkey,origKey,qpKey,jitterKey)
            for row in csv_reader:
                #if i == 6: #for finding indices out of 7sigma
                #    print("index",i,"failed file:",row[0])
                #print(re.search(nBkey,row[0]) , re.search(jitterKey,row[0]) , (re.search(origKey,row[0])))
                if re.search(nBkey,row[0]) and re.search(origKey,row[0]) and re.search(qpKey,row[0]) and re.search(pctKey,row[0]):# and re.search(jitterKey,row[0]):
                    read[i] = float(row[ind])
                    if i == args.samples-1:
                        lenbreak = True
                        #print(row[0])
                        break
                    i+=1
        else: #nominal
            #print("in else")
            if args.physList == "FTFP_BERT_EMZ":
                #print("in FTFP")
                for row in csv_reader:
                    #if i == 0: #for finding indices out of 7sigma
                    #    print("index",i,row[0],args.beamFile,re.search(args.beamFile,row[0]))
                    #print("nominal",re.search(failKey,row[0]))#re.search(nBkey,row[0]) , re.search(pbwKey,row[0]) , (re.search(betaKey,row[0])))
                    if re.search(nBkey,row[0]) and re.search(pbwKey,row[0]) and re.search(betaKey,row[0]) and re.search("FBZ",row[0]) and not re.search("failure",row[0]):
                        #print(args.beamFile,row[0])
                        read[i] = float(row[ind]) #[mm-mrad]
                        if i == args.samples-1:
                            lenbreak = True
                            break
                        i+=1
            else: #nominal
                #print("in physList else")
                for row in csv_reader:
                    #if i == 44: #for finding indices out of 7sigma
                    #    print("index",i,"failed file:",row[0])
                    #print("nominal",re.search(failKey,row[0]))#re.search(nBkey,row[0]) , re.search(pbwKey,row[0]) , (re.search(betaKey,row[0])))
                    if re.search(nBkey,row[0]) and re.search(pbwKey,row[0]) and re.search(betaKey,row[0]) and not re.search("FBZ",row[0]) and not re.search("failure",row[0]) and not re.search("mrad",row[0]):
                    #if re.search(args.beamFile,row[0]):
                        #print(args.beamFile,row[0])
                        read[i] = float(row[ind]) #[mm-mrad]
                        if i == args.samples-1:
                            lenbreak = True
                            break
                        i+=1
    csv_file.close()

    nonEmpty = greater(read,-750) #remove unused elements
    read = read[nonEmpty]

    #if ind == 13:
    #    print(read)
    #pFitLims = (0,[args.samples,100,500])
    print("N",len(read),"guess",paramName, "A",args.samples/4,"Mu",mean(read),"s",std(read),"bins",paramBins,"fitLims",pFitLims)
    #sanity check on values
    #print(mean(read),len(read))
    for i in range(len(read)):
        if read[i] < mean(read)-std(read)*7 or read[i] > mean(read)+std(read)*7:
            print(read[i],", index",i,"is outside 7 sigma!")

                              #gaussian(x,  amplitude,          mu,      sigma)
    mu, sigma,ampl,interval = findFit(read,[args.samples/4,mean(read),std(read)],pFitLims,paramBins)
    #print(mu,sigma,ampl)
    close()
    fig,ax = subplots(dpi=args.dpi,figsize=(6.4,4.8))
    fs=16
    if paramName == "chiSq":#chiSq
        _, bins, _ = ax.hist(read,interval,density=True)
        from scipy.stats import chi2 as chiSq
        from numpy import linspace
        #mean, var, skew, kurt = chiSq.stats(mu, moments='mvsk')
        x = linspace(chiSq.ppf(0.001, mu), chiSq.ppf(0.999, mu), 100)
        ax.plot(x, chiSq.pdf(x, mu),'r--', label='chi2 pdf')
        #rv = chiSq(mu)
        #ax.plot(x, rv.pdf(x), 'r--', lw=2, label='frozen pdf')
        ax.set_ylabel("PDF",fontsize=fs)
    else:
        _, bins, _ = ax.hist(read,interval)
        ax.plot(bins, gaussian(bins,ampl,mu,sigma), "r--", linewidth=2)
        ax.set_ylabel("Counts",fontsize=fs)

    ##Extract jitter parameters from args.twissFile
    nPs = round((2*10+round(0.04*1/14*1e6)/1e3)*args.nP)*args.Nb
    if re.search("Jitter",args.twissFile):
        if int(re.search("(([0-9]*)+(?=Jitter))",args.twissFile)[1]) == 4:
            pct = 1
        elif int(re.search("(([0-9]*)+(?=Jitter))",args.twissFile)[1]) == 3:
            pct = 10
        ax.set_title("{:.0f} Samples of {:.0f}% ".format(len(read),pct)+r"$\beta$"+" Variation Around Nominal\nwith {:.2e} Macro-Particles".format(nPs),fontsize=fs)
        if args.qpNum != "":
            fieldFrac = int(re.search("(([0-9]*)+(?=pctField))",args.twissFile)[1])
            ax.set_title("{:.0f} Samples of QP{:.0f} at {:.0f}% Field\n with {:.0f}% ".format(len(read),int(args.qpNum),fieldFrac,pct)+r"$\beta$"+" Variation Around Nominal \nwith {:.2e} Macro-Particles".format(nPs),fontsize=fs)
    else:#nominal
        ax.set_title("{:.0f} Nominal Samples \nwith {:.2e} Macro-Particles".format(len(read),nPs),fontsize=fs)
        pct=0

    ##Set more plot settings
    if args.failure == 0 and args.qpNum == "" and args.Nb > 99:
        ax.set_xlim(pHistLims)
    elif args.Nb < 99:
        if ind in{7,8,10}:
            ax.set_xlim(pHistLims)
    #elif args.qpNum != "":
    #    if ind in{7,10}:
    #        xlim(pHistLims)
    #    if ind == 8:
    #        xlim(2.85,5.1)

    ax.set_xlabel(paramLabel,fontsize=fs)
    #setp(ax,title=Title,xlabel=xLabel,ylabel=yLabel)

    ##Values are at reasonable positions in plots due to ax.transAxes
    bgdbox=dict(pad=1,fc='w',ec='none')
    if mu < (pHistLims[0]+(pHistLims[1]-pHistLims[0])*0.5):
        propsR = dict(horizontalalignment="right",verticalalignment="top", backgroundcolor = 'w',bbox=bgdbox,fontsize=fs-3,transform=ax.transAxes)
        ax.text(0.99, 0.97, "Beam Twiss at PBW:", propsR)
        ax.text(0.99, 0.91, r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[\mu m]}$",propsR)
        ax.text(0.99, 0.84, r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$", propsR)
        ax.text(0.99, 0.77, r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]), propsR)
        if sigma <1e-2:
            ax.text(0.99, 0.65,r"$\mu$="+"{:.3f}".format(mu)+unit+"\n"+r"$\sigma$="+"{:.2e}".format(sigma)+unit, propsR)
        else:
            ax.text(0.99, 0.65,r"$\mu$="+"{:.3f}".format(mu)+unit+"\n"+r"$\sigma$="+"{:.3f}".format(sigma)+unit, propsR)
        if args.physList == "FTFP_BERT_EMZ":
            ax.text(0.99, 0.50, args.physList, propsR)
    else:
        propsR = dict(horizontalalignment="left",verticalalignment="top", backgroundcolor = 'w',bbox=bgdbox,fontsize=fs-3,transform=ax.transAxes)
        ax.text(0.01, 0.97, "Beam Twiss at PBW:", propsR)
        ax.text(0.01, 0.91, r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[\mu m]}$",propsR)
        ax.text(0.01, 0.84, r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$", propsR)
        ax.text(0.01, 0.77, r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]), propsR)
        if sigma <1e-2:
            ax.text(0.01, 0.65,r"$\mu$="+"{:.3f}".format(mu)+unit+"\n"+r"$\sigma$="+"{:.2e}".format(sigma)+unit, propsR)
        else:
            ax.text(0.01, 0.65,r"$\mu$="+"{:.3f}".format(mu)+unit+"\n"+r"$\sigma$="+"{:.3f}".format(sigma)+unit, propsR)
        if args.physList == "FTFP_BERT_EMZ":
            ax.text(0.01, 0.50, args.physList, propsR)    

    if args.twissFile == "":
        name=statsPWD+"Nb{:.0f}_{:.0f}x".format(args.Nb,len(read))+paramName+"Hist"
    else:
        name=statsPWD+args.twissFile+"Nb{:.0f}_{:.0f}x".format(args.Nb,len(read))+paramName+"Hist"

    if args.failure != 0:
            name +="_failure"+str(int(args.failure))
    if args.physList == "FTFP_BERT_EMZ":
        name += "_FBZ"
    if args.saveSpread:
        tight_layout()
        savefig(name+"."+args.picFormat,dpi=args.dpi)#,pad_inches=0)
    close()
    print(len(read),paramName,ampl,mu,sigma)
    print(name+"."+args.picFormat)
    return mu,sigma,ampl,len(read)












##Controls plotSpread function, sends the plot info which varies with type of run
def spreadHist(args,Twiss,paths,origBX,origBY,beamFile):
    #from os import uname
    #from plotFit import plotSpread,plotSpreadBroad
    from numpy import zeros
    print("statsFile:",paths['statsPWD']+args.statsFile)

    paramName=["jMax","beamPOut","coreJMean","chiSq"]#"rValue","coreArea","coreJMean","centX","centY","chi2"] #len=6 
    paramLabel=[r"Peak Current Density [$\mu$A/cm$^2$]","Beam % Outside Target Area",r"Core Average Current Density [$\mu$A/cm$^2$]",r"$\chi^2$"]
                #r"Core Area [mm$^2$]","Beam Center X [mm]","Beam Center Y [mm]",,r"$r$ Comparison"
                #]#,"R Difference Value"]
    #"Peak Current Density [uA/cm^2]","Beam % Outside Target Area","Core Area [mm^2]",
        #"Core Average Current Density [uA/cm^2]","Beam Center X [mm]","Beam Center Y [mm]","R Divide Value","R Difference Value"
    ind = [7,8,10,14]#,14]#9,11,12,
    unit=[r"$\mu$A/cm$^2$","%",r"$\mu$A/cm$^2$",""]#r"mm$^2$","mm","mm"]#,"",""
                #   ampl, mu, sigma                  ([10,1e3,20],[5e3,1e4,1e3])
                        #jMax             Pout               coreJMean          chi2                                       rVal    
    paramFitLims = [(0,[1e3,500,500]),(0,[1e3,100,500]),(0,[1e3,500,500]),(0,[1e3,1e7,1e5])]#([10,1e1,20],[5e6,1e6,1e6]),(0,[1e3,10,1]),
                        #centX                                  centY                   
                    #([-1e3,-100,-500],[1e3,100,500]),([-1e3,-100,-500],[1e3,100,500]),#,(0,[1e3,10,1])]
    #pHistLimsold = [[35,45],[8.2,8.55],[25,40],[.072,.073]]#,[4.55,4.6]]#nominal#[5000,7000],,[-10,10],[-10,10]
                    #jMax       Pout  coreJMean    chi2    coreArea    centX   centY         rVal 
    pHistLims = [[51.5,52.9],[3.62,3.72],[40,46.4],[1000,11000]]#for sameRefImg, 3800,6300. For individual RefImgs, 3800,5000 #,[5000,8000],[-10,10],[-10,10] ,[.025,.075]
    paramBins = [20,20,20,20]#,20,20,20,20,20]

    mus = zeros(len(paramName))
    sigmas = zeros(len(paramName))
    ampl = zeros(len(paramName))
    lens = zeros(len(paramName))
    for i in range(len(paramName)):
        print(paramName[i],paramFitLims[i],paramBins)
        mus[i], sigmas[i],ampl[i],lens[i] = plotSpread(args,Twiss,paths['statsPWD'],paramName[i],ind[i],unit[i],paramLabel[i],paramFitLims[i],paramBins[i],pHistLims[i],origBX,origBY,beamFile)

    import csv
    found=False
    i=0
    foundRow=0
    with open(paths['statsPWD']+"pOffStats.csv",mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            i+=1
            if row[0] == args.beamFile: #No duplicates
                found = True
                foundRow = i
                print("Values already in row",foundRow)
                break
        csv_file.close()
    if not found:
        with open(paths['statsPWD']+"pOffStats.csv",mode = 'a') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter = ',')
            for i in range(len(paramName)):
                if mus[i] < 5e2:
                    if mus[i] < 10:
                        if sigmas[i] < 1e-2:
                            csv_writer.writerow([args.beamFile,lens[i],paramName[i],"{:.3f} $\pm$ {:.3e}".format(mus[i],sigmas[i])])
                        else:
                            csv_writer.writerow([args.beamFile,lens[i],paramName[i],"{:.3f} $\pm$ {:.3f}".format(mus[i],sigmas[i])])
                    else:
                        csv_writer.writerow([args.beamFile,lens[i],paramName[i],"{:.1f} $\pm$ {:.2f}".format(mus[i],sigmas[i])])
                else: #chi2 
                    csv_writer.writerow([args.beamFile,lens[i],paramName[i],"{:.3e} $\pm$ {:.0f}".format(mus[i],sigmas[i])])

            csv_file.close()
    
        for i in range(len(paramName)):
            print(args.beamFile,lens[i],paramName[i],"{:.3f} +/- {:.3f}".format(mus[i],sigmas[i]))









##Gets Twiss of beam depending on which input arguments are given, loads if already present
def getTwiss(args,sample,paths):
    ##Get Twiss to run, depending on user configuration
    # Twiss= [NemtX,BetaX,AlphX,NemtY,BetaY,AlphY]
    if args.beamClass == 'Yngve': #smallest from Yngve
        Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
    elif args.beamClass == 'ESS': #from my OpenXAL calculation
        Twiss = [0.118980737408497,1085.63306926394,-65.1638312921654,0.123632934174567,136.062409365455,-8.12599512314246] #updated to HEBT-A2T Combo Twiss
    elif args.beamClass == 'pencil': #"pencil" beam of ~0 emittance
        Twiss = [0.0001,0.15,0,0.0001,0.15,0]
    else: #for twiss source, jitter and qpFail cases for an initial definition
        Twiss = [0.118980737408497,1085.63306926394,-65.1638312921654,0.123632934174567,136.062409365455,-8.12599512314246] 

    if args.twiss:
        Twiss = [args.twiss[0],args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5]]
        args.beamClass = "Twiss"

    ##For pulling from an already made CSV
    if args.source == "csv":
        import re
        if re.search("beta",paths['csvPWD']+args.beamFile):
            betaX = float(re.search("(([-+]?[0-9]*\.?[0-9]*)+(?=(,([-+]?[0-9]*\.?[0-9]*)+m_RM)))",args.beamFile)[0])
            betaY = float(re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(m_RM)))",args.beamFile)[0])
            print(betaX,betaY)

    ##The backup of the Twiss
    origBX = Twiss[1]
    origBY = Twiss[4]

    ##Tailored to OpenXAL failure input
    if args.source == "twiss":
        if args.twissFile != "":
            import csv
            if args.beamClass == "jitter": #for getting all samples out of jitter OpenXAL output file
                print("jitter study")
                ##
                ## To Do : Double check this is doing what I want, along with sample appending...
                ##
                #for running all Twiss in OpenXAL Combo Sequence output CSV
                with open(paths['oXALPWD']+args.twissFile+".csv",mode='r') as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    rowNum = 0 #to increment QP number
                    for row in csv_reader:
                        if rowNum == sample: #i #not sure how to go back to QP finding...
                            #print(sample,rowNum,row[0])
                            if float(row[1]) < 1e-4: #Be sure emittance is correctly formatted
                                um = 1e6 #if from OpenXAL 
                            elif float(row[1]) > 1e-4:
                                um = 1
                            #print("um",um)
                            Twiss[0] = float(row[1])*um #[mm-mrad]
                            Twiss[1] = float(row[2])
                            Twiss[2] = float(row[3])
                            Twiss[3] = float(row[4])*um #[mm-mrad]
                            Twiss[4] = float(row[5])
                            Twiss[5] = float(row[6])
                            break
                        rowNum+=1
                csv_file.close()
            elif args.beamClass == "qpFail": #changed to this is the single QP study, since only care about 1-3 per section
                print("In QPs",args.qpNum,args.twissRowN)
                #for running all QPs in OpenXAL Combo Sequence output CSV
                with open(paths['oXALPWD']+args.twissFile+".csv",mode='r') as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    for row in csv_reader:
                        #print(row[0],args.qpNum)
                        if args.qpNum == row[0]: #Must use one process or else does same QP!!
                            if float(row[1]) < 1e-4: #Be sure emittance is correctly formatted
                                um = 1e6 #if from OpenXAL 
                            elif float(row[1]) > 1e-4:
                                um = 1
                            #print("um",um)
                            Twiss[0] = float(row[1])*um #[mm-mrad]
                            Twiss[1] = float(row[2])
                            Twiss[2] = float(row[3])
                            Twiss[3] = float(row[4])*um #[mm-mrad]
                            Twiss[4] = float(row[5])
                            Twiss[5] = float(row[6])
                            break
                csv_file.close()
            else: #not sure what for, since single QP runs takes over above...
                print("Else case in getTwiss from twissFile")
                with open(paths['oXALPWD']+args.twissFile+".csv",mode='r') as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    for row in csv_reader:
                        if row[0] == args.qpNum:
                            if float(row[1]) < 1e-4:
                                um = 1e6 #if from OpenXAL 
                            elif float(row[1]) > 1e-4:
                                um = 1
                            #print("um",um)
                            Twiss[0] = float(row[1])*um #[mm-mrad]
                            Twiss[1] = float(row[2])
                            Twiss[2] = float(row[3])
                            Twiss[3] = float(row[4])*um #[mm-mrad]
                            Twiss[4] = float(row[5])
                            Twiss[5] = float(row[6])
                    print("read Twiss",Twiss)
                csv_file.close()
        else: print("I need a --twissFile argument")


    return Twiss,origBX,origBY


##WIP: Attempt to fit pencil beam with multiple gaussians with python
def fitGaussians(Img,Error,path,lim,axis,gauss):
    import numpy as np
    import matplotlib.pyplot as plt
    print(np.shape(Img))
    if axis in {"y","Y"}:
        lenX = np.shape(Img)[1]
        proj = Img[:,round(lenX/2)-1:round(lenX/2)+1]
        if Error.all() != 1:
            projErr = Error[:,round(lenX/2)-1:round(lenX/2)+1]
        else:
            projErr = None
        label = "Y"
        #print(lenX)
        #print(np.shape(proj))
        proj = (proj[:,0] + proj[:,1]) / 2
        if Error.all() != 1:
            projErr = (projErr[:,0] + projErr[:,1]) / 2
    elif axis in {"x","X"}:
        lenX = np.shape(Img)[0]
        proj = Img[round(lenX/2)-1:round(lenX/2)+1,:]
        if Error.all() != 1:
            projErr = Error[round(lenX/2)-1:round(lenX/2)+1,:]
        else:
            projErr = None
        label = "X"
        #print(lenX)
        #print(np.shape(proj))
        proj = (proj[0,:] + proj[1,:]) / 2
        if Error.all() != 1:
            projErr = (projErr[0,:] + projErr[1,:]) / 2
    #print(np.shape(proj))
    x = np.linspace(-lim,lim,len(proj))

    if Error.all() != 1:
        plt.errorbar(x,proj,yerr=projErr,ecolor="y",elinewidth=1,fmt='.',markersize=2)
    else:
        plt.scatter(x,proj,s=2)
    plt.yscale("log")
    plt.ylim(1e-10,3e-2)
    xlim=20
    plt.xlim(-xlim,xlim)
    plt.savefig(path+"proj"+label+".png",dpi=300)

    def func1(x,a1,s1):
        from numpy import exp as npExp
        return a1 * npExp(-x**2 / (2*s1**2))
    def func2(x,a1,s1,a2,s2):
        from numpy import exp as npExp
        return a1 * npExp(-x**2 / (2*s1**2)) + a2 * npExp(-x**2 / (2*s2**2))
    def func3(x,a1,s1,a2,s2,a3,s3):
        from numpy import exp as npExp
        return a1 * npExp(-x**2 / (2*s1**2)) + a2 * npExp(-x**2 / (2*s2**2)) + a3 * npExp(-x**2 / (2*s3**2))
    def func4(x,a1,s1,a2,s2,a3,s3,a4,s4):
        from numpy import exp as npExp
        return a1 * npExp(-x**2 / (2*s1**2)) + a2 * npExp(-x**2 / (2*s2**2)) + a3 * npExp(-x**2 / (2*s3**2)) + a4 * npExp(-x**2 / (2*s4**2))
    def func5(x,a1,s1,a2,s2,a3,s3,a4,s4,a5,s5):
        from numpy import exp as npExp
        return a1 * npExp(-x**2 / (2*s1**2)) + a2 * npExp(-x**2 / (2*s2**2)) + a3 * npExp(-x**2 / (2*s3**2)) + a4 * npExp(-x**2 / (2*s4**2)) + a5 * npExp(-x**2 / (2*s5**2))

    from scipy.optimize import curve_fit #https://www.datatechnotes.com/2020/09/curve-fitting-with-curve-fit-function-in-python.html
    from numpy import exp as npExp

    params1, covs1 = curve_fit(func1, x, proj,sigma=projErr)
    print("params1: ", params1)
    coeffs = params1
    print("covariance: ", covs1)
    yfit1 = params1[0] * npExp(-x**2 / (2*params1[1]**2))
    plt.plot(x, yfit1, label="1Gaussian")

    if gauss in {2,3,4,5}:
        params2, _ = curve_fit(func2, x, proj,p0=[params1[0],params1[1],1,10],sigma=projErr)#,bounds=((params1[0]-1e-3,params1[0]+1e-3),(params1[1]-1e-3,params1[1]+1e-3),(0,1e5),(0,1e5)))
        print("params2: ", params2)
        coeffs = params2
        #print("covariance: ", covs2)
        yfit2 = params2[0] * npExp(-x**2 / (2*params2[1]**2)) + params2[2] * npExp(-x**2 / (2*params2[3]**2))
        plt.plot(x, yfit2, label="2Gaussian")

    if gauss in {3,4,5}:
        params3, _ = curve_fit(func3, x, proj,p0=(params1[0],params1[1],1,10,1,30),sigma=projErr)
        print("params3: ", params3)
        coeffs = params3
        #print("covariance: ", covs3)
        yfit3 = params3[0] * npExp(-x**2 / (2*params3[1]**2)) + params3[2] * npExp(-x**2 / (2*params3[3]**2)) + params3[4] * npExp(-x**2 / (2*params3[5]**2))
        plt.plot(x, yfit3, label="3Gaussian")

    if gauss in {4,5}:
        #a=18
        #x1=[]
        #y1=[]
        #for i in range(a):
        #    x1.append(x[i])
        #    x1.append(x[-i])
        #    y1.append(proj[i])
        #    y1.append(proj[-i])
        #x1.remove(-76.0)
        #y1.remove(1.254665e-05)
        #print(x1,y1)
        #plt.scatter(x1,y1)
        #params4, _ = curve_fit(func2, x1, y1,p0=(params3[4],params3[5],0.1,100))
        params4, _ = curve_fit(func4, x, proj,p0=(params1[0],params1[1],params3[2],params3[3],1,30,1,50),sigma=projErr,bounds=([0,0,0,0,0,0,0,1],[1,1,1,10,1,60,1,60]))
        print("params4: ", params4)
        coeffs = params4
        #print("covariance: ", covs4)
        #yfit4 = params3[0] * npExp(-x**2 / (2*params3[1]**2)) + params3[2] * npExp(-x**2 / (2*params3[3]**2)) + params3[4] * npExp(-x**2 / (2*params3[5]**2)) + params4[2]*npExp(-x**2 / (2*params4[3]**2)) #p
        yfit4 = params4[0] * npExp(-x**2 / (2*params4[1]**2)) + params4[2] * npExp(-x**2 / (2*params4[3]**2)) + params4[4] * npExp(-x**2 / (2*params4[5]**2)) + params4[6]*npExp(-x**2 / (2*params4[7]**2)) #p
        plt.plot(x, yfit4, label="4Gaussian")

    if gauss == 5:
        params5, _ = curve_fit(func5, x, proj,sigma=projErr)
        print("params5: ", params5)
        coeffs = params5
        #print("covariance: ", covs5)
        yfit5 = params4[0] * npExp(-x**2 / (2*params4[1]**2)) + params4[2] * npExp(-x**2 / (2*params4[3]**2)) + params4[4] * npExp(-x**2 / (2*params4[5]**2)) + params4[6]*npExp(-x**2 / (2*params4[7]**2)) + params5[8] * npExp(-x**2 / (2*params5[9]**2))
        plt.plot(x, yfit5, label="5Gaussian")

    plt.legend()
    plt.grid()
    plt.savefig(path+"proj"+label+".png",dpi=600)
    print(path+"proj"+label+".png")
    plt.close()


    return coeffs