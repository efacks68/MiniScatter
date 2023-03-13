#plotFit.py
#Eric Fackelman
#29 March- Present

#This holds the functions to plot the figures for the analysisScripts
#-plotFit - 

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
    plt.savefig(name,bbox_inches="tight")
    #plt.show()
    plt.close() #be sure to close the plot














def gaussian(x,amplitude,mu,sigma):
    from numpy import exp as npExp
    return amplitude * npExp( - (x - mu) ** 2 / (2 * sigma ** 2))














def findFit(data,guess,lims,nBins):
    #with help from MSeifert in stackoverflow fit a curve to a histogram in python
    #https://stackoverflow.com/questions/35544233/fit-a-curve-to-a-histogram-in-python
    from numpy import diff,linspace
    from matplotlib.pyplot import hist,plot,close,savefig
    from scipy.optimize import curve_fit
    #from plotFit import gaussian

    close()
    #print(data)
    bin_heights, bin_borders, _ = hist(data,bins=nBins,label="histogram")
    pwd="/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    savefig(pwd+"hist.png"); print(pwd+"hist.png")
    #print(len(bin_borders))
    bin_centers = bin_borders[:-1] + diff(bin_borders) / 2

    popt,_ = curve_fit(gaussian, bin_centers, bin_heights,p0=guess,bounds=lims) #p0 should be good start
    x_interval_for_fit = linspace(bin_borders[0],bin_borders[-1],len(bin_borders))
    plot(x_interval_for_fit, gaussian(x_interval_for_fit,*popt),label="fit")
    #plt.legend()
    #plt.xlim([bin_borders[0],bin_borders[-1]])
    #plt.show()
    #pwd="/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    #savefig(pwd+"fit.png")
    close() #be sure to close this or else these show up in the multi plots!
    #print(popt)
    #return the mean, abs(sigma), interval #sigma can sometimes be - so must abs() it.
        #   mu,     sigma,       ampl
    return(popt[1],abs(popt[2]),popt[0],x_interval_for_fit) #for Vac and Air sigma sometimes give - number...













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













def getMoments(Twiss):
    from numpy import sqrt
    mm=1e-3
    #[beta,alpha,gemt,gamma]
    sigma = sqrt(Twiss[2]*Twiss[0]) /mm #sigma x rms = beam size = sqrt(beta*epsG)
    sigmapx = sqrt(Twiss[2]*Twiss[3]) #mrad actually
    #print(sigma,sigmapx,sigpx)
    return [sigma,sigmapx]













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
    mux,sigmax,amplx,xinterval =    findFit(xs, [0.01,0.1,10],(0,[1,10,100]))
    mupx,sigmapx,amply,pxinterval = findFit(pxs,[0.01,0.1,10],(0,[1,10,100])) #this is a least square fit to Gaussian, so less sensitive to outliers.
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
    plt.savefig(name,bbox_inches="tight")
    plt.close()
    print(name)













def toTarget(TwPm,label):
    #Extend the distributions to the Target Location
    import numpy as np
    #Twiss=[beta,alpha,gemt,gamma]
    PBWexitBetaMx = np.array([[TwPm[0],-TwPm[1]],[-TwPm[1],TwPm[3]]])

    d_PBW_Targ = 4.4 #[m] Is this correct?
    drift_PBW_Targ = np.array([ [1, d_PBW_Targ],[ 0, 1]])
    Calc_PBW_Targ = np.linalg.multi_dot([drift_PBW_Targ,PBWexitBetaMx,np.transpose(drift_PBW_Targ)])
    #print(label,"PBWexit:",PBWexitBetaMx,"\n",drift_PBW_Targ,"\n",Calc_PBW_Targ)

    return [Calc_PBW_Targ[0][0],-Calc_PBW_Targ[1][0],TwPm[2],Calc_PBW_Targ[1][1]]













def compareTargets(targx,targy,targTwx,targTwy,fitTwx,fitTwy,fitlabel,savename,mat,PBWTwx,PBWTwy,args):
    #Now compare the MiniScatter Target distribution (targxTwissf) to initTarg, exitTarg, e8Targ and e16Targ PDFs
    import numpy as np
    #from plotFit import findFit, getMoments
    from scipy.stats import norm 
    import matplotlib.pyplot as plt
    mm = 1e-3
    fs = 14

    #Find fit to histogram
    mux,sigmax,amplx,xinterval = findFit(targx,[0.01,0.1,15],(0,[5e10,7e10,7e10]),"auto") #to accomodate 1e7 protons
    muy,sigmay,amply,yinterval = findFit(targy,[0.01,0.1,15],(0,[5e10,7e10,7e10]),"auto")

    #Find range of particles that are outside 3 sigma
    sigx=np.abs(targx)>3*sigmax# and np.abs(xs)<10*sigma)
    sigy=np.abs(targy)>3*sigmay

    #Find % of particles outside 3 sigma
    pOut3sigx = len(targx[sigx])/len(targx)*100
    pOut3sigy = len(targy[sigy])/len(targy)*100
    print("Distribution at Target: {:0.3f}% outside 3 sigma X".format(pOut3sigx))
    print("Distribution at Target: {:0.3f}% outside 3 sigma Y".format(pOut3sigy))

    #Get sigma for position and angle
    Mfx = getMoments(fitTwx)
    Mfy = getMoments(fitTwy)
    Mhx = getMoments(targTwx)
    Mhy = getMoments(targTwy)
    #print(fitlabel,Mhx,Mfx)

    #Sigmas of Histogram vs fit
    print("\nDistribution at Target \t sigma X: {:.2f}mm".format(sigmax))
    print("PDF at Target \t\t sigma X: {:.2f}mm".format(Mhx[0]))
    print(fitlabel,"\t\t sigma X: {:.2f}mm".format(Mfx[0]))

    print("\nDistribution at Target \t sigma Y: {:.2f}mm".format(sigmay))
    print("PDF at Target \t\t sigma Y: {:.2f}mm".format(Mhy[0]))
    print(fitlabel,"\t\t sigma Y: {:.2f}mm".format(Mfy[0]))

    #Create the fig with 2 plots side by side
    plt.clf()
    fig = plt.figure(figsize=(15,6.0))
    plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
    s1 = fig.add_subplot(1,2,1)
    s2 = fig.add_subplot(1,2,2)

    #Make the histogram of the full energy distrubtion for X with findFit intervals
    nx, binsx, patchesx = s1.hist(targx, xinterval, density=True, log=True, facecolor="green", alpha=0.75,label="MiniScatter Target X Histogram")

    #Add the "best fit" line using the findFit mu and sigma for x and display sigma
    y1a = norm.pdf(binsx, mux, sigmax)
    l1a = s1.plot(binsx, y1a, "k--", linewidth=1,label="Histogram Least Square")
    y1b = norm.pdf(binsx, mux, Mhx[0])
    l1b = s1.plot(binsx, y1b, "r--", linewidth=2,label="Target Twiss RMS")
    y1c = norm.pdf(binsx, mux, Mfx[0])
    l1c = s1.plot(binsx, y1c, "b--", linewidth=2,label=fitlabel+" Twiss RMS")

    #Make the histogram of the full energy distrubtion for Y with findFit intervals
    ny, binsy, patchesy = s2.hist(targy, yinterval, density=True, log=True, facecolor="green", alpha=0.75,label="MiniScatter Target Y Histogram")

    #Add the "best fit" line using the findFit mu and sigma for Y and display sigma
    y2a = norm.pdf(binsy, muy, sigmay)
    l2a = s2.plot(binsy, y2a, "k--", linewidth=1,label="Histogram Least Square")
    y2b = norm.pdf(binsy, muy, Mhy[0])
    l2b = s2.plot(binsy, y2b, "r--", linewidth=2,label="Target Twiss RMS")
    y2c = norm.pdf(binsy, muy, Mfy[0])
    l2c = s2.plot(binsy, y2c, "b--", linewidth=2,label=fitlabel+" Twiss RMS")

    #Write texts to display sigmas
    sigmatextx = r"$\sigma_{X}$:"
    sigmatextx +="\nHistogram = "+"{:.2f}".format(sigmax)+"mm"
    sigmatextx +="\nTwiss RMS = "+"{:.2f}".format(Mhx[0])+"mm"
    import re
    if re.search(",",fitlabel): #Not sure what this is for...
        label = re.sub(".+(?<=(,))","",fitlabel)
        label = label.replace(" ","")
    else:
        label = fitlabel
    sigmatextx +="\n"+label+" RMS = "+"{:.2f}".format(Mfx[0])+"mm"
    sigmatextx +="\n\n{:.3f}% outside".format(pOut3sigx)+r" 3$\sigma$"

    sigmatexty = r"$\sigma_{Y}$:"
    sigmatexty +="\nHistogram = "+"{:.2f}".format(sigmay)+"mm"
    sigmatexty +="\nTwiss RMS = "+"{:.2f}".format(Mhy[0])+"mm"
    sigmatexty +="\n"+label+" RMS = "+"{:.2f}".format(Mfy[0])+"mm"
    sigmatexty +="\n\n{:.3f}% outside".format(pOut3sigy)+r" 3$\sigma$"

    #Set various plot variables, often found through trial and error
    xlim = 6*sigmax
    plt.setp(s1,title="X Distribution At Target after "+fitlabel,xlabel="X Position [mm]",
            ylabel="Probability Density",xlim=([-xlim,xlim]),ylim=([1e-6,1])) #if don't set ylim, log goes to e-58
    #s1.set_title("X Distribution At Target after "+fitlabel,fontsize=fs) #+"\n"+rf"$\sigma_D=${{:.1f}}mm".format(sigmax)
    #s1.set_xlabel("X Position [mm]",fontsize=fs)
    #s1.set_ylabel("Probability Density",fontsize=fs)
    #s1.set_xlim([-xlim,xlim])
    ##s1.set_ylim([1e-6,0.1])
    #s1.set_ylim([1e-6,1])  #if don't set, then log goes to e-58
    
    xlim1 = s1.get_xlim()
    ylim1 = s1.get_ylim()
    s1.text(xlim1[0]*0.97,3e-2,sigmatextx,fontsize=fs-2)
    s1.text(xlim1[0]*0.97,1.2e-2,"Beam Twiss at PBW:",fontsize=fs-4)
    s1.text(xlim1[0]*0.97,7e-3,r"$\epsilon_{Nx}$ = "+"{:.3f} [mm*mrad]".format(PBWTwx[2]),fontsize=fs-4)
    s1.text(xlim1[0]*0.97,4e-3,r"$\beta_{x}$ = "+"{:.1f} [m]".format(PBWTwx[0]),fontsize=fs-4)
    s1.text(xlim1[0]*0.97,2.4e-3,r"$\alpha_{x}$ = "+"{:.2f}".format(PBWTwx[1]),fontsize=fs-4)

    handles1, labels1 = plt.gca().get_legend_handles_labels()
    order1=[0,1,2]
    s1.legend([handles1[idx] for idx in order1],[labels1[idx] for idx in order1],fontsize=9)

    #Set s2
    plt.setp(s2,title="Y Distribution At Target after "+fitlabel,xlabel="Y Position [mm]",
            ylabel="Probability Density",xlim=([-xlim,xlim]),ylim=([1e-6,1])) #if don't set ylim, log goes to e-58
    #s2.set_xlim([-xlim,xlim])
    #s2.set_ylim([1e-6,1]) #if don't set, then log goes to e-58
    #s2.set_title("Y Distribution At Target after "+fitlabel,fontsize=fs)
    #s2.set_xlabel("Y Position [mm]",fontsize=fs)
    #s2.set_ylabel("Probability Density",fontsize=fs)

    handles2, labels2 = plt.gca().get_legend_handles_labels()
    order2=[0,1,2]
    s2.legend([handles2[idx] for idx in order2],[labels2[idx] for idx in order2],fontsize=9)
    ylim2=s2.get_ylim() #dynamically get the graph limits
    xlim2=s2.get_xlim()

    s2.text(xlim2[0]*0.97,3e-2,sigmatexty,fontsize=fs-2)
    #PBWTwx = [Ibetax,Ialphx,Inemtx]
    s2.text(xlim2[0]*0.97,1.2e-2,"Beam Twiss at PBW:",fontsize=fs-4)
    s2.text(xlim2[0]*0.97,7e-3,r"$\epsilon_{Ny}$ = "+"{:.3f} [mm*mrad]".format(PBWTwy[2]),fontsize=fs-4)
    s2.text(xlim2[0]*0.97,4e-3,r"$\beta_{y}$ = "+"{:.1f} [m]".format(PBWTwy[0]),fontsize=fs-4)
    s2.text(xlim2[0]*0.97,2.4e-3,r"$\alpha_{y}$ = "+"{:.2f}".format(PBWTwy[1]),fontsize=fs-4)

    #Can date stamp the multi plot for easier tracking of changes, if necessary
    from datetime import datetime
    dt = datetime.now()

    name = savename+"_compTargs_"+dt.strftime("%H-%M-%S")+".pdf"##
    #plt.show()
    if args.savePics:
        plt.tight_layout()
        plt.savefig(name,bbox_inches="tight")
    plt.close()
    print(name)














def printParticles(savename,xinit,pxinit,yinit,pyinit,Einit):
    import csv
    from datetime import datetime
    dt = datetime.now()
    z = -5

    fname = savename+'.csv'
    with open(fname,mode = 'w') as part_file:
        part_writer = csv.writer(part_file,delimiter = ',')
        for i in range(len(Einit)):
            part_writer.writerow(["proton",xinit[i],pxinit[i],yinit[i],pyinit[i],z,Einit[i]])
    part_file.close()
    print(fname)














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
    plt.savefig(name+".pdf",bbox_inches="tight")
    plt.close()
    print(name)













def numLines(filename):
    import csv
    with open(filename+".csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        N = 0
        for row in csv_reader:
            N += 1
        csv_file.close()
    #from datetime import datetime
    #dt = datetime.now()
    #print("There are ",N,"lines "+dt.strftime("%H-%M-%S"))
    return N














def gaussianFit(hist,axis,width,maxim,options,name,y1,y2,saveFits):
    #from Kyrre's doubleGaussian.ipynb example
    import ROOT

    #Project center slice 
    if   axis == "y" or axis == "Y": 
        proj = hist.ProjectionY(axis,hist.GetYaxis().FindBin(-width),hist.GetYaxis().FindBin(width))
    elif axis == "x" or axis == "X":
        proj = hist.ProjectionX(axis,hist.GetXaxis().FindBin(-width),hist.GetXaxis().FindBin(width))
    differenceNG = 100
    total = 100
    differencePG = 100
    coeffsG = [100]

    #Define a gaussian function and fit it to the projection
    f1 = ROOT.TF1('f1','gaus',-maxim,maxim)
    f1_res = proj.Fit(f1, 'RSQ')
    #print(f1_res)
    p0 = f1.GetParameter(0)
    p1 = f1.GetParameter(1)
    p2 = f1.GetParameter(2)
    print("initial Gaussian fit parameters:",p0,p1,p2)

    #Define function of a sum of multiple Gaussians to fit to projection
    r=0.1
    f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3])) + [4] * exp(-x*x/(2*[5]*[5])) + [6] * exp(-x*x/(2*[7]*[7]))',-maxim,maxim)
    #f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3])) + [4] * exp(-x*x/(2*[5]*[5]))',-maxim,maxim)
    #f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1])) + [2] * exp(-x*x/(2*[3]*[3]))',-maxim,maxim)
    #f2 = ROOT.TF1('f2','[0] * exp(-x*x/(2*[1]*[1]))',-maxim,maxim)

    #constrain parameters, trial and error for Nb=500, RM Amplitudes=0
    if axis == "y" or axis == "Y":
        f2.SetParameters(p0*(1-r),p2,p0*r*r,p2*y1,p2*y1,p2*5,p2,p2*y2)
        #print(f2.GetParameter(0),f2.GetParameter(1),f2.GetParameter(2),f2.GetParameter(3),f2.GetParameter(4),f2.GetParameter(5),f2.GetParameter(6),f2.GetParameter(7))
        f2.SetParLimits(0,p0*0.5,p0*y1) #p0=8e4,7e5
        f2.SetParLimits(1,p2*0.25,p0*y1) #
        f2.SetParLimits(2,p1, p0*r*2) #
        f2.SetParLimits(3,p2, p2*y1) #p3=11,22
        f2.SetParLimits(4,p2, p2*y2) #p5=22,300
        f2.SetParLimits(5,p2, p2*y2) #p5=22,300
        #f2.SetParLimits(6,p2,p2*y2) #p5=22,300
        f2.SetParLimits(7,p2, p2*y2*2) #p7 =300,3000
    elif axis == "x" or axis == "X":
        f2.SetParameters(p0*(1-r),p2,p0*r,p2*2,p2*y1,p2*y1,p2,p2*y2)
        #print(9532,16.3,395,36.5,6511,14,2,334)
        #print(f2.GetParameter(0),f2.GetParameter(1),f2.GetParameter(2),f2.GetParameter(3),f2.GetParameter(4),f2.GetParameter(5),f2.GetParameter(6),f2.GetParameter(7))
        f2.SetParLimits(0,p0*r,p0*y1)
        f2.SetParLimits(1,p2*0.25,p2*y1)
        f2.SetParLimits(2,p2, p0*r*2)
        f2.SetParLimits(3,p2, p2*y1)
        f2.SetParLimits(4,p2, p2*y2)
        f2.SetParLimits(5,p2, p2*y2)
        #f2.SetParLimits(6,0.1,p0)
        f2.SetParLimits(7,p2*y2*0.5,p2*y2*2)
    f2_res = proj.Fit(f2, 'RSQ')
    print("Gaussian",f2.GetParameter(0),f2.GetParameter(1),f2.GetParameter(2),f2.GetParameter(3),f2.GetParameter(4),f2.GetParameter(5),f2.GetParameter(6),f2.GetParameter(7))
    #print(f2_res)

    if saveFits: #if uncommented it opens a canvas even if false...?
        c1 = ROOT.TCanvas()
        proj.Draw()
        f2.SetLineColor(ROOT.kRed)
        f2.Draw('same')
        #f1.SetLineColor(ROOT.kRed)
        #f1.Draw('same')
        if axis == "y" or axis == "Y": proj.GetXaxis().SetRangeUser(-maxim,maxim)
        elif axis =="x" or axis == "X": proj.GetXaxis().SetRangeUser(-maxim,maxim)
        c1.SetLogy()
        c1.Print(name+"_"+axis+"GaussFit"+str(maxim)+".pdf")

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

    differenceNG = proj.Integral(proj.GetXaxis().FindBin(-maxim),proj.GetXaxis().FindBin(maxim),'width')  -  f2.Integral(-maxim,maxim)
    total = proj.Integral(proj.GetXaxis().FindBin(-maxim),proj.GetXaxis().FindBin(maxim),'width') #need better name
    differencePG = differenceNG/total*100
    coeffsG = [f2.GetParameter(0),f2.GetParameter(1),f2.GetParameter(2),f2.GetParameter(3),f2.GetParameter(4),f2.GetParameter(5),f2.GetParameter(6),f2.GetParameter(7)]
    print(axis," difference: {:.0f},\t total: {:.0f},\t{:.3f}%".format(differenceNG,total,differencePG))#,"\n")#,coeffs)
    
    
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
    f3 = ROOT.TF1('f3L','[0] * [1]/(x**2 + [1]**2) + [2] * exp(-x*x/(2*[3]*[3]))',-maxim,maxim)
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
    f3_res = proj.Fit(f3, 'RSQ')

    differenceNL = proj.Integral(proj.GetXaxis().FindBin(-maxim),proj.GetXaxis().FindBin(maxim),'width')  -  f3.Integral(-maxim,maxim)
    total = proj.Integral(proj.GetXaxis().FindBin(-maxim),proj.GetXaxis().FindBin(maxim),'width') #need better name
    differencePL = differenceNL/total*100
    coeffsL = [f3.GetParameter(0),f3.GetParameter(1),f3.GetParameter(2),f3.GetParameter(3)]
    print("Voigt",coeffsL)
    print(axis," L difference: {:.0f},\t total: {:.0f},\t{:.3f}%".format(differenceNL,total,differencePL))#,"\n")#,coeffs)

    if saveFits: #if uncommented it opens a canvas even if false...?
        c2 = ROOT.TCanvas()
        proj.Draw()
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
    um = 1e-6
    from datetime import datetime
    #from plotFit import converter

    #print(datetime.now().strftime("%H-%M-%S"))
    (Img, xax, yax) = converter(histogram2D,args.saveHist,name,paths,args.reBin) #convert from TH2D to numpy map
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
    sumCharge = np.sum(ImgJ)
    widths = np.zeros(len(boxes))
    heights = np.zeros(len(boxes))
    pLargers = np.zeros(len(boxes)) #percent larger than initial box(?)
    #coreImax = np.zeros(len(boxes))
    coreMeanI = np.zeros(len(boxes))
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
    jMax,pOutsideBox,rValue,edges,EI,beamArea,centX,centY,rDiff = PMAS(ImgJ,args,parts,xax,yax,name,paths) #,dispY,dispX
    fPMAS = datetime.now()
    #print(EI)#"finished PMAS in",datetime.now()-sPMAS)#.strftime("%H-%M-%S"))

    #Warning print is outside PMAS because it takes a long time
    jMaxLim = 53
    pOutLim = 4
    centYLim = 11
    centXLim = 14
    rValLim = 0.20
    nomArea = 5600
    #if jMax >= jMaxLim or pOut >= pOutLim or rValue >= rValLim or dispY >= dispYLim or dispX >= dispXLim:
    if jMax >= jMaxLim:
        print("Current Density Warning: {:0.1f}".format(jMax),"uA/cm^2!")
    elif pOutsideBox >= pOutLim:
        print("% Outside Box Warning:",pOutsideBox,"%!")
    elif beamArea <= nomArea:
        print("Beam Area Warning:",beamArea,"mm^2")
    elif centY >= centYLim:
        print("Vertical Displacement Warning:",centY - centYLim,"mm!")
    elif centY <= -centYLim:
        print("Vertical Displacement Warning:",centY + centYLim,"mm!")
    elif centX >= centXLim:
        print("Vertical Displacement Warning:",centX - centXLim,"mm!")
    elif centX <= -centXLim:
        print("Horizontal Displacement Warning:",centX + centXLim,"mm!")
    elif rValue >= rValLim:
        print("R Value Warning",rValue)
    else:
        print("No Warnings!")

    #print(jMax,pOutsideBox,rValue,rDiff,beamArea,dispX,dispY)
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
    Protmax = Img.max() #protons/mm^2
    #for i in range(len(Img)):
    #    for j in range(len(Img[i])):
    #        Img[i][j] = Img[i][j] * C #[uA/cm^2]
    #jMax = Img.max()
    jMin = 0.9 * C * parts/3e5 #background (0 hits) will be un-colored, scaled to 1e5 parts
    print(sample,"PMAS:",fPMAS-sPMAS,"s: Jmax {:.1f}, J in core {:.0f}mm^2: {:.2f}uA/cm^2, nominal A: {:.0f}mm^2, RVal {:.3f}, Jmin {:.1f}, % Out TA {:.2f}"\
                .format(jMax,coreArea,coreJMean,128*42,rValue,jMin,pOutsideBoxes[0]))
    print("Beam Top: {:.1f}mm, Bottom: {:.1f}mm, Left: {:.1f}mm, Right: {:.1f}mm".format(edges[0],edges[1],edges[2],edges[3]))

    if args.gaussFit:
                                           #(hist,     axis, width,maxim,options,name,   y1,y2,saveFits)
        diffNy,diffPy,coeffsy, differenceNLy,differencePLy,coeffsLy = gaussianFit(histogram2D,"y",yBinSize,500,options,savename,2,30,args.saveFits)
        diffNx,diffPx,coeffsx, differenceNLx,differencePLx,coeffsLx = gaussianFit(histogram2D,"x",xBinSize,500,options,savename,3,20,args.saveFits)
        #add minimize function for these?

    if args.savePics:
        from matplotlib.pyplot import subplots,pcolormesh,close,tight_layout,savefig,setp
        from matplotlib.patches import Rectangle
        from matplotlib.colors import LogNorm
        X, Y = np.meshgrid(xax,yax) #Set mesh of plot from the axes given from converter function
        close() #make sure no prior plotting messes up

        fig,ax = subplots(dpi=150,figsize=(6.4,4.8),tight_layout=True)
        #print(datetime.now().strftime("%H-%M-%S"))
    
        #Set maximum value depending on the maximum current density
        from math import log10,ceil
        maxim = ceil(jMax / 10) * 10
        if args.maxim != 0: maxim = args.maxim #user provided maximum
        minim = 10**ceil(log10(jMin))
        cbarVals  = [minim,minim*10,minim*100,minim*0.455,maxim] #make array for color bar values
        cbarLabels = ["{:.2f}".format(cbarVals[0]),"{:.1f}".format(cbarVals[1]),"{:.1f}".format(cbarVals[2]),
                    "{:.0f}".format(cbarVals[3]),"{:.0f}".format(cbarVals[4])]#,"{:.1f}".format(cbarVals[5])] #make labels of Value
        cbarLabel = r"$\mu A / cm^2$"
        #print("Max current Density: ",Img.max(),"/",maxim,datetime.now().strftime("%H-%M-%S"))

        #Use pcolor to show density map, use log scale
        c = ax.pcolormesh(X,Y,Img,shading='auto',norm=LogNorm(vmin=minim, vmax=maxim), cmap='viridis') #viridis or magma are perceptually uniform
        lw=1
        col='k'
        fs=14

        #Set Plot Properties
        setp(ax,xlim=([-args.xlim,args.xlim]),ylim=([-args.ylim,args.ylim]),
                title="Proton Beam Distribution at "+position,ylabel="Vertical [mm]",xlabel="Horizontal [mm]")
        cbar = fig.colorbar(c, ax=ax,pad=0.01)
        cbar.set_label(cbarLabel,labelpad=3,fontsize=fs-2)
        cbar.set_ticks(cbarVals)
        cbar.set_ticklabels(cbarLabels)

        if position =="Target":
            #name = name+"_2212only"
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

        #Show 99% box
        if not args.noBox: #for user clarity, call is noBox, include % outside text for this
            ax.add_patch(Rectangle((boxLLxs[0],boxLLys[0]),widths[0],heights[0],linewidth=lw,edgecolor=cols[0],fill=False))
            ax.text(xlim[0]*0.90, ylim[1]*0.85, "{:.2f}".format(pOutsideBox)+r"% Outside 160x64mm$^2$ Target", color=cols[0], fontsize = fs-2, fontweight='bold', backgroundcolor = 'w')
        else: #yes no Boxes
            name+="_noBox"

        if not args.noText: #for user clarity, call is noText
            ax.set_title("Distribution at "+position+"\n{:.3f}% of {:.2e} total protons".format(Pprotons,parts),fontsize=fs)

            #Display beam characteristics
            bgdbox=dict(pad=1,fc='w',ec='none')
            propsR = dict(horizontalalignment="right",verticalalignment="bottom", backgroundcolor = 'w',bbox=bgdbox)
            #do this better, since set by args.x,ylim 
            ax.text(xlim[0]*0.97, ylim[0]*0.60, "Beam Twiss at PBW:", fontsize=fs-4, backgroundcolor = 'w',bbox=bgdbox)
            ax.text(xlim[0]*0.97, ylim[0]*0.71, r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[mm \cdot mrad]}$", fontsize=fs-4, backgroundcolor = 'w',bbox=bgdbox)
            ax.text(xlim[0]*0.97, ylim[0]*0.83, r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$", fontsize=fs-4, backgroundcolor = 'w',bbox=bgdbox)
            ax.text(xlim[0]*0.97, ylim[0]*0.95, r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]), fontsize=fs-4, backgroundcolor = 'w',bbox=bgdbox)
            ax.text(xlim[1]*0.97, ylim[0]*0.60, r"Core <$\bf{J}$>: "+"{:.1f}".format(coreJMean)+r" $\mu$A/cm$^2$", propsR,fontsize=fs-4)
            #ax.text(xlim[1]*0.97, ylim[0]*0.60, "R={:.4f}".format(rValue), propsR,fontsize=fs-2)
            ax.text(xlim[1]*0.97, ylim[0]*0.75, r"Max $\bf{J}$: "+"{:.1f}".format(jMax)+r" $\mu$A/cm$^2$", propsR,fontsize=fs-4)
            ax.text(xlim[1]*0.97, ylim[0]*0.85, "RM Amplitudes: {:.1f}, {:.1f}mm".format(args.aX,args.aY),propsR,fontsize=fs-4)
            ax.text(xlim[1]*0.97, ylim[0]*0.96, options['physList'], propsR,fontsize=fs-4)

            for i in range(1,len(boxes)): #make multiple boxes
                ax.add_patch(Rectangle((boxLLxs[i],boxLLys[i]),widths[i],heights[i],linewidth=lw,edgecolor=cols[i],fill=False))
                ax.text(xlim[0]*0.90, ylim[1]*(0.85-i*0.1), "{:.2f}".format(pOutsideBoxes[i])+"% Outside {:.0f}% Larger Box".format(pLargers[i]*100), 
                                  color=cols[i], fontweight='bold',fontsize = fs-2, backgroundcolor = 'w',bbox=dict(pad=1))#,path_effects=[path_effects.withStroke(linewidth=1, foreground='k')])
        else: #yes no text
            name+="_noText"


        if args.saveEdges:
            edgeCol = 'k'
            ax.hlines(edges[0],edges[2],edges[3],colors=edgeCol,linewidths=1)
            ax.hlines(edges[1],edges[2],edges[3],colors=edgeCol,linewidths=1)
            ax.vlines(edges[2],edges[1],edges[0],colors=edgeCol,linewidths=1)
            ax.vlines(edges[3],edges[1],edges[0],colors=edgeCol,linewidths=1)
            #print("Edges printed!!")
            name+="_Edges"

        dt = datetime.now()
        #from os.path import isfile
        #if isfile(name+"*.png"):
        #  print("already present")
        tight_layout()
        savefig(name+"_"+dt.strftime("%H-%M-%S")+"."+args.picFormat,bbox_inches='tight')
        close(fig)
        close()
        print(name+"_"+dt.strftime("%H-%M-%S")+"."+args.picFormat)
    #dt = datetime.now()
    #print(dt-start)

    #passing lists to decrease length of argument calls
    return [jMax,pOutsideBox,beamArea,coreJMean,centX,centY,rValue,rDiff]













def converter(hIn,saveHist,name,paths,reBin):
    import ROOT
    from numpy import zeros

    if reBin > 1:
        #print("Rebinning")
        oldName=hIn.GetName()
        hIn = hIn.Rebin2D(reBin,reBin,oldName+"_"+str(reBin))
        #print("Shape of old:",hIn.GetXaxis().GetNbins(),"\nShape of New:",hIn.GetXaxis().GetNbins())

    #Get X axis from ROOT
    xax = zeros(hIn.GetXaxis().GetNbins()+1)
    for i in range(len(xax)-1):
        xax[i] = hIn.GetXaxis().GetBinLowEdge(i+1) #Set elements as bin edges
    xax[-1] = hIn.GetXaxis().GetBinUpEdge(hIn.GetXaxis().GetNbins())

    #Get Y axis from ROOT
    yax = zeros(hIn.GetYaxis().GetNbins()+1)
    for i in range(len(yax)-1):
        yax[i] = hIn.GetYaxis().GetBinLowEdge(i+1) #Set elements as bin edges
    yax[-1] = hIn.GetYaxis().GetBinUpEdge(hIn.GetYaxis().GetNbins())
    
    hOut = zeros( (len(yax), len(xax)) ) #2D Image map

    #Fill Image map with 2D histogram values
    for xi in range(hIn.GetXaxis().GetNbins()):
        for yi in range(hIn.GetYaxis().GetNbins()):
            bx = hIn.GetBin(xi+1,yi+1)
            hOut[yi,xi] = hIn.GetBinContent(bx)

    #Must add Overflow options!!!

    if saveHist: #for rValue calculations
        from os import uname
        if uname()[1] in {"tensor.uio.no","heplab01.uio.no", "heplab04.uio.no","heplab03.uio.no"}:
            import csv,re
            #Remove upper directories that may have come with name for appending outname to scratch folder
            if re.search("/PBW_",name):
                #print("\n",outname,"\n")
                name = re.sub(".+(?=(PBW_))","",name) #substitutes "" for all preceeding PBW_
                #print("Histogram-removed",name)
            if re.search("/Vac_",name):
                #print("\n",outname,"\n")
                name = re.sub(".+(?=(Vac_))","",name) #substitutes "" for all preceeding Vac_
              #print("Histogram-removed",name)

        with open(paths['scratchPath']+name+".csv",mode = 'w',newline=None) as hist_file:
            hist_writer = csv.writer(hist_file,delimiter = ',')
            hist_writer.writerows(hOut)
        hist_file.close()
        print(paths['scratchPath']+name+".csv")

    return (hOut,xax,yax)















#add function to check if ROOT file is complete, else delete and make new
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

    #need to add completeness check
    return pwd















def rCompare(Im,Nb,paths,reBin):
    import numpy as np
    from os import uname

    #Find reference files
    if uname()[1] in {"tensor.uio.no","heplab01.uio.no","heplab04.uio.no","heplab03.uio.no"}:
        if Nb == 100: #do fill in from args
            if reBin == 5:
                Iref = np.genfromtxt(open(paths['scratchPath']+"PBW_570MeV_beta1006.80,129.72m_RMamp55,18mm_N2.9e+06_NpB100_runW_QBZ_TargetImage_rB5.csv"),delimiter=",")
            elif reBin == 4:
                Iref = np.genfromtxt(open(paths['scratchPath']+"PBW_570MeV_beta1006.80,129.72m_RMamp55,18mm_N2.9e+06_NpB100_runW_QBZ_TargetImage_rB4.csv"),delimiter=",")
            elif reBin == 1 or reBin == 2:
                Iref = np.genfromtxt(open(paths['scratchPath']+"PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+06_NpB100_NPls1e+03_runW_QBZ_TargetImage.csv"),delimiter=",")
        elif Nb == 500:
              Iref = np.genfromtxt(open(paths['scratchPath']+"Vac_570MeV_beta1007,130m_RMamp55,18mm_N1.4e+07_NpB500_NPls1e+03_run_QBZ_TargetImage.csv"),delimiter=",")
        else: #Nb=10
            Iref = np.genfromtxt(open(paths['scratchPath']+"PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03_4_runW_QBZ_TargetImage.csv"),delimiter=",")
            #Iref = np.genfromtxt(open(paths['scratchPath']+"Vac_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03_run_QBZ_TargetImage.csv"),delimiter=",")
            #Iref = np.genfromtxt(open("/scratch2/ericdf/PBWScatter/PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03_runW_QBZ_TargetImage.csv"),delimiter=",")
        
    elif uname()[1] == "mbef-xps-13-9300":
        Iref = np.genfromtxt(open("/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/Vac_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03_run_QBZ_TargetImage.csv"),delimiter=",")

    lenx = np.shape(Im)[0]
    leny = np.shape(Im)[1]
    diff = np.zeros((leny,lenx))
    div = np.zeros((leny,lenx))
    offset = 0 #1e-10

    for i in range(lenx):
        for j in range(leny):
            #Iref[j,i] += offset
            if Iref[j,i] == 0: continue #not great, but it produces better values 15.1.23
            div[j,i] = ( ( ( (Im[j,i] + offset) / (Iref[j,i] + offset) - 1) ** 2 ) / (leny * lenx) )
    rDiv = np.sqrt(np.sum(div))

    #for i in range(lenx):
    #    for j in range(leny):
    #        #Iref[j,i] += offset
    #        if Iref[j,i] == 0: continue #not great, but it produces better values 15.1.23
    #        diff[j,i] = ( ( ( (Im[j,i] + offset) - (Iref[j,i] + offset) ) ) / (leny * lenx) ) #
    rDiff = 0#np.sqrt(np.sum(diff))

    return rDiv,rDiff


















def findEdges(Img,jMax,graph,savename,xax,yax,args):
    from numpy import array
    from scipy.signal import convolve2d as scipySigConv2d
    #does it need to be normalized?

    gradXX = array([[1,0,-1],[2,0,-2],[1,0,-1]])
    gradYY = array([[1,2,1],[0,0,0],[-1,-2,-1]])

    gradX = scipySigConv2d(Img,gradXX,'same')
    gradY = scipySigConv2d(Img,gradYY,'same')

    #Prominence of gradient to determine if actual edge or not
    promX = 3
    promY = 4

    EI,edges = localExtrema(gradX,promX,gradY,promY,xax,yax)
    #print(edges)

    #Only do check for Top and Left?
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
        plt.savefig(savename+str(promY)+"_GradPics."+args.picFormat,bbox_inches='tight')
        print(savename,promY,"_GradPics.",args.picFormat,sep="")
        plt.close()

    return EI,edges
















def localExtrema(gradX,nx,gradY,n,xax,yax):
    #Finds Edge Indices from Gradient maps
    from numpy import sum as npSum,argpartition,max as npMax,min as npMin,ndenumerate#,ndenumerate,argmax,argmin,amax,amin,
    from math import ceil

    ##Set conspicuous initial values in case not found
    #avgMaxSumGradX = 0
    #avgMaxSumGradY = 0
    #avgMinSumGradX = 0
    #avgMinSumGradY = 0
    
    sumGradX = npSum(gradX,axis=0)
    sumGradY = npSum(gradY,axis=1)
    #print(len(sumGradY),npMin(sumGradY))
    a=1 #decrease core size

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
    #print([topInd,botInd,lefInd,rigInd])
    ##n=6 #this gave pretty even results
    ##argpartition implemented like this selects the maximum n values
    ##  argpartition sorts and puts the passed k in the k-th index. 
    ##    I call the final n indices, meaning the n greater than k, 
    #maxSumGradXInds = argpartition(sumGradX,-n)[-n:]
    #minSumGradXInds = argpartition(sumGradX,n)[n:]
    #maxSumGradYInds = argpartition(sumGradY,-n)[-n:]
    #minSumGradYInds = argpartition(sumGradY,n)[n:]
    #print("len",len(minSumGradYInds),minSumGradYInds[-n:])
    ##then those max are summed 
    #for i in range(n):
    #    avgMaxSumGradX += maxSumGradXInds[i]
    #    avgMaxSumGradY += maxSumGradYInds[i]
    #    avgMinSumGradX += minSumGradXInds[i]
    #    avgMinSumGradY += minSumGradYInds[i]
    #    print("ind",minSumGradYInds[i])
    ##print("gradCalcs",round(avgMaxSumGradY/n),round(avgMinSumGradY/n),round(avgMaxSumGradX/n),round(avgMinSumGradX/n))

    ##and averaged by the number of points. Ideally this would find the peak of the gradient
    #topInd = ceil(avgMinSumGradY/n) #is ceil necessary?
    #botInd = round(avgMaxSumGradY/n)
    #lefInd = round(avgMaxSumGradX/n)
    #rigInd = ceil(avgMinSumGradX/n)

    lim=550
    if topInd > lim:
        topInd = 0
        #print(avgMinSumGradY,avgMinSumGradY/n,avgMaxSumGradY)
        print("\n\nError in Gradient calculation of top Index!\n\n")
    elif botInd > lim:
        botInd = 0
        #print(maxSumGradYInds,avgMaxSumGradY,avgMaxSumGradY/n)
        print("\n\nError in Gradient calculation of bottom Index!\n\n")
    elif lefInd > lim:
        lefInd = 0
        print("\n\nError in Gradient calculation of left Index!\n\n")
    elif rigInd > lim: 
        rigInd = 0
        print("\n\nError in Gradient calculation of right Index!\n\n")

    #check the sum of the rows?
    def indCheck(ind,arr):
        if arr[ind+1] < arr[ind]*0.5:
            return False
        elif arr[ind-1] < arr[ind]*0.5:
            return False
        else:
            return True
    #Add a check in localExtrema that the gradient of the row/column on either side of the max is at least half(?) of the max. 
        #Could also check if the index is high (<100,>900) and tell it to recalculate.
    #ideally this would then delete the row for this reconing if False is returned and find the next highest...

    lMaxX = xax[lefInd]
    lMinX = xax[rigInd] 
    lMaxY = yax[botInd]
    lMinY = yax[topInd]

    return [topInd,botInd,lefInd,rigInd], [lMinY,lMaxY,lMaxX,lMinX] 













#Detection Algorithm!!!
def PMAS(Img,args,parts,xax,yax,name,paths):
    import numpy as np
    printEdges = True

    xBinSize = xax[round(len(xax)*0.5)+1]-xax[round(len(xax)*0.5)]
    yBinSize = yax[round(len(xax)*0.5)+1]-yax[round(len(xax)*0.5)]

    #99% box
    width = 160
    height = 64
    boxLLx = round(-width/2 / 2) * 2 #-90
    boxLLy = round(-height/2 / 2) * 2 #-80
    #print(len(xax),boxLLx,boxLLy,width,height)
    rdVal=args.reBin*2
    for idx, val in np.ndenumerate(xax):
        if int(val) == round(boxLLx/rdVal)*rdVal:
            boxLInd = idx[0] #bc it is rounding
        if int(val) == round(boxLLx/rdVal)*rdVal+round(width/rdVal)*rdVal:
            boxRInd = idx[0] #545
    for idx, val in np.ndenumerate(yax):
        if int(val) == round(boxLLy/rdVal)*rdVal:
            boxBInd = idx[0] #bc it is rounding
        if int(val) == round(boxLLy/rdVal)*rdVal+round(height/rdVal)*rdVal:
            boxTInd = idx[0] #540
    #Find % outside the 99% Box area
    core = Img[boxBInd:boxTInd,boxLInd:boxRInd]
    sumTot = np.sum(Img)+1
    sumCore = np.sum(core)+1
    pOut = (sumTot-sumCore)/sumTot*100
    #print(pOut)

    #Normalize to full current, see 
    I_pulse = 62.5*1e3 #[uA]
    #ImgJ = np.zeros_like(Img) #for redefinition
    C = I_pulse / parts / (xBinSize * yBinSize * 1e-2) * 0.04 #[uA/cm^2]: /mm^2->/cm^2 = /1e-2, A->uA = 1e6, 4% duty cycle
    for i in range(len(Img)):
        for j in range(len(Img[i])):
            Img[i][j] = Img[i][j] * C #[uA/cm^2] #redefinition was messing with proton sum, but not necessary for final algorithm
    jMax = Img.max()
    #sumCharge = np.sum(Img)+1
    #print("PMAS Sum",sumTot)

    #R value for algorithm. Works when use Current Density, not Nprotons
    rValue=0;rDiff=0
    rValue,rDiff = rCompare(Img,args.Nb,paths,args.reBin)
    #print("R = ",rValue,rDiff)
    #print("Converted in",datetime.now() - start)

#    #Doesn't work, 12 Feb
#    from cv2 import Sobel,COLOR_BGR2GRAY
#    nomEdges = [20,-22,-66,68]
#    edgeImgX = Sobel(Img,-1,1,0)
#    edgeImgY = Sobel(Img,-1,0,1)
#    printEdges=True
#    if printEdges:
#        import csv
#        G = np.sqrt(edgeImgX**2 + edgeImgY**2)
#        with open(name+"_edges.csv",mode = 'w',newline=None) as edge_file:
#            edge_writer = csv.writer(edge_file,delimiter = ',')
#            edge_writer.writerows(G)
#        edge_file.close()
#        print(name+"_edges.csv")
#    for idx, val in np.ndenumerate(edgeImgX):
#        if val == 1:
#            lefInd = idx
#        if val == -1:
#            rigInd = idx
#    for idx, val in np.ndenumerate(edgeImgY):
#        if val == 1:
#            topInd = idx
#        if val == -1:
#            botInd = idx
#    EI = [topInd,botInd,lefInd,rigInd]
#    lEdge = xax[lefInd]
#    rEdge = xax[rigInd] 
#    bEdge = yax[botInd]
#    tEdge = yax[topInd]
#    dispT = np.abs(tEdge-nomEdges[0])
#    dispB = np.abs(bEdge-nomEdges[1])
#    dispL = np.abs(lEdge-nomEdges[2])
#    dispR = np.abs(rEdge-nomEdges[3])
#    print("Beam moved ",dispT,",",dispB,"vertical and",dispL,",",dispR,"horizontally")
#    dispY = (dispT+dispB)/2
#    dispX = (dispL+dispR)/2

    #Edges
    #Not great, may need to reshape and average. Not sure how to do that...
    #from plotFit import findEdges
    EI,edges = findEdges(Img,jMax,args.saveGrads,name,xax,yax,args) #[topInd,botInd,lefInd,rigInd]
    nomEdges = [20,-22,-66,68]
    nomArea = (nomEdges[0]-nomEdges[1])*(nomEdges[3]-nomEdges[2])
    dispT = np.abs(edges[0]-nomEdges[0])
    dispB = np.abs(edges[1]-nomEdges[1])
    dispL = np.abs(edges[2]-nomEdges[2])
    dispR = np.abs(edges[3]-nomEdges[3])
    dispY = (dispT+dispB)/2
    dispX = (dispL+dispR)/2
    beamArea = (np.abs(edges[0])+np.abs(edges[1]))*(np.abs(edges[3])+np.abs(edges[2]))
    #Need to displat beam Center index in X,Y, not displacement 
    centX = edges[2] + edges[3]
    centY = edges[0] + edges[1]
    #print("Beam Center at ",centX,",",centY)

    return jMax,pOut,rValue,edges,EI,beamArea,centX,centY,rDiff #make this a list!












def saveStats(statsPWD,Twiss,rasterBeamFile,jMax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue,rDiff,fileName,reBin):
    #save values in Stats CSV
    import csv
    from os import path as osPath
    from re import sub
    rasterBeamFile = sub(".+(?<=(/CSVs/))","",rasterBeamFile)
    if fileName != "":  
            statsName = fileName+".csv"
    else:
        statsName = statsPWD+"EvalStats12MarHEBTA2T.csv" #don't forget to include statsPWD!!

    #print("Is file present?",osPath.isfile(statsPWD+statsName))
    writeMore = False
    if not osPath.isfile(statsName): #does this work?
        with open(statsName,mode = 'w') as csv_file:
            csv_writer = csv.writer(csv_file,delimiter = ',')
            csv_writer.writerow(["BeamFile","Emittance X [mm-mrad]","Beta X [m]","Alpha X","Emittance Y [mm-mrad]","Beta Y [m]",
                            "Alpha Y","Peak Current Desnity [uA/cm2]","% Outside Boxes","Core Area","Core J Mean","Center X","Center Y","R Value"])
            #csv_writer.writerow([rasterBeamFile,Twiss[0],Twiss[1],Twiss[2],Twiss[3],Twiss[4],Twiss[5],jMax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue])
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
                if row[0] == rasterBeamFile: #No duplicates
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
            csv_writer.writerow([rasterBeamFile,Twiss[0],Twiss[1],Twiss[2],Twiss[3],Twiss[4],Twiss[5],jMax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue])
        csv_file.close()
        print("saved",statsName,"{:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.0f}, {:.0f}, {:.3f}".format(jMax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue))#,rDiff)
        










def plotSpread(args,Twiss,statsPWD,paramName,ind,unit,paramLabel,pFitLims,paramBins,pHistLims,origBX,origBY,beamFile):
    import csv,re
    from matplotlib.pyplot import hist,savefig,close,plot,xlim,ylim,text,title,xlabel,ylabel,tight_layout,xticks
    from numpy import greater,zeros,mean,std
    #from plotFit import findFit, gaussian
    from math import floor,log10

    read = zeros(args.samples)
    read.fill(-750) #allow better filter than nonzero?
    #print(Twiss,beamFile)
    i=0
    lenbreak=False
    #How to select which entries? Make keys to search for
    pctKey = "sampleIn"+"{:.0f}Pct".format(args.betaSpread) #"(PBW_570MeV)"
    pbwKey = "PBW_{:.0f}MeV".format(args.energy)
    nBkey = "{:.0f}_NpB{:.0f}".format(floor(log10(2.88e4*args.Nb)),args.Nb)
    betaKey = "beta{:.2f},{:.2f}m".format(Twiss[1],Twiss[4])
    origKey = "OrigbX{:.2f},bY{:.2f}".format(origBX,origBY)
    failKey = "failure{:.0f}-{:.0f}f".format(args.failure,args.magFails)
    #print(nBkey,pctKey,origKey,betaKey)
    #print(re.search(nBkey,beamFile) , re.search(pctKey,beamFile) , re.search(origKey,beamFile), re.search(betaKey,beamFile) ,"\n\n")
    with open(statsPWD+"EvalStats12MarHEBTA2T.csv",mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        #requires differentiating between type of run we're looking at. Can't have both in one line
        if args.betaSpread != 0: #sampleIn files
            for row in csv_reader:
                #if i == 6: #for finding indices out of 7sigma
                #    print("index",i,"failed file:",row[0])
                #print(re.search(nBkey,row[0]) , re.search(pctKey,row[0]) , (re.search(origKey,row[0])))
                if re.search(nBkey,row[0]) and re.search(pctKey,row[0]) and re.search(origKey,row[0]):
                    read[i] = float(row[ind])
                    if i == args.samples-1:
                        lenbreak = True
                        break
                    i+=1
                if lenbreak:
                    break
        else: #PBW files:
            if args.failure != 0: #RM fail
                for row in csv_reader:
                    if re.search(nBkey,row[0]) and re.search(failKey,row[0]) and re.search(pbwKey,row[0]) and re.search(betaKey,row[0]) :
                        read[i] = float(row[ind]) #[mm-mrad]
                        if i == args.samples-1:
                            lenbreak = True
                            break
                        i+=1
                    if lenbreak:
                        break
            else:
                if args.qpNum != "":
                    for row in csv_reader:
                        #if i == 6: #for finding indices out of 7sigma
                        #    print("index",i,"failed file:",row[0])
                        #print(re.search(nBkey,row[0]) , re.search(pctKey,row[0]) , (re.search(origKey,row[0])))
                        if re.search(nBkey,row[0]) and re.search(pctKey,row[0]) and re.search(origKey,row[0]):
                            read[i] = float(row[ind])
                            if i == args.samples-1:
                                lenbreak = True
                                break
                            i+=1
                        if lenbreak:
                            break
                else: #nominal
                    for row in csv_reader:
                        #if i == 44: #for finding indices out of 7sigma
                        #    print("index",i,"failed file:",row[0])
                        #print(re.search(nBkey,row[0]) , re.search(pbwKey,row[0]) , (re.search(betaKey,row[0])))
                        if re.search(nBkey,row[0]) and re.search(pbwKey,row[0]) and re.search(betaKey,row[0]):
                            read[i] = float(row[ind]) #[mm-mrad]
                            if i == args.samples-1:
                                lenbreak = True
                                break
                            i+=1
                        if lenbreak:
                            break
    csv_file.close()

    nonEmpty = greater(read,-750) #remove unused elements
    read = read[nonEmpty]

    #if ind == 13:
    #    print(read)
    pFitLims = (0,[args.samples,100,500])
    print(len(read),"guess",paramName, args.samples/4,mean(read),std(read),paramBins,len(read),pFitLims)
    #sanity check on values
    #print(read)
    for i in range(len(read)):
        if read[i] < mean(read)-std(read)*7 or read[i] > mean(read)+std(read)*7:
            print(read[i],", index",i,"is outside 7 sigma!")

                              #gaussian(x,  amplitude,          mu,      sigma)
    mu, sigma,ampl,interval = findFit(read,[args.samples/4,mean(read),std(read)],pFitLims,paramBins)
    #print(mu,sigma,ampl)

    _, bins, _ = hist(read,interval)#,range=pHistLims)
    plot(bins, gaussian(bins,ampl,mu,sigma), "r--", linewidth=2)

    nPs = round((2*10+round(0.04*1/14*1e6)/1e3)*args.nP)*args.Nb
    if args.qpNum != "":
        #if ind == 7:
            #xlim(pHistLims)
        title("{:.0f} Samples of {:.0f}% Jitter for QP{:.0f}\nwith {:.2e} Protons".format(len(read),args.betaSpread,int(args.qpNum),nPs))
    else:
        if args.betaSpread != 0:
            title("{:.0f} Samples of {:.0f}% Jitter Around Nominal\nwith {:.2e} Protons".format(len(read),args.betaSpread,nPs))
            #xlim(pHistLims)
        else:
            title("{:.0f} Samples with {:.2e} Protons".format(len(read),nPs))
            #xlim(pHistLims)

    xlabel(paramLabel)
    ylabel("Counts")
    #setp(ax,title=Title,xlabel=xLabel,ylabel=yLabel)

    bgdbox=dict(pad=1,fc='w',ec='none')
    fs=14
    if mu > 25:
        delta = 3e-4
    elif mu > 8.5 and mu <= 25:
        delta = 2e-4
    elif mu >2 and  mu <= 8.5:
        delta = 1.2e-4
    else: #mu<1
        delta=1e-6
        xticks(rotation=45)
    #left=xlim()[0] #=7
    #midPoint=left+(xlim()[1]-left)*0.5 #7+(11-7=4/2=2)=9,
    #print(mu,left,midPoint)
    #if mu < midPoint:
    #    xlim(left=mu - (midPoint-left)) # mu=8.3<9->left=8.3-(9-7=2)=6.3 
    propsR = dict(horizontalalignment="right",verticalalignment="top", backgroundcolor = 'w',bbox=bgdbox,fontsize=fs-2)
    text(xlim()[1]*(1-delta), ylim()[1]*0.97, "Beam Twiss at PBW:", propsR)
    text(xlim()[1]*(1-delta), ylim()[1]*0.90, r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[\mu m]}$",propsR)
    text(xlim()[1]*(1-delta), ylim()[1]*0.83, r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$", propsR)
    text(xlim()[1]*(1-delta), ylim()[1]*0.75, r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]), propsR)

    if sigma <1e-2:
        text(xlim()[1]*(1-delta), ylim()[1]*0.65,r"$\mu$="+"{:.3f}".format(mu)+unit+"\n"+r"$\sigma$="+"{:.2e}".format(sigma)+unit, propsR)
    else:
        text(xlim()[1]*(1-delta), ylim()[1]*0.65,r"$\mu$="+"{:.3f}".format(mu)+unit+"\n"+r"$\sigma$="+"{:.3f}".format(sigma)+unit, propsR)
    if args.twissFile == "":
        name=statsPWD+"Nb{:.0f}_{:.0f}x{:.0f}betaSpreadNom".format(args.Nb,len(read),args.betaSpread)+paramName+"Hist"
    else:
        name=statsPWD+args.twissFile+args.qpNum+"Nb{:.0f}_{:.0f}x{:.0f}betaSpreadNom".format(args.Nb,len(read),args.betaSpread)+paramName+"Hist"
    if args.saveSpread:
        tight_layout()
        savefig(name+"."+args.picFormat,bbox_inches='tight')#,pad_inches=0)
    close()
    print(len(read),paramName,ampl,mu,sigma)
    print(name)
    return mu,sigma,ampl,len(read)













def spreadHist(args,Twiss,paths,origBX,origBY,beamFile):
    #from os import uname
    #from plotFit import plotSpread,plotSpreadBroad
    from numpy import zeros

    paramName=["jMax","beamPOut","coreJMean","rValue"]#,"rDiff"] #len=6 #"coreArea",,"centX","centY"
    paramLabel=[r"Peak Current Density [$\mu$A/cm$^2$]","Beam % Outside Target Area",#,r"Core Area [mm$^2$]",
                r"Core Average Current Density [$\mu$A/cm$^2$]",#,"Beam Center X [mm]","Beam Center Y [mm]",
                "R Value"]#,"R Difference Value"]
    #"Peak Current Density [uA/cm^2]","Beam % Outside Target Area","Core Area [mm^2]","Core Average Current Density [uA/cm^2]","Beam Center X [mm]","Beam Center Y [mm]","R Divide Value","R Difference Value"
    ind = [7,8,10,13]#,14]#9#,11,12,
    unit=[r"$\mu$A/cm$^2$","%",r"$\mu$A/cm$^2$",""]#,""]#r"mm$^2$",,"mm","mm"
                #   ampl, mu, sigma                  ([10,1e3,20],[5e3,1e4,1e3])
    paramFitLims = [(0,[1e3,100,500]),(0,[1e3,100,500]),(0,[1e3,100,500]),#([10,1e1,20],[5e6,1e6,1e6]),([-1e3,-100,-500],[1e3,100,500]),([-1e3,-100,-500],[1e3,100,500]),#
                    (0,[1e3,10,1])]#,(0,[1e3,10,1])]
    pHistLims = [[35,45],[8.2,8.55],[25,40],[.072,.073]]#,[4.55,4.6]]#nominal#[5000,7000],,[-10,10],[-10,10]
    pHistLims = [[35,45],[8.0,8.80],[25,40],[.072,.073]]#,[4.55,4.6]] #[5000,8000],,[-10,10],[-10,10]#if range is outside this, the binning is off...
    paramBins = 20#[30,35,40,30,25,25,25,25]
    #bR = round(args.samples/8)
    #paramBins = [bR,bR,round(args.samples/5),bR,bR,bR,bR,bR] #the coreArea hist goes bonkers with auto--shrugs--
    mus = zeros(len(paramName))
    sigmas = zeros(len(paramName))
    ampl = zeros(len(paramName))
    lens = zeros(len(paramName))
    for i in range(len(paramName)):
        print(paramName[i],paramFitLims[i],paramBins)
        mus[i], sigmas[i],ampl[i],lens[i] = plotSpread(args,Twiss,paths['statsPWD'],paramName[i],ind[i],unit[i],paramLabel[i],paramFitLims[i],paramBins,pHistLims[i],origBX,origBY,beamFile)
        #mus[i], sigmas[i],ampl[i],lens[i] = plotSpreadBroad(args,Twiss,paths['statsPWD'],paramName[i],ind[i],unit[i],paramLabel[i],pFitLims[i],paramBins[i])

    import csv
    found=False
    i=0
    foundRow=0
    with open(paths['statsPWD']+"pOffStats.csv",mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            i+=1
            if float(row[0]) == args.betaSpread: #No duplicates
                found = True
                foundRow = i
                print("Values already in row",foundRow)
                break
        csv_file.close()
    if not found:
        with open(paths['statsPWD']+"pOffStats.csv",mode = 'a') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter = ',')
            for i in range(len(paramName)):
                if sigmas[i] < 1e-3:
                    csv_writer.writerow([args.betaSpread,lens[i],paramName[i],"{:.3f} $\pm$ {:.3e}".format(mus[i],sigmas[i])]) #add JCore
                else:
                    csv_writer.writerow([args.betaSpread,lens[i],paramName[i],"{:.3f} $\pm$ {:.3f}".format(mus[i],sigmas[i])]) #add JCore
            csv_file.close()
    
        for i in range(len(paramName)):
            print(args.betaSpread,lens[i],paramName[i],"{:.3f} +/- {:.3e}".format(mus[i],sigmas[i]))












def plotSpreadBroad(args,Twiss,statsPWD,paramName,ind,unit,paramLabel,pFitLims,paramBins,pHistLims):
    import csv,re
    from matplotlib.pyplot import hist,savefig,close,plot,xlim,ylim,text,title,xlabel,ylabel
    from numpy import greater as npGreater,zeros,mean,std
    #from plotFit import findFit, gaussian

    read = zeros(args.samples)
    read.fill(-10)
    pFitLims[1][2] = 1e3 #increase
    #for i in range(len(Twiss)):
    #    Twiss[i] = round(Twiss[i]*1e5)/1e5
    #print(Twiss)
    i=0
    lenbreak=False
    key = "(N2.9e+05protons)" #"(PBW_570MeV)"
    with open(statsPWD+"EvalStats12MarHEBTA2T.csv",mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            print(re.search(row[0],key))
            if re.search(row[0],key):
            #if row[0]  
            #if float(row[1]) == Twiss[0] and float(row[2]) == Twiss[1] and float(row[3]) == Twiss[2] and float(row[4]) == Twiss[3] and float(row[5]) == Twiss[4] and float(row[6]) == Twiss[5]:
                #print("found Twiss")
                read[i] = float(row[ind]) #[mm-mrad]
                #print(i,read[i])
                if i+1 == args.samples:
                    lenbreak = True
                    break
                i+=1
                #print(lenbreak)
            if lenbreak:
                break
    csv_file.close()

    i#print(paramLabel,len(read))
    nonzero = npGreater(read,-10)
    read = read[nonzero]
    print("guess",paramName,len(read),args.samples/40,mean(read),std(read))
    mu, sigma, ampl,interval = findFit(read,[args.samples/40,mean(read),std(read)],pFitLims,paramBins)
    print(mu,sigma,ampl)

    _, bins, _ = hist(read,interval) #ax =
    #y1 = norm.pdf(bins, mu, sigma) #need to pass it ampl to get the scale right...
    plot(bins, gaussian(bins,ampl,mu,sigma), "r--", linewidth=2)
    title("{:.0f}".format(len(read))+r" Samples $\pm$50% Around Nominal")
    xlabel(paramLabel)
    ylabel("Counts")
    #setp(ax,title=Title,xlabel=xLabel,ylabel=yLabel)
    
    bgdbox=dict(pad=1,fc='w',ec='none')
    fs=14
    delta = 1e-3
    propsR = dict(horizontalalignment="right",verticalalignment="top", backgroundcolor = 'w',bbox=bgdbox,fontsize=fs-4)
    text(xlim()[1]*(1-delta), ylim()[1]*0.95,r"$\mu$="+"{:.3f}".format(mu)+unit+"\n"+r"$\sigma$="+"{:.3f}".format(sigma)+unit, propsR)
    name=statsPWD+paramName+"BroadHist_{:.0f}x{:.0fpOffNom".format(len(read),args.betaSpread)
    if args.saveSpread:
        plt.tight_layout()
        savefig(name+"."+args.picFormat)
    close()
    print(len(read),paramName,mu,sigma)
    print(name)
    return mu,sigma,ampl,len(read)













def getTwiss(args,sample,paths):
    #Get Twiss to run, depending on user configuration
    # Twiss= [NemtX,BetaX,AlphX,NemtY,BetaY,AlphY]
    if args.beamClass == 'Yngve': #smallest from Yngve
        Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
    elif args.beamClass == 'ESS': #from my OpenXAL calculation
        Twiss = [0.118980737408497,1085.63306926394,-65.1638312921654,0.123632934174567,136.062409365455,-8.12599512314246] #updated to HEBT-A2T Combo Twiss
    elif args.beamClass == 'pencil': #"pencil" beam of ~0 emittance
        Twiss = [0.0001,0.15,0,0.0001,0.15,0]

    if args.twiss:
        Twiss = [args.twiss[0],args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5]]
        args.beamClass = "Twiss"

    #For pulling from an already made CSV
    if args.source == "csv":
        import re
        if re.search("beta",paths['csvPWD']+args.beamFile):
            betaX = float(re.search("(([-+]?[0-9]*\.?[0-9]*)+(?=(,([-+]?[0-9]*\.?[0-9]*)+m)))",args.beamFile)[0])
            betaY = float(re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(m_N)))",args.beamFile)[0])
            print(betaX,betaY)

    #Tailored to OpenXAL failure input
    if args.source == "twiss":
        if args.twissFile != "":
            import csv
            twissPWD = paths['homePWD']+"Documents/UiO/Forske/ESSProjects/OpenXAL/OXALNotebooks/failureTwiss/"
            if args.qpNum == "all": #for getting all samples out of betaSample OpenXAL output file
                ##
                ## To Do : Double check this is doing what I want, along with sample appending...
                ##
                #for running all Twiss in OpenXAL Combo Sequence output CSV
                with open(twissPWD+args.twissFile+".csv",mode='r') as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    rowNum = 0 #to increment QP number
                    for row in csv_reader:
                        if rowNum == sample: #i #not sure how to go back to QP finding...
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
                            #Adjusts names so the files are in order when it is a QP failure
                            if len(row[0]) == 2:
                                args.qpNum = "0"
                            else: 
                                args.qpNum = ""
                            args.qpNum += row[0]
                            #print(rowNum,args.qpNum,Twiss)
                            break
                        rowNum+=1
                csv_file.close()
            elif args.qpNum == "qps": #for getting all QPs out of OpenXAL Failure file
                print("In QPs")
                #for running all QPs in OpenXAL Combo Sequence output CSV
                with open(twissPWD+args.twissFile+".csv",mode='r') as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    rowNum = 99 #to increment QP number
                    ###
                    ### To Do 13 Mar: Figure out how to read in each QP in FailureHEBT-A2T_80pctField_1e-04Jitter rather than individually
                    ###
                    for row in csv_reader:
                        if rowNum == int(row[0]): #i #not sure how to go back to QP finding...
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
                            #Adjusts names so the files are in order when it is a QP failure
                            if len(row[0]) == 2:
                                args.qpNum = "0"
                            else: 
                                args.qpNum = ""
                            args.qpNum += row[0]
                            print(sample,row[0],args.qpNum,Twiss)
                            break
                        rowNum+=1
                csv_file.close()
            else:
                #for single QP fail runs
                with open(twissPWD+args.twissFile+".csv",mode='r') as csv_file:
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

    
    #The backup of the Twiss
    origBX = Twiss[1]
    origBY = Twiss[4]
    #If want to do this one set at a time, do this.
    #if args.betaSpread != 0:
    #    origBX = Twiss[1]
    #    origBY = Twiss[4]
    #    print(args.betaSpread,"% Spread around BetaX ({:.1f}m) ranges: {:.1f}-{:.1f}".format(origBX,origBX-origBX*args.betaSpread/100,origBX+origBX*args.betaSpread/100)) #requires 3sigma
    #    print(args.betaSpread,"% Spread around BetaY ({:.1f}m) ranges: {:.1f}-{:.1f}".format(origBY,origBY-origBY*args.betaSpread/100,origBY+origBY*args.betaSpread/100))
    #    #Assuming the "scale" in the following function is the %
    #    from numpy.random import default_rng
    #    rng = default_rng()
    #    Twiss[1] = rng.normal(loc=origBX,scale=origBX*args.betaSpread/100/3) #%/3 so 3sigma range = range
    #    Twiss[4] = rng.normal(loc=origBY,scale=origBY*args.betaSpread/100/3)
    #    print("New BetaX: {:.1f}m, BetaY: {:.1f}m".format(Twiss[1],Twiss[4]))
    #If want to write sample to CSV and read in, follow Carl's way.
    if args.betaSpread != 0:
        #requires that a Twiss is defined above
        origBX = Twiss[1]
        origBY = Twiss[4]
        #print("old",Twiss)

        writeMore = False
        #check if the file is already made
        from os.path import isfile
        name = paths['statsPWD']+"betaSpread_betaX{:.0f},{:.0f}m_{:.0f}Pct".format(origBX,origBY,args.betaSpread)
        print(name)

        #Don't have to manually make the file for a new %
        if not isfile(name+".csv"):
            #start a new file
            import csv
            with open(name+".csv",mode = 'w') as csv_file:
                csv_writer = csv.writer(csv_file,delimiter = ',')
                csv_writer.writerow(["i","Beta X [m]","Beta Y [m]"])
                csv_file.close()
            writeMore = True
        else:
            #If file present, read in values for the samples requested
            #Get length of file, in case user requested more samples than previously made
            #from plotFit import numLines
            iters = numLines(name)-1
            #print("File has",iters,"lines")

            #read in what is there
            import csv
            with open(name+".csv",mode='r') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                next(csv_reader) #skips header
                rowNum = 0
                for row in csv_reader:
                    #print(sample,rowNum,row[0],row[1],row[2])
                    #for the sample passed in, get that row values
                    if int(row[0]) == int(sample):
                        #print("Reading",rowNum, sample,row[0],row[1],row[2])
                        #immediately put in value as the Twiss
                        Twiss[1] = float(row[1])
                        Twiss[4] = float(row[2])
                        print("Read in Betas:",Twiss[1],Twiss[4])
                        break

                    #if need to make more samples than in file, flag and exit before EOF error
                    if int(sample) >= int(iters):
                        print("User requested more Twiss than currently written")
                        writeMore = True
                        break

                    rowNum += 1
                csv_file.close()
        #only if need to write more 
        #print("Writing more?",writeMore)
        if writeMore:
            import csv, numpy as np
            from numpy.random import default_rng
            #save original values
            origBX = Twiss[1]
            origBY = Twiss[4]

            #get a random value for betaX and betaY and saves into Twiss array
            Twiss[1] = default_rng().normal(loc=origBX,scale=origBX*args.betaSpread/100)
            Twiss[4] = default_rng().normal(loc=origBY,scale=origBY*args.betaSpread/100)
            if Twiss[4] <= 50 or Twiss [1] <= 50: #make sure not too small
                Twiss[1] = default_rng().normal(loc=origBX,scale=origBX*args.betaSpread/100)
                Twiss[4] = default_rng().normal(loc=origBY,scale=origBY*args.betaSpread/100)
                if Twiss[4] <= 50 or Twiss [1] <= 50:
                    Twiss[1] = default_rng().normal(loc=origBX,scale=origBX*args.betaSpread/100)
                    Twiss[4] = default_rng().normal(loc=origBY,scale=origBY*args.betaSpread/100)
                    if Twiss[4] <= 50 or Twiss [1] <= 50:
                        Twiss[1] = default_rng().normal(loc=origBX,scale=origBX*args.betaSpread/100)
                        Twiss[4] = default_rng().normal(loc=origBY,scale=origBY*args.betaSpread/100)
                        if Twiss[4] <= 50 or Twiss [1] <= 50:
                            Exception("Get Twiss, a Beta has been <50 four times!")
                        #can't happen 4 times in a row that the value is 2 sigma off
            print(origBX,origBY,"Newly written Betas:",Twiss[1],Twiss[4])

            #append new value to the file so don't have to redo it
            with open(name+".csv",mode = 'a') as csv_file:
                csv_writer = csv.writer(csv_file,delimiter = ',')
                csv_writer.writerow([sample,Twiss[1],Twiss[4]]) #save the betas for this sample number, to be used next time.
                csv_file.close()

        ##Yngve: calculate changed alpha as well as beta:
        ##Erik - no need for now
        #from numpy import sqrt
        #origAX = Twiss[2]
        #b1b0X = Twiss[1]/origBX
        #Twiss[2] = origAX + 0.5*(-sqrt(4*(origAX**2) * b1b0X + 8*origAX*b1b0X + 4*b1b0X - 4 / (origBX**2))) - origAX
        #print("{:.2f}\t\t{:.2f}\t\t{:.3f}\t\t{:.3f}".format(origBX,Twiss[1],origAX,Twiss[2]))
        #origAY = Twiss[5]
        #b1b0Y = Twiss[4]/origBY
        #Twiss[5] = origAY + 0.5*(-sqrt(4*(origAY**2) * b1b0Y + 8*origAY*b1b0Y + 4*b1b0Y - 4 / (origBY**2))) - origAY
        #print("{:.3f}\t\t{:.3f}\t\t{:.3f}\t\t{:.3f}".format(origBY,Twiss[4],origAY,Twiss[5]))
    #print(origBX,origBY,"New Betas:",Twiss[1],Twiss[4],"new Alphas:",Twiss[2],Twiss[5])

    return Twiss,origBX,origBY
