#plotFit.py
#Eric Fackelman
#29 March- 27 Sep 2022

#This holds the functions to plot the figures for the analysisScripts
#-plotFit - 


def plotFit(xs,ys,savename,xlim,ylim,material,thick): 
    import numpy as np
    import matplotlib.pyplot as plt
    from plotFit import findFit
    #For creating the PDF of particle distribution from findFit function optimized mu and sigma.
    from scipy.stats import norm 

    fs = 14 #set the axis label font size early

    #Use Scipy.optimize.curve_fit in my findFit function to get mus and sigmas:
    mux, sigmax, xinterval = findFit(xs) #dynamically gets parameters AND histogram intervals!
    muy, sigmay, yinterval = findFit(ys)

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
    s1.set_title("At {:.2f}mm ".format(thick)+mat+" PBW Exit, X Distribution \n"+rf"$\sigma = {{:.3f}}$,".format(sigmax)+
            " {:.3f}% outside".format(pOut3sigx)+r" 3$\sigma$",fontsize=fs)
    s1.set_xlabel("X Position [mm]",fontsize=fs)
    s1.set_ylabel("Probability Density",fontsize=fs)

    #s2.set_title("After "+mat+", Fitted Y Distribution \n"+rf"$\mu= {{:.1e}}, \sigma = {{:.3f}}$,".format(muy,sigmay)+
    #        " {:.3f}% outisde".format(pOut3sigy)+r" 3$\sigma$",fontsize=fs)
    s2.set_title("At {:.2f}mm ".format(thick)+mat+" PBW Exit, Y Distribution \n"+rf"$\sigma = {{:.3f}}$,".format(sigmay)+
            " {:.3f}% outside".format(pOut3sigy)+r" 3$\sigma$",fontsize=fs)
    s2.set_xlabel("Y Position [mm]",fontsize=fs)
    s2.set_ylabel("Probability Density",fontsize=fs)

    #Can date stamp the multi plot for easier tracking of changes, if necessary
    #from datetime import datetime
    #dt = datetime.now()

    #Update Name
    #name = savename+"_"+dt.strftime("%H-%M-%S")+".png"  
    name = savename+".png" #I'd prefer to have it overwrite old files for now
    print(name) #show so it is easy to find the file
    plt.savefig(name,bbox_inches="tight")
    #plt.show()
    plt.close() #be sure to close the plot

def findFit(data):
    #with help from MSeifert in stackoverflow fit a curve to a histogram in python
    #https://stackoverflow.com/questions/35544233/fit-a-curve-to-a-histogram-in-python
    #import gaussian
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    def gaussian(x,mu,amplitude,sigma):
        #import numpy.exp as exp
        return amplitude * np.exp( - (x - mu) ** 2 / (2*sigma ** 2))

    plt.close()
    bin_heights, bin_borders, _ = plt.hist(data,bins="auto",label="histogram")
    #print(len(bin_borders))
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    popt,_ = curve_fit(gaussian, bin_centers, bin_heights,p0=[0.01,1.,1.],bounds=(0,[1e2,1e5,1e2])) #p0 should be good start, though not what is recommended
    x_interval_for_fit = np.linspace(bin_borders[0],bin_borders[-1],len(bin_borders))
    #plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit,*popt),label="fit")
    #plt.legend()
    #plt.xlim([bin_borders[0],bin_borders[-1]])
    #plt.show()
    plt.close() #be sure to close this or else these show up in the multi plots!
    #print(popt)
    #return the mean, abs(sigma), interval #sigma can sometimes be - so must abs() it.
    return(popt[0],abs(popt[2]),x_interval_for_fit) #for Vac and Air sigma sometimes give - number...

def calcTwiss(labelx,labely,x,y):
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
    #16.8.22 No longer useful, as it doesn't use layer by layer
    #For plotting particle Distribution at Target Exit with PDF from Twiss parameters.
    import numpy as np
    import matplotlib.pyplot as plt
    from plotFit import findFit, calcTwiss, getMoments
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
    mux,sigmax,xinterval = findFit(xs)
    mupx,sigmapx,pxinterval = findFit(pxs) #this is a least square fit to Gaussian, so less sensitive to outliers.
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

    name = savename+"_Twiss"+axis+"'"+dt.strftime("%H-%M-%S")  +".png"##
    #plt.show()
    plt.savefig(name,bbox_inches="tight")
    plt.close()
    print(name)

def toTarget(TwPm,label):
    #Extend the distributions to the Target Location
    import numpy as np
    #Twiss=[beta,alpha,gemt,gamma]
    PBWexitBetaMx = np.array([[TwPm[0],-TwPm[1]],[-TwPm[1],TwPm[3]]])

    d_PBW_Targ = 5
    drift_PBW_Targ = np.array([ [1, d_PBW_Targ],[ 0, 1]])
    Calc_PBW_Targ = np.linalg.multi_dot([drift_PBW_Targ,PBWexitBetaMx,np.transpose(drift_PBW_Targ)])
    #print(label,"PBWexit:",PBWexitBetaMx,"\n",drift_PBW_Targ,"\n",Calc_PBW_Targ)

    return [Calc_PBW_Targ[0][0],-Calc_PBW_Targ[1][0],TwPm[2],Calc_PBW_Targ[1][1]]

def compareTargets(targx,targy,targTwx,targTwy,fitTwx,fitTwy,fitlabel,savename,mat,PBWTwx,PBWTwy):
    import numpy as np
    from plotFit import findFit, getMoments
    from scipy.stats import norm 
    import matplotlib.pyplot as plt
    mm = 1e-3
    fs = 14

    #Find fit to histogram
    mux,sigmax,xinterval = findFit(targx)
    muy,sigmay,yinterval = findFit(targy)

    #Find range of particles that are outside 3 sigma
    sigx=np.abs(targx)>3*sigmax# and np.abs(xs)<10*sigma)
    sigy=np.abs(targy)>3*sigmay

    #Find % of particles outside 3 sigma
    pOut3sigx = len(targx[sigx])/len(targx)*100
    pOut3sigy = len(targy[sigy])/len(targy)*100
    print(fitlabel,pOut3sigx,"% outside 3 sigma X")
    print(fitlabel,pOut3sigy,"% outside 3 sigma Y")

    #Get sigma for position and angle
    Mfx = getMoments(fitTwx)
    Mfy = getMoments(fitTwy)
    Mhx = getMoments(targTwx)
    Mhy = getMoments(targTwy)
    #print(fitlabel,Mhx,Mfx)

    print("Histogram",fitlabel,"sigma X: {:.2f}mm".format(Mhx[0]))
    print("Histogram",fitlabel,"sigma Y: {:.2f}mm".format(Mhy[0]))
    print("Fit      ",fitlabel,"sigma X: {:.2f}mm".format(Mfx[0]))
    print("Fit      ",fitlabel,"sigma Y: {:.2f}mm".format(Mfy[0]))

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
    if re.search(",",fitlabel):
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
    s1.set_title("X Distribution At Target after "+fitlabel,fontsize=fs) #+"\n"+rf"$\sigma_D=${{:.1f}}mm".format(sigmax)
    s1.set_xlabel("X Position [mm]",fontsize=fs)
    s1.set_ylabel("Probability Density",fontsize=fs)
    xlim = 4*sigmax
    s1.set_xlim([-xlim,xlim])
    #s1.set_ylim([1e-6,0.1])
    s1.set_ylim([1e-6,1])  #if don't set, then log goes to e-58
    xlim1 = s1.get_xlim()
    ylim1 = s1.get_ylim()
    s1.text(xlim1[0]*0.97,3e-2,sigmatextx,fontsize=fs-2)
    #PBWTwx = [Ibetax,Ialphx,Inemtx]
    s1.text(xlim1[0]*0.97,1.2e-2,"Beam Twiss at PBW:",fontsize=fs-4)
    s1.text(xlim1[0]*0.97,7e-3,r"$\epsilon_{Nx}$ = "+"{:.3f} [mm*mrad]".format(PBWTwx[2]),fontsize=fs-4)
    s1.text(xlim1[0]*0.97,4e-3,r"$\beta_{x}$ = "+"{:.1f} [m]".format(PBWTwx[0]),fontsize=fs-4)
    s1.text(xlim1[0]*0.97,2.4e-3,r"$\alpha_{x}$ = "+"{:.2f}".format(PBWTwx[1]),fontsize=fs-4)

    handles1, labels1 = plt.gca().get_legend_handles_labels()
    order1=[0,1,2]
    s1.legend([handles1[idx] for idx in order1],[labels1[idx] for idx in order1],fontsize=9)

    s2.set_xlim([-xlim,xlim])
    s2.set_ylim([1e-6,1]) #if don't set, then log goes to e-58
    s2.set_title("Y Distribution At Target after "+fitlabel,fontsize=fs)
    s2.set_xlabel("Y Position [mm]",fontsize=fs)
    s2.set_ylabel("Probability Density",fontsize=fs)

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

    name = savename+"_log_"+dt.strftime("%H-%M-%S")+".pdf"##
    #plt.show()
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
    from plotFit import findFit, getMoments
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
    mux,sigmax,xinterval = findFit(targx)
    muy,sigmay,yinterval = findFit(targy)

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

    s1.set_title("X Raster Distribution At "+position,fontsize=fs) #+"\n"+rf"$\sigma_D=${{:.1f}}mm".format(sigmax)
    s1.set_xlabel("X Position [mm]",fontsize=fs)
    s1.set_ylabel("Probability Density",fontsize=fs)
    s1.set_xlim([-xlim1,xlim1])
    s1.set_ylim(ylim1) #if don't set, then log goes to e-58

    s2.set_title("Y Raster Distribution At "+position,fontsize=fs)
    s2.set_xlabel("Y Position [mm]",fontsize=fs)
    s2.set_ylabel("Probability Density",fontsize=fs)
    s2.set_xlim([-xlim2,xlim2])
    s2.set_ylim(ylim2) #if don't set, then log goes to e-100

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
    from datetime import datetime
    dt = datetime.now()
    #print("There are ",N,"lines "+dt.strftime("%H-%M-%S"))
    return N


def rasterImage(savename,position,histogram2D,parts,savePics,Twiss,rasterXAmplitude,rasterYAmplitude,options,boxes):
    import numpy as np
    name=savename+"_"+position+"Image"
    um = 1e-6
    from datetime import datetime
    start= datetime.now()
    #print(start.strftime("%H-%M-%S"))
    from plotFit import converter
  
    #print(datetime.now().strftime("%H-%M-%S"))
    (Img, xax, yax) = converter(histogram2D,options['saveHist'],name) #convert from TH2D to numpy map
    #print(datetime.now().strftime("%H-%M-%S"))
    xBinSize = xax[501]-xax[500]
    yBinSize = yax[501]-yax[500]

    #99% box #for multiple boxes, use arrays
    #From Dmitriy Work on https://stackoverflow.com/questions/432112/is-there-a-numpy-function-to-return-the-first-index-of-something-in-an-array
    #  np.ndenumerate scales the best of the methods compared to 10^4 elements.
    #Find % outside the 99% Box area
    sumTot = np.sum(Img)+1 #so no 'divide by 0 errors'
    Pprotons2 = sumTot / parts * 100
    #print(sumTot,parts,sumCore,pOutsideBox,Img.max())
    #Img[boxBInd:boxTInd,boxLInd:boxRInd] = 0
    widths = np.zeros(len(boxes))
    heights = np.zeros(len(boxes))
    pLargers = np.zeros(len(boxes)) #percent larger than initial box(?)
    coreImax = np.zeros(len(boxes))
    coreMeanI = np.zeros(len(boxes))
    widths[0] = 160
    heights[0] = 64
    boxLLxs = np.zeros(len(boxes))
    boxLLys = np.zeros(len(boxes))
    pOutsideBoxes = np.zeros(len(boxes))
    Pprotons = sumTot / parts * 100 #this is constant
    cols = ["r","cyan","gold","lime","k"]
    for i in range(len(boxes)): 
        pLargers[i] = boxes[i]
        widths[i] = round(widths[0]*(1+pLargers[i]) / 2) * 2
        heights[i] = round(heights[0]*(1+pLargers[i]) / 2) * 2
        boxLLxs[i] = round(-widths[i]/2 / 2) * 2 #-90
        boxLLys[i] = round(-heights[i]/2 / 2) * 2 #-80
        #print(len(xax),boxLLxs[i],boxLLys[i],widths[i],heights[i])

        for idx, val in np.ndenumerate(xax):
            if int(val) == boxLLxs[i]:
                boxLInd = idx[0] #455
            if int(val) == boxLLxs[i]+widths[i]:
                boxRInd = idx[0] #545
        for idx, val in np.ndenumerate(yax):
            if int(val) == boxLLys[i]:
                boxBInd = idx[0] #460
            if int(val) == boxLLys[i]+heights[i]:
                boxTInd = idx[0] #540
        if i == 0: #for base box case
            boxBIndC = boxBInd #set to use later
            boxTIndC = boxTInd
            boxLIndC = boxLInd
            boxRIndC = boxRInd
        #Find % outside the 99% Box area
        core = Img[boxBInd:boxTInd,boxLInd:boxRInd]
        sumCore = np.sum(core)+1
        pOutsideBoxes[i] = (sumTot-sumCore)/sumTot*100
        #print(sumCore,pOutsideBoxes[i])
        coreN = Img[boxBInd:boxTInd,boxLInd:boxRInd]
        coreImax[i] = coreN.max() #[uA/cm^2]
        coreMeanI[i] = np.mean(coreN) #[uA/cm^2]

    #Normalize to full current, see 
    I_pulse = 62.5*1e3 #[uA]
    C = I_pulse / parts / (xBinSize * yBinSize * 1e-2) * 0.04 #[uA/cm^2]: /mm^2->/cm^2 = /1e-2, A->uA = 1e6, 4% duty cycle
    Protmax = Img.max() #protons/mm^2
    for i in range(len(Img)):
        for j in range(len(Img[i])):
            Img[i][j] = Img[i][j] * C #[uA/cm^2]
    Imax = Img.max()
    Imin = 0.9 * C #background (0 hits) will be un-colored
    print(coreMeanI*C,coreImax*C)
    print("C {:.3f}, Proton max {:.0f}, Imax {:.1f}, Imin {:.3f}, coreMeanI {:.1f}, pOutsideBox".format(C,Protmax,Imax,Imin,coreMeanI[0]*C),pOutsideBoxes)

    #R value for algorithm. Works when use Current Density, not Nprotons
    rValue = rCompare(Img,options['Nb'])
    print("R = ",rValue)
    #print("Converted in",datetime.now() - start)

    #Flat Top Current density calculations
    top=0
    Itop = 40
    idxMinX = 1000
    idxMaxX = 1
    idxMinY = 1000
    idxMaxY = 1
    for idx,val in np.ndenumerate(Img):
        if val >= Itop:
            top += 1
            if idx[0] < idxMinX: idxMinX = idx[0]
            if idx[0] > idxMaxX: idxMaxX = idx[0]
            if idx[1] < idxMinY: idxMinY = idx[1]
            if idx[1] > idxMaxY: idxMaxY = idx[1]
    print("Current above",Itop,"uA/cm^2 in",top,"mm^2",idxMaxX-idxMinX,"x",idxMaxY-idxMinY,"mm^2,","{:.2f} uA/cm^2 average".format(np.sum(Img[idxMinY:idxMaxY,idxMinX:idxMaxX])/top))

    if savePics:
        from matplotlib.pyplot import subplots,pcolor,close,tight_layout,savefig
        from matplotlib.patches import Rectangle
        from matplotlib.colors import LogNorm
        X, Y = np.meshgrid(xax,yax) #Set mesh of plot from the axes given from converter function
        close() #make sure no prior plotting messes up

        fig,ax = subplots(dpi=150,figsize=(6.4,4.8),tight_layout=True)
        ax.set_xlim([-150,150]) #for close up of beam
        ax.set_ylim([-75,75])
        #print(datetime.now().strftime("%H-%M-%S"))
    
        #Set maximum value depending on the maximum current density
        from math import log10,ceil
        maxim = ceil(Imax / 10) * 10
        if options['maxim'] != 0: maxim = options['maxim'] #user provided maximum
        minim = 10**ceil(log10(Imin))
        cbarVals  = [minim,minim*10,minim*100,minim*0.455,maxim] #make array for color bar values
        cbarLabels = ["{:.0f}".format(cbarVals[0]),"{:.0f}".format(cbarVals[1]),"{:.0f}".format(cbarVals[2]),
                    "{:.0f}".format(cbarVals[3]),"{:.0f}".format(cbarVals[4])]#,"{:.1f}".format(cbarVals[5])] #make labels of Value
        cbarLabel = r"$\mu A / cm^2$"
        #print("Max current Density: ",Img.max(),"/",maxim,datetime.now().strftime("%H-%M-%S"))

        #Use pcolor to show density map, use log scale
        c = ax.pcolor(X,Y,Img,shading='auto',norm=LogNorm(vmin=minim, vmax=maxim), cmap='viridis') #viridis or magma are perceptually uniform
        lw=1
        col='k'
        fs=14

        #Set Plot Properties
        ax.set_xlim([-options['xlim'],options['xlim']])
        ax.set_ylim([-options['ylim'],options['ylim']])
        ax.set_title("Proton Beam Distribution at "+position,fontsize=fs)
        ax.set_xlabel("Horizontal [mm]")
        ax.set_ylabel("Vertical [mm]")
        cbar = fig.colorbar(c, ax=ax,pad=0.01)
        cbar.set_label(cbarLabel,labelpad=2,fontsize=fs-2)
        cbar.set_ticks(cbarVals)
        cbar.set_ticklabels(cbarLabels)

        if position =="Target":
            #name = name+"_2212only"
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            
        #Show 99% box
        if not options['noBox']: #for user clarity, call is noBox
            ax.add_patch(Rectangle((boxLLxs[0],boxLLys[0]),widths[0],heights[0],linewidth=lw,edgecolor=cols[0],fill=False))

        if not options['noText']: #for user clarity, call is noText
            ax.set_title("Distribution at "+position+"\n{:.3f}% of {:.2e} total protons".format(Pprotons,parts),fontsize=fs)

            #Display beam characteristics
            bgdbox=dict(pad=2,fc='w',ec='none')
            propsR = dict(horizontalalignment="right",verticalalignment="bottom", backgroundcolor = 'w',bbox=bgdbox)
            ax.text(xlim[0]*0.90, ylim[1]*0.85, "{:.2f}".format(pOutsideBoxes[0])+r"% Outside 160x64mm$^2$ Box", color=cols[0], fontsize = fs-2, fontweight='bold', backgroundcolor = 'w')
            ax.text(xlim[0]*0.97, ylim[0]*0.60, "Beam Twiss at PBW:", fontsize=fs-4, backgroundcolor = 'w',bbox=bgdbox)
            ax.text(xlim[0]*0.97, ylim[0]*0.71, r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[mm \cdot mrad]}$", fontsize=fs-4, backgroundcolor = 'w',bbox=bgdbox)
            ax.text(xlim[0]*0.97, ylim[0]*0.83, r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$", fontsize=fs-4, backgroundcolor = 'w',bbox=bgdbox)
            ax.text(xlim[0]*0.97, ylim[0]*0.95, r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]), fontsize=fs-4, backgroundcolor = 'w',bbox=bgdbox)
            #ax.text(xlim[1]*0.97, ylim[0]*0.57, r"Box <$\bf{J}$>: "+"{:.1f}".format(coreMeanI)+r" $\mu$A/cm$^2$", propsR,fontsize=fs-4)
            ax.text(xlim[1]*0.97, ylim[0]*0.60, "R={:.4f}".format(rValue), propsR,fontsize=fs-2)
            ax.text(xlim[1]*0.97, ylim[0]*0.75, r"Max $\bf{J}$: "+"{:.1f}".format(Imax)+r" $\mu$A/cm$^2$", propsR,fontsize=fs-4)
            ax.text(xlim[1]*0.97, ylim[0]*0.85, "RM Amplitudes: {:.1f}, {:.1f}mm".format(rasterXAmplitude,rasterYAmplitude),propsR,fontsize=fs-4)
            ax.text(xlim[1]*0.97, ylim[0]*0.96, options['physList'], propsR,fontsize=fs-4)

            for i in range(1,len(boxes)): #make multiple boxes
                ax.add_patch(Rectangle((boxLLxs[i],boxLLys[i]),widths[i],heights[i],linewidth=lw,edgecolor=cols[i],fill=False))
                ax.text(xlim[0]*0.90, ylim[1]*(0.85-i*0.1), "{:.2f}".format(pOutsideBoxes[i])+"% Outside {:.0f}% Larger Box".format(pLargers[i]*100), 
                                  color=cols[i], fontweight='bold',fontsize = fs-2, backgroundcolor = 'w',bbox=dict(pad=1))#,path_effects=[path_effects.withStroke(linewidth=1, foreground='k')])
    
        dt = datetime.now()
        #from os.path import isfile
        #if isfile(name+"*.png"):
        #  print("already present")
        savefig(name+"_"+dt.strftime("%H-%M-%S")+".png")
        close(fig)
        close()
        print(name+"_"+dt.strftime("%H-%M-%S")+".png")

    #dt = datetime.now()
    #print(dt-start)

    from plotFit import findEdges
    findEdges(Img,Imax,savename,xax,yax)

    return pOutsideBoxes[0], Imax, coreMeanI

def converter(hIn,saveHist,name):
    import ROOT
    import numpy as np

    #Get X axis from ROOT
    xax = np.zeros(hIn.GetXaxis().GetNbins()+1)
    for i in range(len(xax)-1):
        xax[i] = hIn.GetXaxis().GetBinLowEdge(i+1) #Set elements as bin edges
    xax[-1] = hIn.GetXaxis().GetBinUpEdge(hIn.GetXaxis().GetNbins())

    #Get Y axis from ROOT
    yax = np.zeros(hIn.GetYaxis().GetNbins()+1)
    for i in range(len(yax)-1):
        yax[i] = hIn.GetYaxis().GetBinLowEdge(i+1) #Set elements as bin edges
    yax[-1] = hIn.GetYaxis().GetBinUpEdge(hIn.GetYaxis().GetNbins())
    
    hOut = np.zeros( (len(yax), len(xax)) ) #2D Image map

    #Fill Image map with 2D histogram values
    for xi in range(hIn.GetXaxis().GetNbins()):
        for yi in range(hIn.GetYaxis().GetNbins()):
            bx = hIn.GetBin(xi+1,yi+1)
            hOut[yi,xi] = hIn.GetBinContent(bx)
    
    #Must add Overflow options!!!

    if saveHist:
        import os
        if os.uname()[1] == "tensor.uio.no":
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

        with open("/scratch2/ericdf/PBWScatter/"+name+".csv",mode = 'w',newline=None) as hist_file:
            hist_writer = csv.writer(hist_file,delimiter = ',')
            hist_writer.writerows(hOut)
        hist_file.close()
        print(name+".csv")

    return (hOut,xax,yax)


#add function to check if ROOT file is complete, else delete and make new



def rCompare(Im,Nb):
    import numpy as np, os

    #Find reference files
    if os.uname()[1] == "tensor.uio.no":
    if Nb == 10:
          Iref = np.genfromtxt(open("/scratch2/ericdf/PBWScatter/Vac_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03_run_QBZ_TargetImage.csv"),delimiter=",")
          #Iref = np.genfromtxt(open("/scratch2/ericdf/PBWScatter/PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03_runW_QBZ_TargetImage.csv"),delimiter=",")
    elif Nb == 100:
          Iref = np.genfromtxt(open("/scratch2/ericdf/PBWScatter/Vac_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+06_NpB100_NPls1e+03_run_QBZ_TargetImage.csv"),delimiter=",")
          #Iref = np.genfromtxt(open("/scratch2/ericdf/PBWScatter/PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+06_NpB100_NPls1e+03_runW_QBZ_TargetImage.csv"),delimiter=",")

    lenx = np.shape(Im)[0]
    leny = np.shape(Im)[1]
    diff = np.zeros((leny,lenx))

    for i in range(lenx):
        for j in range(leny):
            if Iref[j,i] == 0: continue #not great, but it produces better values 15.1.23
            diff[j,i] = ( ( ( Im[j,i] / Iref[j,i] - 1) ** 2 ) / (leny * lenx) )

    return np.sqrt(np.sum(diff))

def findEdges(Img,Imax,savename,xax,yax):
    import numpy as np
    from scipy.signal import convolve2d as sciSigConv2d
    #does it need to be normalized?

    Gradxx = np.array([[1,0,-1],[2,0,-2],[1,0,-1]])
    Gradyy = np.array([[1,2,1],[0,0,0],[-1,-2,-1]])

    Gradx = sciSigConv2d(Img,Gradxx,'same')
    Grady = sciSigConv2d(Img,Gradyy,'same')

    Gradient = np.sqrt(Gradx ** 2 + Grady ** 2)

    import matplotlib.pyplot as plt
    plt.close()
    X,Y = np.meshgrid(xax,yax)
    fig = plt.figure(figsize=(15,5))
    plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    ax1.pcolor(X,Y,Gradx,shading='auto')
    ax2.pcolor(X,Y,Grady,shading='auto')
    ax3.pcolor(X,Y,Gradient,shading='auto')
    plt.setp(ax1,xlim=(-100,100),ylim=(-100,100))
    plt.setp(ax2,xlim=(-100,100),ylim=(-100,100))
    plt.setp(ax3,xlim=(-100,100),ylim=(-100,100))
    #plt.savefig(savename+"_GradPics.png")
    print(savename,"_GradPics.png",sep="")
    plt.close()

def fitGaussians(Img):
    import numpy as np, ROOT

