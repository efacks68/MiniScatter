#plotFit.py
#Eric Fackelman
#29 March- 5 April 2022

#This is to plot the SimpledemoPBWscript histograms 
# Fits are made in findFit function using scipy.optimize.curve_fit found parameters
# And passed back into plotFit.
#import stuff here

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

  #Find range of particles that are within 3 sigma
  sigx=np.abs(xs)>3*sigmax# and np.abs(xs)<10*sigma)
  sigy=np.abs(ys)>3*sigmay

  #Find % of particles within 3 sigma
  pOut3sigx = len(xs[sigx])/len(xs)*100
  pOut3sigy = len(ys[sigy])/len(ys)*100
  #print(pOut3sigx,"% outisde 3 sigma")
  #print(pOut3sigy,"% outisde 3 sigma")

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
  print("init",Twiss)
  #Twiss=[beta,alpha,gemt,gamma]
  #Calculations from Eq 7 and 8 from Twiss MCS Formalism Calculations from https://cds.cern.ch/record/499590
  e8dgem = 0.5 * Twiss[0]*m * thetasq * thetasq #[m*rad^2]
  e8alph = Twiss[2] * Twiss[1] / (Twiss[2] + e8dgem)
  e8beta = Twiss[2] * Twiss[0]*m / (Twiss[2] + e8dgem) #[m]
  e8gamma = (Twiss[2] * Twiss[3] + thetasq * thetasq ) / (Twiss[2] + e8dgem) #[m^-1]
  e8gemt = Twiss[2] + e8dgem
  #e8nemt = e8nemt/(beta_rel*gamma_rel)
  print("e8",thetasq,e8beta,e8alph,e8gemt,e8gamma)
  #28.7-supposed to have thetasq*thetasq in gamma. Previously did NOT have it! Now numbers are great!
  return [e8beta,e8alph,e8gemt,e8gamma]

def calcEq16(thetasq,Twiss,thick,beta_rel,gamma_rel):
  #import numpy as np
  m=1
  um=1e-6
  mm=1e-3

  #Twiss=[beta,alpha,gemt,gamma]
  #Calculations from Eq 15 and 16 from Twiss MCS Formalism Calculations from https://cds.cern.ch/record/499590
  e16dgem = 0.5 * thetasq * thetasq * (Twiss[0]*m + thick*mm * Twiss[1] + thick*mm*thick*mm/3 * Twiss[3]) #[m*rad^2]
  e16alph = (Twiss[2] * Twiss[1] - thick*mm * 0.5 * thetasq* thetasq ) / (Twiss[2] + e16dgem)
  e16beta = (Twiss[2] * Twiss[0]*m + thick*mm * thick*mm / 3 * thetasq* thetasq ) / (Twiss[2] + e16dgem) #[m]
  e16gamma = (Twiss[2] * Twiss[3] + thetasq * thetasq ) / (Twiss[2] + e16dgem) #[m^-1]
  e16gemt = Twiss[2]*um + e16dgem
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
    y1c = norm.pdf(binsx, mux, Me8[0])
    l1c = s1.plot(binsx, y1c, "b.", linewidth=1,label=PDFlabele8+" "+axis+" Twiss")
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
  #finish setting plot characteristics
  s1.set_title(titlestart+", Twiss Fit "+axis+" Distribution",fontsize=fs) #+"\n"+rf"$\sigma_D=${{:.1f}}mm".format(sigmax)
  s1.set_xlabel(axis+" Position [mm]",fontsize=fs)
  s1.set_ylabel("Probability Density",fontsize=fs)
  handles1, labels1 = plt.gca().get_legend_handles_labels()
  #print(labels1)
  s1.legend([handles1[idx] for idx in order],[labels1[idx] for idx in order],fontsize=9)
  ylim1=s1.get_ylim() #dynamically get the ylimits
  xlim1=s1.get_xlim()
  #text="Original:\n"+rf"$ \epsilon_N: {{:.1e}}$m, $\sigma_F:{{:.1f}}$mm".format(orig,sigmaTwx)+"\n"+"Modified:\n"+rf"$\epsilon_N:{{:.1e}}$m, $\sigma_F:${{:.1f}}mm".format(calcxTwf[2],csigmaTwx)
  if axis=="X":
    sigmatextx = r"$\sigma_{X}$:"
    sigmatextpx = r"$\sigma_{X'}$:"
  elif axis=="Y":
    sigmatextx = r"$\sigma_{Y}$:"
    sigmatextpx = r"$\sigma_{Y'}$:"
  sigmatextx +="\nLeast Square = "+"{:.2f}".format(sigmax)+"mm"
  sigmatextx +="\nTwiss RMS = "+"{:.2f}".format(Mcal[0])+"mm"
  sigmatextx +="\n"+PDFlabele8+" RMS = "+"{:.2f}".format(Me8[0])+"mm"
  sigmatextx +="\n"+PDFlabele16+" RMS = "+"{:.2f}".format(Me16[0])+"mm"
  sigmatextpx +="\nLeast Square = "+"{:.2f}".format(sigmapx/mm)+"mrad"
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
  #print("xlim2: {:.1e}-{:.1e}, ylim2: {:.1e}-{:.1e}".format(xlim2[0],xlim2[1],ylim2[0],ylim2[1]))
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
  #Now extend the distributions to the Target Location
  #But I should also extend the Eq Twiss Params to Target separately from the MiniScatter Output.
  import numpy as np
  #Twiss=[beta,alpha,gemt,gamma]
  PBWexitBetaMx = np.array([[TwPm[0],-TwPm[1]],[-TwPm[1],TwPm[3]]])

  d_PBW_Targ = 5
  drift_PBW_Targ = np.array([ [1, d_PBW_Targ],[ 0, 1]])
  Calc_PBW_Targ = np.linalg.multi_dot([drift_PBW_Targ,PBWexitBetaMx,np.transpose(drift_PBW_Targ)])
  #print(label,"PBWexit:",PBWexitBetaMx,"\n",drift_PBW_Targ,"\n",Calc_PBW_Targ)
  
  return [Calc_PBW_Targ[0][0],-Calc_PBW_Targ[1][0],TwPm[2],Calc_PBW_Targ[1][1]]

def compareTargets(targx,targy,targTwx,targTwy,fitTwx,fitTwy,fitlabel,savename,mat):
  import numpy as np
  from plotFit import findFit, getMoments
  from scipy.stats import norm 
  import matplotlib.pyplot as plt
  mm = 1e-3
  fs = 14

  mux,sigmax,xinterval = findFit(targx)
  muy,sigmay,yinterval = findFit(targy)

  Mfx = getMoments(fitTwx)
  Mfy = getMoments(fitTwy)
  Mtx = getMoments(targTwx)
  Mty = getMoments(targTwy)
  print(fitlabel,Mtx,Mfx)

  #Create the fig with 2 plots side by side
  plt.clf()
  fig = plt.figure(figsize=(15,6.0))
  plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
  s1 = fig.add_subplot(1,2,1)
  s2 = fig.add_subplot(1,2,2)

  #Make the histogram of the full energy distrubtion for X with findFit intervals
  nx, binsx, patchesx = s1.hist(targx, xinterval, density=True, facecolor="green", alpha=0.75,label="MiniScatter Target X Histogram")

  #Add the "best fit" line using the findFit mu and sigma for x and display sigma
  y1a = norm.pdf(binsx, mux, sigmax)
  l1a = s1.plot(binsx, y1a, "k--", linewidth=1,label="Histogram Least Square")
  y1b = norm.pdf(binsx, mux, Mtx[0])
  l1b = s1.plot(binsx, y1b, "r--", linewidth=2,label="Target Twiss RMS")
  y1c = norm.pdf(binsx, mux, Mfx[0])
  l1c = s1.plot(binsx, y1c, "b--", linewidth=2,label=fitlabel+" Twiss RMS")

  #Make the histogram of the full energy distrubtion for Y with findFit intervals
  ny, binsy, patchesy = s2.hist(targy, yinterval, density=True, facecolor="green", alpha=0.75,label="MiniScatter Target Y Histogram")

  #Add the "best fit" line using the findFit mu and sigma for Y and display sigma
  y2a = norm.pdf(binsy, muy, sigmay)
  l2a = s2.plot(binsy, y2a, "k--", linewidth=1,label="Histogram Least Square")
  y2b = norm.pdf(binsy, muy, Mty[0])
  l2b = s2.plot(binsy, y2b, "r--", linewidth=2,label="Target Twiss RMS")
  y2c = norm.pdf(binsy, muy, Mfy[0])
  l2c = s2.plot(binsy, y2c, "b--", linewidth=2,label=fitlabel+" Twiss RMS")

  sigmatextx = r"$\sigma_{X}$:"
  sigmatextx +="\nLeast Square = "+"{:.2f}".format(sigmax)+"mm"
  sigmatextx +="\nTwiss RMS = "+"{:.2f}".format(Mtx[0])+"mm"
  sigmatextx +="\n"+fitlabel+" RMS = "+"{:.2f}".format(Mfx[0])+"mm"
  sigmatexty = r"$\sigma_{Y}$:"
  sigmatexty +="\nLeast Square = "+"{:.2f}".format(sigmay)+"mm"
  sigmatexty +="\nTwiss RMS = "+"{:.2f}".format(Mty[0])+"mm"
  sigmatexty +="\n"+fitlabel+" RMS = "+"{:.2f}".format(Mfy[0])+"mm"

  s1.set_title(fitlabel+" at Target after "+mat+" PBW, X Distribution",fontsize=fs) #+"\n"+rf"$\sigma_D=${{:.1f}}mm".format(sigmax)
  s1.set_xlabel("X Position [mm]",fontsize=fs)
  s1.set_ylabel("Probability Density",fontsize=fs)
  xlim = 40
  s1.set_xlim([-xlim,xlim])
  xlim1 = s1.get_xlim()
  ylim1 = s1.get_ylim()
  s1.text(xlim1[0]*0.97,ylim1[1]*0.7,sigmatextx,fontsize=fs-2)
  handles1, labels1 = plt.gca().get_legend_handles_labels()
  #print(labels1)
  order1=[0,1,2]
  #order1 = [2,0,1]#for no Eq fits for initial Hist plot
  s1.legend([handles1[idx] for idx in order1],[labels1[idx] for idx in order1],fontsize=9)

  s2.set_xlim([-xlim,xlim])
  s2.set_title(fitlabel+" at Target after "+mat+" PBW, Y Distribution",fontsize=fs)
  s2.set_xlabel("Y Position [mm]",fontsize=fs)
  s2.set_ylabel("Probability Density",fontsize=fs)
  handles2, labels2 = plt.gca().get_legend_handles_labels()
  #print(labels2)
  order2=[0,1,2]
  #order2 = [2,0,1] #for no Eq fits for initial Hist plot
  s2.legend([handles2[idx] for idx in order2],[labels2[idx] for idx in order2],fontsize=9)
  ylim2=s2.get_ylim() #dynamically get the graph limits
  xlim2=s2.get_xlim()
  #print("xlim2: {:.1e}-{:.1e}, ylim2: {:.1e}-{:.1e}".format(xlim2[0],xlim2[1],ylim2[0],ylim2[1]))
  s2.text(xlim2[0]*0.97,ylim2[1]*0.7,sigmatexty,fontsize=fs-2)

  #Can date stamp the multi plot for easier tracking of changes, if necessary
  from datetime import datetime
  dt = datetime.now()

  name = savename+fitlabel.replace(' ','')+"Twiss_Target_"+dt.strftime("%H-%M-%S")+".png"##
  #plt.show()
  plt.savefig(name,bbox_inches="tight")
  plt.close()
  print(name)