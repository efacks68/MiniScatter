#plotFit.py
#Eric Fackelman
#29 March- 5 April 2022

#This is to plot the SimpledemoPBWscript histograms 
# Fits are made in findFit function using scipy.optimize.curve_fit found parameters
# And passed back into plotFit.

def plotFit(xs,ys,savename,xlim,ylim,material): 
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
  s1.set_title("At "+mat+" PBW Exit, Fitted X Distribution \n"+rf"$\mu= {{:.1e}}, \sigma = {{:.3f}}$,".format(mux,sigmax)+
          " {:.3f}% outisde".format(pOut3sigx)+r" 3$\sigma$",fontsize=fs)
  s1.set_xlabel("X Position [mm]",fontsize=fs)
  s1.set_ylabel("Probability Density",fontsize=fs)

  #s2.set_title("After "+mat+", Fitted Y Distribution \n"+rf"$\mu= {{:.1e}}, \sigma = {{:.3f}}$,".format(muy,sigmay)+
  #        " {:.3f}% outisde".format(pOut3sigy)+r" 3$\sigma$",fontsize=fs)
  s2.set_title("At "+mat+" PBW Exit, Fitted Y Distribution \n"+rf"$\mu= {{:.1e}}, \sigma = {{:.3f}}$,".format(muy,sigmay)+
          " {:.3f}% outisde".format(pOut3sigy)+r" 3$\sigma$",fontsize=fs)
  s2.set_xlabel("Y Position [mm]",fontsize=fs)
  s2.set_ylabel("Probability Density",fontsize=fs)

  #Can date stamp the multi plot for easier tracking of changes, if necessary
  #from datetime import datetime
  #dt = datetime.now()

  #Update Name
  #name = savename+"_"+dt.strftime("%H-%M-%S")+".png"  
  name = savename+".png" #I'd prefer to have it overwrite old files for now
  print(name) #show so it is easy to find the file
  plt.savefig(name)
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

def plotTwissFit(Twissx,Twissy,xs,pxs,ys,pys,savename,mat,PDFlabel,Twiss2x,Twiss2y,titledescr):
  #For plotting particle Distribution at Target Exit with PDF from Twiss parameters.
  import numpy as np
  import matplotlib.pyplot as plt
  from plotFit import findFit, calcTwiss
  from scipy.stats import norm 

  mm = 1e-3 #[m]
  #Recalculate Twiss
  calcxTwf = calcTwiss("Recheck X","Recheck X'",xs*mm,pxs)
  calcyTwf = calcTwiss("Recheck Y","Recheck Y'",ys*mm,pys)

  #Define material short name for plot labeling
  if mat == "G4_Galactic":
    mat = "Vacuum"
  elif mat == "G4_Al":
    mat ="Al"
  elif mat == "G4_AIR":
    mat ="Air"
  elif mat == "G4_Au":
    mat = "Au"

  mm = 1e-3 #the position numbers are in mm so must convert the PDF parameters to mm!
  fs = 14 #set the axis label font size early
  sigmatextx = r"$\sigma_{Fx}$:"
  sigmatexty = r"$\sigma_{Fy}$:"
  mutextx = r"$\mu_{Fx}$:"
  mutexty = r"$\mu_{Fy}$:"

  #Twiss1x = [betax,alphax,epsx,EepsGx]
  print("Sent Twiss")
  sigmaTwx = np.sqrt(Twissx[0]*Twissx[2]) /mm #sigma rms = beam size = sqrt(beta*epsN)
  muTwx = Twissx[0]*Twissx[3] /mm            #posVar = beta/epsG 
  sigmaTwy = np.sqrt(Twissy[0]*Twissy[2]) /mm
  muTwy = Twissy[0]*Twissy[3] /mm
  print("Betax:{:.2f}m, alphax: {:.2f}, Norm Emitt x: {:.3e}m, Geo Emitt x: {:.3e}m".format(Twissx[0],Twissx[1],Twissx[2],Twissx[3]))
  print("Betay:{:.2f}m, alphay: {:.2f}, Norm Emitt y: {:.3e}m, Geo Emitt y: {:.3e}m".format(Twissy[0],Twissy[1],Twissy[2],Twissy[3]))
  print("muTwx: {:.2e}mm, muTwy: {:.2e}mm, sigmaTwx; {:.2f}mm, sigmaTwy: {:.2f}mm".format(muTwx,muTwy,sigmaTwx,sigmaTwy))

  #Calculate for calculated Twiss and plot to try to find why different. 10.5
  print("Recalculated Twiss")
  csigmaTwx = np.sqrt(calcxTwf[0]*calcxTwf[2]) /mm #sigma rms = beam size = sqrt(beta*epsN)
  cmuTwx = calcxTwf[0]*calcxTwf[3] /mm            #posVar = beta/epsG 
  csigmaTwy = np.sqrt(calcyTwf[0]*calcyTwf[2]) /mm
  cmuTwy = calcyTwf[0]*calcyTwf[3] /mm
  print("Betax:{:.2f}m, alphax: {:.2f}, Norm Emitt x: {:.3e}m, Geo Emitt x: {:.3e}m".format(calcxTwf[0],calcxTwf[1],calcxTwf[2],calcxTwf[3]))
  print("Betay:{:.2f}m, alphay: {:.2f}, Norm Emitt y: {:.3e}m, Geo Emitt y: {:.3e}m".format(calcyTwf[0],calcyTwf[1],calcyTwf[2],calcyTwf[3]))
  print("muTwx: {:.2e}mm, muTwy: {:.2e}mm, sigmaTwx; {:.2f}mm, sigmaTwy: {:.2f}mm".format(muTwx,muTwy,sigmaTwx,sigmaTwy))

  #Get intervals and info about distribution
  mux,sigmax,xinterval = findFit(xs)
  muy,sigmay,yinterval = findFit(ys)
  #print(mux,muy,sigmax,sigmay)

  #Create the fig with 2 plots side by side
  plt.clf()
  fig = plt.figure(figsize=(15,6.0))
  plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
  s1 = fig.add_subplot(1,2,1)
  s2 = fig.add_subplot(1,2,2)
  #Make the histogram of the full energy distrubtion for X with findFit intervals
  nx, binsx, patchesx = s1.hist(xs, xinterval, density=True, facecolor="green", alpha=0.75,label="MiniScatter Distribution")

  #Add the "best fit" line using the findFit mu and sigma for X and display mu and sigma
  y1a = norm.pdf(binsx, cmuTwx, csigmaTwx)
  l1a = s1.plot(binsx, y1a, "r--", linewidth=2,label=PDFlabel+" Twiss PDF") #important to keep it as l# or it won't work
  y1b = norm.pdf(binsx, mux, sigmax)
  l1b = s1.plot(binsx, y1b, "k--", linewidth=1,label=PDFlabel+" Dist PDF")
  mutextx +="\n"+"Distribution = "+"{:.2f}".format(mux)+"mm"
  sigmatextx +="\n"+"Distribution = "+"{:.2f}".format(sigmax)+"mm"
  mutextx +="\n"+PDFlabel+" calc = "+"{:.2f}".format(cmuTwx)+"mm"
  sigmatextx +="\n"+PDFlabel+" calc = "+"{:.2f}".format(csigmaTwx)+"mm"

  #Make the histogram of the full energy distrubtion for Y with findFit intervals
  ny, binsy, patchesy = s2.hist(ys, yinterval, density=True, facecolor="green", alpha=0.75,label="MiniScatter Distribution")

  #Add the "best fit" line using the findFit mu and sigma for Y and display mu and sigma
  y2 = norm.pdf(binsy, cmuTwy, csigmaTwy)
  l2 = s2.plot(binsy, y2, "r--", linewidth=2,label=PDFlabel+" Twiss PDF")
  y2b = norm.pdf(binsx, muy, sigmay)
  l2b = s2.plot(binsx, y2b, "k--", linewidth=1,label=PDFlabel+" Dist PDF")
  mutexty +="\n"+"Distribution = "+"{:.2f}".format(muy)+"mm"
  sigmatexty +="\n"+"Distribution = "+"{:.2f}".format(sigmay)+"mm"
  mutexty +="\n"+PDFlabel+" calc = "+"{:.2f}".format(cmuTwy)+"mm"
  sigmatexty +="\n"+PDFlabel+" calc = "+"{:.2f}".format(csigmaTwy)+"mm"

  #if the 2nd Twiss set != first Twiss set, plot the 2nd Twiss set PDF, as 2nd is the real PBW Exit Twiss
  if Twissx[0] != Twiss2x[0]:
    sigmaTw2x = np.sqrt(Twiss2x[0]*Twiss2x[2]) /mm #sigma rms = beam size = sqrt(beta*epsN)
    muTw2x = Twiss2x[0]*Twiss2x[3] /mm            #posVar = beta/epsG 
    sigmaTw2y = np.sqrt(Twiss2y[0]*Twiss2y[2]) /mm
    muTw2y = Twiss2y[0]*Twiss2y[3] /mm

    y12 = norm.pdf(binsx,muTw2x,sigmaTw2x)
    l12 = s1.plot(binsx, y12, "b--", linewidth=1,label="MiniScatter Twiss PDF")
    mutextx +="\n"+"MiniScatter = "+"{:.2f}".format(muTw2x)+"mm"
    sigmatextx +="\n"+"MiniScatter = "+"{:.2f}".format(sigmaTw2x)+"mm"

    y22 = norm.pdf(binsy, muTw2y, sigmaTw2y)
    l22 = s2.plot(binsy, y22, "b--", linewidth=1,label="MiniScatter Twiss PDF")
    mutexty +="\n"+"MiniScatter = "+"{:.2f}".format(muTw2y)+"mm"
    sigmatexty +="\n"+"MiniScatter = "+"{:.2f}".format(sigmaTw2y)+"mm"

  #Set Plot characterists
  #ylim=5/(10**(5+0)) #change if not N=10^5
  #s1.set_ylim([0,ylim])
  #s2.set_ylim([0,ylim])
  s1.set_xlim([-25,25])
  s2.set_xlim([-25,25])

  xlim1=s1.get_xlim()
  ylim1=s1.get_ylim() #dynamically get the ylimits
  s1.text(xlim1[0]*0.97,ylim1[1]*0.85,mutextx)
  s1.text(xlim1[0]*0.97,ylim1[1]*0.65,sigmatextx,fontsize=10) #set the texts at fraction of ylim

  xlim2=s2.get_xlim()
  ylim2=s2.get_ylim() #dynamically get the ylimits
  s2.text(xlim2[0]*0.97,ylim2[1]*0.85,mutexty)
  s2.text(xlim2[0]*0.97,ylim2[1]*0.65,sigmatexty,fontsize=10)

  #If this is the instance before PBW, change title and display Twiss of PBW Exit
  import re
  x = re.search("init",savename)
  if x:
    titlestart = "Before & After "+ mat + " PBW"
    print("\nPBW Exit Twiss")
    print("Betax:{:.2f}m, alphax: {:.2f}, Norm Emitt x: {:.3e}m, Geo Emitt x: {:.3e}m".format(Twiss2x[0],Twiss2x[1],Twiss2x[2],Twiss2x[3]))
    print("Betay:{:.2f}m, alphay: {:.2f}, Norm Emitt y: {:.3e}m, Geo Emitt y: {:.3e}m".format(Twiss2y[0],Twiss2y[1],Twiss2y[2],Twiss2y[3]))
    print("muTwx: {:.2e}mm, muTwy: {:.2e}mm, sigmaTwx; {:.2f}mm, sigmaTwy: {:.2f}mm".format(muTw2x,muTw2y,sigmaTw2x,sigmaTw2y))
  else:
    titlestart = "At "+ mat + " " + titledescr

  #finish setting plot characteristics
  s1.set_title(titlestart+", Twiss Fit X Distribution",fontsize=fs)
  s1.set_xlabel("X Position [mm]",fontsize=fs)
  s1.set_ylabel("Probability Density",fontsize=fs)
  s1.legend(fontsize=9)

  s2.set_title(titlestart+", Twiss Fit Y Distribution",fontsize=fs)
  s2.set_xlabel("Y Position [mm]",fontsize=fs)
  s2.set_ylabel("Probability Density",fontsize=fs)
  s2.legend(fontsize=9)

  #Display formulas for mu and sigma calculations
  twisstext = r"$\mu_F =  \epsilon_G \beta$"+"\n"+r"$\sigma_F = \sqrt{\epsilon_N \beta}$"#+"\n from Eq 8"
  ylim1=s1.get_ylim() #dynamically get the ylimits
  xlim1=s1.get_xlim()
  s1.text(xlim1[1]*1.015,ylim1[1]*1.01,twisstext,fontsize=fs-2) #set the texts at 3/4 and 1/2 of ylim

  #Can date stamp the multi plot for easier tracking of changes, if necessary
  from datetime import datetime
  dt = datetime.now()

  name = savename+"_Twiss"+dt.strftime("%H-%M-%S")+".png"  #
  #plt.show()
  plt.savefig(name)
  plt.close()
  print(name)

def calcTwiss(labelx,labely,x,y):
  import numpy as np
  beta_rel = 0.948
  gamma_rel = 3.132
  #get Twiss from np.cov
  xcov = np.cov(x,y)
  #calculated Twiss parameters
  cemtg = np.sqrt(np.linalg.det(xcov)) 
  cbeta = xcov[0,0]/cemtg 
  calpha = -xcov[0,1]/cemtg
  cemtn = cemtg*beta_rel*gamma_rel
  #print(labelx,labely,":", xcov,"\nGeometric Emittance: {:.3e}m, Normalized Emittance: {:.3e}m, Beta: {:.3e}m, Alpha: {:.3e}\n".format(cemtg,cemtn,cbeta,calpha))
  return [cbeta,calpha,cemtn,cemtg]

def calcEq8(thetasq,Twiss,thick,beta_rel,gamma_rel):
  import numpy as np
  m=1
  um=1
  mm=1e-3

  #Twiss=[cbeta,calpha,cemtn,cemtg]
  #Calculations from Eq 7 and 8
  e8dnem = 0.5 * Twiss[0]*m * thetasq * thetasq #[m*rad^2]
  e8alph = Twiss[2] * Twiss[1] / (Twiss[2] + e8dnem)
  e8beta = Twiss[2] * Twiss[0]*m / (Twiss[2] + e8dnem)
  e8nemt = Twiss[2] + e8dnem
  e8gemt = e8nemt/(beta_rel*gamma_rel)

  return [e8beta,e8alph,e8nemt,e8gemt]

def calcEq16(thetasq,Twiss,thick,beta_rel,gamma_rel):
  import numpy as np
  m=1
  um=1
  mm=1e-3

  #Twiss=[cbeta,calpha,cemtn,cemtg]
  #Calculations from Eq 15 and 16
  gamma0 = (1 + Twiss[1] * Twiss[1] ) / Twiss[0]*m
  e16dnem = 0.5 * thetasq * thetasq * (Twiss[0]*m + thick*mm * Twiss[1] + thick*mm*thick*mm/3 * gamma0) #[m*rad^2]
  e16alph = (Twiss[2] * Twiss[1] - thick*mm * 0.5 * thetasq) / (Twiss[2] + e16dnem)
  e16beta = (Twiss[2] * Twiss[0]*m + thick*mm * thick*mm / 3 * thetasq) / (Twiss[2] + e16dnem)
  e16nemt = Twiss[2] + e16dnem
  e16gemt = e16nemt/(beta_rel*gamma_rel)

  return [e16beta,e16alph,e16nemt,e16gemt]

def getMoments(Twiss):
  from numpy import sqrt as sq
  mm=1e-3

  sigma = sq(Twiss[0]*Twiss[3]) /mm #sigma rms = beam size = sqrt(beta*epsG)
  mu = Twiss[0]*Twiss[3] /mm #posVar = beta/epsG 
  sigmapx = sq(Twiss[3]*(1+Twiss[1]*Twiss[1])/Twiss[0]) #mrad actually
  #print(mu,sigma,sigmapx)
  return [mu,sigma,sigmapx]

def plotTwissFit2(Twiss,PDFlabel,xs,pxs,savename,mat,Twiss2,PDFlabel2,titledescr,axis,thick,thetasq,beta_rel,gamma_rel):
#def plotTwissFit2(Twissx,Twissy,xs,pxs,ys,pys,savename,mat,PDFlabel,Twiss2x,Twiss2y,titledescr):
  #For plotting particle Distribution at Target Exit with PDF from Twiss parameters.
  import numpy as np
  import matplotlib.pyplot as plt
  from plotFit import findFit, calcTwiss, getMoments
  from scipy.stats import norm 
  mm = 1e-3 #the position numbers are in mm so must convert the PDF parameters to mm!
  urad = 1e-6
  fs = 14

  #Twiss1x = [betax,alphax,epsx,EepsGx]
  #M = [mu,sigma,sigmapx]
  M1 = getMoments(Twiss)
  M2 = getMoments(Twiss2)

  #Recalculate Twiss = [cbeta,calpha,cemtn,cemtg]
  calcTwf = calcTwiss("Recheck "+axis,"Recheck "+axis+"'",xs*mm,pxs)
  #Calculate for calculated Twiss
  Mc = getMoments(calcTwf)
  Me8 = getMoments(calcEq8(thetasq,calcTwf,thick,beta_rel,gamma_rel))
  Me16 = getMoments(calcEq16(thetasq,calcTwf,thick,beta_rel,gamma_rel))
  print(Mc,Me8,Me16)

  print("Sent Twiss",PDFlabel,"Betax:{:.2f}m, alphax: {:.2f}, Norm Emitt x: {:.3e}m, Geo Emitt x: {:.3e}m".format(Twiss[0],Twiss[1],Twiss[2],Twiss[3]))
  print("sigmaTwx; {:.2f}mm, sigmaTwX': {:.2f}mrad".format(M1[1],M1[2])) 
  print("Recalculated Twiss","Betax:{:.2f}m, alphax: {:.2f}, Norm Emitt x: {:.3e}m, Geo Emitt x: {:.3e}m".format(calcTwf[0],calcTwf[1],calcTwf[2],calcTwf[3]))
  print("PDF")
  print("Calculated sigmaTwx: {:.2f}mm, sigmaTwX': {:.2e}rad".format(Mc[1],Mc[2]))
  print("Equation 8 sigmaTwx: {:.2f}mm, sigmaTwX': {:.2e}rad".format(Me8[1],Me8[2]))
  print("Equation 16 sigmaTwx: {:.2f}mm, sigmaTwX': {:.2e}rad".format(Me16[1],Me16[2]))

  #Get intervals and info about distribution
  mux,sigmax,xinterval = findFit(xs)
  mupx,sigmapx,pxinterval = findFit(pxs)
  print("Distribution sigmaTwx: {:.2f}mm, sigmaTwX': {:.2e}rad".format(sigmax,sigmapx))

  #Create the fig with 2 plots side by side
  plt.clf()
  fig = plt.figure(figsize=(15,6.0))
  plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
  s1 = fig.add_subplot(1,2,1)
  s2 = fig.add_subplot(1,2,2)
  #Make the histogram of the full energy distrubtion for X with findFit intervals
  nx, binsx, patchesx = s1.hist(xs, xinterval, density=True, facecolor="green", alpha=0.75,label="MiniScatter "+axis+" Distribution")

  #Add the "best fit" line using the findFit mu and sigma for x and display mu and sigma
  y1a = norm.pdf(binsx, mux, sigmax)
  l1a = s1.plot(binsx, y1a, "k--", linewidth=1,label="Filtered "+axis+" Dist PDF")
  y1b = norm.pdf(binsx, Mc[0], Mc[1])
  l1b = s1.plot(binsx, y1b, "r--", linewidth=2,label="Filtered "+axis+" Twiss PDF")
  y1c = norm.pdf(binsx, Me8[0], Me8[1])
  l1c = s1.plot(binsx, y1c, "bo", linewidth=2,label=PDFlabel+" "+axis+" Twiss PDF")
  y1d = norm.pdf(binsx, Me16[0], Me8[1])
  l1d = s1.plot(binsx, y1d, "y*", linewidth=1,label=PDFlabel2+" "+axis+" Twiss PDF")

  nx,binspx, patchespx = s2.hist(pxs,pxinterval,density=True, facecolor="green", alpha=0.75,label="MiniScatter "+axis+"' Distribution")
  y2a = norm.pdf(binspx, mupx, sigmapx)
  l2a = s2.plot(binspx, y2a, "k--", linewidth=1,label="Filtered "+axis+"' Dist PDF")
  y2b = norm.pdf(binspx, mupx, Mc[2])
  l2b = s2.plot(binspx, y2b, "r--", linewidth=2,label="Filtered "+axis+"' Twiss PDF")
  y2c = norm.pdf(binspx, mupx, Me8[2])
  l2c = s2.plot(binspx, y2c, "bo", linewidth=2,label=PDFlabel+" "+axis+"' Twiss PDF")
  y2d = norm.pdf(binspx, mupx, Me16[2])
  l2d = s2.plot(binspx, y2d, "y*", linewidth=1,label=PDFlabel2+" "+axis+"' Twiss PDF")

  titlestart = "At "+ mat + " " + titledescr
  #finish setting plot characteristics
  s1.set_title(titlestart+", Twiss Fit "+axis+" Distribution",fontsize=fs) #+"\n"+rf"$\sigma_D=${{:.1f}}mm".format(sigmax)
  s1.set_xlabel(axis+" Position [mm]",fontsize=fs)
  s1.set_ylabel("Probability Density",fontsize=fs)
  s1.legend(fontsize=9)
  ylim1=s1.get_ylim() #dynamically get the ylimits
  xlim1=s1.get_xlim()
  #text="Original:\n"+rf"$ \epsilon_N: {{:.1e}}$m, $\sigma_F:{{:.1f}}$mm".format(orig,sigmaTwx)+"\n"+"Modified:\n"+rf"$\epsilon_N:{{:.1e}}$m, $\sigma_F:${{:.1f}}mm".format(calcxTwf[2],csigmaTwx)
  if axis=="X":
    sigmatextx = r"$\sigma_{X}$:"
    sigmatextpx = r"$\sigma_{X'}$:"
  elif axis=="Y":
    sigmatextx = r"$\sigma_{Y}$:"
    sigmatextpx = r"$\sigma_{Y'}$:"
  sigmatextx +="\nDistribution = "+"{:.2f}".format(sigmax)+"mm"
  sigmatextx +="\nFiltered calc = "+"{:.2f}".format(Mc[1])+"mm"
  sigmatextx +="\n"+PDFlabel+" = "+"{:.2f}".format(Me8[1])+"mm"
  sigmatextx +="\n"+PDFlabel2+" = "+"{:.2f}".format(Me16[1])+"mm"
  sigmatextpx +="\nDistribution = "+"{:.2f}".format(sigmapx/mm)+"mrad"
  sigmatextpx +="\nFiltered calc = "+"{:.2f}".format(Mc[2]/mm)+"mrad"
  sigmatextpx +="\n"+PDFlabel+" = "+"{:.2f}".format(Me8[2]/mm)+"mrad"
  sigmatextpx +="\n"+PDFlabel2+" = "+"{:.2f}".format(Me16[2]/mm)+"mrad"
  s1.text(xlim1[0]*0.97,ylim1[1]*0.7,sigmatextx,fontsize=fs-2)

  ylim=5/(10**(5+0)) #change if not N=10^5
  #s2.set_ylim([0,ylim])
  s2.set_xlim(-sigmapx*8,sigmapx*8)
  s2.set_title(titlestart+", Twiss Fit  "+axis+"' Distribution",fontsize=fs)
  s2.set_xlabel(axis+"' Angle [rad]",fontsize=fs)
  s2.set_ylabel("Probability Density",fontsize=fs)
  s2.legend(fontsize=9)
  ylim2=s2.get_ylim() #dynamically get the ylimits
  xlim2=s2.get_xlim()
  s2.text(xlim2[0]*0.97,ylim2[1]*0.7,sigmatextpx,fontsize=fs-2)

  #Can date stamp the multi plot for easier tracking of changes, if necessary
  from datetime import datetime
  dt = datetime.now()

  name = savename+"_Twiss"+axis+"'"+dt.strftime("%H-%M-%S")  +".png"##
  #plt.show()
  plt.savefig(name)
  plt.close()
  print(name)



    #sigmaTwy = np.sqrt(Twissy[0]*Twissy[2]) /mm
  #muTwy = Twissy[0]*Twissy[3] /mm
  #sigmaTwpy = np.sqrt(Twissy[3]*(1+Twissy[1]*Twissy[1])/Twissy[0]) /mm
  #print("Betay:{:.2f}m, alphay: {:.2f}, Norm Emitt y: {:.3e}m, Geo Emitt y: {:.3e}m".format(Twissy[0],Twissy[1],Twissy[2],Twissy[3]))
  #calcyTwf = calcTwiss("Recheck Y","Recheck Y'",ys*mm,pys)
  #csigmaTwpy = np.sqrt(calcyTwf[3]*(1+calcyTwf[1]*calcyTwf[1])/calcyTwf[0])
  #print("Betay:{:.2f}m, alphay: {:.2f}, Norm Emitt y: {:.3e}m, Geo Emitt y: {:.3e}m".format(calcyTwf[0],calcyTwf[1],calcyTwf[2],calcyTwf[3]))
  #csigmaTwy = np.sqrt(calcyTwf[0]*calcyTwf[2]) /mm
  #cmuTwy = calcyTwf[0]*(calcyTwf[3]) /mm
  #muy,sigmay,yinterval = findFit(ys)
  #mupy,sigmapy,pyinterval = findFit(pys)
  #print(mux,muy,sigmax,sigmay)
  #print(sigmay,sigmapy)
  #, sigmaTwY: {:.2f}mm, sigmaTwY': {:.2f}mrad ,sigmaTwy,sigmaTwpy
