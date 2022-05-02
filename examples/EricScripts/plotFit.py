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
  nx, binsx, patchesx = s1.hist(xs, xinterval, density=True, facecolor='green', alpha=0.75)
  
  #Add the 'best fit' line using the earlier findFit mus and sigmas for X
  y1 = norm.pdf(binsx, mux, sigmax)
  l1 = s1.plot(binsx, y1, 'r--', linewidth=2) #important to keep it as l# or it won't work

  #Make the histogram of the full energy distrubtion for Y with findFit intervals
  ny, binsy, patchesy = s2.hist(ys, yinterval, density=True, facecolor='green', alpha=0.75)

  #Add the 'best fit' line using the earlier findFit mus and sigmas for Y
  y2 = norm.pdf(binsy, muy, sigmay)
  l2 = s2.plot(binsy, y2, 'r--', linewidth=2)

  #change limits based on received parameters
  if ylim !=0: #for halo zoom in
    s1.set_ylim([0,ylim])
    s2.set_ylim([0,ylim])
    savename=savename+'_halo'
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
  s1.set_title("After "+mat+", Fitted X Distribution \n"+rf'$\mu= {{:.1e}}, \sigma = {{:.3f}}$,'.format(mux,sigmax)+
          " {:.3f}% outisde".format(pOut3sigx)+r" 3$\sigma$",fontsize=fs)
  s1.set_xlabel("X Position [mm]",fontsize=fs)
  s1.set_ylabel("Probability Density",fontsize=fs)

  s2.set_title("After "+mat+", Fitted Y Distribution \n"+rf'$\mu= {{:.1e}}, \sigma = {{:.3f}}$,'.format(muy,sigmay)+
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
  #import gaussian
  import numpy as np
  import matplotlib.pyplot as plt
  from scipy.optimize import curve_fit

  def gaussian(x,mu,amplitude,sigma):
    #import numpy.exp as exp
    return amplitude * np.exp( - (x - mu) ** 2 / (2*sigma ** 2))

  plt.close()
  bin_heights, bin_borders, _ = plt.hist(data,bins='auto',label='histogram')
  #print(len(bin_borders))
  bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
  popt,_ = curve_fit(gaussian, bin_centers, bin_heights,p0=[0.1,1.,1.]) #p0 should be good start, though not what is recommended
  x_interval_for_fit = np.linspace(bin_borders[0],bin_borders[-1],len(bin_borders))
  #plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit,*popt),label='fit')
  #plt.legend()
  #plt.xlim([bin_borders[0],bin_borders[-1]])
  #plt.show()
  plt.close() #be sure to close this or else these show up in the multi plots!
  #print(popt)
  #return the mean, abs(sigma), interval #sigma can sometimes be - so must abs() it.
  return(popt[0],abs(popt[2]),x_interval_for_fit) #for Vac and Air sigma sometimes give - number...

def plotTwissFit(Twissx,Twissy,xs,ys,savename,mat,PDFlabel,Twiss2x,Twiss2y,titledescr):
  #For plotting particle Distribution at Target Exit with PDF from Formalism Twiss parameters.
  import numpy as np
  import matplotlib.pyplot as plt
  from plotFit import findFit
  from scipy.stats import norm 

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

  sigmaTwx = np.sqrt(Twissx[0]*Twissx[2]) /mm #sigma rms = beam size = sqrt(beta*epsN)
  muTwx = Twissx[0]*Twissx[3] /mm            #posVar = beta/epsG 
  sigmaTwy = np.sqrt(Twissy[0]*Twissy[2]) /mm
  muTwy = Twissy[0]*Twissy[3] /mm
  print('Betax:{:.2f}m, alphax: {:.2f}, Norm Emitt x: {:.3e}m, Geo Emitt x: {:.3e}m'.format(Twissx[0],Twissx[1],Twissx[2],Twissx[3]))
  print('Betay:{:.2f}m, alphay: {:.2f}, Norm Emitt y: {:.3e}m, Geo Emitt y: {:.3e}m'.format(Twissy[0],Twissy[1],Twissy[2],Twissy[3]))
  print('muTwx: {:.2e}mm, muTwy: {:.2e}mm, sigmaTwx; {:.2f}mm, sigmaTwy: {:.2f}mm'.format(muTwx,muTwy,sigmaTwx,sigmaTwy))

  mux,sigmax,xinterval = findFit(xs)
  muy,sigmay,yinterval = findFit(ys)
  #print(mux,muy,sigmax,sigmay)

  #Create the fig with 2 plots side by side
  plt.clf()
  fig = plt.figure(figsize=(15,6.5))
  plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
  s1 = fig.add_subplot(1,2,1)
  s2 = fig.add_subplot(1,2,2)
  #Make the histogram of the full energy distrubtion for X with findFit intervals
  nx, binsx, patchesx = s1.hist(xs, xinterval, density=True, facecolor='green', alpha=0.75,label='MiniScatter Distribution')

  #Add the 'best fit' line using the earlier findFit mus and sigmas for X
  y1 = norm.pdf(binsx, muTwx, sigmaTwx)
  l1 = s1.plot(binsx, y1, 'r--', linewidth=2,label=PDFlabel+' Twiss PDF') #important to keep it as l# or it won't work
  mutextx +="\n"+PDFlabel+" = "+"{:.2f}".format(muTwx)+"mm"
  sigmatextx +="\n"+PDFlabel+" = "+"{:.2f}".format(sigmaTwx)+"mm"

  #Make the histogram of the full energy distrubtion for Y with findFit intervals
  ny, binsy, patchesy = s2.hist(ys, yinterval, density=True, facecolor='green', alpha=0.75,label='MiniScatter Distribution')

  #Add the 'best fit' line using the earlier findFit mus and sigmas for Y
  y2 = norm.pdf(binsy, muTwy, sigmaTwy)
  l2 = s2.plot(binsy, y2, 'r--', linewidth=2,label=PDFlabel+' Twiss PDF')
  mutexty +="\n"+PDFlabel+" = "+"{:.2f}".format(muTwy)+"mm"
  sigmatexty +="\n"+PDFlabel+" = "+"{:.2f}".format(sigmaTwy)+"mm"

  #if the 2nd Twiss set != first Twiss set, plot the 2nd Twiss set PDF, as 2nd is the real PBW Exit Twiss
  if Twissx[0] != Twiss2x[0]:
    sigmaTw2x = np.sqrt(Twiss2x[0]*Twiss2x[2]) /mm #sigma rms = beam size = sqrt(beta*epsN)
    muTw2x = Twiss2x[0]*Twiss2x[3] /mm            #posVar = beta/epsG 
    sigmaTw2y = np.sqrt(Twiss2y[0]*Twiss2y[2]) /mm
    muTw2y = Twiss2y[0]*Twiss2y[3] /mm


    y12 = norm.pdf(binsx,muTw2x,sigmaTw2x)
    l12 = s1.plot(binsx, y12, 'b--', linewidth=1,label='MiniScatter Twiss PDF')
    mutextx +="\n"+"MiniScatter = "+"{:.2f}".format(muTw2x)+"mm"
    sigmatextx +="\n"+"MiniScatter = "+"{:.2f}".format(sigmaTw2x)+"mm"

    y22 = norm.pdf(binsy, muTw2y, sigmaTw2y)
    l22 = s2.plot(binsy, y22, 'b--', linewidth=1,label='MiniScatter Twiss PDF')
    mutexty +="\n"+"MiniScatter = "+"{:.2f}".format(muTw2y)+"mm"
    sigmatexty +="\n"+"MiniScatter = "+"{:.2f}".format(sigmaTw2y)+"mm"

  #Set Plot characterists
  s1.set_xlim([-25,25])
  s2.set_xlim([-25,25])

  xlim1=s1.get_xlim()
  ylim1=s1.get_ylim() #dynamically get the ylimits
  s1.text(xlim1[0]*0.97,ylim1[1]*0.85,mutextx)
  s1.text(xlim1[0]*0.97,ylim1[1]*0.70,sigmatextx,fontsize=10) #set the texts at 3/4 and 1/2 of ylim

  xlim2=s2.get_xlim()
  ylim2=s2.get_ylim() #dynamically get the ylimits
  s2.text(xlim2[0]*0.97,ylim2[1]*0.85,mutexty)
  s2.text(xlim2[0]*0.97,ylim2[1]*0.70,sigmatexty,fontsize=10) #set the texts at 3/4 and 1/2 of ylim

  import re
  x = re.search('init',savename)
  if x:
    titlestart = 'Before & After '+ mat + ' PBW'
    print('\nPBW Exit Twiss')
    print('Betax:{:.2f}m, alphax: {:.2f}, Norm Emitt x: {:.3e}m, Geo Emitt x: {:.3e}m'.format(Twiss2x[0],Twiss2x[1],Twiss2x[2],Twiss2x[3]))
    print('Betay:{:.2f}m, alphay: {:.2f}, Norm Emitt y: {:.3e}m, Geo Emitt y: {:.3e}m'.format(Twiss2y[0],Twiss2y[1],Twiss2y[2],Twiss2y[3]))
    print('muTwx: {:.2e}mm, muTwy: {:.2e}mm, sigmaTwx; {:.2f}mm, sigmaTwy: {:.2f}mm'.format(muTw2x,muTw2y,sigmaTw2x,sigmaTw2y))
  else:
    titlestart = 'At '+ mat + ' ' + titledescr
  s1.set_title(titlestart+" Twiss Fit X Distribution \n"+rf'$\mu_{{}}= {{:.0e}}$mm, $\sigma_{{}}= {{:.4f}}$mm'.format('D',mux,'D',sigmax),fontsize=fs)
  s1.set_xlabel("X Position [mm]",fontsize=fs)
  s1.set_ylabel("Probability Density",fontsize=fs)
  s1.legend(fontsize=9)

  s2.set_title(titlestart+" Twiss Fit Y Distribution \n"+rf'$\mu_{{}}= {{:.0e}}$mm, $\sigma_{{}}= {{:.4f}}$mm'.format('D',muy,'D',sigmay),fontsize=fs)
  s2.set_xlabel("Y Position [mm]",fontsize=fs)
  s2.set_ylabel("Probability Density",fontsize=fs)
  s2.legend(fontsize=9)

  twisstext = r"$\mu_F =  \epsilon_G \beta$"+'\n'+r"$\sigma_F = \sqrt{\epsilon_N \beta}$"#+"\n from Eq 8"
  ylim1=s1.get_ylim() #dynamically get the ylimits
  xlim1=s1.get_xlim()
  s1.text(xlim1[1]*1.015,ylim1[1]*1.01,twisstext,fontsize=fs-1) #set the texts at 3/4 and 1/2 of ylim

  #Can date stamp the multi plot for easier tracking of changes, if necessary
  from datetime import datetime
  dt = datetime.now()

  name = savename+"_Twiss"+dt.strftime("%H-%M-%S")+".png"  #
  #plt.show()
  plt.savefig(name)
  plt.close()
  print(name)