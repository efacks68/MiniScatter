#plotFit.py
#Eric Fackelman
#29 March 2022

#This is to plot the histograms with fits based on scipy.optimize.curve_fit parameters

def plotFit(xs,ys,savename,xlim,ylim,material): 
  #1.4.22 - removed cuts since using the findfit
  #For creating the fit to histogram of particle distribution
  from scipy.stats import norm
  import numpy as np
  import matplotlib.pyplot as plt
  from plotFit import findFit

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

  #Make the histogram of the full energy distrubtion for X
  nx, binsx, patchesx = s1.hist(xs, xinterval, density=True, facecolor='green', alpha=0.75)
  
  #Add the 'best fit' line using the earlier norm.fit mus and sigmas for X
  y1 = norm.pdf(binsx, mux, sigmax)
  l1 = s1.plot(binsx, y1, 'r--', linewidth=2) #important to keep it as l# or it won't work

  #Make the histogram of the full energy distrubtion for Y
  ny, binsy, patchesy = s2.hist(ys, yinterval, density=True, facecolor='green', alpha=0.75)

  #Add the 'best fit' line using the earlier norm.fit mus and sigmas for Y
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
  s1.set_title("After "+mat+", Fitted X Distribution \n"+rf'$\mu= {{:.3e}}, \sigma = {{:.3f}}$,'.format(mux,sigmax)+
          " {:.3f}% outisde".format(pOut3sigx)+r" 3$\sigma$",fontsize=fs)
  s1.set_xlabel("X Position [mm]",fontsize=fs)
  s1.set_ylabel("Probability Density",fontsize=fs)

  s2.set_title("After "+mat+", Fitted Y Distribution \n"+rf'$\mu= {{:.3e}}, \sigma = {{:.3f}}$,'.format(mux,sigmay)+
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
  print(len(bin_borders))
  bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
  popt,_ = curve_fit(gaussian, bin_centers, bin_heights,p0=[0.1,1.,1.]) #p0 should be good start, though not what is recommended
  x_interval_for_fit = np.linspace(bin_borders[0],bin_borders[-1],len(bin_borders))
  #plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit,*popt),label='fit')
  #plt.legend()
  #plt.xlim([-30,30])
  #plt.show()
  plt.close()
  #print(popt)
  return(popt[0],abs(popt[2]),x_interval_for_fit) #for Vac and Air sigma sometimes give - number...