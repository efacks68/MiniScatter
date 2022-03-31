#plotFit.py
#Eric Fackelman
#29 March 2022

#This is to plot the histograms with fits based on sent parameters

def plotFit(xs,ys,savename,bins,cutx,cuty,xlim,ylim,material):
  #For creating the fit to histogram of particle distribution
  from scipy.stats import norm
  import numpy as np
  import matplotlib.pyplot as plt

  fs = 14 #set the axis label font size early

  #Based on the range that fits for the material, mask those particles
  maskx=np.abs(xs)<cutx
  masky=np.abs(ys)<cuty
  #print(len(xs[maskx]))

  #Use the norm.fit function to get mus and sigmas:
  (mux, sigmax) = norm.fit(xs[maskx])
  (muy, sigmay) = norm.fit(ys[masky])

  #Find range of particles that are within 3 sigma
  sigx=np.abs(xs)>3*sigmax# and np.abs(xs)<10*sigma)
  sigy=np.abs(ys)>3*sigmay

  #Find % of particles within 3 sigma
  pOut3sigx = len(xs[sigx])/len(xs)*100
  pOut3sigy = len(ys[sigy])/len(ys)*100
  #print(pOut3sigx,"% outisde 3 sigma")
  #print(pOut3sigy,"% outisde 3 sigma")

  #Create the fig with 2 plots side by side
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
  nx, binsx, patchesx = s1.hist(xs, bins, density=True, facecolor='green', alpha=0.75)
  
  #Add the 'best fit' line using the earlier norm.fit mus and sigmas for X
  y1 = norm.pdf(binsx, mux, sigmax)
  l1 = s1.plot(binsx, y1, 'r--', linewidth=2) #important to keep it as l# or it won't work

  #Make the histogram of the full energy distrubtion for Y
  ny, binsy, patchesy = s2.hist(ys, bins, density=True, facecolor='green', alpha=0.75)

  #Add the 'best fit' line using the earlier norm.fit mus and sigmas for Y
  y2 = norm.pdf(binsy, muy, sigmay)
  l2 = s2.plot(binsy, y2, 'r--', linewidth=2)

  #change limits based on received parameters
  if xlim <= 10:
    s1.set_xlim([-xlim*sigmax,xlim*sigmax])
    s2.set_xlim([-xlim*sigmax,xlim*sigmax])
    savename = savename + "_halo_" + str(xlim) + "sigma"
  elif xlim > 100: #allows dynamic, expected case of ~4-500
    s1.set_xlim([-xlim,xlim])
    s2.set_xlim([-xlim,xlim])
    savename = savename + "_halo_extended"
  else: #for standard case
    s1.set_xlim([-xlim,xlim])
    s2.set_xlim([-xlim,xlim])
  if ylim !=0: #for halo zoom in
    s1.set_ylim([0,ylim])
    s2.set_ylim([0,ylim])

  #Set Plot characterists
  s1.set_title("After "+mat+r", Fitted X Distribution within $\pm$ {:d}mm".format(cutx)+
          "\n"+r'$\mu= $'+"{:.3e}, ".format(mux)+r'$\sigma = $'+"{:.3f},".format(sigmax)+
          " {:.3f}% outisde".format(pOut3sigx)+r" 3$\sigma$",fontsize=fs)
  s1.set_xlabel("X Position [mm]",fontsize=fs)
  s1.set_ylabel("Probability Density",fontsize=fs)

  s2.set_title("After "+mat+r", Fitted Y Distribution within $\pm$ {:d}mm".format(cuty)+
          "\n"+r'$\mu= $'+"{:.3e}, ".format(muy)+r'$\sigma = $'+"{:.3f},".format(sigmay)+
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
