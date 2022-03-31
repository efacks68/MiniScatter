#SimpleDemoPBWscript.py
#Eric Fackelman
#28-31 Mar 2022
#Uses the MiniScatter Python interface to model the ESS beam interaction with the PBW

# ## Code setup
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
from plotFit import plotFit
from simulation import simulation

# Simulation setup
#simulation(N,material,epsx,epsy,alphax,alphay,betax,betay,energy(float),zoff(float))
materials = ["G4_AIR","G4_Al","G4_Au"] #"G4_Galactic"
N=5e4 #number of particles
n=input("How many particles would you like to run?")
N=float(n)
mag=math.floor(math.log10(N)) #for dynamically scaling the halo plots
#print(mag)
ifplot=True #for plotting the 3 graphs per material

#Create the fig before the loop with 2 plots side by side
fig = plt.figure(figsize=(15,5))
plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
s1 = fig.add_subplot(1,2,1)
s2 = fig.add_subplot(1,2,2)
#plt.show()
fs = 14 #setting the axis label font size early

#Create dictionaries for % outside of 3 Sigma for each material used
pOut3sigx = {}
pOut3sigy = {}

#Start the text boxes for displaying info in the graphs
percenttextx = r"% outside 3$\sigma_x$:"
percenttexty = r"% outside 3$\sigma_y$:"
sigmatextx = r"$\sigma_x$:"
sigmatexty = r"$\sigma_y$:"
for material in materials:
  #function for preparing the run and running miniScatterDriver functions
  #savename,xexit,yexit= simulation(  N,material,epsx ,epsy ,alpx,alpy,betax,betay,energy,zoff,Engcut,engplot):
  savename,xexit,yexit = simulation( N,material,0.113,0.121, -50,  -7, 1000,  200,2000.0,-10.0, 0.95,   True)
  #returns the savename and the x and y distributions of particle positions 
  #These returned arrays are from the MiniScatter detector, 5m after the PBW, ~where the ESS Target is.  
  
  #Now plot the distributions with various views depending on the material
  if material == "G4_Galactic" or material == "G4_AIR":
      bins = 50
      cutx = 27 #since these distributions are still Gaussian, the range matches Al
      cuty = 15
      if ifplot:
        print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xexit),np.max(yexit)))
        #For plotting x,y plots for core view and halo zoom-ins
        #plotFit(xs,    ys, savename,bins,cutx,cuty,xlim,   ylim,material)
        #xlim==0 => -30,30; other=> [-xlim*sigma,xlim*sigma] 
        #ylim==0 => with plot; !=0 => [0,ylim] (0.00005 is for halo)
        plotFit(xexit,yexit,savename,bins,cutx,cuty,  30,      0,material) #standard core
        plotFit(xexit,yexit,savename,bins,cutx,cuty, 400,5/(10**(mag+1)),material) #full range halo, with dynamic halo zoom
        plotFit(xexit,yexit,savename,bins,cutx,cuty,  10,5/(10**(mag+1)),material) #10 sigma range halo
  elif material == "G4_Al":
      bins = 1000 #wider range requires more bins
      cutx = 27 #Distribution isn't Gaussian. This range seems to be best fit the spread
      cuty = 15
      if ifplot:
        print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xexit),np.max(yexit)))
        #plotFit(xs,    ys, savename,bins,cutx,cuty,xlim,   ylim,material)
        plotFit(xexit,yexit,savename,bins,cutx,cuty,  30,      0,material) #standard core
        plotFit(xexit,yexit,savename,bins,cutx,cuty, 500,5/(10**(mag+1)),material) #full range halo
        plotFit(xexit,yexit,savename,bins,cutx,cuty,  10,10/(10**(mag+1)),material) #10 sigma range halo
  elif material == "G4_Au":
      bins = 1000 
      cutx = 50 #Distribution isn't Gaussian. This range seems to be best fit the spread
      cuty = 40
      if ifplot:
        print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xexit),np.max(yexit)))
        #plotFit(xs,    ys, savename,bins,cutx,cuty,xlim,   ylim,material) 
        plotFit(xexit,yexit,savename,bins,cutx,cuty,  75,      0,material) #standard core
        plotFit(xexit,yexit,savename,bins,cutx,cuty, 520,15/(10**(mag+1)),material) #full range halo
        plotFit(xexit,yexit,savename,bins,cutx,cuty,  10,15/(10**(mag+1)),material) #10 sigma range halo
  
  #Multi-Material Plot section, continues the material loop to plot data.
  #For creating the fit to histogram of particle distribution
  from scipy.stats import norm
  
  #Based on the range that fits for the material, mask those particles
  maskx=np.abs(xexit)<cutx
  masky=np.abs(yexit)<cuty
  #print(len(xexit[maskx]))

  #Use the norm.fit function to get mus and sigmas:
  (mux, sigmax) = norm.fit(xexit[maskx])
  (muy, sigmay) = norm.fit(yexit[masky])

  #Find range of particles that are within 3 sigma
  sigx=np.abs(xexit)>3*sigmax# and np.abs(xexit)<10*sigma)
  sigy=np.abs(yexit)>3*sigmay

  #Utilize dictionary to find % of particles within 3 sigma 
  pOut3sigx[material] = len(xexit[sigx])/len(xexit)*100 
  pOut3sigy[material] = len(yexit[sigy])/len(yexit)*100
  #print(material," gives ",pOut3sigx[material],"% outisde 3 sigma in X")
  #print(material," gives ",pOut3sigy[material],"% outisde 3 sigma in Y")

  #Update the texts to include this materials % outside 3 sigma
  percenttextx +="\n"+material+" = "+"{:.2f}".format(pOut3sigx[material])+"%"
  percenttexty +="\n"+material+" = "+"{:.2f}".format(pOut3sigy[material])+"%"
  sigmatextx +="\n"+material+" = "+"{:.2f}".format(sigmax)+"mm"
  sigmatexty +="\n"+material+" = "+"{:.2f}".format(sigmay)+"mm"

  #Define parameters for plotting multiplot depending on material
  if material == "G4_Galactic":
    mat = "Vacuum"
    color = 'black'
    tsp = 0.2 #transparency
  elif material == "G4_Al":
    mat ="Al"
    color = 'blue'
    dash = 'b--'
    tsp = 0.5
  elif material == "G4_AIR":
    mat ="Air"
    color = 'green'
    tsp = 0.2
    dash = 'g--'
  elif material == "G4_Au":
    mat = "Au"
    color = 'gold'
    tsp = 0.5
    dash = 'r--'

  #Make the histogram of the full energy distrubtion for X
  nx, binsx, patchesx = s1.hist(xexit, bins, density=True, facecolor=color, alpha=tsp,label=mat)
  #print("the bins are",binsx)
  
  #Add the 'best fit' line using the earlier norm.fit mus and sigmas for X
  y1 = norm.pdf(binsx, mux, sigmax)
  l1 = s1.plot(binsx, y1, dash, linewidth=1,label=mat) #important to keep it as l# or it won't work

  #Make the histogram of the full energy distrubtion for Y
  ny, binsy, patchesy = s2.hist(yexit, bins, density=True, facecolor=color, alpha=tsp,label=mat)
  
  #Add the 'best fit' line using the earlier norm.fit mus and sigmas for Y
  y2 = norm.pdf(binsy, muy, sigmay)
  l2 = s2.plot(binsy, y2, dash, linewidth=1,label=mat) #labeled by material short-name

#Set Plot characterists
xlim = 40 #makes it easier to change range
s1.set_xlim([-xlim,xlim]) #seems to show all distributions, 30.3.22
s1.set_title("Fitted X Distributions",fontsize=fs)
s1.set_xlabel("X Position [mm]",fontsize=fs)
s1.set_ylabel("Probability Density",fontsize=fs)
s1.legend()
ylim1=s1.get_ylim() #dynamically get the ylimits
s1.text(-xlim+2,ylim1[1]*3/4,sigmatextx) #set the texts at 3/4 and 1/2 of ylim
s1.text(-xlim+2,ylim1[1]/2,percenttextx) #xlim is fixed

s2.set_xlim([-xlim,xlim])
s2.set_title("Fitted Y Distributions",fontsize=fs)
s2.set_xlabel("Y Position [mm]",fontsize=fs)
s2.set_ylabel("Probability Density",fontsize=fs)
s2.legend()
ylim2=s2.get_ylim() #dynamically get the ylimits
s2.text(-xlim+2,ylim2[1]*3/4,sigmatexty) #set the texts at 3/4 and 1/2 of ylim
s2.text(-xlim+2,ylim2[1]/2,percenttexty) #x position is fixed

#Can date stamp the multi plot for easier tracking of changes, if necessary
from datetime import datetime
dt = datetime.now()
name = savename+"_"+dt.strftime("%H-%M-%S")+"_multi.png" #update the savename to not overwrite others
print(name) #show so it is easy to find the file
fig.savefig(name)
#plt.show()
plt.close() #be sure to close the plot

