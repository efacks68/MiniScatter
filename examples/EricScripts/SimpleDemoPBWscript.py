#SimpleDemoPBWscript.py
#Eric Fackelman
#28-31 Mar 2022
#Uses the MiniScatter Python interface to model the ESS beam interaction with the PBW

# ## Code setup
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
from plotFit import plotFit,findFit
from simulation import simulation

# Simulation setup
#simulation(N,material,epsx,epsy,alphax,alphay,betax,betay,energy(float),zoff(float))

#Define something preliminarily
#materials = ["G4_Galactic","G4_AIR","G4_Al","G4_Au"] #full list as of 1.4.22
materials = ["G4_Al","G4_Au"] 
N=1e5 #number of particles
epsx = 0.113 #[um-mrad] #ESS beam at PBW: 0.113
betax = 1000 #[m] 1000
alphx = -50 #[mrad] -50
epsy = 0.121 #[um-mrad] 0.121
betay = 200 #[m] 200
alphy = -7 #[mrad] -7
ifplot=False #for plotting the 3 graphs per material
engplot = False

if len(sys.argv) < 2: #if no extra inputs, ask for them
  N = float(input("How many particles would you like to run? "))
  epsx = float(input("What is the Emittance in X? "))
  betax = float(input("What is the Beta in X? "))
  alphx = float(input("What is the Alpha in X? "))
  epsy = float(input("What is the Emittance in Y? "))
  betay = float(input("What is the Beta in Y? "))
  alphy = float(input("What is the Alpha in Y? "))
  if input("Would you like to graph everything? Yes=y, No=Enter "):
    ifplot=True
elif sys.argv[1] == 'ESS': #auto profiles
  #materials = ["G4_Galactic","G4_AIR","G4_Al","G4_Au"]
  materials = ["G4_Al","G4_Au"]
  N=1e5 #number of particles
  epsx = 0.113 #[um-mrad]
  betax = 1000 #[m]
  alphx = -50 #[mrad]
  epsy = 0.121 #[um-mrad]
  betay = 200 #[m]
  alphy = -7 #[mrad]
  ifplot=False #for plotting the 3 graphs per material
  if len(sys.argv) == 3 :
    N = float(sys.argv[2])
elif sys.argv[1] == 'pencil': #for a reasonable pencil beam
  materials = ["G4_Al","G4_Au"] #only use solids as Vac or Air does nothing
  N=1e5 #number of particles
  epsx = 1e-4 #[um-mrad] 
  betax = 1e-2 #[m] 
  alphx = 0 #[mrad] 
  epsy = 1e-4 #[um-mrad] 
  betay = 1e-2 #[m] 
  alphy = 0 #[mrad] 
  ifplot=False 
elif isfloat(sys.argv[1]) and isfloat(sys.argv[7]): #input variables are inputs
  N = sys.argv[1]
  epsx = sys.argv[2] #[um-mrad] 
  betax = sys.argv[3] #[m] 
  alphx = sys.argv[4] #[mrad] 
  epsy = sys.argv[5] #[um-mrad] 
  betay = sys.argv[6] #[m] 
  alphy = sys.argv[7] #[mrad] 
  if sys.arg[8]:
    ifplot = True #remember, sys.argv[8] anything produces plotting!

mag=math.floor(math.log10(N)) #for dynamically scaling the halo plots
#print(mag)

if ifplot:
  engplot=True

#Create the fig before the loop with 2 plots side by side
fig = plt.figure(figsize=(15,8))
plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
s1 = fig.add_subplot(1,2,1)
s2 = fig.add_subplot(1,2,2)
#plt.show()
fs = 14 #setting the axis label font size early

#Create dictionaries for sigma and % outside of 3 Sigma for each material used
sigmax = {}
sigmay = {}
pOut3sigx = {}
pOut3sigy = {}

#Start the text boxes for displaying info in the graphs
percenttextx = r"% outside 3$\sigma_x$:"
percenttexty = r"% outside 3$\sigma_y$:"
sigmatextx = r"$\sigma_x$:"
sigmatexty = r"$\sigma_y$:"

for material in materials:
  #function for preparing the run and running miniScatterDriver functions
  #savename,xexit,yexit= simulation( N,material,beam,thick,epsx ,epsy ,alphx,alphy,betax,betay,energy,zoff,Engcut,engplot):
  savename,xexit,yexit = simulation( N,material,"proton",    1,epsx ,epsy ,alphx,alphy,betax,betay,2000.0,-10.0, 0.95,engplot)
  #returns the savename and the x and y distributions of particle positions 
  #These returned arrays are from the MiniScatter detector, 5m after the PBW, ~where the ESS Target is.  
  
  #try to get rid of this:
  #Now plot the distributions with various views depending on the material
  if material == "G4_Galactic" or material == "G4_AIR":
    if ifplot:
      print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xexit),np.max(yexit)))
      xmax = math.ceil(np.max(xexit)/100)*100
      #For plotting x,y plots for core view and halo zoom-ins
      #plotFit(xs,    ys, savename,xlim,   ylim,material)
      #xlim<=10 => [-xlim*sigma,xlim*sigma]; xlim>10 => [-xlim,xlim]
      #ylim==0 => with plot; ylim!=0 => [0,ylim] (5 particles /(mag of N) is for halo)
      plotFit(xexit,yexit,savename,  3,      0,material) #3 sigma core
      plotFit(xexit,yexit,savename, xmax,5/(10**(mag+0)),material) #full range halo, with dynamic halo zoom
  elif material == "G4_Al":
    if ifplot:
      print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xexit),np.max(yexit)))
      xmax = math.ceil(np.max(xexit)/100)*100
      print(xmax)
      #plotFit(xs,    ys, savename,xlim,   ylim,material)
      plotFit(xexit,yexit,savename,  3,      0,material) #3 sigma core
      plotFit(xexit,yexit,savename, xmax,5/(10**(mag+0)),material) #full range halo
      plotFit(xexit,yexit,savename,  10,10/(10**(mag+0)),material) #10 sigma range halo
  elif material == "G4_Au":
    if ifplot:
      print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xexit),np.max(yexit)))
      xmax = math.ceil(np.max(xexit)/100)*100
      #plotFit(xs,    ys, savename,xlim,   ylim,material) 
      plotFit(xexit,yexit,savename,  3,      0,material) #3 sigma core
      plotFit(xexit,yexit,savename, xmax,15/(10**(mag+0)),material) #full range halo
      plotFit(xexit,yexit,savename,  10,15/(10**(mag+0)),material) #10 sigma range halo
  
  #Multi-Material Plot section, continues the material loop to plot data.
  #For creating the PDF of particle distribution from findFit function optimized mu and sigma.
  from scipy.stats import norm 
  
  #Define plotting characteristics depending on material
  if material == "G4_Galactic":
    mat = "Vacuum"
    color = 'magenta'
    dash = 'm--'
    tsp = 0.2 #transparency
  elif material == "G4_Al":
    mat ="Al"
    color = 'blue'
    dash = 'b--'
    tsp = 0.5
  elif material == "G4_AIR":
    mat ="Air"
    color = 'cyan'
    tsp = 0.2
    dash = 'c--'
  elif material == "G4_Au":
    mat = "Au"
    color = 'gold'
    tsp = 0.5
    dash = 'r--'

  #Use Scipy.optimize.curve_fit in my findFit function to get mus and sigmas:
  mux, sigmax[material], xinterval = findFit(xexit) #dynamically gets parameters AND histogram intervals!
  muy, sigmay[material], yinterval = findFit(yexit)

  #Find range of particles that are within 3 sigma for each material
  sigx=np.abs(xexit)>3*sigmax[material]# and np.abs(xexit)<10*sigma)
  sigy=np.abs(yexit)>3*sigmay[material]

  #Utilize dictionary to find % of particles within 3 sigma 
  pOut3sigx[material] = len(xexit[sigx])/len(xexit)*100 
  pOut3sigy[material] = len(yexit[sigy])/len(yexit)*100
  #print(material," gives ",pOut3sigx[material],"% outisde 3 sigma in X")
  #print(material," gives ",pOut3sigy[material],"% outisde 3 sigma in Y")

  #Update the texts to include this materials sigma and % outside 3 sigma
  percenttextx +="\n"+mat+" = "+"{:.2f}".format(pOut3sigx[material])+"%"
  percenttexty +="\n"+mat+" = "+"{:.2f}".format(pOut3sigy[material])+"%"
  sigmatextx +="\n"+mat+" = "+"{:.2f}".format(sigmax[material])+"mm"
  sigmatexty +="\n"+mat+" = "+"{:.2f}".format(sigmay[material])+"mm"

  #Make the histogram of the full energy distrubtion for X
  nx, binsx, patchesx = s1.hist(xexit, xinterval, density=True, facecolor=color, alpha=tsp,label=mat)
  
  #Add the 'best fit' line using the earlier norm.fit mus and sigmas for X
  y1 = norm.pdf(binsx, mux, sigmax[material])
  l1 = s1.plot(binsx, y1, dash, linewidth=1,label=mat) #important to keep it as l# or it won't work

  #Make the histogram of the full energy distrubtion for Y
  ny, binsy, patchesy = s2.hist(yexit, yinterval, density=True, facecolor=color, alpha=tsp,label=mat)
  
  #Add the 'best fit' line using the earlier norm.fit mus and sigmas for Y
  y2 = norm.pdf(binsy, muy, sigmay[material])
  l2 = s2.plot(binsy, y2, dash, linewidth=1,label=mat) #labeled by material short-name

#Set Plot characterists
s=3 # number of sigma width to plot
sigvals = sigmax.values() #get values of dictionary
xlim1 = s*max(sigvals) #use s*max as the xlim
s1.set_xlim([-xlim1,xlim1]) #seems to show all distributions, 30.3.22
s1.set_title("Fitted X Distributions",fontsize=fs)
s1.set_xlabel("X Position [mm]",fontsize=fs)
s1.set_ylabel("Probability Density",fontsize=fs)
s1.legend()
ylim1=s1.get_ylim() #dynamically get the ylimits
s1.text(-xlim1*0.95,ylim1[1]*0.75,sigmatextx) #set the texts at 3/4 and 1/2 of ylim
s1.text(-xlim1*0.95,ylim1[1]*0.5,percenttextx) #xlim is fixed
#s1.text(-xlim*0.95,ylim1[1]*0.99,)
#s1.text("Add Twiss Parameters! and have in savename as well!")

sigvals = sigmay.values() #get values of dictionary
xlim2 = s*max(sigvals) #use s*max as the xlim
s2.set_xlim([-xlim2,xlim2])
s2.set_title("Fitted Y Distributions",fontsize=fs)
s2.set_xlabel("Y Position [mm]",fontsize=fs)
s2.set_ylabel("Probability Density",fontsize=fs)
s2.legend()
ylim2=s2.get_ylim() #dynamically get the ylimits
s2.text(-xlim2*0.95,ylim2[1]*0.75,sigmatexty) #set the texts at 3/4 and 1/2 of ylim
s2.text(-xlim2*0.95,ylim2[1]*0.5,percenttexty) #x position is fixed
fig.suptitle(rf"Distributions at ESS Target of 10$^{{:d}}$ Protons".format(mag)+" Through Various Material PBWs\n"+
        rf"Initial beam of $\epsilon_x={{:.1e}}, \beta_x={{:.0e}}, \alpha_x={{:.1f}}$, ".format(epsx,betax,alphx) +
    rf" $\epsilon_y={{:.1e}}, \beta_y={{:.0e}}, \alpha_y={{:.1f}}$ ".format(epsy,betay,alphy),fontsize=18,y=0.99) #fontweight="bold",
#suptitle 2 sets of {{}} fix from "linuxtut.com" blog post

#Can date stamp the multi plot for easier tracking of changes, if necessary
from datetime import datetime
dt = datetime.now()
name = savename+"_"+dt.strftime("%H-%M-%S")+"_multi.png" #update the savename to not overwrite others
print(name) #show so it is easy to find the file
fig.savefig(name)
#plt.show()
plt.close() #be sure to close the plot
print("\n")
