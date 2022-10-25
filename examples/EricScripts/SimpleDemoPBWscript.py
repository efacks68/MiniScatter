#SimpleDemoPBWscript.py
#Eric Fackelman
#28-31 Mar 2022
#Uses the MiniScatter Python interface to model the ESS beam interaction with the PBW

# ## Code setup
import numpy as np
import matplotlib.pyplot as plt
import argparse,sys
import math
from plotFit import plotFit,findFit
from simulation import simulation

# Simulation setup
#simulation(N,material,nemtx,nemty,alphax,alphay,betax,betay,energy(float),zoff(float))

#Define something preliminarily
#materials = ["G4_Galactic","G4_AIR","G4_Al","G4_Au"] #full list as of 1.4.22
materials = ["G4_Al"] 
thick=0 #new default
N=1e5 #number of particles
betax = 941.25 #[m]
alphx = -58.81
nemtx = 0.11315 #[um-mrad]
betay = 120.03 #[m]
alphy = -7.46
nemty = 0.12155 #[um-mrad]
ifplot=False #for plotting the 3 graphs per material
engplot = False
energy = 570 #[MeV]
engcut = 95.0
zoff = "*-50" #[mm] with preappended * to keep covar defined at z=0
beamAngle = 0 #[mm] #0 for no angle ~? is nominal
loadParts = False
beamFile=""

#Set up Argument Parsing
parser = argparse.ArgumentParser()
beam = parser.add_mutually_exclusive_group()
beam.add_argument("--ESS",action="store_true")
beam.add_argument("--pencil",action="store_true")
parser.add_argument("--l",type=str,help="Load Particles or not")
parser.add_argument("--t",type=float,help="PBW Thickness, 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al X [mm] thick")
parser.add_argument("--twiss",type=float,nargs=9,help="N particles,BetaX,AlphX,NemtX,BetaY,AlphY,NemtY,Thickness")
parser.add_argument("--angle",type=float,default=0,help="Beam Angle, as calculated and used in aRasterMaker.py")
args = parser.parse_args()
print(args)

if args.l:
  loadParts=True
if args.l != "":
  beamFile = args.l

thick = args.t
if thick == 0:
  materials = ["G4_Al"]
elif thick == 0.1:
  materials = ["G4_Galactic"]
if args.angle != 0:
  beamAngle = angle

if len(sys.argv) < 2: #if no extra inputs, ask for them
  N = float(input("How many particles would you like to run? "))
  betax = float(input("What is the Beta in X? "))
  alphx = float(input("What is the Alpha in X? "))
  nemtx = float(input("What is the Emittance in X? "))
  betay = float(input("What is the Beta in Y? "))
  alphy = float(input("What is the Alpha in Y? "))
  nemty = float(input("What is the Emittance in Y? "))
  if input("Would you like to graph everything? Yes=y, No=Enter "):
    ifplot=True
elif args.ESS: #auto profiles
  #materials = ["G4_Galactic","G4_Al","G4_Au"]
  #materials = ["G4_Al","G4_Galactic"]#,"G4_Au"]
  #materials = ["G4_Galactic"]
  ifplot=False #for plotting the 3 graphs per material
  #if len(sys.argv) == 4 :
  #  N = float(sys.argv[3])
  N=1e5 #number of particles
  #updated Twiss on 13 June 2022 from OpenXAL script
  betax = 941.25 #[m] original:1000m
  alphx = -58.81 #original: -50
  nemtx = 0.11315 #[mm-mrad]
  betay = 120.03 #[m] original: 200m
  alphy = -7.46 #original: -7
  nemty = 0.12155 #[mm-mrad]
elif args.pencil: #for a reasonable pencil beam
  N=1e5 #number of particles
  betax = 1e-2 #[m]
  alphx = 0
  nemtx = 1e-4 #[mm-mrad]
  betay = 1e-2 #[m]
  alphy = 0
  nemty = 1e-4 #[mm-mrad]
  ifplot=False
elif args.twiss: #input variables are inputs
  N = float(sys.argv[1])
  betax = float(sys.argv[2]) #[m]
  alphx = float(sys.argv[3])
  nemtx = float(sys.argv[4]) #[um-mrad]
  betay = float(sys.argv[5]) #[m]
  alphy = float(sys.argv[6])
  nemty = float(sys.argv[7]) #[um-mrad]
  thick = float(sys.argv[8]) #[mm]

print("You've entered: {:f}mm thick".format(thick),materials,", {:.0e} protons, ".format(N),betax,alphx,nemtx,betay,alphy,nemty)

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
  #savename,xtarg,ytarg= simulation( N,material,    beam,thick,nemtx ,nemty ,alphx,alphy,betax,betay,energy,zoff,Engcut,engplot,loadParts,beamAngle,beamFile):
  savename,xtarg,ytarg = simulation( N,material,"proton",thick,nemtx ,nemty ,alphx,alphy,betax,betay,energy,zoff,engcut,engplot,loadParts,beamAngle,beamFile)
  #returns the savename and the x and y distributions of particle positions 
  #These returned arrays are from the MiniScatter detector, 5m after the PBW, ~where the ESS Target is.  
  
  #Now plot the distributions with various views depending on the material
  if material == "G4_Galactic" or material == "G4_AIR":
    if ifplot:
      print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xtarg),np.max(ytarg)))
      xmax = math.ceil(np.max(xtarg)/100)*100
      #For plotting x,y plots for core view and halo zoom-ins
      #plotFit(xs,    ys, savename,xlim,   ylim,material)
      #xlim<=10 => [-xlim*sigma,xlim*sigma]; xlim>10 => [-xlim,xlim]
      #ylim==0 => with plot; ylim!=0 => [0,ylim] (5 particles /(mag of N) is for halo)
      #plotFit(xtarg,ytarg,savename,  3,      0,material,thick) #3 sigma core
      plotFit(xtarg,ytarg,savename, xmax,5/(10**(mag+0)),material,thick) #full range halo, with dynamic halo zoom
  elif material == "G4_Al":
    if ifplot:
      print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xtarg),np.max(ytarg)))
      xmax = math.ceil(np.max(xtarg)/100)*100
      print(xmax)
      #plotFit(xs,    ys, savename,xlim,   ylim,material)
      #plotFit(xtarg,ytarg,savename,  3,      0,material,thick) #3 sigma core
      plotFit(xtarg,ytarg,savename, xmax,5/(10**(mag+0)),material,thick) #full range halo
      #plotFit(xtarg,ytarg,savename,  10,10/(10**(mag+0)),material,thick) #10 sigma range halo
  elif material == "G4_Au":
    if ifplot:
      print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xtarg),np.max(ytarg)))
      xmax = math.ceil(np.max(xtarg)/100)*100
      #plotFit(xs,    ys, savename,xlim,   ylim,material) 
      #plotFit(xtarg,ytarg,savename,  3,      0,material,thick) #3 sigma core
      plotFit(xtarg,ytarg,savename, xmax,15/(10**(mag+0)),material,thick) #full range halo
      #plotFit(xtarg,ytarg,savename,  10,15/(10**(mag+0)),material,thick) #10 sigma range halo
  
  #Multi-Material Plot section, continues the material loop to plot data.
  #For creating the PDF of particle distribution from findFit function optimized mu and sigma.
  from scipy.stats import norm 
  
  #Define plotting characteristics depending on material
  if material == "G4_Galactic":
    mat = "Vacuum"
    color = "magenta"
    dash = "m--"
    tsp = 0.2 #transparency
  elif material == "G4_Al":
    mat ="Al"
    color = "blue"
    dash = "b--"
    tsp = 0.5
  elif material == "G4_AIR":
    mat ="Air"
    color = "cyan"
    tsp = 0.2
    dash = "c--"
  elif material == "G4_Au":
    mat = "Au"
    color = "gold"
    tsp = 0.5
    dash = "r--"

  #Use Scipy.optimize.curve_fit in my findFit function to get mus and sigmas:
  mux, sigmax[material], xinterval = findFit(xtarg) #dynamically gets parameters AND histogram intervals!
  muy, sigmay[material], yinterval = findFit(ytarg)

  #Find range of particles that are within 3 sigma for each material
  sigx=np.abs(xtarg)>3*sigmax[material]# and np.abs(xtarg)<10*sigma)
  sigy=np.abs(ytarg)>3*sigmay[material]

  #Utilize dictionary to find % of particles within 3 sigma 
  pOut3sigx[material] = len(xtarg[sigx])/len(xtarg)*100 
  pOut3sigy[material] = len(ytarg[sigy])/len(ytarg)*100
  #print(material," gives ",pOut3sigx[material],"% outisde 3 sigma in X")
  #print(material," gives ",pOut3sigy[material],"% outisde 3 sigma in Y")

  #Update the texts to include this materials sigma and % outside 3 sigma
  percenttextx +="\n"+mat+" = "+"{:.2f}".format(pOut3sigx[material])+"%"
  percenttexty +="\n"+mat+" = "+"{:.2f}".format(pOut3sigy[material])+"%"
  sigmatextx +="\n"+mat+" = "+"{:.2f}".format(sigmax[material])+"mm"
  sigmatexty +="\n"+mat+" = "+"{:.2f}".format(sigmay[material])+"mm"

  #Make the histogram of the full energy distrubtion for X
  nx, binsx, patchesx = s1.hist(xtarg, xinterval, density=True, facecolor=color, alpha=tsp,label=mat)
  
  #Add the "best fit" line using the earlier norm.fit mus and sigmas for X
  y1 = norm.pdf(binsx, mux, sigmax[material])
  l1 = s1.plot(binsx, y1, dash, linewidth=1,label=mat) #important to keep it as l# or it won't work

  #Make the histogram of the full energy distrubtion for Y
  ny, binsy, patchesy = s2.hist(ytarg, yinterval, density=True, facecolor=color, alpha=tsp,label=mat)
  
  #Add the "best fit" line using the earlier norm.fit mus and sigmas for Y
  y2 = norm.pdf(binsy, muy, sigmay[material])
  l2 = s2.plot(binsy, y2, dash, linewidth=1,label=mat) #labeled by material short-name

#Set Plot characterists
s=10 # number of sigma width to plot
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
        rf"Initial beam of $\epsilon_x={{:.1e}}, \beta_x={{:.0e}}, \alpha_x={{:.1f}}$, ".format(nemtx,betax,alphx) +
    rf" $\epsilon_y={{:.1e}}, \beta_y={{:.0e}}, \alpha_y={{:.1f}}$ ".format(nemty,betay,alphy),fontsize=18,y=0.99) #fontweight="bold",
#suptitle can have values inside math if use 2 sets of {{}} - fix from "linuxtut.com" blog post

#Can date stamp the multi plot for easier tracking of changes, if necessary
from datetime import datetime
dt = datetime.now()
name = savename+"_"+dt.strftime("%H-%M-%S")+"_multi.png" #update the savename to not overwrite others
#print(name) #show so it is easy to find the file
#fig.savefig(name,bbox_inches="tight")
#plt.show()
plt.close() #be sure to close the plot
print(dt.strftime("%H-%M-%S"),"\n")
