#simulation.py
#Eric Fackelman
#29 March 2022

#This is to produce the simulation for beam through PBW and output the arrays at Target location

def simulation(N,material,epsx,epsy,alphax,alphay,betax,betay,energy,zoff,Engcut,engplot):
  import numpy as np
  import ROOT
  import os
  from datetime import datetime
  import sys

  #Setup MiniScatter -- modify the path to where you built MiniScatter!
  MiniScatter_path="../MiniScatter/build/."
  #sys.path.append(MiniScatter_path) #uncomment this if this is your first time running this.
  #print(os.getcwd())
  if os.getcwd() != "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/":
    os.chdir('/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/')
    #print(os.getcwd())

  import miniScatterDriver
  #import miniScatterScanner
  #import miniScatterPlots

  ### Basic simulation parameters ###

  QUIET   = True #Reduced output, doesn't show events
  TRYLOAD = True  #Try to load already existing data instead of recomputing, only if using getData_TryLoad function.
  NUM_THREADS = 8 #Number of parallel threads to use for scans
  #Where to store temporary data for scans (a fast file system, NOT EOS/AFS)
  TMPFOLDER = "/tmp/miniScatter/SimpleDemo_thicknessScan"

  #When making multiple scans, it's nice to first create a `baseSimSetup` and then modify it for each simulation
  # Note: each argument here corresponds roughly to a command line argument.
  # Look inside miniScatterDriver.runScatter() to see how.

  baseSimSetup = {}
  #baseSimSetup["PHYS"] = "QGSP_BERT__SS" #Use the __SS physics lists for thin foils due to checking each atom cross section
  baseSimSetup["PHYS"]  = "QGSP_BERT_EMZ" #better for scattering through 1mm sheets

  baseSimSetup["ENERGY"] = energy #2000.0 #[MeV] #ESS beam energy

  #Use a distribution defined by Twiss parameters for ESS beam ~where PBW is
  # 3 variables = symmetric, 6 variables = asymetric
  EPSX   = epsx#0.113#e3 #[um]
  BETAX  = betax#1000 #[m]
  ALPHAX = alphax#-50
  EPSY   = epsy#0.121#e3 #[um]
  BETAY  = betay#200 #[m]
  ALPHAY = alphay#-7
  baseSimSetup["COVAR"] = (EPSX,BETAX,ALPHAX,EPSY,BETAY,ALPHAY) 

  #Use a flat distribution or cut the tails of the Gaussian?
  #baseSimSetup["BEAM_RCUT"] = 3.0

  #Where to start the beam [mm]
  #baseSimSetup["ZOFFSET_BACKTRACK"] = True
  baseSimSetup["ZOFFSET"]   = zoff#-10.0 #Auto = 0
  #anything behind Target is NEGATIVE!

  #Beam particle type
  baseSimSetup["BEAM"]    = "proton"

  baseSimSetup["WORLDSIZE"] = 1000.0 #Make the world wider

  #Target is 1 mm of aluminium
  baseSimSetup["THICK"] = 1
  baseSimSetup["MAT"] = material
  #Valid choices: G4_Al, G4_Au, G4_C, G4_Cu, G4_Pb, G4_Ti, G4_Si, G4_W, G4_U, G4_Fe, G4_MYLAR, G4_KAPTON,
  #G4_STAINLESS-STEEL, G4_WATER,G4_SODIUM_IODIDE, G4_Galactic, G4_AIR, Sapphire, ChromoxPure, ChromoxScreen

  #Detector distance from target center [mm] Default is 50mm behind Target
  #baseSimSetup["DIST"] = 500.0 
  #For multiple detector locations, make a list
  baseSimSetup["DIST"] = [5000] #only at ESS Target location #[-5,5,5000]

  #Some output settings
  baseSimSetup["QUICKMODE"] = False #Include slow plots
  baseSimSetup["MINIROOT"]  = False #Skip TTRees in the .root files

  baseSimSetup["EDEP_DZ"]   = 1.0 
  baseSimSetup["CUTOFF_RADIUS"] = 100.0 #Larger radial cutoff

  #Store the .root files in a subfolder from where this script is running,
  # normally MiniScatter/examples, in order to keep things together
  baseSimSetup["OUTFOLDER"] = os.path.join("/scratch/ericdf/Scratch/PBWScatter/") 
  #put in Scratch of HepLab0# for faster processing, as per Kyrre

  baseSimSetup["N"]     = N #Just a few events here! Remember that thicker targets are slower

  baseSimSetup["POSLIM"] = 100 #XY histogram Position Limit for a few, check RootFileWriter.cc
  #print(baseSimSetup)

  #Run the simulation
  #copy so it is not messed up(?)
  simSetup_simple1 = baseSimSetup.copy()

  #Define material nickname
  if baseSimSetup["MAT"] == "G4_Galactic":
    mat = "Vac"
  elif baseSimSetup["MAT"] == "G4_Al":
    mat = "Al"
  elif baseSimSetup["MAT"] == "G4_AIR":
    mat = "Air"
  elif material == "G4_Au":
    mat = "Au"

  #Give the .root file a name
  outname = "simplePBW_"+str(baseSimSetup["THICK"])+"mm"+mat+"_{:.0e}_DistribTest".format(baseSimSetup["N"])
  print(outname)
  simSetup_simple1["OUTNAME"] = outname #set the disctionary value as the outname.
    #Variables for automation
  savepath = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/" #Eric's files location
  savename=savepath+outname #base savename for plots downstream, brings directly to my directory
  #print(savename)
  print("before run",os.getcwd())

  #29.3 - not doing runs, just opening files for now
  #miniScatterDriver.runScatter(simSetup_simple1, quiet=QUIET) #this was Kyrre's, but it wasn't trying to load old runs...
  miniScatterDriver.getData_tryLoad(simSetup_simple1,quiet=QUIET)
  
  #If one wants to use the initial spread for some reason, as I did initially xD
  #myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
  #myTree= myFile.Get('InitParts') #because the name is capitalized!
  #print(myTree)
  #print(myTree.())
  #xinit = np.zeros(myTree.GetEntries())
  #pxinit = np.zeros(myTree.GetEntries())
  #yinit = np.zeros(myTree.GetEntries())
  #pyinit = np.zeros(myTree.GetEntries())
  #Einit = np.zeros(myTree.GetEntries())
  #print(len(xinit))
  #for i in range(myTree.GetEntries()):
  #    myTree.GetEntry(i)
  #    #print(myTree.x,myTree.y,myTree.px,myTree.py,myTree.E,myTree.PDG,myTree.charge,myTree.eventID)
  #    xinit[i] = myTree.x
  #    pxinit[i] = myTree.px
  #    yinit[i] = myTree.y
  #    pyinit[i] = myTree.py
  #    Einit[i] = myTree.E
  #myFile.Close()

  #For creating arrays with particle data
  myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
  myTree= myFile.Get('TrackerHits') #TrackerHits has all trackers, be sure to only have 1!
  #print(myTree)

  #print(myTree.GetEntries())
  xexit = np.zeros(myTree.GetEntries()) #dynamic length arrays
  pxexit = np.zeros(myTree.GetEntries())
  yexit = np.zeros(myTree.GetEntries())
  pyexit = np.zeros(myTree.GetEntries())
  #zexit = np.zeros(myTree.GetEntries())
  #pzexit = np.zeros(myTree.GetEntries())
  Eexit = np.zeros(myTree.GetEntries())
  #PDGexit = np.zeros(myTree.GetEntries())
  #print("The length of the arrays are ",len(xexit))

  for i in range(myTree.GetEntries()): #put data in arrays
      #print(myTree.E)
      myTree.GetEntry(i)
      xexit[i] = myTree.x
      pxexit[i] = myTree.px
      yexit[i] = myTree.y
      pyexit[i] = myTree.py
      #zexit[i] = myTree.z
      #pzexit[i] = myTree.pz
      Eexit[i] = myTree.E
      #PDGexit[i] = myTree.PDG
      #break
  myFile.Close() 

  #For plotting Energy distriution
  if engplot:
    engname=savename+"_EnergyPlot" #same for each plot
    titl = "Energy Distribution after "+mat+" PBW"

    if baseSimSetup["MAT"] == "G4_Galactic":
      print("Vacuum Energy Plot not working now, plots empty histogram. *Shrugs*")
    else: #simple histogram plot
      import matplotlib.pyplot as plt
      plt.close()
      plt.hist(Eexit,100,log=True)
      plt.xlabel("Energy [MeV]")
      plt.ylabel("Counts")
      plt.xlim([0,2005.0])
      plt.title(titl)
      print(titl,engname)
      plt.savefig(engname)
      plt.close()

  #Make Energy cut after plotting Energy distribution
  Ecut = 2000 * Engcut #[MeV]
  #print(np.min(Eexit))
  mask = np.abs(Eexit) > Ecut #mask to exclude low energy particles
  Eexit = Eexit[mask]
  xexit = xexit[mask]
  pxexit = pxexit[mask]
  yexit = yexit[mask]
  pyexit = pyexit[mask]
  #print('g',len(Eexit))

  #Display Full Energy distribution results
  print("Full Energy distribution of {:d} particles with minimum Energy {:.3f}MeV".format(len(Eexit),np.min(Eexit)))
 
  return savename, xexit, yexit