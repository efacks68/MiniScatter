#simulation.py
#Eric Fackelman
#29 March - 5 April 2022

#This is to run the simulation of a beam through a PBW of input material 
# and output the position X and Y arrays at ESS Target location.
# Also can plot the Energy Distribution if requested.

def simulation(N,material,beam,thick,epsx,epsy,alphax,alphay,betax,betay,energy,zoff,Engcut,engplot):
  import numpy as np
  import ROOT
  import os
  import sys

  #Setup MiniScatter -- modify the path to where you built MiniScatter!
  MiniScatter_path="../../MiniScatter/build/."
  sys.path.append(MiniScatter_path) #uncomment this if this is your first time running this.
  #print(os.getcwd())
  if os.getcwd() != "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/":
    os.chdir('/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/')
    #print(os.getcwd())

  import miniScatterDriver
  #import miniScatterScanner
  #import miniScatterPlots

  ### Basic simulation parameters ###
  QUIET   = False #Reduced output, doesn't show events
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
  EPSX   = epsx #[um]
  BETAX  = betax #[m]
  ALPHAX = alphax #[um-mrad]
  EPSY   = epsy #[um]
  BETAY  = betay #[m]
  ALPHAY = alphay #[um-mrad]
  baseSimSetup["COVAR"] = (EPSX,BETAX,ALPHAX,EPSY,BETAY,ALPHAY) 

  #Use a flat distribution or cut the tails of the Gaussian?
  #baseSimSetup["BEAM_RCUT"] = 3.0

  #Where to start the beam [mm]
  #baseSimSetup["ZOFFSET_BACKTRACK"] = True
  baseSimSetup["ZOFFSET"]   = zoff#-10.0 #Auto = 0
  #anything behind Target is NEGATIVE!

  #Beam particle type
  baseSimSetup["BEAM"]    = beam

  baseSimSetup["WORLDSIZE"] = 1000.0 #Make the world wider for seeing where particles go

  #Target is 1 mm of aluminium
  baseSimSetup["THICK"] = thick

  baseSimSetup["MAT"] = material
  #Valid choices: G4_Al, G4_Au, G4_C, G4_Cu, G4_Pb, G4_Ti, G4_Si, G4_W, G4_U, G4_Fe, G4_MYLAR, G4_KAPTON,
  #G4_STAINLESS-STEEL, G4_WATER,G4_SODIUM_IODIDE, G4_Galactic, G4_AIR, Sapphire, ChromoxPure, ChromoxScreen

  #Detector distance from target center [mm] Default is 50mm behind Target
  #For multiple detector locations, make a list, e.g. [-5,5,5000]
  baseSimSetup["DIST"] = [5000] #only at ESS Target location 

  #Some output settings
  baseSimSetup["QUICKMODE"] = False #Include slow plots
  baseSimSetup["MINIROOT"]  = False #Skip TTRees in the .root files
  baseSimSetup["ANASCATTER"] = True #don't do Analytical Scatter Angle Test

  baseSimSetup["EDEP_DZ"]   = 1.0 #Z bin width for energy deposit histogram
  baseSimSetup["CUTOFF_RADIUS"] = 100.0 #Larger radial cutoff
  baseSimSetup["CUTOFF_ENERGYFRACTION"] = Engcut #Minimum percent of full Energy to use in cutoff calculations
  #0.95 is default, but let's be explicit

  #Store the .root files in a subfolder from where this script is running,
  # normally MiniScatter/examples, in order to keep things together
  baseSimSetup["OUTFOLDER"] = os.path.join("/scratch/ericdf/Scratch/PBWScatter/") 
  #put in Scratch of HepLab0# for faster processing, as per Kyrre

  baseSimSetup["N"]     = N #Just a few events here! Remember that thicker targets are slower

  baseSimSetup["POSLIM"] = 100 #XY histogram Position Limit for a few, check RootFileWriter.cc
  #print(baseSimSetup)

  #Run the simulation
  #copy so it is if running multiple scans in a Jupyter notebook
  simSetup_simple1 = baseSimSetup.copy()

  #Define material nickname
  #radiation lengths are from https://pdg.lbl.gov/2019/AtomicNuclearProperties/
  if baseSimSetup["MAT"] == "G4_Galactic":
    mat = "Vac"
    radLen = 1e9 #[mm] basically infinity
    z = 1 #not real
  elif baseSimSetup["MAT"] == "G4_Al":
    mat = "Al"
    radLen = 88.97 #[mm] #0.2401 #[g/mm^2] # #24.01 [g/cm^2]
  elif baseSimSetup["MAT"] == "G4_AIR":
    mat = "Air"
    radLen = 3.122e5 #[mm] -taken from 80% N2 gas(3.260e5mm) and 20% O2 gas (2.571e5mm)
  elif material == "G4_Au":
    mat = "Au"
    radLen = 3.344

  #Give the .root file a dynamic name
  outname = "simplePBW_"+str(baseSimSetup["THICK"])+"mm"+mat+"_N{:.0e}_b{:.0e},a{:.0f},e{:.0e}".format(baseSimSetup["N"],betax,alphax,epsx)
  print(outname)
  simSetup_simple1["OUTNAME"] = outname

  #Variables for automation
  savepath = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/" #Eric's files location
  savename=savepath+outname #base savename for plots downstream, brings directly to my directory
  savedfile=os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root"
  #print(savename,'\n',savedfile)
  #print("before run",os.getcwd())

  #Run simulation or load old simulation root file!
  #miniScatterDriver.runScatter(simSetup_simple1, quiet=QUIET) #this was Kyrre's, but it wasn't even trying to load old runs
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

  #Now get the "exit" distributions that you actually want to plot.
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
      #zexit[i] = myTree.z #not used, but available.
      #pzexit[i] = myTree.pz
      Eexit[i] = myTree.E
      #PDGexit[i] = myTree.PDG
      #break
  myFile.Close() 

  #For plotting Energy distribution of all species
  if engplot:
    import math
    mag=math.floor(math.log10(N)) #for dynamic title name
    engname=savename+"_EnergyPlot.png" #same for each plot
    titl = "Energy Distribution at ESS Target Through PBW of "+mat+"\n For All Particle Species"

    if baseSimSetup["MAT"] == "G4_Galactic":
      print("Vacuum Energy Plot not working now, plots empty histogram. *Shrugs*")
    else: #simple histogram plot
      import matplotlib.pyplot as plt
      plt.close()
      plt.hist(Eexit,100,log=True)
      #future addition - make other PDG entries be another color?
      plt.xlabel("Energy [MeV]")
      plt.ylabel("Counts")
      plt.xlim([0,2005.0])
      plt.title(rf"Energy Distribution at ESS Target of 10$^{{:d}}$ Protons".format(mag)+
        "\nThrough PBW of "+mat+" For All Particle Species")
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
  print("Full Energy distribution of {:d} particles with minimum Energy {:.3f}MeV through ".format(len(Eexit),np.min(Eexit)),
      mat," PBW")

  twiss, numPart,objects = miniScatterDriver.getData(savedfile)
  Pepsx = twiss['tracker_cutoffPDG2212']['x']['eps']
  Pbetax = twiss['tracker_cutoffPDG2212']['x']['beta']
  Palphx = twiss['tracker_cutoffPDG2212']['x']['alpha']
  Pepsy = twiss['tracker_cutoffPDG2212']['y']['eps']
  Pbetay = twiss['tracker_cutoffPDG2212']['y']['beta']
  Palphy = twiss['tracker_cutoffPDG2212']['y']['alpha']


  #print("And use Energy cutoff!!! as per Kyrre morning of 7 April") #already using that by default, but passed it explicitly too
  print("The initial Twiss are: \nX: {:.3f}um, {:.3f}m, {:.3f}um-mrad\nY: {:.3f}um, {:.3f}m, {:.3f}um-mrad".format(epsx,betax,alphax,epsy,betay,alphay))
  print("The Twiss at Target are: \nX: {:.3f}um, {:.3f}m, {:.3f}um-mrad\nY: {:.3f}um, {:.3f}m, {:.3f}um-mrad".format(Pepsx,Pbetax,Palphx,Pepsy,Pbetay,Palphy))

  depsx = Pepsx - epsx
  depsy = Pepsy - epsy
  if beam == "proton":
    partmass = 938.27209 #[MeV]
    z=1
  elif beam == "electron":
    partmass = 0.511 #[MeV]
    z=1
  gamma = (energy - partmass)/partmass 
  beta = np.sqrt(1 - 1/(gamma*gamma))
  p=partmass*partmass/energy
  thetasq = 13.6 * z / (energy-p) * np.sqrt(thick/radLen) * (1 + 0.038 * np.log(thick*z*z/(radLen*(1-p/energy)))) #from MiniScatter
  #thetasq = thetasq / 4
  
  print("gamma: {:.3f}, beta: {:.3f}, theta^2: {:.3e} radians".format(gamma,beta,thetasq))

  delepsx = 0.5 * betax * thetasq * thetasq #[m*rad^2]
  delepsy = 0.5 * betay * thetasq * thetasq
  delalpx = epsx * alphax / (epsx + delepsx)
  delalpy = epsy * alphay / (epsy + delepsy)
  delbetx = epsx * betax / (epsx + delepsx)
  delbety = epsy * betay / (epsy + delepsy)
  print("The change in emittances calculated by GEANT4 are: {:.3f}, {:.3f}um.".format(depsx,depsy))
  print("The expected changes in emittances from CERN paper Eq 7: {:.3f}, {:.3f}um".format(delepsx*1e6,delepsy*1e6))
  print("The expected changes in beta from CERN paper Eq 8: {:.3f}, {:.3f}m".format(delbetx,delbety))
  print("The expected changes in alpha from CERN paper Eq 8: {:.3f}, {:.3f}um-mrad".format(delalpx,delalpy))

  print("\n")
  gammax = 1 + alphax * alphax / betax
  gammay = 1 + alphay * alphay / betay
  delepsx = 0.5 * thetasq * thetasq * (betax + thick * alphax + thick*thick/3 * gammax) #[m*rad^2]
  delepsy = 0.5 * thetasq * thetasq * (betay + thick * alphay + thick*thick/3 * gammay) #[m*rad^2]
  delalpx = (epsx * alphax - thick * 0.5 * thetasq) / (epsx + delepsx)
  delalpy = (epsy * alphay - thick * 0.5 * thetasq) / (epsy + delepsy)
  delbetx = (epsx * betax + thick * thick / 3 * thetasq) / (epsx + delepsx)
  delbety = (epsy * betay + thick * thick / 3 * thetasq) / (epsy + delepsy)
  print("The change in emittances calculated by GEANT4 are: {:.3f}, {:.3f}um.".format(depsx,depsy))
  print("The expected changes in emittances from CERN paper Eq 15: {:.3f}, {:.3f}um".format(delepsx*1e6,delepsy*1e6))
  print("The expected changes in beta from CERN paper Eq 16: {:.3f}, {:.3f}m".format(delbetx,delbety))
  print("The expected changes in alpha from CERN paper Eq 16: {:.3f}, {:.3f}um-mrad".format(delalpx,delalpy))


  #want to add Analytical Scatter Angle Options to MiniScatterDriver. Making New branch 11.4.22
  return savename, xexit, yexit