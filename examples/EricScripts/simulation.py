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

  #Run the simulation
  #copy so it is if running multiple scans in a Jupyter notebook
  simSetup_simple1 = baseSimSetup.copy()

  #Give the .root file a dynamic name
  outname = "simplePBW_"+str(baseSimSetup["THICK"])+"mm"+mat+"_N{:.0e}_b{:.0e},a{:.0f},e{:.0e}_{:.2f}Engcut".format(baseSimSetup["N"],betax,alphax,epsx,Engcut)
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

  #If one wants to use the initial spread
  myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
  myTree= myFile.Get('InitParts')
  #print(myTree)
  xinit = np.zeros(myTree.GetEntries())
  pxinit = np.zeros(myTree.GetEntries())
  yinit = np.zeros(myTree.GetEntries())
  pyinit = np.zeros(myTree.GetEntries())
  #Einit = np.zeros(myTree.GetEntries())
  #print(len(xinit))
  for i in range(myTree.GetEntries()):
      myTree.GetEntry(i)
      #print(myTree.x,myTree.y,myTree.px,myTree.py,myTree.E,myTree.PDG,myTree.charge,myTree.eventID)
      xinit[i] = myTree.x
      pxinit[i] = myTree.px
      yinit[i] = myTree.y
      pyinit[i] = myTree.py
      #Einit[i] = myTree.E
  myFile.Close()

  #27.4 just added target-exit pull in and changed previous xexit arrays to target arrays!
  #Now get the "target-exit" distributions for plotting with the Formalism distribution below. 
  #These are not at the ESS Target location, but are at the far side of the PBW
  myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
  myTree= myFile.Get('TargetExit') #TrackerHits has all trackers, be sure to only have 1!
  #print(myTree)

  xexit = np.zeros(myTree.GetEntries()) #dynamic length arrays
  pxexit = np.zeros(myTree.GetEntries())
  yexit = np.zeros(myTree.GetEntries())
  pyexit = np.zeros(myTree.GetEntries())
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
      Eexit[i] = myTree.E
      #PDGexit[i] = myTree.PDG
      #break
  myFile.Close() 

  #Below are the distributions taken from the ESS Target location (the detector located 5m from PBW)
  myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
  myTree= myFile.Get('TrackerHits') #TrackerHits has all trackers, be sure to only have 1!

  #Now get the "targ" distributions
  xtarg = np.zeros(myTree.GetEntries()) #dynamic length arrays
  pxtarg = np.zeros(myTree.GetEntries())
  ytarg = np.zeros(myTree.GetEntries())
  pytarg = np.zeros(myTree.GetEntries())
  Etarg = np.zeros(myTree.GetEntries())
  #PDGtarg = np.zeros(myTree.GetEntries())
  #print("The length of the arrays are ",len(xtarg))

  for i in range(myTree.GetEntries()): #put data in arrays
      myTree.GetEntry(i)
      xtarg[i] = myTree.x
      pxtarg[i] = myTree.px
      ytarg[i] = myTree.y
      pytarg[i] = myTree.py
      Etarg[i] = myTree.E
      #PDGtarg[i] = myTree.PDG
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
      plt.hist(Etarg,100,log=True)
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
  #Be sure this Energy cutting is for the ESS TARGet arrays 
  Ecut = energy * Engcut #[MeV]
  mask = np.abs(Etarg) > Ecut #mask to exclude low energy particles
  Etarg = Etarg[mask]
  xtarg = xtarg[mask]
  pxtarg = pxtarg[mask]
  ytarg = ytarg[mask]
  pytarg = pytarg[mask]

  #Display Full Energy distribution results
  print("Full Energy distribution of {:d} particles with minimum Energy {:.3f}MeV through ".format(len(Etarg),np.min(Etarg)),
      mat," PBW")

  #Twiss MCS Formalism Calculations from https://cds.cern.ch/record/499590
  #constants for below use
  c = 2.99792e8 #[m/s]
  MeV = 1e6*1.602e-19 
  um = 1e-6 #[m] #need to convert to real units as the equations use real units.
  m = 1 #[m]
  mm = 1e-3 #[m]
  ummrad = um*1e-3

  #get Twiss parameters from tracker at ESS Target location:
  twiss, numPart,objects = miniScatterDriver.getData(savedfile)
  Tepsx = twiss['tracker_cutoffPDG2212']['x']['eps'] *um #[um]
  Tbetax = twiss['tracker_cutoffPDG2212']['x']['beta'] *m #[m]
  Talphx = twiss['tracker_cutoffPDG2212']['x']['alpha'] #[um-mrad]
  Tepsy = twiss['tracker_cutoffPDG2212']['y']['eps'] *um #[um]
  Tbetay = twiss['tracker_cutoffPDG2212']['y']['beta'] *m
  Talphy = twiss['tracker_cutoffPDG2212']['y']['alpha'] #[um-mrad]

  # Get Twiss for Target Exit
  Eepsx = twiss['target_exit_cutoff']['x']['eps'] *um #[um]
  Ebetax = twiss['target_exit_cutoff']['x']['beta'] *m #[m]
  Ealphx = twiss['target_exit_cutoff']['x']['alpha'] #[um-mrad]
  Eepsy = twiss['target_exit_cutoff']['y']['eps'] *um #[um]
  Ebetay = twiss['target_exit_cutoff']['y']['beta'] *m
  Ealphy = twiss['target_exit_cutoff']['y']['alpha'] #[um-mrad]
  depsx = Eepsx - epsx*um #[um] #passed #s aren't in units system
  depsy = Eepsy - epsy*um #[um]

  #particle characteristic values
  if beam == "proton":
    partmass = 938.27209 #[MeV]
    z=1
  elif beam == "electron":
    partmass = 0.511 #[MeV]
    z=1

  #beam calculations
  gamma = 1 + energy/partmass #from PrintTwissParameters
  beta = np.sqrt(gamma*gamma -1 )/gamma
  Gepsx = epsx*um/(beta*gamma)
  Gepsy = epsy*um/(beta*gamma)

  q=partmass*partmass/energy
  p = np.sqrt((energy*MeV*energy*MeV-partmass*MeV*partmass*MeV)/ (c*c))/(MeV/c)
  betap = beta*p
  #betap = energy-q
  #thetasq = 13.6 * z / betap * np.sqrt(thick/radLen) * (1 + 0.038 * np.log(thick/radLen)) #from Eq 5
  thetasq = 13.6 * z / betap * np.sqrt(thick/radLen) * (1 + 0.038 * np.log(thick*z*z/(radLen*(1-q/energy)))) #from MiniScatter
  print("\nbetap: {:.3e}, gamma: {:.3f}, beta: {:.3f}, radLen: {:.3f}mm, theta^2: {:.3e} radians".format(betap,gamma,beta,radLen,thetasq))
  #using beta*q or energy-q produces similar #s (within 7%)

  #Calculations from Eq 7 and 8
  delepsx = 0.5 * betax*m * thetasq * thetasq #[m*rad^2]
  delepsy = 0.5 * betay*m * thetasq * thetasq
  delalpx = epsx*um * alphax / (epsx*um + delepsx)
  delalpy = epsy*um * alphay / (epsy*um + delepsy)
  delbetx = epsx*um * betax*m / (epsx*um + delepsx)
  delbety = epsy*um * betay*m / (epsy*um + delepsy)

  #print('Norm. Emittance, Beta, Alpha')
  #print("The initial Twiss are: \nX: {:.3f}um, {:.3f}m, {:.3f}\nY: {:.3f}um, {:.3f}m, {:.3f}".format(epsx,betax,alphax,epsy,betay,alphay)) #passed #s aren't in system
  #print("The Twiss at PBW Exit are: \nX: {:.3f}um, {:.3f}m, {:.3f}\nY: {:.3f}um, {:.3f}m, {:.3f}".format(Eepsx/um,Ebetax,Ealphx,Eepsy/um,Ebetay,Ealphy))
  #print("The calculated Twiss at PBW Exit are: \nX: {:.3f}um, {:.3f}m, {:.3f}\nY: {:.3f}um, {:.3f}m, {:.3f}".format((epsx*um+delepsx)/um,delbetx,delalpx,(epsy*um+delepsy)/um,delbety,delalpy))
  #print("The Twiss at Target are: \nX: {:.3f}um, {:.3f}m, {:.3f}\nY: {:.3f}um, {:.3f}m, {:.3f}".format(Tepsx/um,Tbetax,Talphx,Tepsy/um,Tbetay,Talphy))

  #Calculations from Eq 15 and 16
  gammax = (1 + alphax * alphax ) / betax*m #[pm-mrad-mrad]??
  gammay = (1 + alphay * alphay ) / betay*m
  deleps2x = 0.5 * thetasq * thetasq * (betax*m + thick*mm * alphax + thick*mm*thick*mm/3 * gammax) #[m*rad^2]
  deleps2y = 0.5 * thetasq * thetasq * (betay*m + thick*mm * alphay + thick*mm*thick*mm/3 * gammay) #[m*rad^2]
  delalp2x = (epsx*um * alphax - thick*mm * 0.5 * thetasq) / (epsx*um + delepsx)
  delalp2y = (epsy*um * alphay - thick*mm * 0.5 * thetasq) / (epsy*um + delepsy)
  delbet2x = (epsx*um * betax*m + thick*mm * thick*mm / 3 * thetasq) / (epsx*um + delepsx)
  delbet2y = (epsy*um * betay*m + thick*mm * thick*mm / 3 * thetasq) / (epsy*um + delepsy)
  print("\nThe change in emittances calculated by GEANT4 are: {:.3f}, {:.3f}um.".format(depsx/um,depsy/um))
  print("The expected changes in emittances from CERN paper Eq 7: {:.3f}, {:.3f}um".format(delepsx/um,delepsy/um))
  print("The expected changes in emittances from CERN paper Eq 15: {:.3f}, {:.3f}um".format(deleps2x/um,deleps2y/um))

  #Plotting PBW Exit distribution vs the PDF produced from the Formalism Equations
  from plotFit import plotTwissFit
  from plotFit import plotFit 
  #plotFit(xs,    ys, savename,xlim,   ylim,material)
  #plotFit(xinit,yinit,savename+'init',  3,      0,material)
  #plotFit(xexit,yexit,savename+'texit',  3,      0,material)

  #Use list of Twiss values for simple passing of data 
  Twissix = [betax,alphax,epsx*um,Gepsx] #initial Twiss
  Twissiy = [betay,alphay,epsy*um,Gepsy]
  Twisstex = [Ebetax,Ealphx,Eepsx,Eepsx/(beta*gamma)]#target exit Twiss
  Twisstey = [Ebetay,Ealphy,Eepsy,Eepsy/(beta*gamma)]
  Twisscx = [delbetx,delalpx,epsx*um+delepsx,(epsx*um+delepsx)/(beta*gamma)] #calculated target exit Twiss
  Twisscy = [delbety,delalpy,epsy*um+delepsy,(epsy*um+delepsy)/(beta*gamma)]
  Twissc2x = [delbet2x,delalp2x,epsx*um+deleps2x,(epsx*um+deleps2x)/(beta*gamma)] #calculated 2 target exit Twiss
  Twissc2y = [delbet2y,delalp2y,epsy*um+deleps2y,(epsy*um+deleps2y)/(beta*gamma)]
  
  #Displays Twiss values and calculated mu and sigma
  print('\nInitialTwiss')
  plotTwissFit(Twissix,Twissiy,xinit,yinit,savename+'init',material,'MiniScatter Initial',Twisstex,Twisstey,'')
  #print('\nPBW Exit Twiss')
  #plotTwissFit(Twisstex,Twisstey,xexit,yexit,savename+'texit',material,'MiniScatter',Twisstex,Twisstey,'PBW Exit')
  print('\nEq8 PBW Exit Twiss')
  plotTwissFit(Twisscx,Twisscy,xexit,yexit,savename+'texitEq8',material,'Formalism Eq 8',Twisstex,Twisstey,'PBW Exit')
  print('\nEq16 PBW Exit Twiss')
  plotTwissFit(Twissc2x,Twissc2y,xexit,yexit,savename+'texitEq16',material,'Formalism Eq 16',Twisstex,Twisstey,'PBW Exit')

  return savename, xtarg, ytarg