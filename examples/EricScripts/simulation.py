#simulation.py
#Eric Fackelman
#29 March - 1 June 2022

#This is to run the simulation of a beam through a PBW of input material 
# and output the position X and Y arrays at ESS Target location.
# Also can plot the Energy Distribution if requested.

def simulation(N,material,beam,thick,Inemx,Inemy,Ialphx,Ialphy,Ibetax,Ibetay,energy,zoff,Engcut,engplot,obj):
  import numpy as np
  import ROOT
  import os
  import sys
  from plotFit import calcTwiss

  #constants for below use
  #particle characteristic values
  if beam == "proton":
    partmass = 938.27209 #[MeV]
    z=1
  elif beam == "electron":
    partmass = 0.511 #[MeV]
    z=1
  c = 2.99792e8 #[m/s]
  MeV = 1e6*1.602e-19 
  um = 1e-6 #[m] #need to convert to real units as the equations use real units.
  m = 1 #[m]
  mm = 1e-3 #[m]
  ummrad = um*1e-3
  gamma_rel = 1 + energy/partmass #from PrintTwissParameters
  beta_rel = np.sqrt(gamma_rel*gamma_rel -1 )/gamma_rel

  #Get rid of Normalized Emittance!
  Igemx = Inemx/(beta_rel*gamma_rel)
  Igemy = Inemy/(beta_rel*gamma_rel)

  #Setup MiniScatter -- modify the path to where you built MiniScatter!
  MiniScatter_path="../../MiniScatter/build/."
  sys.path.append(MiniScatter_path) #uncomment this if this is your first time running this.
  #print(os.getcwd())
  if os.getcwd() != "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/":
    os.chdir("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/")
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
  EPSX   = Inemx #[um]
  BETAX  = Ibetax #[m]
  ALPHAX = Ialphx #[um-mrad]
  EPSY   = Inemy #[um]
  BETAY  = Ibetay #[m]
  ALPHAY = Ialphy #[um-mrad]
  baseSimSetup["COVAR"] = (EPSX,BETAX,ALPHAX,EPSY,BETAY,ALPHAY) 

  #Use a flat distribution or cut the tails of the Gaussian?
  #baseSimSetup["BEAM_RCUT"] = 3.0

  #Where to start the beam [mm]
  #baseSimSetup["ZOFFSET_BACKTRACK"] = True
  baseSimSetup["ZOFFSET"]   = zoff#-10.0 #Auto = 0
  #anything behind Target is NEGATIVE!

  #Beam particle type
  baseSimSetup["BEAM"]    = beam
  baseSimSetup["WORLDSIZE"] = 500.0 #Make the world wider for seeing where particles go

  if obj==0:
    #Target is 1 mm of aluminium
    baseSimSetup["THICK"] = thick
    baseSimSetup["MAT"] = material
    #Valid choices: G4_Al, G4_Au, G4_C, G4_Cu, G4_Pb, G4_Ti, G4_Si, G4_W, G4_U, G4_Fe, G4_MYLAR, G4_KAPTON,
    #G4_STAINLESS-STEEL, G4_WATER,G4_SODIUM_IODIDE, G4_Galactic, G4_AIR, Sapphire, ChromoxPure, ChromoxScreen

    #Detector distance from target center [mm] Default is 50mm behind Target
    #For multiple detector locations, make a list, e.g. [-5,5,5000]
    baseSimSetup["DIST"] = [5000] #Detector location. only at ESS Target location 
  else:
    baseSimSetup["THICK"] = 0.0
    baseSimSetup["DIST"] = [4.25] #Detector locations. At PBW Exit and ESS Target location 
    baseSimSetup["MAGNET"] = []
    #How to construct a magnet for miniScatterDriver, as per kyrsjo/MiniScatter/blob/master/examples/SiRi DeltaE-E detector.ipynb
    #Al first piece
    m1 = {}
    m1["pos"]      = 0.6250
    m1["type"]     = "TARGET"
    m1["length"]   = 1.2500 #[mm]
    m1["gradient"] = 0.0
    m1["keyval"] = {}
    m1["keyval"]["material"] = "G4_Al"
    m1["keyval"]["width"]    = 200.0 #[mm]
    m1["keyval"]["height"]   = 100.0 #[mm]
    baseSimSetup["MAGNET"].append(m1)
    #Water
    m2 = {}
    m2["pos"]      = 2.2501
    m2["type"]     = "TARGET"
    m2["length"]   = 2.0000 #[mm]
    m2["gradient"] = 0.0
    m2["keyval"] = {}
    m2["keyval"]["material"] = "G4_WATER"
    m2["keyval"]["width"]    = 200.0 #[mm]
    m2["keyval"]["height"]   = 100.0 #[mm]
    baseSimSetup["MAGNET"].append(m2)
    #Al Back
    m3 = {}
    m3["pos"]      = 3.7502
    m3["type"]     = "TARGET"
    m3["length"]   = 1.000 #[mm]
    m3["gradient"] = 0.0
    m3["keyval"] = {}
    m3["keyval"]["material"] = "G4_Al"
    m3["keyval"]["width"]    = 200.0 #[mm]
    m3["keyval"]["height"]   = 100.0 #[mm]
    baseSimSetup["MAGNET"].append(m3)

    mat = "Real"
    radLen = 64.08 #mm average of 2.25mm Al and 2.0mm Light Water
  #Some output settings
  baseSimSetup["QUICKMODE"] = False #Include slow plots
  baseSimSetup["MINIROOT"]  = False #Skip TTRees in the .root files
  baseSimSetup["ANASCATTER"] = True #don't do Analytical Scatter Angle Test

  Rcut = 50.0
  baseSimSetup["EDEP_DZ"]   = 1.0 #Z bin width for energy deposit histogram
  baseSimSetup["CUTOFF_RADIUS"] = Rcut #Larger radial cutoff #Decreased 10 May
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
  if material == "G4_Galactic":
    mat = "Vac"
    radLen = 1e9 #[mm] basically infinity
  elif material == "G4_Al":
    mat = "Al"
    radLen = 88.97 #[mm] #0.2401 #[g/mm^2] # #24.01 [g/cm^2]
  elif material == "G4_AIR":
    mat = "Air"
    radLen = 3.122e5 #[mm] -taken from 80% N2 gas(3.260e5mm) and 20% O2 gas (2.571e5mm)
  elif material == "G4_Au":
    mat = "Au"
    radLen = 3.344

  #Run the simulation
  #copy so it is if running multiple scans in a Jupyter notebook
  simSetup_simple1 = baseSimSetup.copy()

  #Give the .root file a dynamic name
  #outname = "simplePBW_"+str(baseSimSetup["THICK"])+"mm"+mat+"_N{:.0e}_b{:.0e},a{:.0f},e{:.0e}_{:.2f}Rcut".format(baseSimSetup["N"],Ibetax,Ialphx,Inemx,Rcut)
  outname = "simplePBW_"+str(baseSimSetup["THICK"])+"mm"+mat+"_N{:.0e}_ESSbeam_real".format(baseSimSetup["N"])
  #outname = "simplePBW_"+str(baseSimSetup["THICK"])+"mm"+mat+"_N{:.0e}_Pencilbeam".format(baseSimSetup["N"])
  print(outname)
  simSetup_simple1["OUTNAME"] = outname

  #Variables for automation
  savepath = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/" #Eric's files location
  savename=savepath+outname #base savename for plots downstream, brings directly to my directory
  savedfile=os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root"
  #print(savename,"\n",savedfile)
  #print("before run",os.getcwd())

  #Run simulation or load old simulation root file!
  #miniScatterDriver.runScatter(simSetup_simple1, quiet=QUIET) #this was Kyrre's, but it wasn't even trying to load old runs
  miniScatterDriver.getData_tryLoad(simSetup_simple1,quiet=QUIET)

  #If one wants to use the initial spread
  myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
  myTree= myFile.Get("InitParts")
  #print(myTree)
  xinit = np.zeros(myTree.GetEntries())
  pxinit = np.zeros(myTree.GetEntries())
  yinit = np.zeros(myTree.GetEntries())
  pyinit = np.zeros(myTree.GetEntries())
  Einit = np.zeros(myTree.GetEntries())
  #print(len(xinit))
  for i in range(myTree.GetEntries()):
      myTree.GetEntry(i)
      #print(myTree.x,myTree.y,myTree.px,myTree.py,myTree.E,myTree.PDG,myTree.charge,myTree.eventID)
      xinit[i] = myTree.x *mm
      pxinit[i] = myTree.px
      yinit[i] = myTree.y *mm
      pyinit[i] = myTree.py
      Einit[i] = myTree.E
  myFile.Close()

  initxTw = calcTwiss("Initial X","Initial X'",xinit,pxinit)
  inityTw = calcTwiss("Initial Y","Initial Y'",yinit,pyinit)

  #Filter with Energy for now as PDG isn't working yet
  Einit_filter = np.greater(Einit,energy*Engcut)
  xinit_filtered = xinit[Einit_filter]
  pxinit_filtered = pxinit[Einit_filter]
  yinit_filtered = yinit[Einit_filter]
  pyinit_filtered = pyinit[Einit_filter]

  initxTwf = calcTwiss("Initial X Filtered","Initial X' Filtered",xinit_filtered,pxinit_filtered)
  inityTwf = calcTwiss("Initial Y Filtered","Initial Y' Filtered",yinit_filtered,pyinit_filtered)

  if obj == 0:
    #27.4 just added target-exit pull in and changed previous xexit arrays to target arrays!
    #Now get the "target-exit" distributions for plotting with the Formalism distribution below. 
    #These are not at the ESS Target location, but are at the far side of the PBW
    myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
    myTree= myFile.Get("TargetExit") #TrackerHits has all trackers, be sure to only have 1!
    #print(myTree)
  else:
    #For now, if using several magnets, just get the TrackerHits 7.6.22
    myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
    myTree= myFile.Get("TrackerHits")

  xexit = np.zeros(myTree.GetEntries()) #dynamic length arrays
  pxexit = np.zeros(myTree.GetEntries())
  yexit = np.zeros(myTree.GetEntries())
  pyexit = np.zeros(myTree.GetEntries())
  pzexit = np.zeros(myTree.GetEntries())
  Eexit = np.zeros(myTree.GetEntries())
  PDGexit = np.zeros(myTree.GetEntries())
  #print("The length of the arrays are ",len(xexit))

  #import warnings
  #warnings.filterwarnings("error")

  for i in range(myTree.GetEntries()): #put data in arrays
      #print(myTree.E)
      myTree.GetEntry(i)
      #print(i, myTree.pz, myTree.px,myTree.PDG)
      pzexit[i] = myTree.pz
      xexit[i] = myTree.x *mm #m
      if pzexit[i] == 0.0: continue #11.5.22 recommended by Kyrre
      pxexit[i] = myTree.px/pzexit[i] #from Kyrre 5.5.22 to make it true angle!
      yexit[i] = myTree.y *mm #m
      pyexit[i] = myTree.py/pzexit[i]
      Eexit[i] = myTree.E
      PDGexit[i] = myTree.PDG
  myFile.Close() 

  exitxTw = calcTwiss("Exit X","Exit X'",xexit,pxexit)
  exityTw = calcTwiss("Exit Y","Exit Y'",yexit,pyexit)

  #Filter the relevant distributions
  posmax=5e-3 #mm (?)
  #print(len(Eexit),len(PDGexit),len(xexit))
  PDGexit_filter = np.equal(PDGexit,2212) #first filter for proton PDG
  Eexit_filtered = Eexit[PDGexit_filter]
  Eexit_filter = np.greater(Eexit_filtered,energy*Engcut) #then filter Eng with PDG and create Eng_filter that is filtered
  Eexit_filtered = Eexit_filtered[Eexit_filter]
  #XLfilter = np.less(Eexit_filtered,posmax) #10.May adding to try to match Twiss?
  #XGfilter = np.less(Eexit_filtered,-posmax)
  #YLfilter = np.less(Eexit_filtered,posmax)
  #YGfilter = np.less(Eexit_filtered,-posmax)

  xexit_filtered = xexit[PDGexit_filter]
  xexit_filtered = xexit_filtered[Eexit_filter]
  #xexit_filtered = xexit_filtered[XLfilter][XGfilter]

  pxexit_filtered = pxexit[PDGexit_filter]
  pxexit_filtered = pxexit_filtered[Eexit_filter]
  pxfilterL = np.less(pxexit_filtered,posmax) #[rad]
  pxexit_filtered = pxexit_filtered[pxfilterL]
  xexit_filtered = xexit_filtered[pxfilterL]
  pxfilterG = np.greater(pxexit_filtered,-5e-3)
  pxexit_filtered = pxexit_filtered[pxfilterG]
  xexit_filtered = xexit_filtered[pxfilterG]

  yexit_filtered = yexit[PDGexit_filter]
  yexit_filtered = yexit_filtered[Eexit_filter]
  #yexit_filtered = yexit_filtered[YLfilter][YGfilter]

  pyexit_filtered = pyexit[PDGexit_filter]
  pyexit_filtered = pyexit_filtered[Eexit_filter]
  #print(np.mean(pyexit_filtered),np.std(pyexit_filtered))
  lim=5e-3
  pyfilterL = np.less(pyexit_filtered,lim) #[rad]
  pyexit_filtered = pyexit_filtered[pyfilterL]
  yexit_filtered = yexit_filtered[pyfilterL]
  pyfilterG = np.greater(pyexit_filtered,-lim) #[rad]
  pyexit_filtered = pyexit_filtered[pyfilterG]
  yexit_filtered = yexit_filtered[pyfilterG]

  exitxTwf = calcTwiss("Exit X Filtered","Exit X' Filtered",xexit_filtered,pxexit_filtered)
  exityTwf = calcTwiss("Exit Y Filtered","Exit Y' Filtered",yexit_filtered,pyexit_filtered)

  #Below are the distributions taken from the ESS Target location (the detector located 5m from PBW)
  myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
  myTree= myFile.Get("TrackerHits") #TrackerHits has all trackers, be sure to only have 1!

  #Now get the "targ" distributions
  xtarg = np.zeros(myTree.GetEntries()) #dynamic length arrays
  pxtarg = np.zeros(myTree.GetEntries())
  ytarg = np.zeros(myTree.GetEntries())
  pytarg = np.zeros(myTree.GetEntries())
  pztarg = np.zeros(myTree.GetEntries())
  Etarg = np.zeros(myTree.GetEntries())
  PDGtarg = np.zeros(myTree.GetEntries())
  #print("The length of the arrays are ",len(xtarg))

  for i in range(myTree.GetEntries()): #put data in arrays
      myTree.GetEntry(i)
      pztarg[i] = myTree.pz
      if pztarg[i] == 0.0: continue
      xtarg[i] = myTree.x *mm
      pxtarg[i] = myTree.px/pztarg[i]
      ytarg[i] = myTree.y *mm
      pytarg[i] = myTree.py/pztarg[i]
      Etarg[i] = myTree.E
      PDGtarg[i] = myTree.PDG
      #break
  myFile.Close()

  #targxTwiss = calcTwiss("Target X","Target X'",xtarg,pxtarg)
  #targyTwiss = calcTwiss("Target Y","Target Y'",ytarg,pytarg)

  #For plotting Energy distribution of all species
  if engplot:
    import math
    mag=math.floor(math.log10(N)) #for dynamic title name
    engname=savename+"_EnergyPlot.png" #same for each plot
    titl = "Energy Distribution at ESS Target Through PBW of "+mat+"\n For All Particle Species"

    if material == "G4_Galactic":
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

  #Filter the relevant distributions
  #print(len(Etarg),len(PDGtarg),len(xtarg))
  PDGtarg_filter = np.equal(PDGtarg,2212) #first filter for proton PDG
  Etarg_filter = np.greater(Etarg[PDGtarg_filter],energy*Engcut) #then filter Eng with PDG and create Eng_filter that is filtered
  xtarg_filtered = xtarg[PDGtarg_filter]
  xtarg_filtered = xtarg_filtered[Etarg_filter]
  pxtarg_filtered = pxtarg[PDGtarg_filter]
  pxtarg_filtered = pxtarg_filtered[Etarg_filter]
  ytarg_filtered = ytarg[PDGtarg_filter]
  ytarg_filtered = ytarg_filtered[Etarg_filter]
  pytarg_filtered = pytarg[PDGtarg_filter]
  pytarg_filtered = pytarg_filtered[Etarg_filter]

  #targxTwissf = calcTwiss("Target X Filtered","Target X' Filtered",xtarg_filtered,pxtarg_filtered)
  #targyTwissf = calcTwiss("Target Y Filtered","Target Y' Filtered",ytarg_filtered,pytarg_filtered)

  #Display Full Energy distribution results
  print("Full Energy distribution of {:d} particles with minimum Energy {:.3f}MeV through ".format(len(Eexit),np.min(Eexit_filtered)),mat," PBW")

  #get Twiss parameters from tracker at ESS Target location:
  twiss, numPart,objects = miniScatterDriver.getData(savedfile)
  #Tnemx = twiss["tracker_cutoffPDG2212"]["x"]["eps"] *um #[um]
  #Tbetax = twiss["tracker_cutoffPDG2212"]["x"]["beta"] *m #[m]
  #Talphx = twiss["tracker_cutoffPDG2212"]["x"]["alpha"] #[um-mrad]
  #Tnemy = twiss["tracker_cutoffPDG2212"]["y"]["eps"] *um #[um]
  #Tbetay = twiss["tracker_cutoffPDG2212"]["y"]["beta"] *m
  #Talphy = twiss["tracker_cutoffPDG2212"]["y"]["alpha"] #[um-mrad]

  # Get Twiss for Target Exit
  #Enemtx = twiss["target_exit_cutoff"]["x"]["eps"] *um #[um]
  #Ebetax = twiss["target_exit_cutoff"]["x"]["beta"] *m #[m]
  #Ealphx = twiss["target_exit_cutoff"]["x"]["alpha"] #[um-mrad]
  #Enemty = twiss["target_exit_cutoff"]["y"]["eps"] *um #[um]
  #Ebetay = twiss["target_exit_cutoff"]["y"]["beta"] *m
  #Ealphy = twiss["target_exit_cutoff"]["y"]["alpha"] #[um-mrad]
  #delnemx = Enemtx - Inemx*um #[um] #passed #s aren't in units system
  #delnemy = Enemty - Inemy*um #[um]

  q=partmass*partmass/energy
  p = np.sqrt((energy*MeV*energy*MeV-partmass*MeV*partmass*MeV)/ (c*c))/(MeV/c)
  betap = beta_rel*p
  #betap = energy-q
  #thetasq = 13.6 * z / betap * np.sqrt(thick/radLen) * (1 + 0.038 * np.log(thick/radLen)) #from Eq 5
  thetasq = 13.6 * z / betap * np.sqrt(thick/radLen) * (1 + 0.038 * np.log(thick*z*z/(radLen*(1-q/energy)))) #from MiniScatter
  print("\nbetap: {:.3e}, gamma: {:.3f}, beta: {:.3f}, radLen: {:.3f}mm, theta^2: {:.3e} radians".format(betap,gamma_rel,beta_rel,radLen,thetasq))
    #using beta*q or energy-q produces similar #s (within 7%)

  #Plotting PBW Exit distribution vs the PDF produced from the Formalism Equations
  from plotFit import plotTwissFit,plotFit,calcEq8,calcEq16
  #plotFit(xs,    ys, savename,xlim,   ylim,material)
  #plotFit(xinit,yinit,savename+"init",   3,      0,material)
  #plotFit(xexit,yexit,savename+"texit",  3,      0,material)

  #Use list of Twiss values for simple passing of data: [beta,alpha,gemt]
  TwissIx   = [Ibetax,Ialphx,Igemx*um] #Initial Twiss
  TwissIy   = [Ibetay,Ialphy,Igemy*um]
  ##TwisstEx  = [Ebetax,Ealphx,Enemtx/(beta_rel*gamma_rel)]#target Exit Twiss
  ##TwisstEy  = [Ebetay,Ealphy,Enemty/(beta_rel*gamma_rel)]
  #Twisse8x  = calcEq8(thetasq, TwissIx,thick,beta_rel,gamma_rel) #calculated target exit Twiss
  #Twisse8y  = calcEq8(thetasq, TwissIy,thick,beta_rel,gamma_rel)
  #Twisse16x = calcEq16(thetasq,TwissIx,thick,beta_rel,gamma_rel) #calculated 2 target exit Twiss
  #Twisse16y = calcEq16(thetasq,TwissIy,thick,beta_rel,gamma_rel)
  #print("Twiss=[beta,alpha,Igemx]")
  #print("Func\nEq8 Twiss X:",calcEq8(thetasq,TwissIx,thick,beta_rel,gamma_rel),"\nEq8 Twiss Y:",calcEq8(thetasq,TwissIy,thick,beta_rel,gamma_rel))
  #print("Eq16 Twiss X:",calcEq16(thetasq,TwissIx,thick,beta_rel,gamma_rel),"\nEq16 Twiss Y:",calcEq16(thetasq,TwissIy,thick,beta_rel,gamma_rel))
  
  #Displays Twiss values and calculated mu and sigma
  print("\nPre PBW Twiss")
  plotTwissFit(xinit_filtered/mm,pxinit_filtered,savename+"init",mat,"Pre PBW","X",thick,thetasq,beta_rel,gamma_rel)
  
  print("\nPBW Exit Twiss Calculated")
  plotTwissFit(xexit_filtered/mm,pxexit_filtered,savename+"texitFiltered",mat,"PBW Exit","X",thick,thetasq,beta_rel,gamma_rel)
  plotTwissFit(yexit_filtered/mm,pyexit_filtered,savename+"texitFiltered",mat,"PBW Exit","Y",thick,thetasq,beta_rel,gamma_rel)
  
  #plotTwissFit2(Twisstex,Twisstey,xexit_filtered/mm,pxexit_filtered,yexit_filtered/mm,pyexit_filtered,savename+"texitFiltered",material,"Filtered",exitxTwf,exityTwf,"PBW Exit")
  #print("\nEq8 PBW Exit Twiss")
  #plotTwissFit(Twisse8x,Twisse8y,xexit_filtered/mm,pxexit_filtered,yexit_filtered/mm,pyexit_filtered,savename+"texitEq8",material,"Formalism Eq 8",Twisstex,Twisstey,"PBW Exit")
  #print("\nEq16 PBW Exit Twiss")
  #plotTwissFit(Twisse16x,Twisse16y,xexit_filtered/mm,pxexit_filtered,yexit_filtered/mm,pyexit_filtered,savename+"texitEq16",material,"Formalism Eq 16",Twisstex,Twisstey,"PBW Exit")

  #return savename, xtarg/mm, ytarg/mm
  return savename, xexit/mm, yexit/mm #bc why?