#simulation.py
#Eric Fackelman
#29 March - 1 June 2022

#This is to run the simulation of a beam through a PBW of input material 
# and output the position X and Y arrays at ESS Target location.
# Also can plot the Energy Distribution if requested.

def simulation(N,material,beam,thick,Inemx,Inemy,Ialphx,Ialphy,Ibetax,Ibetay,energy,zoff,Engcut,engplot):
  import numpy as np
  import ROOT
  import os
  import sys
  from plotFit import calcTwiss

  #constants for below use
  #particle characteristic values
  if beam == "proton":
    partA = 938.27209 #[MeV/c2]
    partZ = 1
  elif beam == "electron":
    partA = 0.511 #[MeV/c2]
    partZ = 1
  c = 2.99792e8 #[m/s]
  MeV = 1e6*1.602e-19 
  um = 1e-6 #[m] #need to convert to real units as the equations use real units.
  m = 1 #[m]
  mm = 1e-3 #[m]
  ummrad = um*1e-3
  gamma_rel = 1 + energy/partA #from PrintTwissParameters
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

  baseSimSetup["N"]     = N #Just a few events here! Remember that thicker targets are slower
  baseSimSetup["POSLIM"] = 100 #XY histogram Position Limit for a few, check RootFileWriter.cc
  #print(baseSimSetup)

  #Define material nickname
  #radiation lengths are from https://pdg.lbl.gov/2019/AtomicNuclearProperties/
  if material == "G4_Galactic":
    mat = "Vac"
    radLen = 1e12 #[mm] basically infinity
    atomZ = 0.1
    atomA = 0.1
  elif material == "G4_Al":
    mat = "Al"
    radLen = 88.97 #[mm] #0.2401 #[g/mm^2] # #24.01 [g/cm^2]
    atomZ = 13
    atomA = 26.981
  elif material == "G4_AIR":
    mat = "Air"
    radLen = 3.122e5 #[mm] -taken from 80% N2 gas(3.260e5mm) and 20% O2 gas (2.571e5mm)
    atomZ = 2 #the lbl page says Z/A = 0.49919, so 2/4 is close enough
    atomA = 4
  elif material == "G4_Au":
    mat = "Au"
    radLen = 3.344
    atomZ = 79
    atomA = 196.966

  if thick!=0.0:
    baseSimSetup["THICK"] = thick
    baseSimSetup["MAT"] = material
    #Valid choices: G4_Al, G4_Au, G4_C, G4_Cu, G4_Pb, G4_Ti, G4_Si, G4_W, G4_U, G4_Fe, G4_MYLAR, G4_KAPTON,
    #G4_STAINLESS-STEEL, G4_WATER,G4_SODIUM_IODIDE, G4_Galactic, G4_AIR, Sapphire, ChromoxPure, ChromoxScreen

    #Detector distance from target center [mm] Default is 50mm behind Target
    #For multiple detector locations, make a list, e.g. [-5,5,5000] but they stack in TTree.
    baseSimSetup["DIST"] = [5000] #Detector location. only at ESS Target location
    outname = "simplePBW_"+str(baseSimSetup["THICK"])+"mm"+mat+"_N{:.0e}_betterESSbeam".format(baseSimSetup["N"])
  else:
    baseSimSetup["THICK"] = 0.0
    baseSimSetup["DIST"] = [5000] #Detector locations. At PBW Exit and ESS Target location 
    baseSimSetup["MAGNET"] = []
    #How to construct a magnet for miniScatterDriver, as per kyrsjo/MiniScatter/blob/master/examples/SiRi DeltaE-E detector.ipynb
    #Al first piece
    m1 = {}
    m1["pos"]      = 0.6250 #pos is at CENTER
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
    m2["pos"]      = 2.2501 #pos is at CENTER
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
    m3["pos"]      = 3.7502 #pos is at CENTER
    m3["type"]     = "TARGET"
    m3["length"]   = 1.000 #[mm]
    m3["gradient"] = 0.0
    m3["keyval"] = {}
    m3["keyval"]["material"] = "G4_Al"
    m3["keyval"]["width"]    = 200.0 #[mm]
    m3["keyval"]["height"]   = 100.0 #[mm]
    baseSimSetup["MAGNET"].append(m3)

    mat = "Real"
    radLenAl = 88.97 #Al
    radLenH2O = 36.08 #liquid Water 
    radLen = 64.08 #mm average of 2.25mm Al and 2.0mm Light Water
    atomZ = 13 + 8 + 2 / 4 #not great, but ok
    atomA = 26.9 + 16 + 2 / 4
    #from https://cds.cern.ch/record/1279627/files/PH-EP-Tech-Note-2010-013.pdf
    #radLen = 4.25/((.225*2.70)/8.897 + (2.0*1.0)/3.608) 

    outname = "simplePBW_"+str(thick)+"mm"+mat+"_N{:.0e}_betterESSbeam_full".format(baseSimSetup["N"])
    #outname = "simplePBW_N{:.0e}_betterESSbeam_real_".format(baseSimSetup["N"])+str(baseSimSetup["PHYS"])

  #Some output settings
  baseSimSetup["QUICKMODE"] = False #Include slow plots
  baseSimSetup["MINIROOT"]  = False #Skip TTRees in the .root files
  baseSimSetup["ANASCATTER"] = True #don't do Analytical Scatter Angle Test

  Rcut = 100.0
  baseSimSetup["EDEP_DZ"]   = 1.0 #Z bin width for energy deposit histogram
  baseSimSetup["CUTOFF_RADIUS"] = Rcut #Larger radial cutoff #Decreased 10 May
  baseSimSetup["CUTOFF_ENERGYFRACTION"] = Engcut #Minimum percent of full Energy to use in cutoff calculations
  #0.95 is default, but let's be explicit

  #Store the .root files in a subfolder from where this script is running,
  # normally MiniScatter/examples, in order to keep things together
  baseSimSetup["OUTFOLDER"] = os.path.join("/scratch/ericdf/Scratch/PBWScatter/") 
  #put in Scratch of HepLab0# for faster processing, as per Kyrre

  #copy so it is if running multiple scans in a Jupyter notebook
  simSetup_simple1 = baseSimSetup.copy()

  #Give the .root file a dynamic name, DONE IN MAGNET IF !!!
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

  #Filter with Energy for now as PDG isn't working yet
  Einit_filter = np.greater(Einit,energy*Engcut)
  xinit_filtered = xinit[Einit_filter]
  pxinit_filtered = pxinit[Einit_filter]
  yinit_filtered = yinit[Einit_filter]
  pyinit_filtered = pyinit[Einit_filter]

  initxTwf = calcTwiss("Initial X Filtered","Initial X' Filtered",xinit_filtered,pxinit_filtered)
  inityTwf = calcTwiss("Initial Y Filtered","Initial Y' Filtered",yinit_filtered,pyinit_filtered)

  ##Get the distributions from the PBW exit face
  if thick != 0:
    #27.4 just added target-exit pull in and changed previous xexit arrays to target arrays!
    #Now get the "target-exit" distributions for plotting with the Formalism distribution below. 
    #These are not at the ESS Target location, but are at the far side of the PBW
    myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
    myTree= myFile.Get("TargetExit") #TrackerHits has all trackers, be sure to only have 1!
  else:
    myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
    myTree = myFile.Get("magnet_3_ExitHits")
    thick=4.25 #[mm] set here for thickness calculations

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
      myTree.GetEntry(i)
      pzexit[i] = myTree.pz
      if pzexit[i] == 0.0:
        print("warning: PZexit[{}]==0".format(i))
        continue #11.5.22 recommended by Kyrre
      xexit[i] = myTree.x *mm #m
      pxexit[i] = myTree.px / pzexit[i] #from Kyrre 5.5.22 to make it true X'!
      yexit[i] = myTree.y *mm #m
      pyexit[i] = myTree.py/pzexit[i]
      Eexit[i] = myTree.E
      PDGexit[i] = myTree.PDG
  myFile.Close() 

  #Filter the relevant distributions
  PDGexit_filter = np.equal(PDGexit,2212) #first filter for proton PDG
  Eexit_filtered = Eexit[PDGexit_filter]
  Eexit_filter = np.greater(Eexit_filtered,energy*Engcut) #then create Eng_filter that is filtered
  Eexit_filtered = Eexit_filtered[Eexit_filter]

  angmax=4e-3 #[rad] one angle filter limit 
  #Apply PDG, Energy Filters
  xexit_filtered = xexit[PDGexit_filter][Eexit_filter]
  yexit_filtered = yexit[PDGexit_filter][Eexit_filter]
  pxexit_filtered = pxexit[PDGexit_filter][Eexit_filter]
  pyexit_filtered = pyexit[PDGexit_filter][Eexit_filter]
  #Apply >,< filters
  #X, <
  pxfilterL = np.less(pxexit_filtered,angmax) #[rad]
  pxexit_filtered = pxexit_filtered[pxfilterL]
  xexit_filtered = xexit_filtered[pxfilterL]
  #Y, <
  pyfilterL = np.less(pyexit_filtered,angmax) #[rad]
  pyexit_filtered = pyexit_filtered[pyfilterL]
  yexit_filtered = yexit_filtered[pyfilterL]
  #X, >
  pxfilterG = np.greater(pxexit_filtered,-angmax) #[rad]
  pxexit_filtered = pxexit_filtered[pxfilterG]
  xexit_filtered = xexit_filtered[pxfilterG]
  #Y, >
  pyfilterG = np.greater(pyexit_filtered,-angmax) #[rad]
  pyexit_filtered = pyexit_filtered[pyfilterG]
  yexit_filtered = yexit_filtered[pyfilterG]

  #Get Twiss for the filtered distributions
  exitxTwf = calcTwiss("Exit X Filtered","Exit X' Filtered",xexit_filtered,pxexit_filtered)
  exityTwf = calcTwiss("Exit Y Filtered","Exit Y' Filtered",yexit_filtered,pyexit_filtered)


  ##Below are the distributions taken from the ESS Target location (the detector located 5m from PBW)
  myFile = ROOT.TFile.Open(os.path.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
  myTree= myFile.Get("TrackerHits") #TrackerHits has all trackers, be sure to only have 1!

  #Get the distributions at the ESS Target location
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
      if pzexit[i] == 0.0:
        print("warning: PZexit[{}]==0".format(i))
        continue #11.5.22 recommended by Kyrre
      xtarg[i] = myTree.x *mm
      pxtarg[i] = myTree.px / pztarg[i] #from Kyrre 5.5.22 to make it true X'!
      ytarg[i] = myTree.y *mm
      pytarg[i] = myTree.py / pztarg[i]
      Etarg[i] = myTree.E
      PDGtarg[i] = myTree.PDG
  myFile.Close()

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

  #Display Full Energy distribution results
  print("Full Energy distribution of {:d} particles with minimum Energy {:.3f}MeV through ".format(len(Eexit),np.min(Eexit_filtered)),mat," PBW")

  #Filter the relevant distributions
  PDGtarg_filter = np.equal(PDGtarg,2212) #first filter for proton PDG
  Etarg_filter = np.greater(Etarg[PDGtarg_filter],energy*Engcut) #then filter Eng with PDG and create Eng_filter that is filtered
  #These proton only arrays are returned to the original script!
  xtarg_filtered_p = xtarg[PDGtarg_filter]
  pxtarg_filtered_p = pxtarg[PDGtarg_filter]
  ytarg_filtered_p = ytarg[PDGtarg_filter]
  pytarg_filtered_p = pytarg[PDGtarg_filter]

  angmax=6e-3 #[rad] one angle filter limit
  #Apply PDG, Energy Filters
  xtarg_filtered = xtarg[PDGtarg_filter][Etarg_filter]
  #xtarg_filtered = xtarg_filtered[Etarg_filter]
  pxtarg_filtered = pxtarg[PDGtarg_filter][Etarg_filter]
  #pxtarg_filtered = pxtarg_filtered[Etarg_filter]
  ytarg_filtered = ytarg[PDGtarg_filter][Etarg_filter]
  #ytarg_filtered = ytarg_filtered[Etarg_filter]
  pytarg_filtered = pytarg[PDGtarg_filter][Etarg_filter]
  #pytarg_filtered = pytarg_filtered[Etarg_filter]
  #Apply >,< filters
  #X, <
  pxfilterL = np.less(pxtarg_filtered,angmax) #[rad]
  pxtarg_filtered = pxtarg_filtered[pxfilterL]
  xtarg_filtered = xtarg_filtered[pxfilterL]
  #Y, <
  pyfilterL = np.less(pytarg_filtered,angmax) #[rad]
  pytarg_filtered = pytarg_filtered[pyfilterL]
  ytarg_filtered = ytarg_filtered[pyfilterL]
  #X, >
  pxfilterG = np.greater(pxtarg_filtered,-angmax) #[rad]
  pxtarg_filtered = pxtarg_filtered[pxfilterG]
  xtarg_filtered = xtarg_filtered[pxfilterG]
  #Y, >
  pyfilterG = np.greater(pytarg_filtered,-angmax) #[rad]
  pytarg_filtered = pytarg_filtered[pyfilterG]
  ytarg_filtered = ytarg_filtered[pyfilterG]

  #Get Twiss for the filtered distributions
  targxTwissf = calcTwiss("Target X Filtered","Target X' Filtered",xtarg_filtered,pxtarg_filtered)
  targyTwissf = calcTwiss("Target Y Filtered","Target Y' Filtered",ytarg_filtered,pytarg_filtered)

  ##Highland Equation Radiation Length Calculation
  p = np.sqrt((energy+partA)**2 - (partA)**2) #[MeV/c] #derived with Kyrre 15.6.22
  betap = beta_rel*p #Eq 5
  thetasq = 13.6 * partZ / betap * np.sqrt(thick/radLen) * (1 + 0.038 * np.log(thick/radLen)) #from Eq 5
  #q=partA*partA/energy #used for mimicing AnaScatter's "BeamBeta2" RadLen Calculatons
  #betap = energy-q 
  #thetasq = 13.6 * partZ / betap * np.sqrt(thick/radLen) * (1 + 0.038 * np.log(thick*partZ**2/(radLen*(1-q/energy)))) #from MiniScatter
  print("\nradLen: {:.2f}, p: {:.3e}, gamma: {:.3f}, beta: {:.3f}, theta^2: {:.3e} radians".format(radLen,p,gamma_rel,beta_rel,thetasq))
  #using beta*q or energy-q produces similar #s (within 7%)

  ##Characteristic and Screening Angle calculations, 14/15 June 22
  # as per https://doi-org.ezproxy.uio.no/10.1016/0168-583X(91)95671-Y
  centerFrac = [0.7,0.9,0.95,0.98,0.99,0.995,0.997]
  characAng = 0.157 * (atomZ * (atomZ+1) * thick/atomA) * (partZ / (p * beta_rel))**2
  screenAng = 2.007e-5 * (atomZ**(2.0/3)) * (1 + 3.34 * (atomZ * partZ * (1/137) / beta_rel)**2 ) / (p*p)
  meanScatters = characAng / screenAng
  nu = 0.5 * meanScatters / (1-centerFrac[3]) 
  sigma = characAng / (1-centerFrac[3]**2) * ((1+nu)/nu * np.log(1+nu) -1)
  print("For fraction {}, characteristic angle {:0.2e}, screening angle {:0.2e}".format(centerFrac[3],characAng,screenAng))
  print("We get {:0.2e} scatters, nu {:0.2e} and sigma {:0.2e}".format(meanScatters,nu,sigma))

  ##If magnet, use multiple scattering layers instead of averaging!
  from plotFit import plotTwissFit,calcEq8,calcEq16
  #Use list of Twiss values for simple passing of data: [beta,alpha,gemt]
  TwissIx   = [Ibetax,Ialphx,Igemx*um,((1+Ialphx*Ialphx)/Ibetax)] #Initial Twiss
  TwissIy   = [Ibetay,Ialphy,Igemy*um,((1+Ialphy*Ialphy)/Ibetay)]
  ##TwisstEx  = [Ebetax,Ealphx,Enemtx/(beta_rel*gamma_rel)]#target Exit Twiss
  ##TwisstEy  = [Ebetay,Ealphy,Enemty/(beta_rel*gamma_rel)]

  if baseSimSetup["THICK"] == 0.0:
    m1Len = baseSimSetup["MAGNET"][0]["length"]
    #Al Front contribution
    thetasqAl1 = 13.6 * partZ / betap * np.sqrt(m1Len/radLenAl) * (1 + 0.038 * np.log(m1Len/radLenAl))
    Twisse8xAl1 = calcEq8(thetasqAl1, TwissIx,m1Len,beta_rel,gamma_rel)
    Twisse8yAl1 = calcEq8(thetasqAl1, TwissIy,m1Len,beta_rel,gamma_rel)

    #H2O contribution
    m2Len = baseSimSetup["MAGNET"][1]["length"]
    thetasqH2O = 13.6 * partZ / betap * np.sqrt(m2Len/radLenH2O) * (1 + 0.038 * np.log(m2Len/radLenH2O))
    Twisse8xH2O = calcEq8(thetasqH2O, Twisse8xAl1,m2Len,beta_rel,gamma_rel)
    Twisse8yH2O = calcEq8(thetasqH2O, Twisse8yAl1,m2Len,beta_rel,gamma_rel)

    #Al Back contribution
    m3Len = baseSimSetup["MAGNET"][2]["length"]
    thetasqAl2 = 13.6 * partZ / betap * np.sqrt(m3Len/radLenAl) * (1 + 0.038 * np.log(m3Len/radLenAl))
    Twisse8xAl2 = calcEq8(thetasqAl2, Twisse8xH2O,m3Len,beta_rel,gamma_rel)
    Twisse8yAl2 = calcEq8(thetasqAl2, Twisse8yH2O,m3Len,beta_rel,gamma_rel)

  elif baseSimSetup["THICK"] != 0.0: #if only one layer (MiniScatter "target")
    Twisse8x  = calcEq8(thetasq, TwissIx,thick,beta_rel,gamma_rel) #calculated target exit Twiss
    Twisse8y  = calcEq8(thetasq, TwissIy,thick,beta_rel,gamma_rel)
    Twisse16x = calcEq16(thetasq,TwissIx,thick,beta_rel,gamma_rel) #calculated 2 target exit Twiss
    Twisse16y = calcEq16(thetasq,TwissIy,thick,beta_rel,gamma_rel)

  ##Plotting PBW Exit distribution vs the PDF produced from the Formalism Equations
    #Displays Twiss values and calculated mu and sigma
  #print("\nPre PBW Twiss")
  #plotTwissFit(xinit_filtered/mm,pxinit_filtered,savename+"init",mat,"Pre PBW","X",thick,thetasq,beta_rel,gamma_rel,TwissIx)
  #plotTwissFit(yinit_filtered/mm,pyinit_filtered,savename+"init",mat,"Pre PBW","Y",thick,thetasq,beta_rel,gamma_rel,TwissIy)
  #print("\nPBW Exit Twiss Calculated")
  plotTwissFit(xexit_filtered/mm,pxexit_filtered,savename+"texitFiltered",mat,"PBW Exit","X",thick,thetasq,beta_rel,gamma_rel,TwissIx)
  #plotTwissFit(yexit_filtered/mm,pyexit_filtered,savename+"texitFiltered",mat,"PBW Exit","Y",thick,thetasq,beta_rel,gamma_rel,TwissIy)

  #Extension to Target. Needed for compareTargets!
  from plotFit import toTarget
  #initTargx = toTarget(TwissIx,"initX")
  #intiTargy = toTarget(TwissIy,"initY")
  #exitTargx = toTarget(exitxTwf,"exitX")
  #exitTargy = toTarget(exityTwf,"exitY")
  #e8Targx = toTarget(Twisse8x,"e8X")
  #e8Targy = toTarget(Twisse8y,"e8Y")
  #e16Targx = toTarget(Twisse16x,"e16X")
  #e16Targy = toTarget(Twisse16y,"e16Y")
  e8TargxReal = toTarget(Twisse8xAl2,"e8XReal")
  e8TargyReal = toTarget(Twisse8yAl2,"e8YReal")

  #Now compare the MiniScatter Target distribution (targxTwissf) to initTarg, exitTarg, e8Targ and e16Targ PDFs
  from plotFit import compareTargets
  #compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,initTargx,intiTargy,"Initial",savename,mat)
  #compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,exitTargx,exitTargy,"PBW Exit",savename,mat)
  #compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,e8Targx,e8Targy,"Eq 8",savename,mat)
  #compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,e16Targx,e16Targy,"Eq 16",savename,mat)
  compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,e8TargxReal,e8TargyReal,"Eq 8",savename,mat)
  
  return savename, xtarg_filtered_p/mm, ytarg_filtered_p/mm #filter by PDG only
  #return savename, xexit/mm, yexit/mm #bc why?