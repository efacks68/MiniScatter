##simulation.py
##29 March 2022 - 27 June 2023
##For holding configuration and TTree functions

##setup - For loading simulation configuration for use in MiniScatter
def setup(args,mat,beamFile,Twiss,options,paths):
    from numpy import cos,sin
    from os import uname,getcwd,chdir,path as osPath
    from sys import path as sysPath
    from datetime import datetime
    print("setup start",datetime.now())

    ##Setup MiniScatter -- modify the path to where you built MiniScatter
    MiniScatter_path="../../MiniScatter/build/."
    sysPath.append(MiniScatter_path)
    #print(getcwd())
    if uname()[1] in {"tensor.uio.no", "heplab01.uio.no", "heplab04.uio.no","heplab03.uio.no"}:
        if getcwd() != "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/":
            chdir("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/")
    elif uname()[1] == "mbef-xps-13-9300":
        if getcwd() != "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/":
            chdir("/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/MiniScatter/build/")
    else: chdir(input("What directory would you like to save files to? "))
    #print(getcwd())

    ##Make Dictionary for holding options for MiniScatter command, joined together in miniScatterDriver
    baseSimSetup = {}
    baseSimSetup["PHYS"]  = options['physList'] #better for scattering through 1mm sheets
    baseSimSetup["ZOFFSET"]   = options["zoff"]
    baseSimSetup["WORLDSIZE"] = args.rCut #[mm] Make the world wider for seeing where particles go
    
    ##Plot Output settings
    baseSimSetup["QUICKMODE"]  = False #Include slow plots
    baseSimSetup["MINIROOT"]   = options['MiniRoot'] #Skip TTRees in the .root files
    baseSimSetup["ANASCATTER"] = True #don't do Analytical Scatter Angle Test
    baseSimSetup["EDEP_DZ"]    = 1.0 #Z bin width for energy deposit histogram
    baseSimSetup["CUTOFF_RADIUS"] = args.rCut #[mm]Larger radial cutoff #Decreased 10 May
    baseSimSetup["CUTOFF_ENERGYFRACTION"] = args.engCut #Minimum percent of full Energy to use in cutoff calculations
    baseSimSetup["POSLIM"] = args.rCut #XY histogram Position Limit for a few, check RootFileWriter.cc
    #Detector distance from target center [mm] Default is 50mm behind Target
    #For multiple detector locations, make a list, e.g. [-5,5,5000] but they stack in TTree.
    baseSimSetup["DIST"] = [3565] #Detector locations. At ESS Target location, update from Synoptic Viewer 15-3-15

    ##For loading particles from external (e.g. with runARasterMaker)
    if args.sim == "raster" or args.sim == "map":
        from plotFit import numLines
        #print("numLines start",datetime.now())
        parts = numLines(beamFile)
        #print("numLines end",datetime.now())
        #print(parts,"in sim setup N")
        baseSimSetup["N"]        = parts #change to match file particles. Used for file name
        baseSimSetup["BEAMFILE"] = beamFile+".csv" # number of particles >= N
        baseSimSetup["ENERGY"]   = args.energy #570 #[MeV] #ESS beam energy update 15.8
    else: ##particles generated now by Geant4 from a provided Twiss
        baseSimSetup["BEAM"]    = args.particle
        baseSimSetup["ENERGY"]  = args.energy #570 #[MeV] #ESS beam energy update 15.8
        baseSimSetup["COVAR"]   = (Twiss[0],Twiss[1],Twiss[2],Twiss[3],Twiss[4],Twiss[5]) #only for without "--beamFile"
        baseSimSetup["N"]       = args.Nbeamlet #N per beamlet when making only one beamlet
    #print(baseSimSetup)

    ##Define material name from nickname
    if mat == "Vac":
        material = "G4_Galactic"
    elif mat == "Al":
        material = "G4_Al"
    elif mat == "Air":
        material = "G4_AIR"
    elif mat == "Au":
        material = "G4_Au"

    ##What type of PBW to use?
    if args.t != 0.0: ##Flat sheet of material
        baseSimSetup["THICK"] = args.t
        baseSimSetup["MAT"] = material
        #Valid choices: G4_Al, G4_Au, G4_C, G4_Cu, G4_Pb, G4_Ti, G4_Si, G4_W, G4_U, G4_Fe, G4_MYLAR, G4_KAPTON,
        #G4_STAINLESS-STEEL, G4_WATER,G4_SODIUM_IODIDE, G4_Galactic, G4_AIR, Sapphire, ChromoxPure, ChromoxScreen

        ##Make file name from provided beamFile name
        if beamFile != "":
            import re
            name = re.sub(".+(PBW)",mat,beamFile) #change file name
            #print(name)
            outname = name + "_run"
        else: ##Make new name for simple PBW
            outname = "simplePBW_"+str(baseSimSetup["THICK"])+"mm"+mat+"_{:.0f}MeV_emtx{:.0f}um".format(baseSimSetup["ENERGY"],Twiss[1]*1e3)
    
    else: ##Real curved PBW
        baseSimSetup["THICK"] = 0.0
        baseSimSetup["MAGNET"] = []
        ##How to construct a magnet for miniScatterDriver, as per kyrsjo/MiniScatter/blob/master/examples/SiRi DeltaE-E detector.ipynb
        ##Specialized PBW magnet! See MiniScatter/src/MagnetPBW.cc for details
        m1 = {}
        arcPhi = 120
        m1["type"]     = "PBW"
        m1["length"]   = 0 #[mm] Must be 0!
        m1["gradient"] = 0.0
        m1["keyval"] = {}
        m1["keyval"]["material"]   = material
        m1["keyval"]["radius"]     = 88.0 #[mm]
        m1["keyval"]["al1Thick"]   = 1.0  #[mm]
        m1["keyval"]["waterThick"] = 2.0  #[mm]
        m1["keyval"]["al2Thick"]   = 1.25 #[mm]
        #m1["keyval"]["width"]      = 200 #[mm] #(?)
        ##For thickness variation simulations, this modifies the actual PBW to be correct thickness.
        if args.sim == "thick":
            m1["keyval"]["al1Thick"]   = options['PBWT']*0.5 #[mm]
            m1["keyval"]["al2Thick"]   = options['PBWT']*0.5 #[mm]

        ##For Analytical Formular calculations(l505-657),
        ## need to take into account total thickness and modify position of PBW
        from math import pi
        totalThickness = m1["keyval"]["al1Thick"] + m1["keyval"]["waterThick"] + m1["keyval"]["al2Thick"]
        m1["pos"]      = (m1["keyval"]["radius"]*cos(arcPhi/2*pi/180)+totalThickness)/2 #[mm] Z0=24.125mm for r=88m,t=4.25,arcPhi=120
        #print("z0=",m1["pos"])
        baseSimSetup["MAGNET"].append(m1)

        ##Start output name for real PBW
        if beamFile != "": #from beamFile name which was made in runARasterMaker, could be for twiss
            #print("t!=0 beamFile",beamFile)
            outname = beamFile+"_runW"
        else: #This would get used for a beamlet with real PBW
            outname = "PBW_{:.0f}MeV_eX{:.2f},eY{:.2f}um_bX{:.2f},bY{:.2f}m_aX{:.2f},aY{:.2f}_N{:.0e}".format(args.energy,Twiss[0],Twiss[3],Twiss[1],Twiss[4],Twiss[2],Twiss[5],baseSimSetup["N"])
            #if options['MiniRoot']:
            #    outname+="_miniR"
        #print(mat)
        if mat != "Al":
            outname=outname+"_"+mat
            print("material not Al! ",mat,outname)

        ##If including PBIP, include in output name and add to MiniScatter command
        if args.PBIP:
            outname = outname + "_PBIP"
            #PBIP Magnet
            m2 = {}
            m2["pos"] = 1874.0 #[mm]
            m2["type"] = "COLLIMATORRECT"
            m2["length"]   = 450.0 #[mm]
            m2["gradient"] = 0.0
            m2["keyval"] = {}
            m2["keyval"]["material"] = "G4_Al"
            m2["keyval"]["holeWidth"]   = 200.0 #[mm]
            m2["keyval"]["holeHeight"]   = 80.0 #[mm]
            m2["keyval"]["absorberWidth"]    = 950.0 #[mm]
            m2["keyval"]["absorberHeight"]   = 950.0 #[mm]
            baseSimSetup["MAGNET"].append(m2)

    ##Append physics List to end of output name
    if args.physList == "QGSP_BERT_EMZ":
        outname = outname + "_QBZ"
    elif args.physList == "FTFP_BERT_EMZ":
        outname = outname + "_FBZ"
    elif args.physList == "QGSP_BIC_EMZ":
        outname = outname + "_QBICZ"
    elif args.physList == "FTF_BIC_EMZ":
        outname = outname + "_FBICZ"
    elif args.physList == "FTFP_BIC_EMZ":
        outname = outname + "_QBBC"
    else:
        outname = outname + options['physList']

    ##Remove upper directories that may have come with beamFile for appending outname to scratch folder
    import re
    if re.search("/CSVs/",outname):
        #print("\n",outname,"\n")
        outname = re.sub(".+(?<=(/CSVs/))","",outname) #removes all preceeding real name
        #print("Directories removed. Now:",outname)

    ##Make name accessible later on
    args.beamFile = outname

    ##Find which folder root file is in and check if complete - not yet working
    #loc = findRoot(savename) #still need to work on this function
    #Where to save/load root and Picture files

    ##Set output folder depending on configuration
    if uname()[1] == "mbef-xps-13-9300": #my laptop
        baseSimSetup["OUTFOLDER"] = osPath.join(paths['scratchPath'])
    elif uname()[1] in {"tensor.uio.no", "heplab01.uio.no", "heplab04.uio.no","heplab03.uio.no"}: #UiO desktop
        if args.source == "twiss":
            baseSimSetup["OUTFOLDER"] = osPath.join(paths['scratchPath']+"failures/")
        elif args.source == "particles":
            baseSimSetup["OUTFOLDER"] = osPath.join(paths['scratchPath']+"ESS/")
        elif Twiss[1] < 1:
            baseSimSetup["OUTFOLDER"] = osPath.join(paths['scratchPath']+"pencil/")
        elif args.source == "csv":
            baseSimSetup["OUTFOLDER"] = osPath.join(paths['scratchPath']+"ESS/")
        elif options['dependence'] == "RA":
            baseSimSetup["OUTFOLDER"] = osPath.join(paths['scratchPath']+"2Dmaps/")
        else:
            baseSimSetup["OUTFOLDER"] = osPath.join(paths['scratchPath'])
    else: print("Help! Unknown build directory!, simulation.py l 243")

    ##Copy the configuration in case it is if running multiple scans in a Jupyter notebook
    simSetup_simple1 = baseSimSetup.copy()
    #print(outname,"\n")
    simSetup_simple1["OUTNAME"] = outname #"PBW_570MeV_pencil_N1e+05"#
    savename=paths['statsPWD']+outname #base savename for plots downstream, brings directly to my directory
    #savedfile=osPath.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root"
    #print("setup end",datetime.now())#savename, savedfile)

    return savename,simSetup_simple1

##simulation - For loading the TTrees for various plots to show scattered beam 
def simulation(args,mat,beamXAngle,beamYAngle,beamFile,Twiss,options,boxes,paths):
    import numpy as np
    import ROOT
    from os import path as osPath
    from sys import path as sysPath
    from plotFit import calcTwiss
    from datetime import datetime
    um = 1e-6 #[m] 
    mm = 1e-3 #[m]

    ##Define material nickname
    ##radiation lengths are from https://pdg.lbl.gov/2019/AtomicNuclearProperties/
        ##also https://cds.cern.ch/record/1279627/files/PH-EP-Tech-Note-2010-013.pdf
    #print(material)
    radLenH2O = 360.8 #[mm] liquid Water 
    if mat == "Vac":
        radLen = 1e24 #[mm] basically infinity
        atomZ = 0.1
        atomA = 0.1
    elif mat == "Al":
        radLen = 88.97 #[mm]
        atomZ = 13
        atomA = 26.981
    elif mat == "Air":
        radLen = 3.122e5 #[mm] -taken from 80% N2 gas(3.260e5mm) and 20% O2 gas (2.571e5mm)
        atomZ = 2 #the lbl page says Z/A = 0.49919, so 2/4 is close enough
        atomA = 4
    elif mat == "Au":
        radLen = 3.344
        atomZ = 79
        atomA = 196.966

    from simulation import setup
    savename,simSetup_simple1 = setup(args,mat,beamFile,Twiss,options,paths)
    #print(savename,simSetup_simple1)
    MiniScatter_path="../../MiniScatter/build/."
    sysPath.append(MiniScatter_path)
    import miniScatterDriver
    #import miniScatterScanner
    #import miniScatterPlots

    ##Constants for loading with MiniScatter
    QUIET     = False #Reduced output, doesn't show events
    ### Basic simulation parameters ###
    TRYLOAD = True  #Try to load already existing data instead of recomputing, only if using getData_TryLoad function.
    NUM_THREADS = 8 #Number of parallel threads to use for scans
    #Where to store temporary data for scans (a fast file system, NOT EOS/AFS)
    TMPFOLDER = "/tmp/miniScatter/SimpleDemo_thicknessScan"

    ##Set variable place holders for below use
    xtarg_filtered_p = np.zeros(10)
    ytarg_filtered_p = np.zeros(10)
    e8SigX = 0
    e8SigY = 0
    targSigX = 0
    targSigY = 0

    #Run simulation or load old simulation root file!

    if options['targetTree'] or options['initTree'] or options['exitTree']:
        #miniScatterDriver.runScatter(simSetup_simple1, quiet=QUIET) #this was Kyrre's, but it wasn't even trying to load old runs
        miniScatterDriver.getData_tryLoad(simSetup_simple1,quiet=QUIET)
        #print("Simulation Finished",datetime.now().strftime("%H-%M-%S"))

    #For accessing the initial spread TTrees
    if options['initTree']:
        myFile = ROOT.TFile.Open(osPath.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
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
            xinit[i] = myTree.x *mm #[m]?
            pxinit[i] = myTree.px
            yinit[i] = myTree.y *mm #[m]?
            pyinit[i] = myTree.py
            Einit[i] = myTree.E
        myFile.Close()

        ##Save Particle spreaf with momentum
        if args.saveParts:
            from plotFit import printParticles
            printParticles(savename,xinit,pxinit,yinit,pyinit,Einit)

    ##Get the particle distributions from the PBW exit face TTrees
    if options['exitTree']: #save some time by only getting when needed
        if args.t != 0:
            myFile = ROOT.TFile.Open(osPath.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
            myTree= myFile.Get("TargetExit") #TrackerHits has all trackers, be sure to only have 1!
        else:
            myFile = ROOT.TFile.Open(osPath.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
            myTree = myFile.Get("magnet_1_ExitHits")

        xexit = np.zeros(myTree.GetEntries()) #dynamic length arrays
        pxexit = np.zeros(myTree.GetEntries())
        yexit = np.zeros(myTree.GetEntries())
        pyexit = np.zeros(myTree.GetEntries())
        pzexit = np.zeros(myTree.GetEntries())
        Eexit = np.zeros(myTree.GetEntries())
        PDGexit = np.zeros(myTree.GetEntries())

        #import warnings
        #warnings.filterwarnings("error")

        ##Read in position and momentum of particles
        for i in range(myTree.GetEntries()): #put data in arrays
            myTree.GetEntry(i)
            pzexit[i] = myTree.pz
            if pzexit[i] == 0.0: #Catch for 0 momentum particles
                print("warning: PZexit[{}]==0".format(i))
                continue #11.5.22 recommended by Kyrre
            xexit[i] = myTree.x *mm #[m] ?
            pxexit[i] = myTree.px / pzexit[i] #from Kyrre 5.5.22 to make it true X'!
            yexit[i] = myTree.y *mm #[m] ?
            pyexit[i] = myTree.py / pzexit[i]
            Eexit[i] = myTree.E #[MeV]
            PDGexit[i] = myTree.PDG
        myFile.Close() 

        ##Filter the relevant distributions for protons above energy cut
        PDGexit_filter = np.equal(PDGexit,2212) #first filter for proton PDG
        Eexit_filtered = Eexit[PDGexit_filter]
        Eexit_filter = np.greater(Eexit_filtered,args.energy*simSetup_simple1["CUTOFF_ENERGYFRACTION"]) #then create Eng_filter that is filtered
        Eexit_filtered = Eexit_filtered[Eexit_filter]

        angmax=4e-3 #[rad] one angle filter limit 
        ##Apply PDG, Energy Filters
        xexit_filtered = xexit[PDGexit_filter][Eexit_filter]
        yexit_filtered = yexit[PDGexit_filter][Eexit_filter]
        pxexit_filtered = pxexit[PDGexit_filter][Eexit_filter]
        pyexit_filtered = pyexit[PDGexit_filter][Eexit_filter]
        ##Apply >,< filters
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

        ##Get Twiss for the filtered distributions
        exitxTwf = calcTwiss("Exit X Filtered","Exit X' Filtered",xexit_filtered,pxexit_filtered)
        exityTwf = calcTwiss("Exit Y Filtered","Exit Y' Filtered",yexit_filtered,pyexit_filtered)
        
        ##Drift the Twiss to the Target
        from plotFit import Drift
        exitTargx = Drift(exitxTwf,"exitX")
        exitTargy = Drift(exityTwf,"exitY")
    ##end of exit Distribution if

    #if args.PBIP:
    #  myFile = ROOT.TFile.Open(osPath.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
    #  myTree= myFile.Get("magnet_2_edeps")
    #Get depositions?

    ##Get TTree of particles at Target plane (BEW plane)
    if options['targetTree']:
        ##Distributions at ESS Target location (the detector located 4.4m from PBW)
        myFile = ROOT.TFile.Open(osPath.join(simSetup_simple1["OUTFOLDER"],simSetup_simple1["OUTNAME"])+".root")
        myTree= myFile.Get("TrackerHits") #TrackerHits has all trackers, be sure to only have 1!

        ##Get the distributions at the ESS Target location
        xtarg = np.zeros(myTree.GetEntries()) #dynamic length arrays
        pxtarg = np.zeros(myTree.GetEntries())
        ytarg = np.zeros(myTree.GetEntries())
        pytarg = np.zeros(myTree.GetEntries())
        pztarg = np.zeros(myTree.GetEntries())
        Etarg = np.zeros(myTree.GetEntries())
        PDGtarg = np.zeros(myTree.GetEntries())
        print("The length of the arrays are ",myTree.GetEntries())#len(pztarg))

        ##Read in position and momentum of particles
        for i in range(myTree.GetEntries()): #put data in arrays
            myTree.GetEntry(i)
            pztarg[i] = myTree.pz
            if pztarg[i] == 0.0:
                print("warning: PZtarg[{}]==0".format(i))
                continue #11.5.22 recommended by Kyrre
            xtarg[i] = myTree.x *mm #[m] ?
            pxtarg[i] = myTree.px / pztarg[i] #from Kyrre 5.5.22 to make it true X'!
            ytarg[i] = myTree.y *mm #[m] ?
            pytarg[i] = myTree.py / pztarg[i]
            Etarg[i] = myTree.E #[MeV]
            PDGtarg[i] = myTree.PDG
        myFile.Close()

        ##Filter the relevant distributions
        protarg_filter = np.equal(PDGtarg,2212) #PDG filter for protons
        eletarg_filter = np.equal(PDGtarg,11) #PDG filter for electrons
        neutarg_filter = np.equal(PDGtarg,2112) #PDG filter for neutrons
        Etarg_filter = np.greater(Etarg[protarg_filter],args.energy*simSetup_simple1["CUTOFF_ENERGYFRACTION"]) #Engcut is % number not decimal
        ##These proton only arrays are returned to the original script!
        Etarg_filtered_p = Etarg[protarg_filter]
        xtarg_filtered_p = xtarg[protarg_filter]
        #pxtarg_filtered_p = pxtarg[protarg_filter]
        ytarg_filtered_p = ytarg[protarg_filter]
        #pytarg_filtered_p = pytarg[protarg_filter]
        protons = Etarg[protarg_filter]
        electrons = Etarg[eletarg_filter]
        neutrons = Etarg[neutarg_filter]

        ##For plotting Energy distribution of all species
        if options['engPlot']:
            from math import floor,log10,ceil
            if args.sim == "raster":
                mag = floor(log10(args.Nb*args.nP)) #for dynamic title name
            elif args.sim == "beamlet":
                mag = floor(log10(args.Nbeamlet))
            engname=savename+"_EnergyPlot" #same for each plot

            if mat == "Vac":
                print("Vacuum Energy Plot not working now, plots empty histogram. *Shrugs*")
            else: #simple histogram plot
                import matplotlib.pyplot as plt
                #Protons Only
                plt.close()
                plt.hist(Etarg_filtered_p,100,log=True)
                plt.xlabel("Energy [MeV]")
                plt.ylabel("Counts")
                plt.xlim([0,args.energy])
                plt.title(rf"Energy Distribution at ESS Target of 10$^{{:d}}$ Protons".format(mag)+
                "\nThrough PBW of "+mat+", Protons Only")
                print("Protons Only",engname,args.picFormat)
                if args.savePics:
                    plt.savefig(engname+"."+args.picFormat)
                plt.close()

                ##Various Species Energy Plot
                #print(len(protons),len(electrons),len(neutrons))
                plt.hist(protons,100,log=True,color='b',label="Protons")
                plt.hist(electrons,25,log=True,color='r',label="Electrons") #bin 25 or else really thin
                plt.hist(neutrons,100,log=True,color='g',label="Neutrons")
                plt.xlabel("Energy [MeV]")
                plt.ylabel("Counts")
                #plt.xlim([0,args.energy])
                plt.legend()
                titl = "Energy Distribution at ESS Target Through PBW of "+mat+"\n For Various Species"
                plt.title(titl)
                if args.savePics:
                    plt.savefig(engname+"_Various."+args.picFormat)
                plt.close()                

        ##Display Full Energy distribution results
        #print("Full Energy distribution of {:d} particles with minimum Energy {:.3f}MeV through ".format(len(Eexit),np.min(Eexit_filtered)),mat," PBW")

        angMax=6e-3 #[rad] one angle filter limit to make the fits to the core of the beam
        posMax=500
        ##Apply PDG, Energy Filters
        xtarg_filtered = xtarg[protarg_filter][Etarg_filter]
        pxtarg_filtered = pxtarg[protarg_filter][Etarg_filter]
        ytarg_filtered = ytarg[protarg_filter][Etarg_filter]
        pytarg_filtered = pytarg[protarg_filter][Etarg_filter]
        ##Apply >,< filters
        #X, <
        pxfilterLa = np.less(pxtarg_filtered,angMax) #[rad]
        pxtarg_filteredA = pxtarg_filtered[pxfilterLa]
        xtarg_filteredA = xtarg_filtered[pxfilterLa]

        pxfilterLp = np.less(xtarg_filtered,posMax) #[mm]
        pxtarg_filteredP = pxtarg_filtered[pxfilterLp]
        xtarg_filteredP = xtarg_filtered[pxfilterLp]
        #Y, <
        pyfilterLa = np.less(pytarg_filtered,angMax) #[rad]
        pytarg_filteredA = pytarg_filtered[pyfilterLa]
        ytarg_filteredA = ytarg_filtered[pyfilterLa]

        pyfilterLp = np.less(ytarg_filtered,posMax) #[mm]
        pytarg_filteredP = pytarg_filtered[pyfilterLp]
        ytarg_filteredP = ytarg_filtered[pyfilterLp]

        #X, >
        pxfilterGa = np.greater(pxtarg_filteredA,-angMax) #[rad]
        pxtarg_filteredA = pxtarg_filteredA[pxfilterGa]
        xtarg_filteredA = xtarg_filteredA[pxfilterGa]

        pxfilterGp = np.greater(xtarg_filteredP,-posMax)
        pxtarg_filteredP = pxtarg_filteredP[pxfilterGp]
        xtarg_filteredP = xtarg_filteredP[pxfilterGp]
        #Y, >
        pyfilterGa = np.greater(pytarg_filteredA,-angMax) #[rad]
        pytarg_filteredA = pytarg_filteredA[pyfilterGa]
        ytarg_filteredA = ytarg_filteredA[pyfilterGa]

        pyfilterGp = np.greater(ytarg_filteredP,-posMax)
        pytarg_filteredP = pytarg_filteredP[pyfilterGp]
        ytarg_filteredP = ytarg_filteredP[pyfilterGp]

        ##Get Twiss for the filtered distributions #why not in /mm?
        targxTwissf = calcTwiss("Target X Filtered","Target X' Filtered",xtarg_filteredA,pxtarg_filteredA) #angle filtered to fit core 4.3.23
        targyTwissf = calcTwiss("Target Y Filtered","Target Y' Filtered",ytarg_filteredA,pytarg_filteredA)

    ##Analytical Formula Calculations
    if options['MCS']:
        ##If magnet, use multiple scattering layers instead of averaging!
        from plotFit import plotTwissFit,calcEq8,calcEq16

        ##Particle characteristic values
        if args.particle == "proton":
            partA = 938.27209 #[MeV/c2]
            partZ = 1
        elif args.particle == "electron":
            partA = 0.511 #[MeV/c2]
            partZ = 1
        gamma_rel = (args.energy+partA)/partA #from PrintTwissParameters
        beta_rel = np.sqrt(gamma_rel*gamma_rel -1 )/gamma_rel
        ##Get rid of Normalized Emittance!
        Igemtx = Twiss[0]/(beta_rel*gamma_rel)
        Igemty = Twiss[3]/(beta_rel*gamma_rel)

        ##Use list of Twiss values for simple passing of data: [beta,alpha,gemt,gamma]
        TwissIx = [Twiss[1],Twiss[2],Igemtx*um,((1+Twiss[2]*Twiss[2])/Twiss[1])] #Initial Twiss
        TwissIy = [Twiss[4],Twiss[5],Igemty*um,((1+Twiss[5]*Twiss[5])/Twiss[4])]
        PBWTwx = [Twiss[1],Twiss[2],Twiss[0]] #Twiss for each for printing on graphs
        PBWTwy = [Twiss[4],Twiss[5],Twiss[3]]


        print("ORIGINAL:", TwissIx)

        ##Highland Equation Radiation Length Calculation
        p = np.sqrt((args.energy+partA)**2 - (partA)**2) #[MeV/c] #derived with Kyrre 15.6.22
        #pAppleby = beta_rel * gamma_rel * partA / c #Appleby2022 Eqn 2.16
        betap = beta_rel*p #Eq 5

        ##For PBW magnet case:
        if simSetup_simple1["THICK"] == 0.0:
            from plotFit import Drift
            m1Len = simSetup_simple1["MAGNET"][0]["keyval"]["al1Thick"] #[mm]?
            print("m1Len =", m1Len)
            if beamXAngle != 0: #account for angle contribution to thickness
                m1Len = m1Len / np.cos(beamXAngle) 
            ##Al Front contribution
            thetasqAl1 = 13.6 * partZ / betap * np.sqrt(m1Len/radLen) * (1 + 0.038 * np.log(m1Len/radLen))
            #             calcEq8(thetasq,      Twiss,thick,beta_rel,gamma_rel
            Twisse8xAl1  = Drift(calcEq8(thetasqAl1,  TwissIx,m1Len,beta_rel,gamma_rel),"driftAl1",m1Len)
            Twisse8yAl1  = Drift(calcEq8(thetasqAl1,  TwissIy,m1Len,beta_rel,gamma_rel),"driftAl1",m1Len)
            Twisse16xAl1 = Drift(calcEq16(thetasqAl1, TwissIx,m1Len,beta_rel,gamma_rel),"driftAl1",m1Len)
            Twisse16yAl1 = Drift(calcEq16(thetasqAl1, TwissIy,m1Len,beta_rel,gamma_rel),"driftAl1",m1Len)
            #MUST ADD DRIFTS HERE, as per Kyrre 3.5.23
            ##H2O contribution
            m2Len = simSetup_simple1["MAGNET"][0]["keyval"]["waterThick"]
            if beamXAngle != 0: #account for angle contribution to thickness
                m2Len = m2Len / np.cos(beamXAngle) 
            thetasqH2O = 13.6 * partZ / betap * np.sqrt(m2Len/radLenH2O) * (1 + 0.038 * np.log(m2Len/radLenH2O))
            #             calcEq8(thetasq,          Twiss,thick,beta_rel,gamma_rel
            Twisse8xH2O  = Drift(calcEq8(thetasqH2O,  Twisse8xAl1, m2Len,beta_rel,gamma_rel),"drifth20",m2Len)
            Twisse8yH2O  = Drift(calcEq8(thetasqH2O,  Twisse8yAl1, m2Len,beta_rel,gamma_rel),"drifth20",m2Len)
            Twisse16xH2O = Drift(calcEq16(thetasqH2O, Twisse16xAl1,m2Len,beta_rel,gamma_rel),"drifth20",m2Len)
            Twisse16yH2O = Drift(calcEq16(thetasqH2O, Twisse16yAl1,m2Len,beta_rel,gamma_rel),"drifth20",m2Len)

            ##Al Back contribution
            m3Len = simSetup_simple1["MAGNET"][0]["keyval"]["al2Thick"]
            if beamXAngle != 0: #account for angle contribution to thickness
                m3Len = m3Len / np.cos(beamXAngle) 
            thetasqAl2 = 13.6 * partZ / betap * np.sqrt(m3Len/radLen) * (1 + 0.038 * np.log(m3Len/radLen))
            #          calcEq8(thetasq,          Twiss,   thick,beta_rel,gamma_rel
            Twisse8x  = Drift(calcEq8(thetasqAl2,  Twisse8xH2O, m3Len,beta_rel,gamma_rel),"driftAl2",m3Len)
            Twisse8y  = Drift(calcEq8(thetasqAl2,  Twisse8yH2O, m3Len,beta_rel,gamma_rel),"driftAl2",m3Len)
            Twisse16x = Drift(calcEq16(thetasqAl2, Twisse16xH2O,m3Len,beta_rel,gamma_rel),"driftAl2",m3Len)
            Twisse16y = Drift(calcEq16(thetasqAl2, Twisse16yH2O,m3Len,beta_rel,gamma_rel),"driftAl2",m3Len)

            totalThickness = m1Len+m2Len+m3Len
            print(totalThickness,"mm SCATTER:", Twisse8x, Twisse16x)

        ##Analytical Calculations for simple flat sheet
        elif simSetup_simple1["THICK"] != 0.0: #if only one layer (MiniScatter "target")
            #need to check if this starts at 0 or is centered on 0!!
            from plotFit import Drift
            ##Highland Equation Radiation Length Calculation for single layer
            if beamXAngle != 0: #account for angle contribution to thickness
                args.t = args.t / np.cos(beamXAngle) 
            thetasq = 13.6 * partZ / betap * np.sqrt(args.t/radLen) * (1 + 0.038 * np.log(args.t/radLen)) #from Eq 5
            #print("\nradLen: {:.2f}, p: {:.3e}, gamma: {:.3f}, beta: {:.3f}, theta^2: {:.3e} radians".format(radLen,p,gamma_rel,beta_rel,thetasq))
            print("THICKNESS:",args.t)
            Twisse8x  = calcEq8(thetasq, TwissIx,args.t,beta_rel,gamma_rel) #calculated target exit Twiss
            Twisse8y  = calcEq8(thetasq, TwissIy,args.t,beta_rel,gamma_rel)
            Twisse16x = calcEq16(thetasq,TwissIx,args.t,beta_rel,gamma_rel) #calculated 2 target exit Twiss
            Twisse16y = calcEq16(thetasq,TwissIy,args.t,beta_rel,gamma_rel)
            totalThickness = args.t
            Twisse8x = Drift(Twisse8x,"driftSlab",args.t)
            Twisse8y = Drift(Twisse8y,"driftSlab",args.t)
            Twisse16x = Drift(Twisse16x,"driftSlab",args.t)
            Twisse16y = Drift(Twisse16y,"driftSlab",args.t)
            print(totalThickness,"mm SCATTER:", Twisse8x, Twisse16x)

        ##Plotting PBW Exit distribution vs the PDF produced from the Formalism Equations
        #Displays Twiss values and calculated mu and sigma
        if options['TwissFits']:
            print("\nPre PBW Twiss")
            plotTwissFit(xinit/mm,pxinit,savename+"init",mat,"Pre PBW","X",args.t,thetasq,beta_rel,gamma_rel,TwissIx)
            plotTwissFit(yinit/mm,pyinit,savename+"init",mat,"Pre PBW","Y",args.t,thetasq,beta_rel,gamma_rel,TwissIy)
            print("\nPBW Exit Twiss Calculated")
            plotTwissFit(xexit/mm,pxexit,savename+"texitHalo",mat,"PBW Exit","X",args.t,thetasq,beta_rel,gamma_rel,TwissIx)
            plotTwissFit(xexit_filtered/mm,pxexit_filtered,savename+"texitFiltered",mat,"PBW Exit","X",args.t,thetasq,beta_rel,gamma_rel,TwissIx)
            plotTwissFit(yexit_filtered/mm,pyexit_filtered,savename+"texitFiltered",mat,"PBW Exit","Y",args.t,thetasq,beta_rel,gamma_rel,TwissIy)

        ##Extension to Target. Needed for compareTargets!
        from plotFit import Drift,compareTargets,getMoments
        if args.t==0:
            z0Diff = simSetup_simple1["MAGNET"][0]["pos"] - 24.125
        else:z0Diff=0
        driftLength = 3565 - totalThickness #[mm]
        e8TargxReal = Drift(Twisse8x,"e8XReal",driftLength)
        e8TargyReal = Drift(Twisse8y,"e8YReal",driftLength)
        e16TargxReal = Drift(Twisse16x,"e8XReal",driftLength)
        e16TargyReal = Drift(Twisse16y,"e8YReal",driftLength)

        print(z0Diff, totalThickness,driftLength,"DRIFT:", e8TargxReal, e16TargxReal)

        ##Now compare the MiniScatter Target distribution (targxTwissf) to initTarg, exitTarg, e8Targ and e16Targ PDFs
        if args.compTargs:
            print("args.compTargs")
            e2t = "Extended to Target"
            ##This is to look at exit distribution, select which to use:
            #compareTargets(xexit/mm,yexit/mm,exitxTwf,exityTwf,TwissIx,TwissIy,"PBW Exit",savename+"PBWExit",mat,PBWTwx,PBWTwy,args)
            #compareTargets(targx distrib,      targy distrib,   targTwx,   targTwy,    fitTwx,   fitTwy,   fitlabel,               savename,   mat,PBWTwx,PBWTwy,args):
            #compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,initTargx,initTargy,"Initial Twiss "+e2t,savename+"NoPBW",mat,PBWTwx,PBWTwy,args)
                    #The above comparison btwn initial Twiss and Target distrib shows that the drift propagation of the initial twiss is != what the beamlet will actually be. 
                            #The actual Y will be 2.3x larger than a Twiss drift and actual X will be 1.25x larger than a Twiss drift.
            #compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,exitTargx,exitTargy,"PBW Exit Twiss "+e2t,savename,mat,PBWTwx,PBWTwy,args)
                    #The above comparison btwn PBW Exit Twiss and Target distrib shows that the beam distrib transforms during the drift from PBW Exit.
                            #That it isn't a mere drift, but the MCS causes a 14% larger beam in Y and 10% larger beam in X.
            targSigX,targSigY = compareTargets(xtarg_filteredP/mm,ytarg_filteredP/mm,targxTwissf,targyTwissf,e8TargxReal,e8TargyReal,"Müller Eqn 8",savename,mat,PBWTwx,PBWTwy,args)
                    #Position filtered to show spread
                    #The above comparison btwn the Twiss at Target with drift resulting from Muellers formalism and actual Target distrib shows that they are in strong agreement.
                            #The Mueller formalism equations with correct radiation lengths, etc produce beam sigmas within 1% in X and Y (-0.92% in X and +0.38% in Y)
            #compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,e16TargxReal,e16TargyReal,"Eq 16 Calculation from Initial Twiss "+e2t,savename,mat,PBWTwx,PBWTwy,args)
                    #This shows that the eq 16 calculations don't work xD
            #print(args.t)
            #print("make plots",datetime.now().strftime("%H-%M-%S"))
            targSigsX = [targSigX,0]
            targSigsY = [targSigY,0]
            e8SigsX = getMoments(e8TargxReal)
            e8SigsY = getMoments(e8TargyReal)

            ##Target comparison plots
            if args.t == 0.1:
                print("Vacuum")
                compareTargets(xtarg_filtered_p/mm,ytarg_filtered_p/mm,targxTwissf,targyTwissf,e8TargxReal,e8TargyReal,"No PBW, Eq 8",savename+"HaloPDGFiltered_Eq8",mat,PBWTwx,PBWTwy,args)
            elif args.t == 4.25:
                print("PBW!")
                compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,initTargx,initTargy,"No PBW",savename+"NoPBW",mat,PBWTwx,PBWTwy,args)
                compareTargets(xtarg_filtered_p/mm,ytarg_filtered_p/mm,targxTwissf,targyTwissf,e8TargxReal,e8TargyReal,"Real PBW, Eq 8",savename+"HaloPDGFiltered_Eq8",mat,PBWTwx,PBWTwy,args)
                #compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,e8TargxReal,e8TargyReal,"Real PBW, Eq 8",savename+"Eq8",mat,PBWTwx,PBWTwy,args)
                #compareTargets(xtarg_filtered/mm,ytarg_filtered/mm,targxTwissf,targyTwissf,e16TargxReal,e16TargyReal,"Eq 16",savename,mat,PBWTwx,PBWTwy,args)

    ##For Raster Images
    if beamFile != "":
        from plotFit import plot1DRaster,rasterImage
        #plot1DRaster(xtarg_filtered_p/mm,ytarg_filtered_p/mm,"Traster",savename,mat,"Target")
        (twiss_PBW, numPart_PBW, objects_PBW) = miniScatterDriver.getData_tryLoad(simSetup_simple1, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212"])
        PMASreturn = rasterImage(savename,"Target",objects_PBW["tracker_cutoff_xy_PDG2212"],simSetup_simple1["N"],args,Twiss,options,boxes,paths,0)
    
        #if options['initTree']:
        #    (twiss_PBW, numPart_PBW, objects_PBW) = miniScatterDriver.getData_tryLoad(simSetup_simple1, tryload=TRYLOAD,getObjects=["init_xy"])
        #    #plot1DRaster(xinit/mm,yinit/mm,"Iraster",savename,mat,"PBW")
        #    initPOutBox = rasterImage(savename,"PBW",objects_PBW["init_xy"],simSetup_simple1["N"],args,Twiss,options,boxes,paths)
    else:
        ##Place holder variables
        jMax = 0
        pOutsideBoxes = 0
        beamArea = 0
        coreJMean = 0
        centX = 0
        centY = 0
        rValue = 0
        rDiff = 0
        PMASreturn = [jMax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue,rDiff]

        ##For fitting scattered beamlets
        if args.gaussFit:
            from plotFit import converter,gaussianFit,fitGaussians
            print("else GaussFit")
            (twiss_PBW, numPart_PBW, objects_PBW) = miniScatterDriver.getData_tryLoad(simSetup_simple1, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212"])
            (Img, xax, yax) = converter(objects_PBW["tracker_cutoff_xy_PDG2212"],args.saveHist,savename,paths,args,False) #convert from TH2D to numpy map
            maxim = 500
            xBinSize = xax[maxim+1]-xax[maxim]
            yBinSize = yax[maxim+1]-yax[maxim]
            maxim=500
            yBinSize = 500
            xBinSize = 500
            Lorentz = True
            Error = np.ones((2,2))
            #fitGaussians(Img,Error,paths["statsPWD"],maxim,"y",5)
            diffNy,diffPy,coeffsy, differenceNLy,differencePLy,coeffsLy = gaussianFit(objects_PBW["tracker_cutoff_xy_PDG2212"],"y",yBinSize,maxim,savename,2,25,args.saveFits,True,args.gauss2Fit,Lorentz)
            diffNx,diffPx,coeffsx, differenceNLx,differencePLx,coeffsLx = gaussianFit(objects_PBW["tracker_cutoff_xy_PDG2212"],"x",xBinSize,maxim,savename,3,10,args.saveFits,True,args.gauss2Fit,Lorentz)
            #print(diffNx,diffPx,diffNy,diffPy,"\n",coeffsx,"\n",coeffsy)
            if args.savePics:
                from plotFit import rasterImage
                PMASreturn = rasterImage(savename,"Target",objects_PBW["tracker_cutoff_xy_PDG2212"],simSetup_simple1["N"],args,Twiss,options,boxes,paths,0)

    ##For thickness variation simulations, return Analytical Mueller calculated beam sizes for thicknesses
    if args.sim == "thick":
        from plotFit import getMoments,compareTargets,findFit
        e8SigsX = getMoments(e8TargxReal)
        e8SigsY = getMoments(e8TargyReal)
        targSigsX = getMoments(targxTwissf)
        targSigsY = getMoments(targyTwissf)
        print("getting targx",np.shape(xtarg_filteredP))
        #_,targSigX,_,_ = findFit(xtarg_filteredP/mm,[0.01,0.1,10],(0,[1e8,1,20]),"auto")
        #targSigX,targSigY = compareTargets(xtarg_filteredP/mm,ytarg_filteredP/mm,targxTwissf,targyTwissf,e8TargxReal,e8TargyReal,"Müller Eqn 8",savename,mat,PBWTwx,PBWTwy,args)
        targSigsX = [targSigX,0]
        print("getting targy")
        #_,targSigY,_,_ = findFit(ytarg_filteredP/mm,[0.01,0.1,10],(0,[1e8,1,20]),"auto")
        targSigsY = [targSigY,0]
        print(e8SigX,e8SigY,targSigX,targSigY)
    else:
        e8SigsX = [0,0];e8SigsY = [0,0]; targSigsX = [0,0]; targSigsY = [0,0]

    ##Passing lists to decrease length of argument calls
    return savename, xtarg_filtered_p/mm, ytarg_filtered_p/mm, PMASreturn,e8SigsX,e8SigsY,targSigsX,targSigsY #filter by PDG only