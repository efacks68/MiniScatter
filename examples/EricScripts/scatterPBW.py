##scatterPBW.py
##Script for directing the rastered beam simulations

##Comes from rasterScatter
def scatterPBW(args,Twiss,sample,paths,origBx,origBY):
    from datetime import datetime
    origin = datetime.now()
    from runARasterMaker import runARasterMaker

    ##Constants for running MiniScatter
    zoff = "*-1" #[mm] with preappended * to keep covar defined at z=0
    options     = {'physList':args.physList, 'dependence':"Twiss", 'zoff':zoff, 'initTree':False,
                    'exitTree':False, 'targetTree':False, 'MCS':False, 'engPlot':False,
                    'mat3Plot':False,'MiniRoot':True }

    ##Custom settings, if t==0, want standard PBW. If t=0.1, want no PBW, make it Vacuum
    if args.t == 0:
        args.material = "Al" #overwrites potential user input. Needs work
    elif args.t == 0.1:
        args.material = "Vac"
    boxes = [0]#,-.25,-.375,-.45,-.50]#,0.125,0.25,0.375] #make an args for 24.11.22

    ##Get/Create Rastered Beam file. runARasterMaker checks if the CSV is already present
    rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,paths['csvPWD'],options,sample,origBx,origBY)

    ##Get run info from setup function
    from simulation import setup
    savename,simSetup_simple1 = setup(args,args.material,rasterBeamFile,Twiss,options,paths)

    ##Use MiniScatter to load ROOT file
    import miniScatterDriver
    TRYLOAD = True  #Try to load already existing data instead of recomputing, only if using getData_TryLoad function.
    (twiss_PBW, numPart_PBW, objects_PBW) = miniScatterDriver.getData_tryLoad(simSetup_simple1, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212","init_xy"])
    
    ##Make an image of the rastered beam and return beam info
    from plotFit import rasterImage
    [Jmax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue,rDiff] = rasterImage(savename,"Target",objects_PBW["tracker_cutoff_xy_PDG2212"],simSetup_simple1["N"],args,Twiss,options,boxes,paths,sample)

    ##Save the information in a statistics csv for quick loading
    from plotFit import saveStats
    saveStats(paths['statsPWD'],Twiss,rasterBeamFile,Jmax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue,rDiff,"",args.reBin,args)

    print("Single Simulation took ",datetime.now()-origin,"s long",sep="")