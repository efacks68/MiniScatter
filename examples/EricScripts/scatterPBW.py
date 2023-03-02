#Script for single plots of 
# optimization value: % outside 99% box
#typical command: (no PBW, no text, no boxes, save pic of 200x200mm^2 window)
#python3 scatterPBW.py --savePics --xlim 200 --ylim 200 --t 0.1 --noText --noBox

def scatterPBW(args,Twiss,iteration,paths):
    from datetime import datetime
    origin = datetime.now()
    print(origin)
    from runARasterMaker import runARasterMaker
    from simulation import setup
    from plotFit import rasterImage

    #Constants for running scripts
    physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ" or "QGSP_BERT__SS"
    zoff = "*-1" #[mm] with preappended * to keep covar defined at z=0
    options     = {'physList':physList, 'dependence':"Twiss", 'zoff':zoff, 'initTree':False,
                    'exitTree':False, 'targetTree':False, 'MCS':False, 'engPlot':False,
                    'mat3Plot':False }

    #Important things
    if args.t == 0:
        args.material = "Al" #overwrites potential user input. Needs work
    elif args.t == 0.1:
        args.material = "Vac"
    boxes = [0]#,-.25,-.375,-.45,-.50]#,0.125,0.25,0.375] #make an args for 24.11.22

    #Create Rastered Beam file, runARasterMaker checks if the CSV is already present
    rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,paths['csvPWD'],options,iteration)

    #def                        setup(args,material,           beamFile,Twiss,options,paths):
    savename,simSetup_simple1 = setup(args,args.material,rasterBeamFile,Twiss,options,paths)

    import miniScatterDriver
    TRYLOAD = True  #Try to load already existing data instead of recomputing, only if using getData_TryLoad function.
    (twiss_PBW, numPart_PBW, objects_PBW) = miniScatterDriver.getData_tryLoad(simSetup_simple1, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212","init_xy"])
    Jmax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue,rDiff = rasterImage(savename,"Target",objects_PBW["tracker_cutoff_xy_PDG2212"],simSetup_simple1["N"],args,Twiss,options,boxes,paths)

    from plotFit import saveStats
    saveStats(paths['statsPWD'],Twiss,rasterBeamFile,Jmax,pOutsideBoxes,beamArea,coreJMean,centX,centY,rValue,rDiff)

    print("Simulation took ",datetime.now()-origin,"s long",sep="")