#Script for single plots of 
# optimization value: % outside 99% box
#typical command: (no PBW, no text, no boxes, save pic of 200x200mm^2 window)
#python3 scatterPBW.py --savePics --xlim 200 --ylim 200 --t 0.1 --noText --noBox

def scatterPBW(args,Twiss,iteration):
    from datetime import datetime
    origin = datetime.now()
    print(origin)
    import numpy as np
    from os import uname
    from runARasterMaker import runARasterMaker
    from simulation import simulation

    #Constants for running scripts
    physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ" or "QGSP_BERT__SS"
    zoff = "*-10" #[mm] with preappended * to keep covar defined at z=0
    options     = {'physList':physList, 'dependence':"Twiss", 'zoff':zoff, 'initTree':False,
                    'exitTree':False, 'targetTree':False, 'MCS':False, 'engPlot':False,
                    'mat3Plot':False }

    #Important things
    if args.t == 0:
        args.material = "Al" #overwrites potential user input. Needs work
    elif args.t == 0.1:
        args.material = "Vac"
    boxes = [0]#,-.25,-.375,-.45,-.50]#,0.125,0.25,0.375] #make an args for 24.11.22

    if uname()[1] == "tensor.uio.no":
        csvPWD = "/scratch2/ericdf/PBWScatter/CSVs/"
        statsPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    elif uname()[1] == "mbef-xps-13-9300":
        csvPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/"
        statsPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    else:
        csvPWD = input("Path from home to direction you like to save root files to: ")
        statsPWD = "."

    #Create Rastered Beam file, runARasterMaker checks if the CSV is already present
    rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,csvPWD,options,iteration)
    ##print(rasterBeamFile,beamXAngle,beamYAngle)
    #Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
    #                                       simulation(args,material,   engplot,beamXAngle,beamYAngle,rasterBeamFile,Twiss,options,boxes):
    #_,_,_,Jmax,pOutsideBoxes,dispY,dispX,rValue = simulation(args,args.material,False,beamXAngle,beamYAngle,rasterBeamFile,Twiss,options,boxes)
    #Returns filtered particle histograms when doing beamlet simulation, not needed here!
    
    from simulation import setup
    #def                        setup(args,material,           beamFile,Twiss,options):
    savename,simSetup_simple1 = setup(args,args.material,rasterBeamFile,Twiss,options)

    import miniScatterDriver
    from plotFit import rasterImage
    TRYLOAD = True  #Try to load already existing data instead of recomputing, only if using getData_TryLoad function.
    (twiss_PBW, numPart_PBW, objects_PBW) = miniScatterDriver.getData_tryLoad(simSetup_simple1, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212","init_xy"])
    Jmax,pOutsideBoxes,dispY,dispX,rValue,rDiff = rasterImage(savename,"Target",objects_PBW["tracker_cutoff_xy_PDG2212"],simSetup_simple1["N"],args,Twiss,options,boxes)

    from plotFit import saveStats
    saveStats(statsPWD,Twiss,rasterBeamFile,Jmax,pOutsideBoxes,dispY,dispX,rValue,rDiff)

    print("Simulation took ",datetime.now()-origin,"s long",sep="")