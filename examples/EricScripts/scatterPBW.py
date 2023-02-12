#Script for single plots of 
# optimization value: % outside 99% box
#typical command: (no PBW, no text, no boxes, save pic of 200x200mm^2 window)
#python3 scatterPBW.py --savePics --text --xlim 200 --ylim 200 --box --t 0.1 



#Command Line arguments for save control
#parser = ArgumentParser()
#parser.add_argument("--beamClass", type=str,    default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
#parser.add_argument("--l",         type=str,    help="Load Particles or not",   default="PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03")
#parser.add_argument("--twiss",     type=float,  nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
#parser.add_argument("--t",         type=float,  default=0,     help="PBW Thickness [mm], 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al Xmm thick. Default=0")
#parser.add_argument("--energy",    type=float,  default=570,   help="Beam Energy [MeV]. Default=570")
#parser.add_argument("--Nb",        type=int,    default=10,    help="Number of macroparticles per beamlet. Default=10")
#parser.add_argument("--nP",        type=float,  default=1e3,   help="Numper of beamlets in pulse. Default=1e3")
#parser.add_argument("--rX",        type=float,  default=0,     help="X distance from beam axis [mm]. Default=0")
#parser.add_argument("--rY",        type=float,  default=0,     help="Y distance from beam axis [mm]. Default=0")
#parser.add_argument("--aX",        type=float,  default=54.65, help="RM X Amplitude [mm]. Default=54.65")
#parser.add_argument("--aY",        type=float,  default=18.37, help="RM Y Amplitude [mm]. Default=18.37")
#parser.add_argument("--failure",   type=float,  default=0,     choices = range(0,5),  help="Which RM Failure case, 0-4. Default=0")
#parser.add_argument("--magFails",  type=int,    default=2,     choices = range(0,5),  help="Number of Raster Magnets that fail, 1-4. Default=2")
#parser.add_argument("--xlim",      type=float,  default=150,   help="+/- value for horizontal axis of output rastered image [mm]. Default=150")
#parser.add_argument("--ylim",      type=float,  default=150,   help="+/- value for vertical axis of output rastered image [mm]. Default=150")
#parser.add_argument("--maxim",     type=float,  default=0  ,   help="Maximum current density value for output rastered imagem[uA/cm^2]. Default=0")
#parser.add_argument("--edgeRaster",action="store_true",  help="Only populate edges of raster. Default=False")
#parser.add_argument("--PBIP",      action="store_true",  default=False,   help="Is PBIP present? Default=False")
#parser.add_argument("--noText",    action="store_true",  default=False,   help="Turns off printed text when called. Default=False")
#parser.add_argument("--noBox",     action="store_true",  default=False,   help="Turns off printed box when called. Default=False")
#parser.add_argument("--savePics",  action="store_true",  default=False,   help="Saves Rastered Image. Default=False")
#parser.add_argument("--saveGrads", action="store_true",  default=False,   help="Plots gradients of beam at Target. Default=False")
#parser.add_argument("--saveEdges", action="store_true",  default=False,   help="Plots edges on Raster Image. Default=False")
#parser.add_argument("--gaussFit",  action="store_true",  default=False,   help="Computes sum of Gaussian fits for central axis projection. Default=False")
#parser.add_argument("--saveFits",  action="store_true",  default=False,   help="Saves plots of Gaussian Fitting. Default=False")
#parser.add_argument("--saveHist",  action="store_true",  default=False,   help="Saves Histogram of proton density at target. Default=False")
#parser.add_argument("--saveRaster",action="store_true",  default=False,   help="Saves plot of rastered beam. Default=False")
#args = parser.parse_args()

#to do: split the arg parse into the script to run which you can call scatter or map or others from?
def scatterPBW(args,Twiss):
    from datetime import datetime
    origin = datetime.now()
    print(origin)
    import numpy as np
    #import matplotlib.pyplot as plt
    from os import uname
    #from sys import path as sysPath
    #from argparse import ArgumentParser
    from runARasterMaker import runARasterMaker
    from simulation import simulation

    #Constants for running scripts
    physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ" or "QGSP_BERT__SS"
    zoff = "*-10" #[mm] with preappended * to keep covar defined at z=0

    options     = {'physList':physList, 'dependence':"Twiss", 'zoff':zoff }
                        #'xlim':args.xlim, 'ylim':args.ylim, 'maxim':args.maxim, 'saveHist':args.saveHist,
                        #'PBIP':args.PBIP, 'beamClass':args.beamClass, 'Nb':args.Nb, 'failure':args.failure,
                        #'magFails':args.magFails, 'saveRaster':args.saveRaster, 'saveFits':args.saveFits,
                        #'saveGrads':args.saveGrads, 'saveEdges':args.saveEdges, 'gaussFit':args.gaussFit,
                        #'material':args.material,

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
    rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,csvPWD,options)
    ##print(rasterBeamFile,beamXAngle,beamYAngle)
    #Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
    #                                       simulation(args,material,   engplot,beamXAngle,beamYAngle,rasterBeamFile,Twiss,options,boxes):
    #_,_,_,Jmax,pOutsideBoxes,dispY,dispX,rValue = simulation(args,args.material,False,beamXAngle,beamYAngle,rasterBeamFile,Twiss,options,boxes)
    #Returns filtered particle histograms when doing beamlet simulation, not needed here!
    
    from simulation import setup
    #def                                  setup(args,material,     beamXAngle,beamYAngle,      beamFile,Twiss,options):
    savename,simSetup_simple1 = setup(args,args.material,beamXAngle,beamYAngle,rasterBeamFile,Twiss,options)
    from sys import path as sysPath
    MiniScatter_path="../../MiniScatter/build/."
    sysPath.append(MiniScatter_path)
    import miniScatterDriver
    TRYLOAD = True  #Try to load already existing data instead of recomputing, only if using getData_TryLoad function.
    
    from plotFit import rasterImage
    (twiss_PBW, numPart_PBW, objects_PBW) = miniScatterDriver.getData_tryLoad(simSetup_simple1, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212","init_xy"])
    Jmax,pOutsideBoxes,dispY,dispX,rValue = rasterImage(savename,"Target",objects_PBW["tracker_cutoff_xy_PDG2212"],simSetup_simple1["N"],args,Twiss,options,boxes)

    from plotFit import saveStats
    saveStats(statsPWD,rasterBeamFile,Jmax,pOutsideBoxes,dispY,dispX,rValue)

    print("Simulation took ",datetime.now()-origin,"s long",sep="")