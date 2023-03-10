#rasterScatter.py
#Takes config for scattering simulations and runs
#Includes setup and output preferences

#possible commands:
#Vary N particles for nominal convergence study
    #python3 rasterScatter.py --Nb 100 --samples 200 --saveSpread
#Jitter study:
    #python3 rasterScatter.py --samples 200 --saveSpread --betaSpread 10
#Fit Gaussian and Voigt to a beamlet:
    #python3 rasterScatter.py --sim beamlet --gaussFit --saveFits --Nbeamlet 1e7
#Failure QP 139 (close to nominal, but slightly 1/4 smaller X) spread approx:
    #python3 rasterScatter.py --source twiss --twissFile FailureHEBT-A2T --qpNum 139 --betaSpread 10 --samples 5 --saveSpread
#Map RA Dependence:
    #python3 rasterScatter.py --sim map --ampl map --Nstep 7
#thickness dependence plot:
    #python3 rasterScatter.py --sim thick --stepThick 0.25 --maxThick 3
#raster magnet failures:
    #python3 rasterScatter.py --failure 1 --magFails 4 --savePics

#load from CSV file name:
    #python3 rasterScatter.py --source csv --beamFile sampleIn50Pct_OrigbX1006.80,bY129.72m_beta1006.80,129.72m_N2.9e+06_NpB100 --saveGrads --saveEdges


#Not sure it's useful:Various Twiss ranges with samples:
    #python3 rasterScatter.py --source twiss --twissFile TwissRange0,12mm-mrad_3Bx-3By-2Ax-2Ay --qpNum 0 --samples 10 --betaSpread 10


from datetime import datetime
origin = datetime.now()
from argparse import ArgumentParser
from os import uname
from plotFit import getTwiss
#look into argument groups!

#Command Line arguments for save control
parser = ArgumentParser()
parser.add_argument("--sim",       type=str,   default="raster",   choices=("raster","map","beamlet","thick"), help="Type of simulation to perform")
parser.add_argument("--source",    type=str,   default="particles",choices=("particles","twiss","csv"), help="From load particles or Twiss or CSV (put name in --beamFile)")
parser.add_argument("--twissFile", type=str,   default="",    help="Load file with Twiss, auto look in OpenXAL folder")
parser.add_argument("--beamClass", type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
parser.add_argument("--particle",  type=str,   default="proton",choices=("proton","electron"), help="Which particle to simulate?")
parser.add_argument("--beamFile",  type=str,   help="Load Particles or not", default="PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03")
parser.add_argument("--twiss",     type=float, nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
parser.add_argument("--qpNum",     type=str,   default="",    help="Either a number between 099 and 148, or all")
parser.add_argument("--betaSpread",type=float, default=0,     help="What % around provided Beta should we sample from")
parser.add_argument("--samples",   type=int,   default=1,     help="How many times to sample this setting")
parser.add_argument("--csvFile",   type=str,   default="",    help="Load Beam of already made csv")
#General Beam Setup Options
parser.add_argument("--t",         type=float, default=0,     help="PBW Thickness [mm], 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al Xmm thick. Default=0")
parser.add_argument("--energy",    type=float, default=570,   help="Beam Energy [MeV]. Default=570")
parser.add_argument("--Nb",        type=int,   default=100,    help="Number of macroparticles per beamlet. Default=10")
parser.add_argument("--nP",        type=float, default=1e3,   help="Numper of beamlets in pulse. Default=1e3")
parser.add_argument("--rX",        type=float, default=0,     help="X distance from beam axis at BPM94 [mm]. Default=0")
parser.add_argument("--rY",        type=float, default=0,     help="Y distance from beam axis at BPM94 [mm]. Default=0")
parser.add_argument("--aX",        type=float, default=54.65, help="RM X Amplitude [mm]. Default=54.65")
parser.add_argument("--aY",        type=float, default=18.37, help="RM Y Amplitude [mm]. Default=18.37")
parser.add_argument("--failure",   type=float, default=0,     choices = range(0,5),  help="Which RM Failure case, 0-4. Default=0")
parser.add_argument("--magFails",  type=int,   default=2,     choices = range(0,5),  help="Number of Raster Magnets that fail, 1-4. Default=2")
parser.add_argument("--xlim",      type=float, default=150,   help="+/- value for horizontal axis of output rastered image [mm]. Default=150")
parser.add_argument("--ylim",      type=float, default=100,   help="+/- value for vertical axis of output rastered image [mm]. Default=150")
parser.add_argument("--maxim",     type=float, default=0  ,   help="Maximum current density value for output rastered imagem[uA/cm^2]. Default=0")
parser.add_argument("--edgeRaster",action="store_true",  help="Only populate edges of raster. Default=False")
parser.add_argument("--PBIP",      action="store_true",  default=False,   help="Is PBIP present? Default=False")
parser.add_argument("--material",  type=str,   default="Al",  choices=("Al","Au","C","Vac"),  help="What material PBW?")
parser.add_argument("--Nbeamlet",  type=float, default=1e5,   help="For beamlet simulation, how many protons?")
#Output Options                         Lowering the rCut increases the jMax...
parser.add_argument("--rCut",      type=float, default=1e3,  help="Radial cut, defines worldSize and histLims")
parser.add_argument("--engCut",    type=float, default=0.9,  help="Energy cut, see MiniScatter description")
parser.add_argument("--noText",    action="store_true",  default=False,   help="Turns off printed text when called. Default=False")
parser.add_argument("--noBox",     action="store_true",  default=False,   help="Turns off printed box when called. Default=False")
parser.add_argument("--savePics",  action="store_true",  default=False,   help="Saves Rastered Image. Default=False")
parser.add_argument("--saveGrads", action="store_true",  default=False,   help="Plots gradients of beam at Target. Default=False")
parser.add_argument("--saveEdges", action="store_true",  default=False,   help="Plots edges on Raster Image. Default=False")
parser.add_argument("--gaussFit",  action="store_true",  default=False,   help="Computes sum of Gaussian fits for central axis projection. Default=False")
parser.add_argument("--saveFits",  action="store_true",  default=False,   help="Saves plots of Gaussian Fitting. Default=False")
parser.add_argument("--saveHist",  action="store_true",  default=False,   help="Saves Histogram of proton density at target. Default=False")
parser.add_argument("--saveRaster",action="store_true",  default=False,   help="Saves plot of rastered beam. Default=False")
parser.add_argument("--picFormat", choices=("png","svg","pdf"), type=str, default="png",  help="Whic file format extension?")
parser.add_argument("--matPlots",  action="store_true",  default=False,   help="Whether to do various material plots for beamlets")
parser.add_argument("--saveSpread",action="store_true",  default=False,   help="Saves PMAS parameter spread histograms. Default=False")
parser.add_argument("--compTargs", action="store_true",  default=False,   help="Whether to compare Mueller formula with Target beamlet")
parser.add_argument("--reBin",     type=int, default=5,  help="Number of bins to make into 1 in 2D histogram for smoothing")
parser.add_argument("--processes", type=int, default=4,  help="Number of processes to use in multiProcessing of raster sampling")
#Maps options:
parser.add_argument("--ampl",   type=str,     default='map', help="Range of amplitudes: map(x by y), short(nominal-10%) or large(nominal-70%)")
parser.add_argument("--eX",     type=int,     default=60,    help="End ampl X")
parser.add_argument("--eY",     type=int,     default=25,    help="End ampl Y")
parser.add_argument("--startX", type=int,     default=30,    help="Start ampl for X")
parser.add_argument("--startY", type=int,     default=10,    help="Start ampl for Y")
parser.add_argument("--Nstep", type=int,     default=7,     help="N steps")
#parser.add_argument("--NstepX", type=int,     default=7,     help="N steps for X")
#parser.add_argument("--NstepY", type=int,     default=7,     help="N steps for Y")
#Thickness Dependence option:
parser.add_argument("--minThick",type=float,  default=0.5,   help="Minimum Thickness")
parser.add_argument("--maxThick",type=float,  default=3,     help="Maximum Thickness")
parser.add_argument("--stepThick",type=float, default=0.5,   help="Step of Thicknesses")
parser.add_argument("--thickInd",type=int, default=0, choices=(0,1),   help="Position (0) or Angle (1) dependence on thickness?")
args = parser.parse_args()

#Where to save CSVs and statistics
if uname()[1] == "tensor.uio.no":
    scratchPath = "/scratch2/ericdf/PBWScatter/"
elif uname()[1] in {"heplab01.uio.no", "heplab04.uio.no","heplab03.uio.no"}:
    scratchPath = "/scratch/ericdf/Scratch/PBWScatter/"
    print(scratchPath)

if uname()[1] in {"tensor.uio.no", "heplab01.uio.no", "heplab04.uio.no","heplab03.uio.no"}:
    csvPWD = scratchPath+"CSVs/"
    homePWD = "/uio/hume/student-u52/ericdf/"
    statsPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
elif uname()[1] == "mbef-xps-13-9300":
    scratchPath = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/"
    csvPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/"
    homePWD = "/home/efackelman/"
    statsPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
else:
    csvPWD = input("Path from home to direction you like to save root files to: ")
    statsPWD = "."
paths = {'scratchPath':scratchPath, 'csvPWD':csvPWD, 'statsPWD':statsPWD, 'homePWD':homePWD}

def loopSamples(args,paths,i):
    #Get Twiss and run; if no samples argment defined, will run through once.
    print("\n\nSample",i)
    Twiss,origBX,origBY = getTwiss(args,i,paths)

    #Get full rastered for this one Twiss
    #if args.sim == "raster":
    from scatterPBW import scatterPBW
    scatterPBW(args,Twiss,i,paths,origBX,origBY)

    return Twiss,origBX,origBY

#attempting multiProcessing with Python from https://www.digitalocean.com/community/tutorials/python-multiprocessing-example
from multiprocessing import Lock, Process, Queue, current_process, Manager
import time
import queue #imported for using queue.Empty exception
print("Number of processes to use:",args.processes)

def do_job(tasks_to_accomplish, outLock,outTwiss,outOrigBX,outOrigBY):
    while True:
        try:
            '''
                try to get task from the queue. get_nowait() function will 
                raise queue.Empty exception if the queue is empty. 
                queue(False) function would do the same task also.
            '''
            task = tasks_to_accomplish.get_nowait()
        except queue.Empty:
            break
        else:
            '''
                if no exception has been raised, add the task completion 
                message to task_that_are_done queue
            '''
            #print(task)
            Twiss,origBX,origBY = loopSamples(args,paths,task)
            #print(task,Twiss)

        with outLock:
            #print("in outLock, task=",task, "process=",current_process().name)
            outTwiss[task] = Twiss #says this is undefined before it runs the first sampling, ???
            outOrigBX[task] = origBX
            outOrigBY[task] = origBY
            #print(outTwiss)

    return True


def multiLoop(args,paths,outLock,outTwiss,outOrigBX,outOrigBY): #from outputTest.txt, this is running double everything. Why?
    number_of_task = args.samples
    number_of_processes = args.processes #Number of cpus :  8 on tensor, 32 on heplab04
    tasks_to_accomplish = Queue()
    processes = []

    for i in range(number_of_task):
        tasks_to_accomplish.put(i)

    # creating processes
    if number_of_processes > 0:
        for w in range(number_of_processes):
            p = Process(target=do_job, args=(tasks_to_accomplish, outLock,outTwiss,outOrigBX,outOrigBY))
            processes.append(p)
            p.start()

        # completing process
        for p in processes:
            p.join()
    else:
        raise ValueError("number of processes < 1")

#import pdb; pdb.set_trace()

#Twiss,origBX,origBY = multiLoop(args,paths)
if args.sim =="raster":
    Twiss = {}
    origBX = {}
    origBY = {}
    with Manager() as multiManager:
        outTwiss  = multiManager.dict()
        outOrigBX = multiManager.dict()
        outOrigBY = multiManager.dict()
        outLock = Lock()

        multiLoop(args,paths, outLock,outTwiss,outOrigBX,outOrigBY)

        Twiss  = outTwiss.copy()
        origBX = outOrigBX.copy()
        origBY = outOrigBY.copy()
    print()
    print("final:")
    print(len(Twiss))

    #import pdb; pdb.set_trace()

    i=0
else:
    i=0
    Twiss,origBX,origBY = getTwiss(args,i,paths)
    #Examine individual beamlet of Twiss
    if args.sim == "beamlet":
        from beamletScatter import beamletScatter
        beamletScatter(args,Twiss,i,paths)

    #Get graph of % Outside Box and Current Density on Target for the thickness range specified
    if args.sim == "thick": #for sending only 1 Twiss.
        from thicknessDependence import thicknessDependence
        from os import uname
        if uname()[1] == "mbef-xps-13-9300": args.nP = 1e1
        print("it works!")
        thicknessDependence(args,Twiss,i,paths)

    #Get map of % Outside Box and Current Density on Target for the RMA range specified
    if args.sim == "map": #for sending only 1 Twiss.
        from mapRADependence import mapRADependence
        from os import uname
        if uname()[1] == "mbef-xps-13-9300": args.nP = 1e1
        print("it works!")
        mapRADependence(args,Twiss,i,paths,origBX,origBY)

print("Simulation took ",datetime.now()-origin,"s long\n",sep="")

if args.samples >= 2:
    from plotFit import spreadHist
    spreadHist(args,Twiss[0],paths,origBX[0],origBX[0],args.beamFile) #need to think about this, maybe pass whole discitonary down, but not sure...
