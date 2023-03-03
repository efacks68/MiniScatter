#rasterScatter.py
#Takes config for scattering simulations and runs
#Includes setup and output preferences

#possible commands:
#python3 rasterScatter.py --source twiss --twissFile TwissRange0,12mm-mrad_3Bx-3By-2Ax-2Ay --qpNum 0 --iterations 10 --betaSpread 10
    #10x smallest Twiss -10%
#python3 rasterScatter.py --iterations 50 --betaSpread 1
    #50x nominal Twiss - 5%

from datetime import datetime
origin = datetime.now()
from argparse import ArgumentParser
from os import uname
from plotFit import getTwiss
#look into argument groups!

#Command Line arguments for save control
parser = ArgumentParser()
parser.add_argument("--sim",       type=str,    default="raster",   choices=("raster","map","beamlet"), help="Type of simulation to perform")
parser.add_argument("--source",    type=str,    default="particles",choices=("particles","twiss","csv"), help="From load particles or Twiss or CSV (put name in --beamFile)")
parser.add_argument("--twissFile", type=str,   default="", help="Load file with Twiss, auto look in OpenXAL folder")
parser.add_argument("--beamClass", type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
parser.add_argument("--particle",  type=str,   default="proton", choices=("proton","electron"), help="Which particle to simulate?")
parser.add_argument("--beamFile",  type=str,   help="Load Particles or not",   default="PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03")
parser.add_argument("--twiss",     type=float, nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
parser.add_argument("--qpNum",     type=str,   default="138", help="Either a number between 099 and 148, or all")
parser.add_argument("--betaSpread",type=float, default=0,     help="What % around provided Beta should we sample from")
parser.add_argument("--iterations",type=int,   default=1,     help="How many times to iterate this setting")
parser.add_argument("--csvFile",   type=str,   default="",    help="Load Beam of already made csv")
#General Beam Setup Options
parser.add_argument("--t",         type=float, default=0,     help="PBW Thickness [mm], 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al Xmm thick. Default=0")
parser.add_argument("--energy",    type=float, default=570,   help="Beam Energy [MeV]. Default=570")
parser.add_argument("--Nb",        type=int,   default=10,    help="Number of macroparticles per beamlet. Default=10")
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
parser.add_argument("--picFormat", type=str,   default="png",  choices=("png","svg","pdf"),help="Whic file format extension?")
parser.add_argument("--matPlots",  action="store_true",  default=False,   help="Whether to do various material plots for beamlets")
parser.add_argument("--saveSpread",action="store_true",  default=False,   help="Saves PMAS parameter spread histograms. Default=False")
#Maps options:
parser.add_argument("--ampl",   type=str,     default='map', help="Range of amplitudes: map(x by y), short(nominal-10%) or large(nominal-70%)")
parser.add_argument("--eX",     type=int,     default=55,    help="End ampl X")
parser.add_argument("--eY",     type=int,     default=20,    help="End ampl Y")
parser.add_argument("--startX", type=int,     default=30,    help="Start ampl for X")
parser.add_argument("--startY", type=int,     default=10,    help="Start ampl for Y")
parser.add_argument("--NstepX", type=int,     default=6,     help="N steps for X")
parser.add_argument("--NstepY", type=int,     default=6,     help="N steps for Y")
args = parser.parse_args()

#Where to save CSVs and statistics
if uname()[1] == "tensor.uio.no":
    scratchPath = "/scratch2/ericdf/PBWScatter/"
elif uname()[1] in {"heplab01.uio.no", "heplab04.uio.no"}:
    scratchPath = "/scratch/ericdf/Scratch/PBWScatter/"
    #print(scratchPath)

if uname()[1] in {"tensor.uio.no", "heplab01.uio.no", "heplab04.uio.no"}:
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
paths = {'scratchPath':scratchPath, 'csvPWD':csvPWD, 'statsPWD':statsPWD}

#get Twiss and run; if no iterations argment defined, will run through once.
for i in range(args.iterations):
    print("\n\nIteration",i)
    Twiss = getTwiss(args,i,paths)

    #Get full rastered for this one Twiss
    if args.sim == "raster":
        from scatterPBW import scatterPBW
        scatterPBW(args,Twiss,i,paths)

    #Get map of % Outside Box and Current Density on Target for the range specified
    if args.sim == "map":
        from mapRADependence import mapRADependence
        from os import uname
        if uname()[1] == "mbef-xps-13-9300": args.nP = 1e1
        print("it works!")
        mapRADependence(args,Twiss,i,paths)

    #Examine individual beamlet of Twiss
    if args.sim == "beamlet":
        from beamletScatter import beamletScatter
        beamletScatter(args,Twiss,i,paths)

print("Simulation took ",datetime.now()-origin,"s long",sep="")

if args.iterations >= 2:
    from plotFit import spreadHist
    spreadHist(args,Twiss,paths)
