#plot spread of observables from spread of Twiss
#python3 spreadTwiss.py --betaSpread 50 --samples 200 --saveSpread

from argparse import ArgumentParser
from plotFit import spreadHist,getTwiss
from runARasterMaker import runARasterMaker
from os import uname

#Command Line arguments for save control
parser = ArgumentParser()
parser.add_argument("--sim",       type=str,   default="raster",   choices=("raster","map","beamlet","thick"), help="Type of simulation to perform")
parser.add_argument("--source",    type=str,   default="particles",choices=("particles","twiss","csv"), help="From load particles or Twiss or CSV (put name in --beamFile)")
parser.add_argument("--twissFile", type=str,   default="",    help="Load file with Twiss, auto look in OpenXAL folder")
parser.add_argument("--beamClass", type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
parser.add_argument("--particle",  type=str,   default="proton",choices=("proton","electron"), help="Which particle to simulate?")
parser.add_argument("--beamFile",  type=str,   help="Load Particles or not", default="PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03")
parser.add_argument("--twiss",     type=float, nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
parser.add_argument("--qpNum",     type=str,   default="138", help="Either a number between 099 and 148, or all")
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
args = parser.parse_args()

#Where to save CSVs and statistics
if uname()[1] == "tensor.uio.no":
    scratchPath = "/scratch2/ericdf/PBWScatter/"
elif uname()[1] in {"heplab01.uio.no", "heplab04.uio.no","heplab03.uio.no"}:
    scratchPath = "/scratch/ericdf/Scratch/PBWScatter/"
    #print(scratchPath)

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


i = 0
Twiss,origBX,origBY = getTwiss(args,i,paths)

rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,paths['csvPWD'],options,i,origBX,origBY)


spreadHist(args,Twiss,paths,origBX,origBY,rasterBeamFile)
