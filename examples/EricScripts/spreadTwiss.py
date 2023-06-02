#plot spread of observables from spread of Twiss
    #python3 spreadTwiss.py --samples 200 --saveSpread --Nb 100
#Jitter:
    #python3 spreadTwiss.py --source twiss --twissFile HEBT-A2T_100pctField_1.0e-03Jitter_200x --beamClass jitter --saveSpread --samples 200
#QP Fail:
    #python3 spreadTwiss.py --source twiss --twissFile HEBT-A2T_QP100_80pctField_1.0e-04Jitter_200x  --beamClass qpFail --qpNum 100 --saveSpread --samples 200
#Twiss param vs pOut
    #python3 spreadTwiss.py --source twiss --twissFile HEBT-A2T_100pctField_1.0e-03Jitter_200x --beamClass jitter --samples 200 --saveFits

from argparse import ArgumentParser
from plotFit import spreadHist,getTwiss
from runARasterMaker import runARasterMaker
from os import uname


#Command Line arguments for save control
parser = ArgumentParser()
parser.add_argument("--sim",       type=str,   default="raster",   choices=("raster","map","beamlet","thick"), help="Type of simulation to perform")
parser.add_argument("--source",    type=str,   default="particles",choices=("particles","twiss","csv"), help="From load particles, Twiss, or CSV (put name in --beamFile)")
parser.add_argument("--twissFile", type=str,   default="",    help="Load file with Twiss, auto look in OpenXAL folder")
parser.add_argument("--twissRowN", type=int,   default=99,    help="What rowNum to use in getting QP failure Twiss")
parser.add_argument("--beamClass", type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', 'pencil','qpFail'(expects a qpNum), or 'jitter'. If other, just do --twiss. Default='ESS'")
parser.add_argument("--particle",  type=str,   default="proton",choices=("proton","electron"), help="Which particle to simulate?")
parser.add_argument("--beamFile",  type=str,   help="Load Particles or not", default="PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100")
parser.add_argument("--twiss",     type=float, nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
parser.add_argument("--qpNum",     type=str,   default="",    help="Either a number between 099 and 148, qps, or all, see getTwiss function")
parser.add_argument("--betaSpread",type=float, default=0,     help="What % around provided Beta should we sample from")
parser.add_argument("--samples",   type=int,   default=1,     help="How many times to sample this setting")
parser.add_argument("--statsFile", type=str,   default="EvalStats27Apr", help="Load Beam of already made csv")
#General Beam Setup Options
parser.add_argument("--t",         type=float, default=0,     help="PBW Thickness [mm], 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al Xmm thick. Default=0")
parser.add_argument("--energy",    type=float, default=570,   help="Beam Energy [MeV]. Default=570")
parser.add_argument("--Nb",        type=int,   default=100,   help="Number of macroparticles per beamlet. Default=10")
parser.add_argument("--nP",        type=float, default=1e3,   help="Numper of beamlets in pulse. Default=1e3")
parser.add_argument("--rX",        type=float, default=0,     help="X distance from beam axis at BPM94 [mm]. Default=0")
parser.add_argument("--rY",        type=float, default=0,     help="Y distance from beam axis at BPM94 [mm]. Default=0")
parser.add_argument("--aX",        type=float, default=48.7867, help="RM X Amplitude [mm]. Default=48.7867") #14.3.23-new RMA at PBW calculations with Synoptic Viewer distances. old - 54.65
parser.add_argument("--aY",        type=float, default=16.3991, help="RM Y Amplitude [mm]. Default=16.3991") #old - 18.37
parser.add_argument("--failure",   type=float, default=0,     choices = range(0,5),  help="Which RM Failure case, 0-4. Default=0")
parser.add_argument("--magFails",  type=int,   default=2,     choices = range(0,5),  help="Number of Raster Magnets that fail, 1-4. Default=2")
parser.add_argument("--xlim",      type=float, default=150,   help="+/- value for horizontal axis of output rastered image [mm]. Default=150")
parser.add_argument("--ylim",      type=float, default=100,   help="+/- value for vertical axis of output rastered image [mm]. Default=150")
parser.add_argument("--maxim",     type=float, default=0  ,   help="Maximum current density value for output rastered imagem[uA/cm^2]. Default=0")
parser.add_argument("--edgeRaster",action="store_true",       help="Only populate edges of raster. Default=False")
parser.add_argument("--PBIP",      action="store_true",       default=False,   help="Is PBIP present? Default=False")
parser.add_argument("--material",  type=str,   default="Al",  choices=("Al","Au","C","Vac"),  help="What material PBW?")
parser.add_argument("--Nbeamlet",  type=float, default=1e5,   help="For beamlet simulation, how many protons?")
#Output Options                         Lowering the rCut decreases bin size which increases the jMax...
parser.add_argument("--rCut",      type=float, default=1e3,  help="Radial cut, defines worldSize and histLims")
parser.add_argument("--engCut",    type=float, default=0.9,  help="Energy cut, see MiniScatter description")
parser.add_argument("--noText",    action="store_true",  default=False,   help="Turns off printed text when called. Default=False")
parser.add_argument("--noBox",     action="store_true",  default=False,   help="Turns off printed box when called. Default=False")
parser.add_argument("--savePics",  action="store_true",  default=False,   help="Saves Rastered Image. Default=False")
parser.add_argument("--saveGrads", action="store_true",  default=False,   help="Plots gradients of beam at Target. Default=False")
parser.add_argument("--saveEdges", action="store_true",  default=False,   help="Plots edges on Raster Image. Default=False")
parser.add_argument("--gaussFit",  action="store_true",  default=False,   help="Computes sum of Gaussian fits for central axis projection. Default=False")
parser.add_argument("--saveFits",  action="store_true",  default=False,   help="Saves plots of Gaussian Fitting. Default=False")
parser.add_argument("--saveHist",  action="store_true",  default=False,   help="Saves Histogram of proton density at target for R Compare. Default=False")
parser.add_argument("--saveRaster",action="store_true",  default=False,   help="Saves plot of rastered beam. Default=False")
parser.add_argument("--savePull",  action="store_true",  default=False,   help="Saves Pull Values of Particle Distribution. Default=False")
parser.add_argument("--picFormat", choices=("png","jpeg","svg","pdf"), type=str, default="png",  help="Which picture format extension?")
parser.add_argument("--matPlots",  action="store_true",  default=False,   help="Whether to do various material plots for beamlets")
parser.add_argument("--saveSpread",action="store_true",  default=False,   help="Saves PMAS parameter spread histograms. Default=False")
parser.add_argument("--saveParts", action="store_true",  default=False,   help="Saves initial Particle info in MiniScatter format(x,px,y,py,z,E). Default=False")
parser.add_argument("--compTargs", action="store_true",  default=False,   help="Whether to compare Mueller formula with Target beamlet")
parser.add_argument("--reBin",     type=int, default=4,  help="Number of bins to make into 1 in 2D histogram for smoothing")
parser.add_argument("--processes", type=int, default=4,  help="Number of processes to use in multiProcessing of raster sampling")
parser.add_argument("--dpi",       type=int, default=500,help="DPI for pngs")
parser.add_argument("--physList",  type=str, default="QGSP_BERT_EMZ",help="Physics List, either 'QGSP_BERT_EMZ','FTFP_BERT_EMZ', or 'QGSP_BERT__SS'")
parser.add_argument("--threshold", type=int, default=10,help="threshold of particles/bin to include in Chi2 and pull calculations")
#parser.add_argument("--refImg",    action="store_true",  default=False,   help="Add Img to the Ref Img of this beam")
#Maps options:
parser.add_argument("--ampl",   type=str,     default='map', help="Range of amplitudes: map(x by y), short(nominal-10%) or large(nominal-70%)")
parser.add_argument("--startX", type=float,     default=40,    help="Start ampl for X")
parser.add_argument("--eX",     type=float,     default=50,    help="End ampl X")
parser.add_argument("--startY", type=float,     default=13,    help="Start ampl for Y")
parser.add_argument("--eY",     type=float,     default=18,    help="End ampl Y")
parser.add_argument("--Nstep",  type=int,     default=6,     help="N steps; for defaults gives whole number values for map")
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
    #print(scratchPath)

if uname()[1] in {"tensor.uio.no", "heplab01.uio.no", "heplab04.uio.no","heplab03.uio.no"}:
    csvPWD = scratchPath+"CSVs/"
    homePWD = "/uio/hume/student-u52/ericdf/"
    statsPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    oXALPWD = homePWD+"Documents/UiO/Forske/ESSProjects/OpenXAL/OXALNotebooks/failureTwiss/"
elif uname()[1] == "mbef-xps-13-9300":
    scratchPath = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/"
    csvPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/"
    homePWD = "/home/efackelman/"
    statsPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    oXALPWD = homePWD+"Documents/UiO/Forske/ESSProjects/OpenXAL/OXALNotebooks/failureTwiss/"
else:
    csvPWD = input("Path from home to directory you would like to save root files to: ")
    statsPWD = "."
paths = {'scratchPath':scratchPath, 'csvPWD':csvPWD, 'statsPWD':statsPWD, 'homePWD':homePWD, 'oXALPWD':oXALPWD}

#Constants for running scripts
physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ" or "QGSP_BERT__SS"
zoff = "*-1" #[mm] with preappended * to keep covar defined at z=0
options     = {'physList':physList, 'dependence':"Twiss", 'zoff':zoff, 'initTree':False,
                'exitTree':False, 'targetTree':False, 'MCS':False, 'engPlot':False,
                'mat3Plot':False, 'MiniRoot':True  }

#Important things
if args.t == 0:
    args.material = "Al" #overwrites potential user input. Needs work
elif args.t == 0.1:
    args.material = "Vac"
boxes = [0]#,-.25,-.375,-.45,-.50]#,0.125,0.25,0.375] #make an args for 24.11.22


i = 0
Twiss,origBX,origBY = getTwiss(args,i,paths)

rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,paths['csvPWD'],options,i,origBX,origBY)
from simulation import setup
savename,simSetup_simple1 = setup(args,args.material,rasterBeamFile,Twiss,options,paths)

if not args.saveFits:
    spreadHist(args,Twiss,paths,origBX,origBY,rasterBeamFile)

if args.saveFits:
    from betaPOutFit import betaPOutFit
    from alphaPOutFit import alphaPOutFit
    from emittPOutFit import emittPOutFit
    from betaAlphaFit import betaAlphaFit
    #betaPOutFit(args,Twiss,paths,origBX,origBY,rasterBeamFile,"x")
    #betaPOutFit(args,Twiss,paths,origBX,origBY,rasterBeamFile,"y")
    #alphaPOutFit(args,Twiss,paths,origBX,origBY,rasterBeamFile,"x")
    #alphaPOutFit(args,Twiss,paths,origBX,origBY,rasterBeamFile,"y")
    #emittPOutFit(args,Twiss,paths,origBX,origBY,rasterBeamFile,"x")
    #emittPOutFit(args,Twiss,paths,origBX,origBY,rasterBeamFile,"y")
    betaAlphaFit(args,Twiss,paths,origBX,origBY,rasterBeamFile,"x")
    betaAlphaFit(args,Twiss,paths,origBX,origBY,rasterBeamFile,"y")
