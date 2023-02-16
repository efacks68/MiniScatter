#rasterScatter.py
#Takes config for scattering simulations and runs
#Includes setup and output preferences

from argparse import ArgumentParser
#look into argument groups!

#Command Line arguments for save control
parser = ArgumentParser()
parser.add_argument("--sim",       type=str,    default="raster",choices=("raster","map","beamlet"), help="Type of simulation to perform")
parser.add_argument("--source",    type=str,    default="load",     choices=("load","twiss"), help="From load particles or Twiss?")
#General Beam Setup Options
parser.add_argument("--beamClass", type=str,    default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
parser.add_argument("--particle",  type=str,    default="proton", choices=("proton","electron"), help="Which particle to simulate?")
parser.add_argument("--beamFile",  type=str,    help="Load Particles or not",   default="PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03")
parser.add_argument("--twiss",     type=float,  nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
parser.add_argument("--twissFile", type=str,    default="~/Documents/UiO/Forske/ESSProjects/OpenXAL/OXALNotebooks/failureTwiss/FailureA2T_QP138.csv", help="Load file with Twiss")
parser.add_argument("--t",         type=float,  default=0,     help="PBW Thickness [mm], 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al Xmm thick. Default=0")
parser.add_argument("--energy",    type=float,  default=570,   help="Beam Energy [MeV]. Default=570")
parser.add_argument("--Nb",        type=int,    default=10,    help="Number of macroparticles per beamlet. Default=10")
parser.add_argument("--nP",        type=float,  default=1e3,   help="Numper of beamlets in pulse. Default=1e3")
parser.add_argument("--rX",        type=float,  default=0,     help="X distance from beam axis [mm]. Default=0")
parser.add_argument("--rY",        type=float,  default=0,     help="Y distance from beam axis [mm]. Default=0")
parser.add_argument("--aX",        type=float,  default=54.65, help="RM X Amplitude [mm]. Default=54.65")
parser.add_argument("--aY",        type=float,  default=18.37, help="RM Y Amplitude [mm]. Default=18.37")
parser.add_argument("--failure",   type=float,  default=0,     choices = range(0,5),  help="Which RM Failure case, 0-4. Default=0")
parser.add_argument("--magFails",  type=int,    default=2,     choices = range(0,5),  help="Number of Raster Magnets that fail, 1-4. Default=2")
parser.add_argument("--xlim",      type=float,  default=150,   help="+/- value for horizontal axis of output rastered image [mm]. Default=150")
parser.add_argument("--ylim",      type=float,  default=150,   help="+/- value for vertical axis of output rastered image [mm]. Default=150")
parser.add_argument("--maxim",     type=float,  default=0  ,   help="Maximum current density value for output rastered imagem[uA/cm^2]. Default=0")
parser.add_argument("--edgeRaster",action="store_true",  help="Only populate edges of raster. Default=False")
parser.add_argument("--PBIP",      action="store_true",  default=False,   help="Is PBIP present? Default=False")
parser.add_argument("--material",  type=str,    default="Al",   choices=("Al","Au","C","Vac"),  help="What material PBW?")
#Output Options
parser.add_argument("--rCut",      type=float,     default=1e3,  help="Radial cut, defines worldSize and histLims")
parser.add_argument("--engCut",    type=float,     default=0.9,  help="Energy cut, see MiniScatter description")
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
#Maps options:
parser.add_argument("--ampl",   type=str,     default='map', help="Range of amplitudes: map(x by y), short(nominal-10%) or large(nominal-70%)")
parser.add_argument("--eX",     type=int,     default=55,    help="End ampl X")
parser.add_argument("--eY",     type=int,     default=20,    help="End ampl Y")
parser.add_argument("--startX", type=int,     default=30,    help="Start ampl for X")
parser.add_argument("--startY", type=int,     default=10,    help="Start ampl for Y")
parser.add_argument("--NstepX", type=int,     default=6,     help="N steps for X")
parser.add_argument("--NstepY", type=int,     default=6,     help="N steps for Y")

args = parser.parse_args()

#Get Twiss, put in function in other file? send args, return Twiss?
# Twiss= [NemtX,BetaX,AlphX,NemtY,BetaY,AlphY]
if args.beamClass == 'Yngve': #smallest from Yngve
    Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
elif args.beamClass == 'ESS': #from my OpenXAL calculation
    Twiss = [0.11315,1006.80,-60.44,0.12155,129.72,-7.72]
elif args.beamClass == 'pencil': #"pencil" beam of ~0 emittance
    Twiss = [0.0001,0.15,0,0.0001,0.15,0]

if args.twiss:
    Twiss = [args.twiss[0],args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5]]
    args.beamClass = "Twiss"

if args.source == "twiss":
    if args.twissFile != "":
        import csv
        i=0
        with open(args.twissFile,mode='r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                if i==2:break
                Twiss[0] = float(row[0])
                Twiss[1] = float(row[1])
                Twiss[2] = float(row[2])
                Twiss[3] = float(row[3])
                Twiss[4] = float(row[4])
                Twiss[5] = float(row[5])
                i+=1
            print(Twiss)
        csv_file.close()

if args.sim == "raster":
    from scatterPBW import scatterPBW
    scatterPBW(args,Twiss)

if args.sim == "map":
    from mapRADependence import mapRADependence
    from os import uname
    if uname()[1] == "mbef-xps-13-9300": args.nP = 1e1
    print("it works!")
    mapRADependence(args,Twiss)

if args.sim == "beamlet":
    from beamletScatter import beamletScatter

    beamletScatter(args,Twiss)
    print("Need to figure this out!")