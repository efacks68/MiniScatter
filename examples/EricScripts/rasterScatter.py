#rasterScatter.py
#Takes config for scattering simulations and runs
#Includes setup and output preferences

#possible commands:
#python3 rasterScatter.py --source twiss --twissFile TwissRange0,12mm-mrad_3Bx-3By-2Ax-2Ay --qpNum 0 --iterations 10 --pOff -10
    #10x smallest Twiss -10%
#python3 rasterScatter.py --iterations 50 --pOff -5
    #50x nominal Twiss - 5%


def run(args,iteration):
    #Get Twiss to run, depending on user configuration
    # Twiss= [NemtX,BetaX,AlphX,NemtY,BetaY,AlphY]
    if args.beamClass == 'Yngve': #smallest from Yngve
        Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
    elif args.beamClass == 'ESS': #from my OpenXAL calculation
        Twiss = [0.11315,1006.80,-60.44,0.12155,129.72,-7.72]
        #pOff = 0 #negative to have less than
        Twiss[1] = Twiss[1] * (1 + args.pOff/100)
        Twiss[4] = Twiss[4] * (1 + args.pOff/100)
    elif args.beamClass == 'pencil': #"pencil" beam of ~0 emittance
        Twiss = [0.0001,0.15,0,0.0001,0.15,0]

    if args.twiss:
        Twiss = [args.twiss[0],args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5]]
        args.beamClass = "Twiss"

    #Tailored to OpenXAL failure input
    if args.source == "twiss":
        if args.twissFile != "":
            import csv
            from os import uname
            i=0
            #Find file
            if uname()[1] == "tensor.uio.no":
                twissPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/OpenXAL/OXALNotebooks/failureTwiss/"
            elif uname()[1] == "mbef-xps-13-9300":
                twissPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/OpenXAL/OXALNotebooks/failureTwiss/"
            if args.qpNum == "all":
                #for running all Twiss in OpenXAL Combo Sequence output CSV
                from plotFit import numLines
                n = numLines(twissPWD+args.twissFile)
                for i in range(n):
                    with open(twissPWD+args.twissFile+".csv",mode='r') as csv_file:
                        csv_reader = csv.reader(csv_file, delimiter=',')
                        rowNum = 0 #to increment QP number
                        for row in csv_reader:
                            #if rowNum == 0:
                            #    print("row 0")
                            #    if type(row[1]) == str:
                            #        print("Header skipped")
                            #        continue
                                    #next(csv_reader,None)
                            if rowNum == i: #could be done cleaner?
                                if float(row[1]) < 1e-4:
                                    um = 1e6 #if from OpenXAL 
                                elif float(row[1]) > 1e-4:
                                    um = 1
                                #print("um",um)
                                Twiss[0] = float(row[1])*um #[mm-mrad]
                                Twiss[1] = float(row[2])
                                Twiss[2] = float(row[3])
                                Twiss[3] = float(row[4])*um #[mm-mrad]
                                Twiss[4] = float(row[5])
                                Twiss[5] = float(row[6])
                                print(Twiss)
                                break
                                #Adjusts names so the files are in order
                                if len(row[0]) == 2: 
                                    args.qpNum = "0"
                                else: 
                                    args.qpNum = ""
                                args.qpNum += row[0]
                                #print(i,rowNum,args.qpNum)
                                break
                            rowNum+=1
                    csv_file.close()
                    
                    if args.sim == "raster":
                        from scatterPBW import scatterPBW
                        scatterPBW(args,Twiss,iteration)
            else:
                #for single QP fail runs
                with open(twissPWD+args.twissFile+".csv",mode='r') as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    for row in csv_reader:
                        if row[0] == args.qpNum:
                            if float(row[1]) < 1e-4:
                                um = 1e6 #if from OpenXAL 
                            elif float(row[1]) > 1e-4:
                                um = 1
                            #print("um",um)
                            Twiss[0] = float(row[1])*um #[mm-mrad]
                            Twiss[1] = float(row[2])
                            Twiss[2] = float(row[3])
                            Twiss[3] = float(row[4])*um #[mm-mrad]
                            Twiss[4] = float(row[5])
                            Twiss[5] = float(row[6])
                    Twiss[1] = Twiss[1] * (1 + args.pOff/100)
                    Twiss[4] = Twiss[4] * (1 + args.pOff/100)
                    print(Twiss)
                csv_file.close()

    #Get full rastered for this one Twiss
    if args.sim == "raster":
        from scatterPBW import scatterPBW
        scatterPBW(args,Twiss,iteration)

    #Get map of % Outside Box and Current Density on Target for the range specified
    if args.sim == "map":
        from mapRADependence import mapRADependence
        from os import uname
        if uname()[1] == "mbef-xps-13-9300": args.nP = 1e1
        print("it works!")
        mapRADependence(args,Twiss,iteration)

    #Examine individual beamlet of Twiss
    if args.sim == "beamlet":
        from beamletScatter import beamletScatter
        beamletScatter(args,Twiss,iteration)
    return Twiss

from datetime import datetime
origin = datetime.now()
from argparse import ArgumentParser
#look into argument groups!

#Command Line arguments for save control
parser = ArgumentParser()
parser.add_argument("--sim",       type=str,    default="raster",   choices=("raster","map","beamlet"), help="Type of simulation to perform")
parser.add_argument("--source",    type=str,    default="particles",choices=("particles","twiss"), help="From load particles or Twiss?")
#General Beam Setup Options
parser.add_argument("--beamClass", type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
parser.add_argument("--particle",  type=str,   default="proton", choices=("proton","electron"), help="Which particle to simulate?")
parser.add_argument("--beamFile",  type=str,   help="Load Particles or not",   default="PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03")
parser.add_argument("--twiss",     type=float, nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
parser.add_argument("--twissFile", type=str,   default="FailureA2T", help="Load file with Twiss, auto look in OpenXAL folder")
parser.add_argument("--qpNum",     type=str,   default="138", help="Either a number between 099 and 148, or all")
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
parser.add_argument("--ylim",      type=float, default=150,   help="+/- value for vertical axis of output rastered image [mm]. Default=150")
parser.add_argument("--maxim",     type=float, default=0  ,   help="Maximum current density value for output rastered imagem[uA/cm^2]. Default=0")
parser.add_argument("--edgeRaster",action="store_true",  help="Only populate edges of raster. Default=False")
parser.add_argument("--PBIP",      action="store_true",  default=False,   help="Is PBIP present? Default=False")
parser.add_argument("--material",  type=str,   default="Al",  choices=("Al","Au","C","Vac"),  help="What material PBW?")
parser.add_argument("--Nbeamlet",  type=float, default=1e5,   help="For beamlet simulation, how many protons?")
parser.add_argument("--iterations", type=int,   default=1, help="How many times to iterate this setting")
parser.add_argument("--pOff", type=int,default=0, help="What % off of nominal should the Twiss be?")
#Output Options
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
#Maps options:
parser.add_argument("--ampl",   type=str,     default='map', help="Range of amplitudes: map(x by y), short(nominal-10%) or large(nominal-70%)")
parser.add_argument("--eX",     type=int,     default=55,    help="End ampl X")
parser.add_argument("--eY",     type=int,     default=20,    help="End ampl Y")
parser.add_argument("--startX", type=int,     default=30,    help="Start ampl for X")
parser.add_argument("--startY", type=int,     default=10,    help="Start ampl for Y")
parser.add_argument("--NstepX", type=int,     default=6,     help="N steps for X")
parser.add_argument("--NstepY", type=int,     default=6,     help="N steps for Y")
args = parser.parse_args()

for i in range(args.iterations):
    Twiss = run(args,i)

print("Simulation took ",datetime.now()-origin,"s long",sep="")

if args.iterations >= 2:
    from spreadTwiss import spreadHist
    spreadHist(args,Twiss,args.iterations)
