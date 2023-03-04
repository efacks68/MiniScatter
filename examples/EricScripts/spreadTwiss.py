#plot spread of same Twiss

from argparse import ArgumentParser
from plotFit import spreadHist,getTwiss
from os import uname

parser = ArgumentParser()
parser.add_argument("--beamClass", type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
parser.add_argument("--source",    type=str,    default="particles",choices=("particles","twiss","csv"), help="From load particles or Twiss or CSV (put name in --beamFile)")
parser.add_argument("--twiss",     type=float, nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
parser.add_argument("--betaSpread",type=float, default=0,     help="What % around provided Beta should we sample from")
parser.add_argument("--iterations",type=int,   default=10,    help="How many times to iterate this setting")
parser.add_argument("--twissFile", type=str,   default="",    help="Load file with Twiss, auto look in OpenXAL folder")
parser.add_argument("--qpNum",     type=str,   default="138", help="Either a number between 099 and 148, or all")
parser.add_argument("--Nb",        type=int,   default=10,    help="Number of macroparticles per beamlet. Default=10")
parser.add_argument("--saveSpread",action="store_true",  default=False,   help="Saves PMAS parameter spread histograms. Default=False")
args = parser.parse_args()

#Where to save CSVs and statistics
if uname()[1] == "tensor.uio.no":
    scratchPath = "/scratch2/ericdf/PBWScatter/"
elif uname()[1] in {"heplab01.uio.no", "heplab04.uio.no"}:
    scratchPath = "/scratch/ericdf/Scratch/PBWScatter/"
    #print(scratchPath)

if uname()[1] in {"tensor.uio.no", "heplab01.uio.no", "heplab04.uio.no"}:
    csvPWD = scratchPath+"CSVs/"
    statsPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
elif uname()[1] == "mbef-xps-13-9300":
    csvPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/"
    statsPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
else:
    csvPWD = input("Path from home to direction you like to save root files to: ")
    statsPWD = "."
paths = {'scratchPath':scratchPath, 'csvPWD':csvPWD, 'statsPWD':statsPWD}

i = 0
Twiss = getTwiss(args,i,paths)
#if args.beamClass == 'Yngve': #smallest from Yngve
#    Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
#elif args.beamClass == 'ESS': #from my OpenXAL calculation
#    Twiss = [0.11315,1006.80,-60.44,0.12155,129.72,-7.72]
#    #pOff = 10 #negative to have less than
#    Twiss[1] = Twiss[1] * (1 + args.pOff/100)
#    Twiss[4] = Twiss[4] * (1 + args.pOff/100)
#elif args.beamClass == 'pencil': #"pencil" beam of ~0 emittance
#    Twiss = [0.0001,0.15,0,0.0001,0.15,0]
#if args.twiss:
#        Twiss = [args.twiss[0],args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5]]
#        args.beamClass = "Twiss"

spreadHist(args,Twiss,paths)
