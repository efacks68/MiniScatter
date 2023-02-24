#plot spread of same Twiss

from argparse import ArgumentParser
from plotFit import spreadHist

parser = ArgumentParser()
parser.add_argument("--beamClass", type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
parser.add_argument("--twiss",     type=float, nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
parser.add_argument("--pOff",      type=int,   default=0)
parser.add_argument("--iterations",type=int,   default=10)
args = parser.parse_args()

if args.beamClass == 'Yngve': #smallest from Yngve
    Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
elif args.beamClass == 'ESS': #from my OpenXAL calculation
    Twiss = [0.11315,1006.80,-60.44,0.12155,129.72,-7.72]
    #pOff = 10 #negative to have less than
    Twiss[1] = Twiss[1] * (1 + args.pOff/100)
    Twiss[4] = Twiss[4] * (1 + args.pOff/100)
elif args.beamClass == 'pencil': #"pencil" beam of ~0 emittance
    Twiss = [0.0001,0.15,0,0.0001,0.15,0]
if args.twiss:
        Twiss = [args.twiss[0],args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5]]
        args.beamClass = "Twiss"

spreadHist(args,Twiss)