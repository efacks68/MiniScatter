#Script for single plots of 
# optimization value: % outside 99% box
#typical command: (no PBW, no text, no boxes, save pic of 200x200mm^2 window)
#python3 scatterPBW.py --savePics --text --xlim 200 --ylim 200 --box --t 0.1 

from datetime import datetime
origin = datetime.now()
print(origin,"\n")
import numpy as np
#import matplotlib.pyplot as plt
import os,argparse
from runARasterMaker import runARasterMaker
from runPBW import runPBW

#Command Line arguments for save control
parser = argparse.ArgumentParser()
parser.add_argument("--beamClass",type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss")
parser.add_argument("--l",       type=str,    help="Load Particles or not",   default="PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03")
parser.add_argument("--t",       type=float,  default=0,     help="PBW Thickness, 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al X [mm] thick")
parser.add_argument("--energy",  type=float,  default=570,   help="Beam Energy")
parser.add_argument("--Nb",      type=int,    default=10,    help="Number of macroparticles per beamlet")
parser.add_argument("--nP",      type=float,  default=1e3,   help="Numper of beamlets in pulse")
parser.add_argument("--rX",      type=float,  default=0,     help="X distance from beam axis")
parser.add_argument("--rY",      type=float,  default=0,     help="Y distance from beam axis")
parser.add_argument("--aX",      type=float,  default=54.65, help="RM X Amplitude")
parser.add_argument("--aY",      type=float,  default=18.37, help="RM Y Amplitude")
parser.add_argument("--xlim",    type=float,  default=450,   help="+/- value for horizontal axis of output rastered image")
parser.add_argument("--ylim",    type=float,  default=500,   help="+/- value for vertical axis of output rastered image")
parser.add_argument("--maxim",   type=float,  default=0  ,   help="Maximum current density value for output rastered image")
parser.add_argument("--twiss",   type=float,  nargs=6,       help="Twiss parameters in form: NemtX,BetaX,AlphX,NemtY,BetaY,AlphY")
parser.add_argument("--edges",   action="store_true",  help="Only populate edges of raster?")
parser.add_argument("--PBIP",    action="store_true",  default=False,   help="Is PBIP present?")
parser.add_argument("--text",    action="store_false", default=True,    help="Print texts on images, default is True, calling means no text")
parser.add_argument("--savePics",action="store_true",  default=False,   help="Saves Rastered Image")
parser.add_argument("--box",     action="store_false", default=True,    help="Turns off printed box")
parser.add_argument("--saveHist",action="store_true",  default=False,   help="Saves Histogram of proton density at target")
args = parser.parse_args()

#Constants for running scripts
graph       = True
physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ"
dependence  = "Twiss"
options     = {'text':args.text, 'box':args.box, 'wide':True, 'physList':physList, 'dependence':dependence,
                            'xlim':args.xlim, 'ylim':args.ylim, 'maxim':args.maxim, 'saveHist':args.saveHist,
                            'PBIP':args.PBIP, 'beamClass':args.beamClass}

# Twiss= [NemtX,BetaX,AlphX,NemtY,BetaY,AlphY]
if args.beamClass == 'Yngve': #smallest from Yngve
  Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
elif args.beamClass == 'ESS': #from my OpenXAL calculation
  Twiss = [0.11315,1006.80,-60.44,0.12155,129.72,-7.72]
elif args.beamClass == 'pencil': #"pencil" beam of ~0 emittance
  Twiss = [0.0001,0.15,0,0.0001,0.15,0]
if args.twiss:
  Twiss = [args.twiss[0],args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5]]
  options['beamClass'] = "Twiss"
if os.uname()[1] == "tensor.uio.no":
  csvPWD = "/scratch2/ericdf/PBWScatter/CSVs/"
if os.uname()[1] == "mbarrios-XPS-13-9300":
  csvPWD = "~/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWscripts/"
#Create Rastered Beam file, runARasterMaker checks if the CSV is already present
rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args.energy,graph,args.Nb,args.nP,args.rX,args.rY,args.edges,Twiss,args.aX,args.aY,csvPWD,options)
print(rasterBeamFile,beamXAngle,beamYAngle)
#Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
ImgPOutBox = runPBW(args.energy,rasterBeamFile,args.t,beamXAngle,beamYAngle,args.savePics,Twiss,args.aX,args.aY,options)

print("Simulation took",datetime.now()-origin,"long")