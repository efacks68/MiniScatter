#Script for single plots of 
# optimization value: % outside 99% box

from datetime import datetime
origin = datetime.now()
print(origin,"\n")
import numpy as np
import matplotlib.pyplot as plt
import os,csv,argparse
from runARasterMaker import runARasterMaker
from runPBW import runPBW

#Command Line arguments for save control
parser = argparse.ArgumentParser()
parser.add_argument("--beamType",type=str,default="ESS",help="Determines beam Twiss: 'ESS', 'Yngve', 'pencil', or 'twiss'? 'twiss' expects '--twiss' argument with Twiss in MiniScatter format")
parser.add_argument("--savePics",action="store_true",default=False,help="Save Pictures?")
parser.add_argument("--l",type=str,help="Load Particles or not",default="PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03")
parser.add_argument("--t",type=float,default=0,help="PBW Thickness, 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al X [mm] thick")
parser.add_argument("--twiss",type=float,nargs=8,help="N particles,NemtX,BetaX,AlphX,NemtY,BetaY,AlphY,Thickness")
parser.add_argument("--Nb",type=int,default=10,help="Number of macroparticles per beamlet")
parser.add_argument("--nP",type=float,default=1e3,help="Numper of beamlets in pulse")
parser.add_argument("--rX",type=float,default=0,help="X distance from beam axis")
parser.add_argument("--rY",type=float,default=0,help="Y distance from beam axis")
parser.add_argument("--aX",type=float,default=54.65,help="RM X Amplitude")
parser.add_argument("--aY",type=float,default=18.37,help="RM Y Amplitude")
parser.add_argument("--edges",action="store_true",help="Only populate edges of raster?")
parser.add_argument("--PBIP",action="store_true",default=False,help="Is PBIP present?")
args = parser.parse_args()

#Constants for running scripts
energy      = 570 #[MeV]
graph       = False
physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ"
dependence  = "Twiss"
beamClass   = "ESS" #classification for runARasterMaker function
# Twiss= [NemtX,BetaX,AlphX,NemtY,BetaY,AlphY]
if args.beamType == 'Yngve': #smallest from Yngve
  Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
elif args.beamType == 'ESS': #from my OpenXAL calculation
  Twiss = [0.11315,1006.80,-60.44,0.12155,129.72,-7.72]
elif args.beamType == 'pencil': #"pencil" beam of 0 emittance
  Twiss = [0.0001,0.15,0,0.0001,0.15,0]
  beamClass="pencil" #not ESS class beam
elif args.beamType == 'twiss':
  Twiss = [args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5],args.twiss[6]]

csvPWD = "/scratch2/ericdf/PBWScatter/CSVs/"
#Create Rastered Beam file, runARasterMaker checks if the CSV is already present
rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(beamClass,energy,graph,args.Nb,args.nP,args.rX,args.rY,args.edges,Twiss,args.aX,args.aY,dependence,csvPWD)
print(rasterBeamFile,beamXAngle,beamYAngle)
#Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
ImgPOutBox = runPBW(energy,rasterBeamFile,args.beamType,args.t,beamXAngle,beamYAngle,args.PBIP,args.savePics,physList,Twiss,args.aX,args.aY,dependence)

print("Simulation took",datetime.now()-origin,"long")