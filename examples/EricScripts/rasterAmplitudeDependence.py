#Script for searching Raster Magnet Amplitude space
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
parser.add_argument("--savePics",action="store_true",default=False)
parser.add_argument("--saveCsv",action="store_false",default=True)
parser.add_argument("--Twiss",type=str,help="Which Twiss set to use, the smallest from Yngve ('y') or Eric's OpenXAL calculation,('o')")
parser.add_argument("--ampl",type=str,help="Which range of amplitudes to use, short (nominal - 10% less) or large (nominal to 70% less)")
args = parser.parse_args()

#Constants for running scripts
beamType    = "ESS"
energy      = 570 #[MeV]
graph       = False
NperBunch   = 10
nPulses     = 1e3
envXatBPM94 = 0
envYatBPM94 = 0
edges       = False
thick       = 0
PBIP        = False
physList    = "QGSP_BERT_EMZ" #"FTFP_BERT_EMZ" is alternate
dependence  = "RasterAmplitude"
rasterXAmplitude0 = 54.65 #[mm]
rasterYAmplitude0 = 18.37 #[mm]

#Twiss selection, 'y'=smallest from Yngve, 'o'= OpenXAL calculation
if args.Twiss == 'y':
  Twiss = [144.15027172522036,-8.184063058768368,0.3519001,88.04934327630778,-1.0382192928960423,0.3651098] #smallest from Yngve
elif args.Twiss == 'o':
  Twiss = [1006.80,-60.44,0.11315,129.72,-7.72,0.12155] #from my OpenXAL calculation
elif args.Twiss == 's':
  Twiss = [50,-10,0.5,30,-5,0.5]
elif args.Twiss == 'p':
  Twiss = [0.15,0,0.001,0.15,0,0.001]

#For the input amplitude range selection, 'short' or 'long'
amplRatio = rasterXAmplitude0 / rasterYAmplitude0
if args.ampl == 's':
  rXAmps = np.array([rasterXAmplitude0,49,50,51,52,53,54])
  #Make Y amplitude values based on ratio (necessary?)
  rYAmps = np.zeros(len(rXAmps))
  for i in range(len(rXAmps)):
    rYAmps[i] = rXAmps[i]/amplRatio
  rXRange = rXAmps.max() - rXAmps.min()
  print(rXRange)
  legloc = "center right"
  twTloc = 0.55
  pLloc = 0.95
elif args.ampl == 'l':
  #rXAmps = np.array([rasterXAmplitude0,15,25,30,35,40,45,50])
  rXAmps = np.array([rasterXAmplitude0,0.001,15,25,35])
  #Make Y amplitude values based on ratio (necessary?)
  rYAmps = np.zeros(len(rXAmps))
  for i in range(len(rXAmps)):
    rYAmps[i] = rXAmps[i]/amplRatio
  rXRange = rXAmps.max() - rXAmps.min()
  print(rXRange)
  legloc = "center left"
  twTloc = 0.85
  pLloc = 0.85
elif args.ampl == 'srcy': #Short Range Constant Y amplitude
  rXAmps = np.array([rasterXAmplitude0,49,50,51,52,53,54])
  rYAmps = rasterYAmplitude0*0.9 * np.ones(len(rXAmps)) #90% amplitude constant
  rXRange = rXAmps.max() - rXAmps.min()
  print(rXRange)
  legloc = "center right"
  twTloc = 0.55

POutBoxes = np.zeros(len(rXAmps))
VacPOutBoxes = np.zeros(len(rXAmps))

#Check if there is a CSV with the data already present. Speeds up plot modifications
csvPWD = "/scratch/ericdf/Scratch/PBWScatter/CSVs/"
name = "POutBoxRasterAmplDependence_bX{:.1f}m_{:.1f}mm".format(Twiss[0],rXRange)
if os.path.isfile(csvPWD+name+".csv"):
  print("Found data! Reading in!",name)
  with open(csvPWD+name+".csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    z=1
    for row in csv_reader:
      if line_count == 0: #Header line
        line_count += 1
      else:
        rXAmps[line_count-z] = row[0]
        rYAmps[line_count-z] = row[1]
        POutBoxes[line_count-z] = row[2]
        VacPOutBoxes[line_count-z] = row[3]
        line_count += 1
    csv_file.close()
else:
  print("Run simulations!",name)
  for i in range(len(rXAmps)):
    print("\nline",i,rXAmps[i],rYAmps[i])
    #Create Rastered Beam file, runARasterMaker checks if the CSV is already present
    rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(beamType,energy,graph,NperBunch,nPulses,envXatBPM94,envYatBPM94,edges,Twiss,rXAmps[i],rYAmps[i],dependence)
    #Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
    ImgPOutBox = runPBW(energy,rasterBeamFile,beamType,thick,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss,rXAmps[i],rYAmps[i],dependence)
    #Store % value for plotting
    POutBoxes[i] = ImgPOutBox
    #Run with 0.1 thickness, which is the value that triggers the Vacuum rectangular 'PBW' for reference
    noPBWPOutBox = runPBW(energy,rasterBeamFile,beamType,0.1,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss,rXAmps[i],rYAmps[i],dependence)
    VacPOutBoxes[i] = noPBWPOutBox

  #Save values for future quick use
  if args.saveCsv:
    with open(csvPWD+name+".csv",mode = 'w') as csv_file:
      csv_writer = csv.writer(csv_file,delimiter = ',')
      csv_writer.writerow(["Raster X Amplitude","Raster Y Amplitude","POutBox PBW","POutBox Vacuum"])
      for i in range(len(rXAmps)):
        csv_writer.writerow([rXAmps[i],rYAmps[i],POutBoxes[i],VacPOutBoxes[i]])
      csv_file.close()

for i in range(len(rXAmps)):
  print(rXAmps[i],rYAmps[i],POutBoxes[i],VacPOutBoxes[i])

#Plot for parameter search analysis
fs=14
picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
plt.close()
plt.clf()
plt.scatter(rXAmps,POutBoxes,marker='o',label=r"PBW % Outisde Box")
plt.scatter(rXAmps,VacPOutBoxes,marker='x',label=r"No PBW % Outside Box")
plt.xlabel("Horizontal Raster Amplitude [mm]",fontsize=fs)
plt.ylabel(r"% Outside 99% Box at Target",fontsize=fs)
plt.title("Rastered Beam Halo Growth on Target\nwith Beam Roundness at PBW",fontsize = fs+2)
#Set up texts to include with relevant info
xlim = plt.xlim()
ylim = plt.ylim()
plt.ylim([-0.7,ylim[1]])
ylim = plt.ylim()
plt.text(xlim[1]*0.98,ylim[0]+0.1,physList,fontsize = fs-5,color="k",horizontalalignment="right",verticalalignment="bottom")
plt.text(xlim[0]+0.2,ylim[1]*twTloc,"Beam Twiss at PBW:",fontsize=fs-4) #position Twiss print out depending on plot range
plt.text(xlim[0]+0.2,ylim[1]*(twTloc-0.07*1),r"$\epsilon_{Nx,Ny}$ = "+"{:.3f}, {:.3f}".format(Twiss[2],Twiss[5])+r"$_{[mm \cdot mrad]}$",fontsize=fs-4)
plt.text(xlim[0]+0.2,ylim[1]*(twTloc-0.07*2),r"$\beta_{x,y}$ = "+"{:.0f}, {:.0f}".format(Twiss[0], Twiss[3])+r"$_{[m]}$",fontsize=fs-4)
plt.text(xlim[0]+0.2,ylim[1]*(twTloc-0.07*3),r"$\alpha_{x,y}$ = "+"{:.1f}, {:.1f}".format(Twiss[1],Twiss[4]),fontsize=fs-4)
if args.ampl == 's' or args.ampl == 'l':
  plt.text(2*(xlim[1]-xlim[0])/7+xlim[0],ylim[1]*0.93,"Raster Amplitude H:V = 2.975",fontsize=fs-2)
plt.legend(loc=legloc)
plt.tight_layout()
dt = datetime.now()
plt.savefig(picPWD+name+dt.strftime("%H-%M-%S")+".png")

print(datetime.now()-origin)