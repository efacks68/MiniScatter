#Script for searching Twiss values at PBW parameter space,
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
physList    = "QGSP_BERT_EMZ"
dependence  = "Twiss"
rasterXAmplitude0 = 54.65 #[mm]
rasterYAmplitude0 = 18.37 #[mm]

#Read in CSV file Twiss to simulate
Twiss = np.zeros((17,6)) #known length of Twiss file.
with open('beam_sizes18Aug.csv') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  line_count = 0
  for row in csv_reader:
    if line_count == 0:
      line_count += 1
    else:
      #File has the beam size at column 4 and 7, don't read those values in
      Twiss[line_count-1] = [row[0],row[1],row[2],row[4],row[5],row[6]]
      line_count += 1
  print("Processed ",line_count,"lines")
csv_file.close()

n=line_count-1
POutBoxes = np.zeros(n)
VacPOutBoxes = np.zeros(n)
BeamSRatio = np.zeros(n)

#Check if there is a CSV with the data already present. Speeds up plot modifications
csvPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/TwissDependence/"
if os.path.isfile(csvPWD+"POutBoxTwissDependence.csv"):
  print("Found data!")
  with open(csvPWD+'POutBoxTwissDependence.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    z=1
    for row in csv_reader:
      if line_count == 0: #Header line
        line_count += 1
      else:
        BeamSRatio[line_count-z] = row[0]
        POutBoxes[line_count-z] = row[1]
        VacPOutBoxes[line_count-z] = row[2]
        line_count += 1
    csv_file.close()
else:
  print("Run simulations!")
  for i in range(line_count-1):
    if Twiss[i][0] == 0.0:
      continue
    print("\nline",i,Twiss[i])
    #Create Rastered Beam file, runARasterMaker checks if the CSV is already present
    rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(beamType,energy,graph,NperBunch,nPulses,envXatBPM94,envYatBPM94,edges,Twiss[i],rasterXAmplitude0,rasterYAmplitude0,dependence)
    #Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
    ImgPOutBox = runPBW(energy,rasterBeamFile,beamType,thick,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss[i],rasterXAmplitude0,rasterYAmplitude0,dependence)
    #Store % value for plotting
    POutBoxes[i] = ImgPOutBox
    #Run with 0.1 thickness, which is the value that triggers the Vacuum rectangular 'PBW' for reference
    noPBWPOutBox = runPBW(energy,rasterBeamFile,beamType,0.1,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss[i],rasterXAmplitude0,rasterYAmplitude0,dependence)
    VacPOutBoxes[i] = noPBWPOutBox

    #Make Beam Size Ratio array
    BeamSRatio[i] = Twiss[i][0]/Twiss[i][3]
  
  #Save values for future quick use
  if args.saveCsv:
    with open(csvPWD+'POutBoxTwissDependence.csv',mode = 'w') as csv_file:
      csv_writer = csv.writer(csv_file,delimiter = ',')
      csv_writer.writerow(["Beam Size Ratio","POutBox PBW","POutBox Vacuum"])
      for i in range(len(BeamSRatio)):
        csv_writer.writerow([BeamSRatio[i],POutBoxes[i],VacPOutBoxes[i]])
      csv_file.close()


for i in range(len(BeamSRatio)):
  print(BeamSRatio[i],POutBoxes[i],VacPOutBoxes[i])

#Plot for parameter search analysis
fs=14
picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
plt.close()
plt.clf()
plt.scatter(BeamSRatio,POutBoxes,c="green",label=r"PBW % Outisde Box")
plt.scatter(BeamSRatio,VacPOutBoxes,c="cyan",label=r"No PBW % Outside Box")
plt.xlabel(r"Beam Size Ratio at PBW, $\beta_X / \beta_Y$",fontsize=fs)
plt.ylabel(r"% Outside 99% Box at Target",fontsize=fs)
plt.title("Rastered Beam Halo Growth on Target\nwith Beam Roundness",fontsize = fs+2)
#Set up texts to include with relevant info
xlim = plt.xlim()
ylim = plt.ylim()
plt.xlim([0.25,xlim[1]+0.25])
plt.ylim([-1,ylim[1]+.75])
plt.text(xlim[1]*0.7,ylim[0]+0.1,physList,fontsize = fs-4,color="k")
plt.text(0.35,ylim[1]*0.7,r"$\frac{\beta_x}{\beta_y}$",fontsize=fs+2)
plt.text(0.65,ylim[1]*0.7,r"$\frac{410}{472} _{[\rm m]}$",fontsize=fs)
plt.text(xlim[1]-0.7,ylim[1]*0.5,r"$\frac{1006}{129} _{[\rm m]}$",fontsize=fs)
plt.text(xlim[1]*0.75,ylim[1]*0.95,r"$\epsilon_{Nx} = 0.352_{[\rm mm \cdot \rm mrad]}$",fontsize=fs-4)
plt.text(xlim[1]*0.75,ylim[1]*0.9,r"$\epsilon_{Ny} = 0.365_{[\rm mm \cdot \rm mrad]}$",fontsize=fs-4)
plt.tight_layout()
plt.legend(loc='upper center')
dt = datetime.now()
plt.savefig(picPWD+"POutsideBoxTwissDependence"+dt.strftime("%H-%M-%S")+".png")

print(datetime.now()-origin)