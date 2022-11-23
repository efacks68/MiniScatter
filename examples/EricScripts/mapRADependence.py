#Script for searching Raster Magnet Amplitude space
# making 2D map of % outside box, with input from Carl Lindstrom.

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
parser.add_argument("--eX",type=int,default=50,help="End ampl X")
parser.add_argument("--eY",type=int,default=20,help="End ampl Y")
parser.add_argument("--startX",type=int,default=0,help="Start ampl for X")
parser.add_argument("--startY",type=int,default=0,help="Start ampl for Y")
parser.add_argument("--stepX",type=int,default=11,help="N steps for X")
parser.add_argument("--stepY",type=int,default=11,help="N steps for Y")
args = parser.parse_args()
print(args)

#Constants for running scripts
beamType    = "ESS"
energy      = 570 #[MeV]
graph       = False
NperBunch   = 10
nPulses     = 1e2
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
elif args.Twiss == 's': #smaller beam than Yngve's smallest
  Twiss = [50,-10,0.5,30,-5,0.5]
elif args.Twiss == 'p': #"pencil" beam of 0 emittance
  Twiss = [0.15,0,0.0001,0.15,0,0.0001]

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
elif args.ampl =="map":
  rXAmps = np.linspace(args.startX,args.eX,args.stepX)
  rYAmps = np.linspace(args.startY,args.eY,args.stepY)
  rXRange = args.eX-args.startX
  #print(start,end,step,"\n",rXAmps,"\n",rYAmps)
  print("there are ",args.stepX*args.stepY,"points to plot. Expect that number of minutes.")


POutBoxes = np.zeros([len(rYAmps),len(rXAmps)])
Imaxes = np.zeros([len(rYAmps),len(rXAmps)])
coreMeans = np.zeros([len(rYAmps),len(rXAmps)])

#VacPOutBoxes = np.zeros(len(rXAmps))

#Check if there is a CSV with the data already present. Speeds up plot modifications
csvPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/rAmplDependence/2DMap/"
name = "RasterAmplDependence_POutBox,Imax_bX{:.1f}m_R{:.1f},{:.1f}mm".format(Twiss[0],rXRange,args.eY-args.startY)
if os.path.isfile(csvPWD+name+".csv"):
  print("Found data! Reading in!",name)
  #from plotFit import numLines
  #nLines = numLines(csvPWD+name)
  #print(nLines)
  with open(csvPWD+name+".csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    j = 0
    k = 0
    z=1
    for row in csv_reader:
      if i == 0 and j ==0: #Header line
        i += 1
        j += 1
      else:
        rXAmps[i-z] = row[0]
        rYAmps[j-z] = row[1]
        POutBoxes[j-z][i-z] = row[2]
        Imaxes[j-z][i-z] = row[3]
        coreMeans[j-z][i-z] = row[4]
        #VacPOutBoxes[line_count-z] = row[3]
        #print(i-z,j-z,rXAmps[i-z],rYAmps[j-z],POutBoxes[i-z][j-z])
        j += 1
        if j == len(rXAmps)+1:
          i += 1
          j = 0 + z
    csv_file.close()
else:
  print("Run simulations!",name)
  for i in range(len(rXAmps)):
    for j in range(len(rYAmps)):
      print("\nline [",i,",",j,"]",rXAmps[i],rYAmps[j])
      #Create Rastered Beam file, runARasterMaker checks if the CSV is already present
      rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(beamType,energy,graph,NperBunch,nPulses,envXatBPM94,envYatBPM94,edges,Twiss,rXAmps[i],rYAmps[j],dependence,csvPWD)
      #Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
      POutBoxes[j][i], Imaxes[j][i], coreMeans[j][i] = runPBW(energy,rasterBeamFile,beamType,thick,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss,rXAmps[i],rYAmps[j],dependence)
      #Store % value for plotting
      #Run with 0.1 thickness, which is the value that triggers the Vacuum rectangular 'PBW' for reference
      #noPBWPOutBox = runPBW(energy,rasterBeamFile,beamType,0.1,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss,rXAmps[i],rYAmps[i],dependence)
      #VacPOutBoxes[i] = noPBWPOutBox
    print("time elapsed",datetime.now()-origin)

       #Save values for future quick use
  if args.saveCsv:
    print("Writing CSV")
    with open(csvPWD+name+".csv",mode = 'w') as csv_file:
      csv_writer = csv.writer(csv_file,delimiter = ',')
      csv_writer.writerow(["Raster X Amplitude","Raster Y Amplitude","POutBox PBW","IMaxes [uA/mm2]","CoreIMeans [uA/mm2]"])
      for i in range(len(rXAmps)):
        for j in range(len(rYAmps)):
          csv_writer.writerow([rXAmps[i],rYAmps[j],POutBoxes[j][i],Imaxes[j][i],coreMeans[j][i]])#,VacPOutBoxes[i]])
      csv_file.close()
      print("CSV written",csvPWD+name)
  print("time elapsed",datetime.now()-origin)

for i in range(len(rXAmps)):
  for j in range(len(rYAmps)):
    print(rXAmps[i],rYAmps[j],POutBoxes[j][i],Imaxes[j][i])

#Plot for parameter search analysis
fs=14
minim = POutBoxes.min()+0.01
maxim = POutBoxes.max()*1.1
plotRange = maxim-minim
print(minim,maxim,plotRange)
picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/rAmplDependence/2DMap/"
plt.close()
plt.clf()
X,Y = np.meshgrid(rXAmps,rYAmps)
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(15,6))
plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
from matplotlib.colors import LogNorm
c = ax1.pcolor(X,Y,POutBoxes,norm=LogNorm(vmin=minim, vmax=maxim),shading='auto',cmap='viridis')#,label=r"PBW % Outisde Box")
ax1.set_xlabel("Horizontal Raster Amplitude [mm]",fontsize=fs)
ax1.set_ylabel("Vertical Raster Amplitude [mm]",fontsize=fs)
ax1.set_title("Rastered Beam % Outside Box on Target\nwith Raster Amplitude at PBW",fontsize = fs+2)
cbarVals  = [minim,minim+plotRange*0.1,minim+plotRange*0.4,minim+plotRange*0.7,maxim] #make array for color bar values
cbarLabels = ["{:.1f}".format(cbarVals[0]),"{:.1f}".format(cbarVals[1]),"{:.1f}".format(cbarVals[2]),
              "{:.1f}".format(cbarVals[3]),"{:.1f}".format(cbarVals[4])]#,"{:.1f}".format(cbarVals[5])] #make labels of Value
cbarLabel = "% Outside Box"
cbar = fig.colorbar(c, ax=ax1,pad=0.01,ticks=cbarVals)
cbar.set_label(cbarLabel,labelpad=2,fontsize=fs-2)
#cbar.set_ticks(cbarVals)
cbar.set_ticklabels(cbarLabels)

##Max I Plot
minimMax = Imaxes.min()*0.9999
maximMax = Imaxes.max()*1.0001
plotRangeMax = maximMax-minimMax
print(Imaxes.min(),Imaxes.max())
d = ax2.pcolor(X,Y,Imaxes,norm=LogNorm(vmin=minimMax, vmax=maximMax),shading='auto',cmap='viridis')
ax2.set_xlabel("Horizontal Raster Amplitude [mm]",fontsize=fs)
ax2.set_ylabel("Vertical Raster Amplitude [mm]",fontsize=fs)
ax2.set_title("Peak Current Density on Target\nwith Raster Amplitude at PBW",fontsize = fs+2)
cbarVals2  = [minimMax,minimMax+plotRangeMax*0.25,minimMax+plotRangeMax*0.5,minimMax+plotRangeMax*0.75,maximMax] #make array for color bar values
cbarLabels2 = ["{:.1f}".format(cbarVals2[0]),"{:.1f}".format(cbarVals2[1]),"{:.1f}".format(cbarVals2[2]),
              "{:.1f}".format(cbarVals2[3]),"{:.1f}".format(cbarVals2[4])]
cbarLabel2 = r"Current Density [$\mu$A/cm$^2$]"
cbar2 = fig.colorbar(d, ax=ax2,pad=0.01,ticks=cbarVals2)
cbar2.set_label(cbarLabel2,labelpad=2,fontsize=fs-2)
#cbar2.set_ticks(cbarVals2)
cbar2.set_ticklabels(cbarLabels2)
##Set up texts to include with relevant info
#xlim = plt.xlim()
#ylim = plt.ylim()
#plt.ylim([-0.7,ylim[1]])
#ylim = plt.ylim()
#plt.text(xlim[1]*0.98,ylim[0]+0.1,physList,fontsize = fs-5,color="k",horizontalalignment="right",verticalalignment="bottom")
#plt.text(xlim[0]+0.2,ylim[1]*twTloc,"Beam Twiss at PBW:",fontsize=fs-4) #position Twiss print out depending on plot range
#plt.text(xlim[0]+0.2,ylim[1]*(twTloc-0.07*1),r"$\epsilon_{Nx,Ny}$ = "+"{:.3f}, {:.3f}".format(Twiss[2],Twiss[5])+r"$_{[mm \cdot mrad]}$",fontsize=fs-4)
#plt.text(xlim[0]+0.2,ylim[1]*(twTloc-0.07*2),r"$\beta_{x,y}$ = "+"{:.0f}, {:.0f}".format(Twiss[0], Twiss[3])+r"$_{[m]}$",fontsize=fs-4)
#plt.text(xlim[0]+0.2,ylim[1]*(twTloc-0.07*3),r"$\alpha_{x,y}$ = "+"{:.1f}, {:.1f}".format(Twiss[1],Twiss[4]),fontsize=fs-4)
#if args.ampl == 's' or args.ampl == 'l':
#  plt.text(2*(xlim[1]-xlim[0])/7+xlim[0],ylim[1]*0.93,"Raster Amplitude H:V = 2.975",fontsize=fs-2)
#plt.legend(loc=legloc)
plt.tight_layout()
dt = datetime.now()
plt.savefig(picPWD+name+dt.strftime("%H-%M-%S")+".png")
print("saved",picPWD+name+dt.strftime("%H-%M-%S")+".png")

print(datetime.now()-origin)