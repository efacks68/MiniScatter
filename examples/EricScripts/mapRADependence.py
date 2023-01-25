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
parser.add_argument("--beamClass",type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss")
parser.add_argument("--twiss",   type=float,  nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
parser.add_argument("--t",       type=float,  default=0,     help="PBW Thickness [mm], 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al Xmm thick")
parser.add_argument("--energy",  type=float,  default=570,   help="Beam Energy [MeV]")
parser.add_argument("--Nb",      type=int,    default=10,    help="Number of macroparticles per beamlet")
parser.add_argument("--nP",      type=float,  default=1e2,   help="Numper of beamlets in pulse")
parser.add_argument("--rX",      type=float,  default=0,     help="X distance from beam axis [mm]")
parser.add_argument("--rY",      type=float,  default=0,     help="Y distance from beam axis [mm]")
parser.add_argument("--aX",      type=float,  default=54.65, help="RM X Amplitude [mm]")
parser.add_argument("--aY",      type=float,  default=18.37, help="RM Y Amplitude [mm]")
parser.add_argument("--failure", type=float,  default=0,     choices = range(0,5),  help="Which RM Failure case, 0-4.")
parser.add_argument("--magFails",type=int,    default=2,     choices = range(0,5),  help="Number of Raster Magnets that fail, 1-4.")
parser.add_argument("--xlim",    type=float,  default=450,   help="+/- value for horizontal axis of output rastered image [mm]")
parser.add_argument("--ylim",    type=float,  default=500,   help="+/- value for vertical axis of output rastered image [mm]")
parser.add_argument("--maxim",   type=float,  default=0  ,   help="Maximum current density value for output rastered imagem[uA/cm^2]")
parser.add_argument("--edges",   action="store_true",  help="Only populate edges of raster?")
parser.add_argument("--PBIP",    action="store_true",  default=False,   help="Is PBIP present?")
parser.add_argument("--noText",  action="store_true",  default=False,    help="Turns off printed text when called")
parser.add_argument("--savePics",action="store_true",  default=False,   help="Saves Rastered Image")
parser.add_argument("--noBox",   action="store_true",  default=False,   help="Turns off printed box when called")
parser.add_argument("--saveHist",action="store_true",  default=False,   help="Saves Histogram of proton density at target")
parser.add_argument("--saveRaster",action="store_true",default=False,   help="Saves plot of rastered beam")
parser.add_argument("--saveFits", action="store_true",  default=False,   help="Saves plots of Gaussian Fitting")

parser.add_argument("--ampl",   type=str,     default='map',help="Range of amplitudes to use: short(nominal-10% less) or large(nominal-70% less)")
parser.add_argument("--eX",     type=int,     default=50,    help="End ampl X")
parser.add_argument("--eY",     type=int,     default=20,    help="End ampl Y")
parser.add_argument("--startX", type=int,     default=0,     help="Start ampl for X")
parser.add_argument("--startY", type=int,     default=0,     help="Start ampl for Y")
parser.add_argument("--NstepX",  type=int,     default=6,    help="N steps for X")
parser.add_argument("--NstepY",  type=int,     default=6,    help="N steps for Y")
args = parser.parse_args()
print(args)

#Constants for running scripts
physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ"
options     = {'noText':args.noText, 'noBox':args.noBox, 'wide':True, 'physList':physList, 'dependence':"RA",
                            'xlim':args.xlim, 'ylim':args.ylim, 'maxim':args.maxim, 'saveHist':args.saveHist,
                            'PBIP':args.PBIP, 'beamClass':args.beamClass, 'Nb':args.Nb, 'failure':args.failure,
                            'magFails':args.magFails, 'saveRaster':args.saveRaster, 'saveFits':args.saveFits }

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
elif os.uname()[1] == "mbef-XPS-13-9300":
    csvPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/"
else: print("Help! Unknown build directory!, scatterPBW.py l 61")

#For the input amplitude range selection, 'short' or 'long'
amplRatio = (args.aX + 0.001) / (args.aY + 0.001) #so no /zero
defaultRMAmplX = 54.65
defaultRMAmplY = 18.37
if args.ampl == 's':
    #rXAmps = np.array([args.aX,49,50,51,52,53,54])
    rXAmps = np.array([args.aX,49,51,53,57])
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
    #rXAmps = np.array([args.aX,15,25,30,35,40,45,50])
    rXAmps = np.array([args.aX,0.001,15,25,35])
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
    rXAmps = np.array([args.aX,49,50,51,52,53,54])
    rYAmps = args.aY*0.9 * np.ones(len(rXAmps)) #90% amplitude constant
    rXRange = rXAmps.max() - rXAmps.min()
    print(rXRange)
    legloc = "center right"
    twTloc = 0.55
elif args.ampl =="map":
    rXAmps = np.linspace(args.startX,args.eX,args.NstepX)
    rYAmps = np.linspace(args.startY,args.eY,args.NstepY)
    rXRange = args.eX-args.startX
    #print(start,end,step,"\n",rXAmps,"\n",rYAmps)
    print("there are ",args.NstepX*args.NstepY,"points to plot. Expect that number of minutes.")


POutBoxes = np.zeros([len(rYAmps),len(rXAmps)])
Imaxes = np.zeros([len(rYAmps),len(rXAmps)])
coreMeans = np.zeros([len(rYAmps),len(rXAmps)])
#VacPOutBoxes = np.zeros(len(rXAmps))

#Check if there is a CSV with the data already present. Speeds up plot modifications
mapCsvPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/rAmplDependence/2DMap/"
name = "RasterAmplDependence_POutBox,Imax_bX{:.1f}m_R{:.1f},{:.1f}mm".format(Twiss[0],rXRange,args.eY-args.startY)
if os.path.isfile(mapCsvPWD+name+".csv"):
    print("Found data! Reading in!",name)
    #from plotFit import numLines
    #nLines = numLines(csvPWD+name)
    #print(nLines)
    with open(mapCsvPWD+name+".csv") as csv_file:
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
    print("Run simulations for RA ratio",amplRatio."\t",name)
    for i in range(len(rXAmps)):
        for j in range(len(rYAmps)):
            print("\nline [",i,",",j,"]",rXAmps[i],rYAmps[j])
            #Create Rastered Beam file, runARasterMaker checks if the CSV is already present
            rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args.energy,args.Nb,args.nP,args.rX,args.rY,args.edges,Twiss,rXAmps[i],rYAmps[j],csvPWD,options)
            #Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
            POutBoxes[j][i], Imaxes[j][i], coreMeans[j][i] = runPBW(args.energy,rasterBeamFile,args.t,beamXAngle,beamYAngle,args.savePics,Twiss,rXAmps[i],rYAmps[j],options)
            #Store % value for plotting
            #Run with 0.1 thickness, which is the value that triggers the Vacuum rectangular 'PBW' for reference
            #noPBWPOutBox = runPBW(energy,rasterBeamFile,beamType,0.1,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss,rXAmps[i],rYAmps[i],dependence)
            #VacPOutBoxes[i] = noPBWPOutBox
            print(datetime.now().strftime("%H-%M-%S"))
        print("time elapsed",datetime.now()-origin)

       #Save values for future quick use
    #if args.saveCsv:
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
from math import floor,ceil
xlim1 = ax1.get_xlim()
ylim1 = ax1.get_ylim()
#22Jan-you are adjusting the nominal markers, then runnning several different maps
#try to find an optimal ampl.
#also, you have to adjust the original values in the Gaussian fitting as it keeps throwing 
#  the "lower/upper bounds outside current parameter value." error which is about the initial
#  value being outside the SetParamterLimit bounds

ax1.hlines(defaultRMAmplX,floor(defaultRMAmplX*0.5),ceil(defaultRMAmplX*0.5),color='m')
ax1.vlines(defaultRMAmplX*0.5,defaultRMAmplY*xlim1[1],ceil(defaultRMAmplY),color='m')
ax1.text(defaultRMAmplX*0.5+1,defaultRMAmplY+1,"Nominal Y, Half Nominal X")
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