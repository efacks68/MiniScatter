import numpy as np
from math import pi, asin, sin
import argparse
from datetime import datetime
start = datetime.now()
print(start)
#distance 1e0 unit in this script is mm, so 1e3=1meter

#Use Argument Parser to control
parser = argparse.ArgumentParser()
parser.add_argument("--twiss",   type=float,  nargs=6,       help="Twiss parameters in form: NemtX,BetaX,AlphX,NemtY,BetaY,AlphY")
parser.add_argument("--s",       action='store_true',        default=False)
parser.add_argument("--beamClass",type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', 'pencil', or 'twiss'? 'twiss' expects '--twiss' argument with Twiss in MiniScatter format")
parser.add_argument("--g",       action='store_true',        default=False)
parser.add_argument("--Nb",      type=int,    default=10,    help="Number of macroparticles per beamlet")
parser.add_argument("--nP",      type=float,  default=1e3,   help="Numper of beamlets in pulse")
parser.add_argument("--rX",      type=float,  default=0,     help="X distance from beam axis")
parser.add_argument("--rY",      type=float,  default=0,     help="Y distance from beam axis")
parser.add_argument("--aX",      type=float,  default=54.65, help="RM X Amplitude")
parser.add_argument("--aY",      type=float,  default=18.37, help="RM Y Amplitude")
parser.add_argument("--edges",   action="store_true")
args = parser.parse_args()
NperBunch = args.N
nPulses = args.a
um = 1e-6
mm = 1e-3

# Twiss= [NemtX,BetaX,AlphX,NemtY,BetaY,AlphY]
if args.beamClass == 'Yngve': #smallest from Yngve
  Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
elif args.beamClass == 'ESS': #from my OpenXAL calculation
  Twiss = [0.11315,1006.80,-60.44,0.12155,129.72,-7.72]
elif args.beamClass == 'pencil': #"pencil" beam of 0 emittance
  Twiss = [0.0001,0.15,0,0.0001,0.15,0]
if args.twiss:
  Twiss = [args.twiss[0],args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5]]
  options['beamClass'] = "Twiss"

#Twiss 
nemtX = Twiss[0] #[mm-mrad]
betaX = Twiss[1] #[m]
alphX = Twiss[2]
nemtY = Twiss[3] #[mm-mrad]
betaY = Twiss[4] #[m]
alphY = Twiss[5]
#Calculate Geometric Emittance for Covariance Matrix
partA = 938.27209 #[MeV/c2]
partZ = 1
gamma_rel = 1 + energy / partA #from PrintTwissParameters
beta_rel = np.sqrt(gamma_rel * gamma_rel - 1 ) / gamma_rel
gemtX = nemtX*um / (beta_rel * gamma_rel) #[m]
gemtY = nemtY*um / (beta_rel * gamma_rel) #[m]

t_pulse = round(0.04 * 1/14 * 1e6) # mus
pulse_start = 10
t_end = (2* pulse_start + t_pulse) /1000# - 2#3 ms
time_length = round(t_end * nPulses ) #number of pulses, nominal = 2.86e3

t = np.linspace(0,t_end,time_length) #array of steps of length time_length
N_t = len(t) # number of time samples
n_tii  = 10 
#print("t= ",N_t, time_length,t_end)

totX = np.zeros([N_t*n_tii*NperBunch,2])
totY = np.zeros([N_t*n_tii*NperBunch,2])
#print(len(totX))

##Raster Constants
#Raster Frequency
Fx = 39.953*1e3 #[kHz]
Fy = 28.7051*1e3 #[kHz]
pb = 1
periodX = pb/Fx * np.ones(N_t) #[s] used for beamlet center calculation
periodY = pb/Fy * np.ones(N_t) #[s]
dt =  np.mean(np.diff(t)) #[s]
delta_t = np.linspace(0,dt,n_tii) #[s]
z = -10 #[mm] for z location to generate protons at in MiniScatter

#Assumes Protons
covX = gemtX/mm * np.asarray([[betaX/mm,-alphX],[-alphX,(1+alphX**2)/(betaX/mm)]]) #[mm]
covY = gemtY/mm * np.asarray([[betaY/mm,-alphY],[-alphY,(1+alphY**2)/(betaY/mm)]]) #[mm]

#Calculate Envelope Center Angle
dBPM93to94 = 3031 #[mm] from OpenXAL(?)
dBPM93toPBW = 16822 #[mm] PBW 4400mm upstream of Target: https://gitlab.esss.lu.se/ess-bp/ess-lattice/-/blob/HEBT_RASTER_V29/9.0_HEBT/Beam_Physics/lattice.dat
dBPM93toTarg = 21222 #[mm] from Synoptic Viewer https://confluence.esss.lu.se/pages/viewpage.action?pageId=222397499
dPBWtoTarg = 4400 #[mm] from lattice and Synoptic
envXAngle = envXatBPM94 / dBPM93to94 #x' = distance from beamline axis at BPM94, assume Cross Over at BPM93 / distance BPM 93 to 94
envYAngle = envYatBPM94 / dBPM93to94 #y' not radians, as per Kyrre 2.11.22
beamletXAngle = 0 #default
beamletYAngle = 0 #default

#Envelope Center Offset, i.e. raster centre on BEW
envXCenterOffset = envXAngle * dBPM93toPBW #[mm] Takes into account angular drift since Cross Over
envYCenterOffset =  envYAngle * dBPM93toPBW #[mm]
envXCenOff = envXCenterOffset * np.ones(N_t) #[mm]
envYCenOff = envYCenterOffset * np.ones(N_t) #[mm]

#Raster Amplitude: ax0,ay0 here produces the beamlet displacement AT the Target, as per Cyrille
#rasterXAmplitude0 = 54.65 #[mm] #default
#rasterYAmplitude0 = 18.37 #[mm]
#because generating particles just before PBW, must scale a0 by dPBWtoTarg * envAngle = (1- dPBWtoTarg / dBPM93toTarg)
#amplScale = 1 - dPBWtoTarg / dBPM93toTarg #double check you account for Z before PBW in MiniScatter! beam production plane in GEANT, not exact PBW center!
amplScale = 1 #Cyrille said the RM Amplitude is already scaled
sRasterXAmpl = rasterXAmplitude0 * amplScale * np.ones(N_t) #[mm]
sRasterYAmpl = rasterYAmplitude0 * amplScale * np.ones(N_t) #[mm]

#For weighting edges case (--edges argument)
Right = 50
Top = 17

i=0
#j=0
#k=0

#Pick name based on beam
if options['beamClass'] == "ESS":
  name = "PBW_{:.0f}MeV_beta{:.0f},{:.0f}m_RMamp{:.0f},{:.0f}mm_N{:.1e}_NpB{:.0f}_NPls{:.0e}".format(energy,betaX,betaY,rasterXAmplitude0,rasterYAmplitude0,len(totX[:,0]),NperBunch,nPulses)#+dt.strftime("%H-%M-%S")
elif options['beamClass'] == "pencil":
  name = "PBW_{:.0f}MeV_pencilBeam_RMampl{:.0f},{:.0f}mm_N{:.1e}_NpB{:.0f}_NPls{:.1e}".format(energy,rasterXAmplitude0,rasterYAmplitude0,len(totX[:,0]),NperBunch,nPulses)
elif options['beamClass'] == "Twiss":
  name = "PBW_{:.0f}MeV_eX{:.2f}um,eY{:.2f}um_bX{:.0f}m,bY{:.0f}m_aX{:.0f},aY{:.0f}_N{:.1e}".format(energy,nemtX/um,nemtY/um,betaX,betaY,alphX,alphY,len(totX[:,0]))
if envXatBPM94 != 0:
  name = name + "_X{:.0f}mrad".format(envXAngle*1e3)
if envYatBPM94 != 0:
  name = name + "_Y{:.0f}mrad".format(envYAngle*1e3)

#if file found, don't make again!
#picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
csvPWD = "/scratch2/ericdf/PBWScatter/CSVs/" #put all CSVs in Scratch to save my disk space!

print(name)
centroids = np.zeros([N_t*n_tii,2])
for jj in range(N_t):
  for ii in range(n_tii):
    tjj_ii = t[jj] + delta_t[ii]

    #Calculate the Raster Magnet contribution to Beamlet Center Location relative to beamline axis, as per Cyrille Thomas
    beamletX = envXCenOff[jj] + 2 * sRasterXAmpl[jj] / pi * asin(sin(2 * pi / periodX[jj] * tjj_ii )) #[mm]
    beamletY = envYCenOff[jj] + 2 * sRasterYAmpl[jj] / pi * asin(sin(2 * pi / periodY[jj] * tjj_ii )) #[mm]

    #Calculate the total Beamlet Angle = Envelope Center Angle + the angle given to each beamlet by the Raster Magnets
    beamletXAngle = beamletX / dBPM93toTarg + envXAngle #[mm]
    beamletYAngle = beamletY / dBPM93toTarg + envYAngle #[mm]

    #save total beamlet position
    centroids[i,0] = beamletX
    centroids[i,1] = beamletY
    NperBunch = NperBunch

    #In case of weighting edges of raster 
    if edges:
      if beamletX > -Right and beamletX < Right and beamletY < Top and beamletY > -Top: #set weight depending on position
        NperBunch = 5 #decrease center NperBunch
    #    j += 1
      else:
        NperBunch = NperBunch #Edges get full NperBunch
    #    k += 1

    #Generate beamlet distributions
    rng = np.random.default_rng()
    ptsX = rng.multivariate_normal([beamletX,beamletXAngle],covX,size = NperBunch) #mean is [pos,ang]!
    ptsY = rng.multivariate_normal([beamletY,beamletYAngle],covY,size = NperBunch)
    
    for k in range(NperBunch): #put this beamlet into total. Could just be written, figure that out later.
      totX[NperBunch*i+k,0] = ptsX[k,0]
      totX[NperBunch*i+k,1] = ptsX[k,1]
      totY[NperBunch*i+k,0] = ptsY[k,0]
      totY[NperBunch*i+k,1] = ptsY[k,1]
    i +=1
#print(i,j,k)

#Check on output parameters
print("Centroid X max: {:.2f}mm; Particle X max: {:.2f}mm".format(np.max(centroids[:,0]),np.max(totX[:,0])),"; Shape:",np.shape(totX))
#Remove 0,0 particles, should be none except for weighted edges case
nonzero = np.not_equal(totX[:,0],0) #be careful!
totX = totX[nonzero]
totY = totY[nonzero]
#print(np.shape(totX))
#print(name)

#Check on output parameters
print("Centroid X max: {:.2f}mm; Particle X max: {:.2f}mm".format(np.max(centroids[:,0]),np.max(totX[:,0])),"; Shape:",np.shape(totX))
#Remove 0,0 particles, should be none except for weighted edges case
nonzero = np.not_equal(totX[:,0],0) #be careful!
totX = totX[nonzero]
totY = totY[nonzero]
#print(np.shape(totX))
#print(name)

outname = picPWD+name
if args.s:
  import csv
  with open(outname+".csv",mode = 'w',newline=None) as part_file:
    part_writer = csv.writer(part_file,delimiter = ',')
    for i in range(len(totX)):
      part_writer.writerow(["proton", totX[i,0], totX[i,1], totY[i,0], totY[i,1], z, energy])
  part_file.close()

finish = datetime.now()
print(name, "; started at:",start.strftime("%H-%M-%S"), "finished in: ",finish-start)

if args.g:
  import matplotlib.pyplot as plt
  #found the below method: https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
  import mpl_scatter_density
  fs=14
  plt.close()
  fig = plt.figure()
  s1 = fig.add_subplot(1,1,1,projection="scatter_density")
  x=totX[:,0]
  y=totY[:,0]
  density = s1.scatter_density(x,y,cmap='jet')
  fig.colorbar(density,label=r"Protons/mm^2")
  xlim = plt.xlim()
  ylim = plt.ylim()
  plt.text(xlim[0]*0.95,ylim[0]+Twloc,"Beam Twiss at PBW:",fontsize=fs-4,c='w') #position Twiss print out depending on plot range
  plt.text(xlim[0]*0.95,ylim[0]+Twloc-(a*1),r"$\epsilon_{Nx,Ny}$ = "+"{:.3f}, {:.3f}".format(nemtX/mm,nemtY/mm)+r"$_{[mm \cdot mrad]}$",c='w',fontsize=fs-4)
  plt.text(xlim[0]*0.95,ylim[0]+Twloc-(a*2),r"$\beta_{x,y}$ = "+"{:.0f}, {:.0f}".format(betaX/mm, betaY/mm)+r"$_{[mm]}$",c='w',fontsize=fs-4)
  plt.text(xlim[0]*0.95,ylim[0]+Twloc-(a*3),r"$\alpha_{x,y}$ = "+"{:.1f}, {:.1f}".format(alphX,alphY),c='w',fontsize=fs-4)
  plt.text(xlim[0]*0.95,ylim[0]+Twloc-(a*4),r"$\sigma_{x,y}$ = "+"{:.1f}, {:.1f}".format(np.sqrt(betaX*gemtX)/mm,np.sqrt(betaY*gemtY)/mm)+r"$_{[mm]}$",c='w',fontsize=fs-4)
  s1.set_xlabel("X [mm]")
  s1.set_ylabel("Y [mm]")
  s1.set_title("Rastered Beam Number Density\n{:.1e} protons {:.2f}ms".format(len(totX),time_length*1e-3))
  plt.savefig(outname+".png")
  print(outname+".png")  plt.close()
