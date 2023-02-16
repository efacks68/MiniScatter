import numpy as np
import matplotlib.pyplot as plt
import argparse
from math import pi, asin, sin

#Use Argument Parser to control
parser = argparse.ArgumentParser()
parser.add_argument("--s",action='store_true',default=False)
parser.add_argument("--ESS",action='store_true')
parser.add_argument("--pencil",action='store_true')
parser.add_argument("--g",action='store_true',default=False)
parser.add_argument("--Nb",        type=int,    default=10,    help="Number of macroparticles per beamlet. Default=10")
parser.add_argument("--nP",        type=float,  default=1e3,   help="Numper of beamlets in pulse. Default=1e3")
parser.add_argument("--rX",        type=float,  default=0,     help="X distance from beam axis [mm]. Default=0")
parser.add_argument("--rY",        type=float,  default=0,     help="Y distance from beam axis [mm]. Default=0")
parser.add_argument("--aX",        type=float,  default=54.65, help="RM X Amplitude [mm]. Default=54.65")
parser.add_argument("--aY",        type=float,  default=18.37, help="RM Y Amplitude [mm]. Default=18.37")
parser.add_argument("--edges",action="store_true")
parser.add_argument("--picFormat", type=str,   default="png",  choices=("png","svg","pdf"),help="Whic file format extension?")
args = parser.parse_args()
NperBunch = args.Nb
nPulses = args.nP

t_pulse = round(0.04 * 1/14 * 1e6) # mus
pulse_start = 10
t_end = (2* pulse_start + t_pulse) /1000# - 2#3 ms
time_length = round(t_end * nPulses ) #number of pulses, nominal = 2.86e3

t = np.linspace(0,t_end,time_length) #array of steps of length time_length
N_t = len(t) # number of time samples
n_tii  = 10 
#print("t= ",N_t, time_length,t_end)

##Raster Constants
#Raster Frequency
Fx = 39.953*1e3 #[kHz]39.550
Fy = 28.7051*1e3 #[kHz]28.7051
pb = 1
periodX = pb/Fx * np.ones(N_t) #[s] used for beamlet center calculation
periodY = pb/Fy * np.ones(N_t) #[s]
dt =  np.mean(np.diff(t)) #[s]
delta_t = np.linspace(0,dt,n_tii) #[s]

#Calculate Envelope Center Angle
envXatBPM94 = args.rX
envYatBPM94 = args.rY
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
RasterXAmpl = args.aX * np.ones(N_t) #[mm]
RasterYAmpl = args.aY * np.ones(N_t) #[mm]

x = np.zeros(N_t*n_tii)
y = np.zeros(N_t*n_tii)
time = np.zeros(N_t*n_tii)
print(x.shape)
i=0
for jj in range(N_t):
  for ii in range(n_tii):
    time[i] = t[jj] + delta_t[ii]

    x[i] = envXCenOff[jj] + 2 * RasterXAmpl[jj] / pi * asin(sin(2 * pi / periodX[jj] * time[i] )) #[mm]
    y[i] = envYCenOff[jj] + 2 * RasterYAmpl[jj] / pi * asin(sin(2 * pi / periodY[jj] * time[i] )) #[mm]
    i+=1
    #save total beamlet position
    #centroids[ii,0] = x
    #centroids[ii,1] = y


print(x.shape)
fig,ax = plt.subplots()
from matplotlib.patches import Ellipse
ax.add_patch(Ellipse((0,0),width=10.6/args.aX,height=3.83/args.aY,fill=False,edgecolor="red"))
ax.scatter(x/args.aX,y/args.aY,s=1)
ax.set_xlabel(r"Horizontal Centroid Deflection / a$_X$")
ax.set_ylabel(r"Vertical Centroid Deflection / a$_Y$")
ax.set_xticks([-1,-.5,0,.5,1])
ax.set_yticks([-1,-.5,0,.5,1])
ax.set_title("Lissajous Pattern")
ax.grid(which='major')
plt.tight_layout()
plt.savefig("lissajousPattern."+args.picFormat)
plt.close()
####
n=150
plt.plot(time*1e3,x/args.aX,linewidth=1)
plt.grid(which='major')
plt.gca().set_yticks([-1,-.5,0,.5,1])
plt.xlim([0,periodX[0]*4.4*n*1e3])
plt.xlabel(r"Time [$\mu$s]")
plt.ylabel(r"Horizontal Centroid Deflection / a$_X$")
plt.title("Waveform "+str(n)+" Periods")
plt.tight_layout()
plt.savefig("Sawtooth2.png")