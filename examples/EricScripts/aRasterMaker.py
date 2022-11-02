import numpy as np
from beamConstants import *
import matplotlib.pyplot as plt
import csv, math,sys
import argparse
from datetime import datetime
start = datetime.now()
print(start)
#distance 1e0 unit in this script is mm, so 1e3=1meter

#Use Argument Parser to control
parser = argparse.ArgumentParser()
parser.add_argument("--s",action='store_true',default=False)
parser.add_argument("--ESS",action='store_true')
parser.add_argument("--pencil",action='store_true')
parser.add_argument("--g",action='store_true',default=False)
parser.add_argument("--N",type=int,default=100)
parser.add_argument("--a",type=float,default=1e3)
parser.add_argument("--rX",type=float,default=0)
parser.add_argument("--rY",type=float,default=0)
parser.add_argument("--edges",action="store_true")
args = parser.parse_args()

t_pulse = round(0.04 * 1/14 * 1e6) # mus
pulse_start = 10
t_end = (2* pulse_start + t_pulse) /1000# - 2#3 ms
time_length = round(t_end * args.a ) #number of pulses, nominal = 2.86e3

t = np.linspace(0,t_end,time_length) #array of steps of length time_length
N_t = len(t) # number of time samples
n_tii  = 10 
#print("t= ",N_t, time_length,t_end)

totX = np.zeros([N_t*n_tii*args.N,2])
totY = np.zeros([N_t*n_tii*args.N,2])
print(len(totX))

#Raster Constants
#Raster Frequency
Fx = 39.953*1e3 #[kHz]
Fy = 28.7051*1e3 #[kHz]
pb = 1
px = pb/Fx * np.ones(N_t) #[s] used for beamlet center calculation
py = pb/Fy * np.ones(N_t) #[s]
dt =  np.mean(np.diff(t)) #[s]
delta_t = np.linspace(0,dt,n_tii) #[s]
z = -10 #[mm] for z location to generate protons at in MiniScatter
picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"

if args.ESS: #assume 570MeV
  covx = gemtx/mm * np.asarray([[betax/mm,-alphx],[-alphx,(1+alphx**2)/(betax/mm)]]) #[mm]
  covy = gemty/mm * np.asarray([[betay/mm,-alphy],[-alphy,(1+alphy**2)/(betay/mm)]]) #[mm]

if args.pencil: #assume 570MeV
  covx = 7.945383034919336e-05/mm * np.asarray([[1e-2/mm,-0],[-0,(1+0**2)/1e-2/mm]]) #[mm]
  covy = 7.945383034919336e-05/mm * np.asarray([[1e-2/mm,-0],[-0,(1+0**2)/1e-2/mm]]) #[mm]

#Calculate Envelope Center Angle
dBPM93to94 = 3031 #[mm] from OpenXAL(?)
dBPM93toPBW = 16822 #[mm] PBW 4400mm upstream of Target: https://gitlab.esss.lu.se/ess-bp/ess-lattice/-/blob/HEBT_RASTER_V29/9.0_HEBT/Beam_Physics/lattice.dat
dBPM93toTarg = 21222 #[mm] from Synoptic Viewer https://confluence.esss.lu.se/pages/viewpage.action?pageId=222397499
dPBWtoTarg = 4400 #[mm] from lattice and Synoptic
envXatBPM94 = args.rX #[mm] distance from beamline axis at BPM94, assume Cross Over at BPM93
envYatBPM94 = args.rY #[mm]
envXAngle = envXatBPM94 / dBPM93to94 #x' not radians, as per Kyrre 2.11.22
envYAngle = envYatBPM94 / dBPM93to94 #y'
beamletXAngle = 0 #default
beamletYAngle = 0 #default

#Envelope Center Offset, i.e. raster centre on BEW
cx_r = 0 + envXAngle * dBPM93toPBW #[mm] Takes into account angular drift since Cross Over
cy_r = 0 + envYAngle * dBPM93toPBW #[mm]
cx_r = cx_r * np.ones(N_t) #[mm]
cy_r = cy_r * np.ones(N_t) #[mm]

#Raster Amplitude: ax0,ay0 here produces the beamlet displacement AT the Target, as per Cyrille
ax0 = 54.65 #[mm]
ay0 = 18.37 #[mm]
#because generating particles just before PBW, must scale a0 by dPBWtoTarg * envAngle = (1- dPBWtoTarg / dBPM93toTarg)
amplScale = 1 - dPBWtoTarg / dBPM93toTarg #double check you account for Z before PBW in MiniScatter! beam production plane in GEANT, not exact PBW center!
ax = ax0 * amplScale * np.ones(N_t) #[mm]
ay = ay0 * amplScale * np.ones(N_t) #[mm]

#For weighting edges case (--edges argument)
Right = 50
Top = 17

i=0
j=0
k=0

centroids = np.zeros([N_t*n_tii,2])
for jj in range(N_t):
  for ii in range(n_tii):
    tjj_ii = t[jj] + delta_t[ii]

    #Calculate the Raster Magnet contribution to Beamlet Center Location relative to beamline axis, as per Cyrille Thomas
    beamletX = cx_r[jj] + 2 * ax[jj] / math.pi * math.asin(math.sin(2 * math.pi / px[jj] * tjj_ii )) #[mm]
    beamletY = cy_r[jj] + 2 * ay[jj] / math.pi * math.asin(math.sin(2 * math.pi / py[jj] * tjj_ii )) #[mm]

    #Calculate the total Beamlet Angle = Envelope Center Angle + the angle given to each beamlet by the Raster Magnets
    beamletXAngle = beamletX / dBPM93toTarg + envXAngle #[mm]
    beamletYAngle = beamletY / dBPM93toTarg + envYAngle #[mm]

    #save total beamlet position
    centroids[i,0] = beamletX
    centroids[i,1] = beamletY
    N = args.N

    #In case of weighting edges of raster 
    if args.edges:
      if cx0 > -Right and cx0 < Right and cy0 < Top and cy0 > -Top: #set weight depending on position
        N = 5
        j += 1
      else:
        N = args.N #for later on, write core to csv with difference weight than edge particles for ROOT.
        k += 1

    #Generate beamlet distributions
    rng = np.random.default_rng()
    ptsX = rng.multivariate_normal([beamletX,beamletXAngle],covx,size = N) #mean is [pos,ang]!
    ptsY = rng.multivariate_normal([beamletY,beamletYAngle],covy,size = N)
    
    for k in range(N): #put this beamlet into total. Could just be written, figure that out later.
      totX[N*i+k,0] = ptsX[k,0]
      totX[N*i+k,1] = ptsX[k,1]
      totY[N*i+k,0] = ptsY[k,0]
      totY[N*i+k,1] = ptsY[k,1]
    i +=1
print(i,j,k)

#Check on output parameters
print(np.max(centroids[:,0]),np.max(centroids[:,1]))
print(np.max(totX[:,0]),np.max(totY[:,0]))
#Remove 0,0 particles, should be none except for weighted edges case
print(np.shape(totX))
nonzero = np.not_equal(totX[:,0],0) #be careful!
totX = totX[nonzero]
totY = totY[nonzero]
print(np.shape(totX))

#Pick name based on beam
if args.ESS:
  name = "PBW_{:.0f}MeV_ESSBeam_aX{:.0f}mm_aY{:.0f}mm_N{:.1e}_{:.1e}us_cyrille".format(energy,ax0,ay0,len(totX),time_length)#dt.strftime("%H-%M-%S")
if args.pencil:
  name = "PBW_{:.0f}MeV_pencilBeam_aX{:.0f}mm_aY{:.0f}mm_N{:.1e}_{:.1e}us_cyrille".format(energy,ax0,ay0,len(totX),time_length)#dt.strftime("%H-%M-%S")
if envXatBPM94 != 0:
  name = name + "_X{:.0f}mrad".format(envXAngle*1e3)
if envYatBPM94 != 0:
  name = name + "_Y{:.0f}mrad".format(envYAngle*1e3)
print(name)

outname = picPWD+name
if args.s:
  with open(outname+".csv",mode = 'w',newline=None) as part_file:
    part_writer = csv.writer(part_file,delimiter = ',')
    for i in range(len(totX)):
      part_writer.writerow(["proton", totX[i,0], totX[i,1], totY[i,0], totY[i,1], z, energy])
  part_file.close()

finish = datetime.now()
print(finish-start)

if args.g:
  #found the below method: https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
  import mpl_scatter_density
  plt.close()
  fig = plt.figure()
  s1 = fig.add_subplot(1,1,1,projection="scatter_density")
  x=totX[:,0]
  y=totY[:,0]
  density = s1.scatter_density(x,y,cmap='jet')
  fig.colorbar(density,label=r"Protons/mm^2")
  s1.set_xlabel("X [mm]")
  s1.set_ylabel("Y [mm]")
  s1.set_title("Rastered Beam Number Density\n{:.1e} protons {:.2f}ms".format(len(totX),time_length*1e-3))
  plt.show()
  plt.close()
