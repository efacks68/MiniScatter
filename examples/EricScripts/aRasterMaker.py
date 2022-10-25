import numpy as np
from beamConstants import *
import matplotlib.pyplot as plt
import csv, math,sys
import argparse

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
args = parser.parse_args()

t_pulse = round(0.04 * 1/14 * 1e6) # mus
pulse_start = 10
t_end = (2* pulse_start + t_pulse) /1000# - 2#3 ms
time_length = round(t_end * args.a ) #number of pulses, nominal = 2.86e3

t = np.linspace(0,t_end,time_length)
N_t = len(t) # number of time samples
n_tii  = 10 
#print("t= ",N_t, time_length,t_end)

totX = np.zeros([N_t*n_tii*args.N,2])
totY = np.zeros([N_t*n_tii*args.N,2])
#print(len(totX))

# raster centre on BEW
cx_r = 0 # mm
cy_r = 0 # mm
cx_r = cx_r * np.ones(N_t) # [mm]
cy_r = cy_r * np.ones(N_t) # [mm]
#ax0,ay0 here produce the beam size ON target, not before
b = 1
ax0 = 54.65*b # [mm]
ay0 = 18.37*b # [mm]
ax = ax0 * np.ones(N_t) # [mm]
ay = ay0 * np.ones(N_t) # [mm]
Fx = 39.953*1e3 #[kHz]
Fy = 28.7051*1e3 #[kHz]
pb = 1
px = pb/Fx * np.ones(N_t) #[s]
py = pb/Fy * np.ones(N_t) #[s]
beamAngleX = 0
beamAngleY = 0
dt =  np.mean(np.diff(t)) #[s]
delta_t = np.linspace(0,dt,n_tii) #[s]

if args.ESS: #assume 570MeC
  covx = gemtx/mm * np.asarray([[betax/mm,-alphx],[-alphx,(1+alphx**2)/(betax/mm)]]) #[mm]
  covy = gemty/mm * np.asarray([[betay/mm,-alphy],[-alphy,(1+alphy**2)/(betay/mm)]]) #[mm]
  name = "PBW_{:.0f}MeV_ESSBeam_aX{:.0f}mm_aY{:.0f}mm_N{:.1e}_{:.1e}us_cyrille".format(energy,ax0,ay0,len(totX),time_length)#dt.strftime("%H-%M-%S")

if args.pencil: #assume 570MeV
  covx = 7.945383034919336e-05/mm * np.asarray([[1e-2/mm,-0],[-0,(1+0**2)/1e-2/mm]]) #[mm]
  covy = 7.945383034919336e-05/mm * np.asarray([[1e-2/mm,-0],[-0,(1+0**2)/1e-2/mm]]) #[mm]
  name = "PBW_{:.0f}MeV_pencilBeam_aX{:.0f}mm_aY{:.0f}mm_N{:.1e}_{:.1e}us_cyrille".format(energy,ax0,ay0,len(totX),time_length)#dt.strftime("%H-%M-%S")

#Angle Contribution this is the x', the 2nd part of [cx0,0]
dBPM93to94 = 3031 #[mm]
xAtBPM93 = args.rX #[mm] from axis
yAtBPM93 = args.rY #[mm] from axis
if xAtBPM93 !=0:
  beamAngleX = np.arctan(xAtBPM93/dBPM93to94) #[rad]
  name = name + "_X{:.0f}mrad".format(beamAngleX*1e3)
if yAtBPM93 !=0:
  beamAngleY = np.arctan(yAtBPM93/dBPM93to94) #[rad]
  name = name + "_Y{:.0f}mrad".format(beamAngleY*1e3)
#print("dx",beamAngleX,"dy",beamAngleY)
print(name)

i=0
j=0
#Right = 52*mm
#Top = 17*mm

centroids = np.zeros([N_t*n_tii,2])
for jj in range(N_t):
  for ii in range(n_tii):
    tjj_ii = t[jj] + delta_t[ii]
    cx0 = cx_r[jj] + 2 * ax[jj] / math.pi * math.asin(math.sin(2*math.pi/px[jj] * tjj_ii )) #[mm]
    cy0 = cy_r[jj] + 2 * ay[jj] / math.pi * math.asin(math.sin(2*math.pi/py[jj] * tjj_ii )) #[mm]

    centroids[i,0] = cx0
    centroids[i,1] = cy0

    #if cx0 > -Right and cx0 < Right and cy0 < Top and cy0 > -Top: #set weight depending on position
    #  N = 10
    #  j += 1
    #else:
    #N = args.N #for later on, write core to csv with difference weight than edge particles for ROOT.
    #  k += 1

    rng = np.random.default_rng()
    ptsX = rng.multivariate_normal([cx0,beamAngleX],covx,size = args.N) #mean is [pos,ang]!
    ptsY = rng.multivariate_normal([cy0,beamAngleY],covy,size = args.N)
    
    for k in range(args.N): #put this beamlet into total. Could just be written, figure that out later.
      totX[args.N*i+k,0] = ptsX[k,0]
      totX[args.N*i+k,1] = ptsX[k,1]
      totY[args.N*i+k,0] = ptsY[k,0]
      totY[args.N*i+k,1] = ptsY[k,1]
    i +=1
print(i,j)

print(np.max(centroids[:,0]),np.max(centroids[:,1]))
print(np.min(centroids[:,0]),np.min(centroids[:,1]))
print(np.max(totX[:,0]),np.max(totY[:,0]))
from datetime import datetime
dt = datetime.now()
z = -10
picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
outname = picPWD+name

nonzero = np.not_equal(totX[:,0],0) #be careful!
totX = totX[nonzero]
totY = totY[nonzero]
print(np.shape(totX))

#print(len(totX))
if args.s:
  with open(outname+".csv",mode = 'w',newline=None) as part_file:
    part_writer = csv.writer(part_file,delimiter = ',')
    for i in range(len(totX)):
      part_writer.writerow(["proton", totX[i,0], totX[i,1], totY[i,0], totY[i,1], z, energy])
  part_file.close()

if args.g:
  plt.close()
  fig = plt.figure(figsize=(15,6.0))
  plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
  s1 = fig.add_subplot(1,2,1)
  s2 = fig.add_subplot(1,2,2)
  s1.plot(totX[:,0], totY[:,0], '.', alpha=0.5)
  s2.hist(totX[:,0],bins='auto')
  s1.set_xlabel("X [mm]")
  s2.set_xlabel("X [mm]")
  plt.grid()
  plt.show()
  plt.close()
