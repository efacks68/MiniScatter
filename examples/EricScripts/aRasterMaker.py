import numpy as np
from beamConstants import *
import matplotlib.pyplot as plt
import csv, math,sys
import argparse

Nparts=50
save=False
#Use Argument Parser to control
parser = argparse.ArgumentParser()
parser.add_argument("--s",action='store_true')
args = parser.parse_args()
if args.s:
  save=True

t_pulse = round(0.04 * 1/14 * 1e6) # mus
pulse_start = 10
t_end = (2* pulse_start + t_pulse) /1000# - 2#3 ms
time_length = round(t_end * 1e3 ) #number of pulses
#time_length = int(2.0e4 * 1.0) #200kHz of chopper from paper

t = np.linspace(0,t_end,time_length)
N_t = len(t) # number of time samples
n_tii  = 10 
#print("t= ",N_t, time_length,t_end)

totX = np.zeros([N_t*n_tii*Nparts,2])
totY = np.zeros([N_t*n_tii*Nparts,2])
#print(len(totX))

#check units!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
covx = gemtx * np.asarray([[betax,-alphx],[-alphx,(1+alphx**2)/betax]])
covy = gemty * np.asarray([[betay,-alphy],[-alphy,(1+alphy**2)/betay]])

# raster centre on BEW
cx_r = 0 # mm
cy_r = 0 # mm
cx_r = cx_r * np.ones(N_t) # [mm]
cy_r = cy_r * np.ones(N_t) # [mm]
ax0 = 54.65 # [mm]
ay0 = 18.37 # [mm]
Fx = 39.953 # [kHz]
Fy = 28.7051 # [kHz]
ax = ax0 * np.ones(N_t) # [mm]
ay = ay0 * np.ones(N_t) # [mm]
pb = 1
px = pb/Fx * np.ones(N_t) #[s]
py = pb/Fy * np.ones(N_t) #[s]
#nx = 400 # nx is 'number of columns'
#ny = 350 # ny is 'number of rows'
dt =  np.mean(np.diff(t)) #[s]
delta_t = np.linspace(0,dt,n_tii) #[s]

i=0
j=0
Left = -49
Right = 49
Top = 15
Bottom = -15
centroids = np.zeros([N_t*n_tii,2])
print(np.shape(totX))

rng = np.random.default_rng()
ptsX = rng.multivariate_normal([0,0],covx,size = Nparts)
##plt.scatter(ptsX[:,0],ptsX[:,1])
#plt.hist(ptsX)
#plt.show()
print(np.covar(ptsX))

for jj in range(N_t):
  for ii in range(n_tii):

    tjj_ii = t[jj] + delta_t[ii] 
    #check units!!!!!!!!!!!!!!
    cx0 = cx_r[jj] + 2 * ax[jj] / math.pi * math.asin(math.sin(2*math.pi/px[jj] * tjj_ii )) #[mm]
    cy0 = cy_r[jj] + 2 * ay[jj] / math.pi * math.asin(math.sin(2*math.pi/py[jj] * tjj_ii )) #[mm]

    centroids[i,0] = cx0
    centroids[i,1] = cy0

    #if cx0 < Left or cx0 > Right or cy0 > Top or cy0 < Bottom: #set weight depending on position
    #  N = Nparts
    #  j += 1
    #else:
    N = Nparts#5 #for later on, write core to csv with difference weight than edge particles for ROOT.

    rng = np.random.default_rng()
    #Check units!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 19.10
    ptsX = rng.multivariate_normal([cx0,0],covx,size = N) #mean is [pos,ang]!
    ptsY = rng.multivariate_normal([cy0,0],covy,size = N)
    
    for k in range(N): #put this beamlet into total. Could just be written, figure that out later.
      totX[N*i+k,0] = ptsX[k,0]
      totX[N*i+k,1] = ptsX[k,1]
      totY[N*i+k,0] = ptsY[k,0]
      totY[N*i+k,1] = ptsY[k,1]
    i +=1
print(i,j)

print(np.max(centroids[:,0]),np.max(centroids[:,1]))
print(np.min(centroids[:,0]),np.min(centroids[:,1]))
print(np.max(totX[:,0]),np.max(totY[:,0]))
from datetime import datetime
dt = datetime.now()
z = 0.01
picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
#outname = picPWD+"PBW_{:.0f}MeV_eX{:.0f}um,eY{:.0f}um_bX{:.0f}m,bY{:.0f}m_aX{:.0f},aY{:.0f}_N{:.0e}_m{}_d{:.0f}mm".format(energy,nemtx*1e3,nemty*1e3,betax,betay,alphx,alphy,N,a,d*1e3)
outname = picPWD+"PBW_{:.0f}MeV_eX{:.0f}um,eY{:.0f}um_bX{:.0f}m,bY{:.0f}m_aX{:.0f},aY{:.0f}_N{:.1e}_cyrille".format(energy,nemtx*1e3,nemty*1e3,betax,betay,alphx,alphy,len(totX))#+"20-07-23"#dt.strftime("%H-%M-%S")

#from plotFit import numLines
#N = numLines(outname)

#nonzero = np.not_equal(totX[:,0],0) #be careful!
#totX = totX[nonzero]
#totY = totY[nonzero]
print(np.shape(totX))

#print(len(totX))
if save:
  with open(outname+".csv",mode = 'w',newline=None) as part_file:
    part_writer = csv.writer(part_file,delimiter = ',')
    for i in range(len(totX)):
      part_writer.writerow(["proton", totX[i,0], totX[i,1], totY[i,0], totY[i,1], z, energy])
  part_file.close()
  print(outname)

#plt.plot(totX[:,0], totY[:,0], '.', alpha=0.5)
plt.hist(totX[:,0],bins='auto')
#plt.axis('equal')
#plt.yaxis.set_minor_locator(MultipleLocator(10))
plt.grid()
plt.show()
