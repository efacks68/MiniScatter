import numpy as np
import matplotlib.pyplot as plt

#formThick = np.array([3.35,3.75,4.25,4.75,5.25])
#miniThick = np.array([3.35,4.25,5.25])
formThick = np.array([  0.2,  0.5,  1.0, 1.50, 2.25,  3.0,  3.5, 4.25,    5, 5.75,  6.5, 7.25,    8,   9,   10])
formSizeX = np.array([13.11,13.54,14.29,15.05,16.13,17.18,17.87,18.85,19.80,20.73,21.61,22.48,23.31,24.40,25.44])
formSizeY = np.array([ 6.78, 7.59, 8.86,10.03,11.59,13.02,13.91,15.15,16.32,17.43,18.48,19.48,20.44,21.66,22.84])
miniThick = np.array([  0.2,  0.5,  1.0, 2.25,  3.0, 4.25, 5.75,    8, 10])
miniSizeX = np.array([13.27,13.82,14.56,16.17,17.14,18.13,19.09,19.94,20.34])
miniSizeY = np.array([6.84 , 7.69, 9.10,11.75,13.10,14.75,16.23,17.68,18.61])
noPBW = np.array([0])
origX = np.array([12.07])
origY = np.array([ 4.46])

fig = plt.figure()
plt.tight_layout()
fs=14
plt.scatter(formThick,formSizeX,c='green',label=r"Müller $\sigma_x$",s=65)
plt.scatter(formThick,formSizeY,c='red',  label=r"Müller $\sigma_y$",s=65)
plt.scatter(miniThick,miniSizeX,c='aqua',label=r"GEANT4 $\sigma_x$")
plt.scatter(miniThick,miniSizeY,c='gold',label=r"GEANT4 $\sigma_y$")
plt.scatter(noPBW,origX, c="steelblue",label=r"No PBW $\sigma_x$")
plt.scatter(noPBW,origY,c="darkorange",label=r"No PBW $\sigma_y$")

plt.text(0,23.9,"Beam Twiss at PBW:",fontsize=fs-4)
plt.text(0,22.9,r"$\epsilon_{Nx,Ny} = 0.113, 0.122$ [mm*mrad]",fontsize=fs-4)
plt.text(0,21.9,r"$\beta_{x,y} = 941, 120$ [m]",fontsize=fs-4)
plt.text(0,20.9,r"$\alpha_{x,y} = -59, -7.5$",fontsize=fs-4)
plt.annotate('Planned PBW Al \nThickness=2.25mm',
    xy=(2.4, 11.7), xycoords='data',
    xytext=(40, -29), textcoords='offset points',
    arrowprops=dict(arrowstyle="->",
                    connectionstyle="arc,angleA=90,armA=5,angleB=-30,armB=30,rad=7"))
#plt.annotate('', #down
#    xy=(2.25, 4.5), xycoords='data',
#    xytext=(38.5, 73), textcoords='offset points',
#    arrowprops=dict(arrowstyle="->",
#                    connectionstyle="arc,angleA=-115,armA=40,angleB=40,armB=0,rad=7"))
plt.annotate('', #up
    xy=(2.3, 16), xycoords='data',
    xytext=(41, -66), textcoords='offset points',
    arrowprops=dict(arrowstyle="->",
                    connectionstyle="arc,angleA=140,armA=35,angleB=30,armB=0,rad=7"))

plt.xlabel("PBW Aluminium Thickness [mm]",fontsize=fs)
plt.ylabel(r"$\sigma$ at Target [mm]",fontsize=fs)
plt.title(r"Target $\sigma$ Growth with PBW Thickness",fontsize=fs+2)

#handles1, labels1 = plt.gca().get_legend_handles_labels()
#print(labels1)
#order1=[0,1,2,3]
#([handles1[idx] for idx in order1],[labels1[idx] for idx in order1]
plt.legend(ncol=3,loc='lower right',columnspacing=0.02,handletextpad=0.1,fontsize=fs-4)
from datetime import datetime
savename = "../../../Pictures/AlThicknessDependence_"+datetime.now().strftime("%H-%M-%S")+".pdf"
#savename = "../../../Pictures/AlThicknessDependence.pdf"
fig.savefig(savename)
#plt.show()
plt.close()
