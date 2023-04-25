#createReferenceImage.py

def decode(path,file,lenx,leny):
    import csv
    from numpy import zeros
    Img = zeros((leny,lenx))
    j=0
    if file[-1] != "v":
        print(file)
        file.append(".csv")
        print(file)
    with open(path+file,mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            #print(j)
            for i in range(lenx):
                Img[j,i] = float(row[i])
            j+=1
            if j==1000:
                break
    csv_file.close()
    #print(Img.shape)
    return Img

import csv,numpy as np,os,sys,re
choice=sys.argv[1]
basePath="/scratch2/ericdf/PBWScatter/"
if choice == "nominal":
    path = "/scratch2/ericdf/PBWScatter/ImageCSVs/PBW_Nominal_Nb100/"
    refImgName = "PBW_570MeV_beta1085.63,136.06m_RMamp49,16mm_N2.9e+06_NpB100_runW_QBZ_rB1_refImg"
elif choice == "1":
    path = "/scratch2/ericdf/PBWScatter/ImageCSVs/HEBT-A2T_100pctField_1.0e-04Jitter_QBZ/"
    refImgName = "HEBT-A2T_100pctField_1.0e-04Jitter_250x_570MeV_OrigbX1085.63,bY136.06m_N2.9e+06_NpB100_runW_QBZ_rB1_refImg"
elif choice == "10":
    path = "/scratch2/ericdf/PBWScatter/ImageCSVs/HEBT-A2T_100pctField_1.0e-03Jitter_QBZ/"
    refImgName = "HEBT-A2T_100pctField_1.0e-03Jitter_200x_570MeV_OrigbX1085.63,bY136.06m_N2.9e+06_NpB100_runW_QBZ_rB1_refImg"

if os.uname()[1] == "tensor.uio.no":
    files = os.listdir(path)
    print(path,len(files))
refImg = np.zeros((1000,1000))
lenx = np.shape(refImg)[0]
leny = np.shape(refImg)[1]

if not os.path.isfile(basePath+refImgName+".csv"):
    print("Making new Reference Image!")
    counter = 0
    for n in range(len(files)):
        Img = decode(path,files[n],lenx,leny)
        refImg +=Img
        counter+=1

    refImg /= counter
    print(refImg.max())

    with open(basePath+refImgName+".csv",mode='w') as csv_file:
        csv_writer = csv.writer(csv_file,delimiter = ',')
        csv_writer.writerows(refImg)
        csv_file.close()
    print(path+refImgName+".csv")
else:
    print("Found Reference Image!")
    refImg = decode(basePath,refImgName,lenx,leny)
    print(refImg[500,500])

from matplotlib.pyplot import subplots,pcolormesh,close,tight_layout,savefig,setp,title,xlabel,ylabel
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm
xax = np.arange(-lenx,lenx+1,2)
yax = np.arange(-leny,leny+1,2)
print(xax[-1])
X, Y = np.meshgrid(xax,yax) #Set mesh of plot from the axes given from converter function
close() #make sure no prior plotting messes up

fig,ax = subplots(dpi=300)
#print(datetime.now().strftime("%H-%M-%S"))

#Set maximum value depending on the maximum current density
from math import log10,ceil
maxim = ceil(refImg.max() / 10) * 10 +10
minim = 10**ceil(log10(0.001))
#cbarVals  = [minim,minim*10,minim*100,maxim] #make array for color bar values
#cbarLabels = ["{:.2f}".format(cbarVals[0]),"{:.1f}".format(cbarVals[1]),"{:.1f}".format(cbarVals[2]),"{:.0f}".format(cbarVals[3])]
cbarLabel = "Protons"
#print("Max current Density: ",Img.max(),"/",maxim,datetime.now().strftime("%H-%M-%S"))

#Use pcolor to show density map, use log scale
c = ax.pcolormesh(X,Y,refImg,shading='auto',norm=LogNorm(vmin=minim, vmax=maxim), cmap='viridis',rasterized=True) #viridis or magma are perceptually uniform

lw=1
fs=15

#Set Plot Properties
setp(ax,xlim=([-150,150]),ylim=([-100-10,100]))
ax.set_title("Reference Macro-Particle Distribution at Target",fontsize=fs+2)
ax.set_ylabel("Vertical [mm]",fontsize=fs)
ax.set_xlabel("Horizontal [mm]",fontsize=fs)
cbar = fig.colorbar(c, ax=ax,pad=0.01)
cbar.set_label(cbarLabel,fontsize=fs-1)
#cbar.set_ticks(cbarVals)
#cbar.set_ticklabels(cbarLabels)
picpath=statsPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
tight_layout()
savefig(picpath+refImgName+".png")#
close(fig)
close()
print(picpath+refImgName+".png")