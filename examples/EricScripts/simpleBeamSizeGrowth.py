import numpy as np
import matplotlib.pyplot as plt


formRatio = np.array([0.15,0.5,1.0,1.5,7.84])
miniRatio = np.array([0.15,1,7.84])
miniSizeX = np.array([11.90, 16.27,16.33])
miniSizeY = np.array([16.41,16.44,11.80])
formSizeX = np.array([11.68,13.64,16.13,16.13,16.13])
formSizeY = np.array([16.36,16.36,16.36,14.76,11.59])

fig = plt.figure()
plt.scatter(formRatio,formSizeX,c='r',label="Analytical X")
plt.scatter(formRatio,formSizeY,c='b',label="Analytical Y")
plt.scatter(miniRatio,miniSizeX,c='y',label="GEANT4 X")
plt.scatter(miniRatio,miniSizeY,c='g',label="GEANT4 Y")
plt.text(miniRatio[2],miniSizeX[2]+0.1,"3",fontsize=10)
plt.xlabel(r"Beam Size Ratio at PBW, $\beta_X / \beta_Y$")
plt.ylabel(r"$\sigma$ [mm]")
plt.title(r"Target $\sigma$ Growth with Beam Size")
plt.legend(loc='center right')
#plt.savefig("GrowthBeamSize.svg")
plt.show()
plt.close()
