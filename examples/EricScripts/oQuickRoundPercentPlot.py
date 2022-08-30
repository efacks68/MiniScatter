#read in
import numpy as np
import matplotlib.pyplot as plt
import csv
from datetime import datetime

fs=14
n=16
HistSigmaX = np.zeros(n)
HistSigmaY = np.zeros(n)
EqSigmaX = np.zeros(n)
EqSigmaY = np.zeros(n)
BeamSRatio = np.zeros(n)
OrigX = np.zeros(n)
OrigY = np.zeros(n)
PBWpOut3sigX = np.zeros(n)
PBWpOut3sigY = np.zeros(n)
noPBW = np.zeros(1)
VacpOut3sigX = np.zeros(1)
VacpOut3sigY = np.zeros(1)

with open('roundnessDependence30Aug.csv') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  line_count = 0
  z = 1
  for row in csv_reader:
    if line_count == 0:
      #print(f'Column names are {", ".join(row)}')
      line_count += 1
      z += 1
    elif line_count == 1:
      noPBW[0] = float(row[0])
      VacpOut3sigX = float(row[1])
      VacpOut3sigY = float(row[2])
      line_count += 1
      z += 1
    elif line_count == 2 or line_count == 3 or line_count == 4 or line_count == 5:
      line_count += 1
      z += 1
    else:
      #print(row[0],row[6])
      #print(f"betaX: {row[0]}, betaY: {row[5]}, nemtY: {row[6]}")
      BeamSRatio[line_count-z] = float(row[0])
      HistSigmaX[line_count-z] = float(row[1])
      EqSigmaX[line_count-z] = float(row[2])
      OrigX[line_count-z] = float(row[3])
      PBWpOut3sigX[line_count-z] = float(row[4])
      HistSigmaY[line_count-z] = float(row[5])
      EqSigmaY[line_count-z] = float(row[6])
      OrigY[line_count-z] = float(row[7])
      PBWpOut3sigY[line_count-z] = float(row[8])
      line_count +=1
  #print("Processed ",line_count,"lines")
#print(datetime.now().strftime("%H-%M-%S"),"\n")

#for i in range(len(BeamSRatio)):
#  print(BeamSRatio[i],HistSigmaX[i],EqSigmaX[i],OrigX[i],HistSigmaY[i],EqSigmaY[i],OrigY[i])

#filter out zeros
nonzero = np.greater(BeamSRatio,0)
BeamSRatio = BeamSRatio[nonzero]
HistSigmaX = HistSigmaX[nonzero]
HistSigmaY = HistSigmaY[nonzero]
EqSigmaX = EqSigmaX[nonzero]
EqSigmaY = EqSigmaY[nonzero]
OrigX = OrigX[nonzero]
OrigY = OrigY[nonzero]
PBWpOut3sigX = PBWpOut3sigX[nonzero]
PBWpOut3sigY = PBWpOut3sigY[nonzero]
#VacpOut3sigX = VacpOut3sigX[nonzero]
#VacpOut3sigY = VacpOut3sigY[nonzero]

#Create the fig with 2 plots side by side
plt.clf()
fig = plt.figure(figsize=(15,6.0))
plt.subplots_adjust(wspace=0.17) #increase width space to not overlap
s1 = fig.add_subplot(1,2,1)
s2 = fig.add_subplot(1,2,2)

s1.scatter(BeamSRatio,EqSigmaX,c="green",label=r"Müller $\sigma_x$",s=65)
s1.scatter(BeamSRatio,EqSigmaY,c="red",label=r"Müller $\sigma_y$",s=65)
s1.scatter(BeamSRatio,HistSigmaX,c="aqua",label=r"GEANT4 $\sigma_x$")
s1.scatter(BeamSRatio,HistSigmaY,c="gold",label=r"GEANT4 $\sigma_y$")
s1.scatter(BeamSRatio,OrigX,c="steelblue",label=r"No PBW $\sigma_x$")
s1.scatter(BeamSRatio,OrigY,c="darkorange",label=r"No PBW $\sigma_y$")

s2.scatter(BeamSRatio,PBWpOut3sigX,c="aqua",label=r"PBW % >3$\sigma_x$")
s2.scatter(BeamSRatio,PBWpOut3sigY,c="gold",label=r"PBW % >3$\sigma_y$")
s2.scatter(noPBW,VacpOut3sigX,c="steelblue",label=r"No PBW % >3$\sigma_x$")
s2.scatter(noPBW,VacpOut3sigY,c="darkorange",label=r"No PBW % >3$\sigma_y$")

s1.set_ylim([6,25.5])
s1.text(0.7,21,r"$\frac{\beta_x}{\beta_y}$",fontsize=fs+1)
s1.text(1.1,21,r"$\frac{410}{472}_{[\rm m]}$",fontsize=fs)
s1.text(6.7,19.5,r"$\frac{1055}{142}_{[\rm m]}$",fontsize=fs)
s1.text(0.7,24.4,r"$\epsilon_{Nx} = 0.352$ [mm*mrad]",fontsize=fs-4)
s1.text(0.7,23.5,r"$\epsilon_{Ny} = 0.365$ [mm*mrad]",fontsize=fs-4)

s1.set_xlabel(r"Beam Size Ratio at PBW, $\beta_X / \beta_Y$",fontsize=fs)
s1.set_ylabel(r"$\sigma$ at Target [mm]",fontsize=fs)
s1.set_title(r"Target $\sigma$ Growth"+"\nwith Beam Roundness at PBW",fontsize=fs+2)
s1.legend(ncol=3,loc="lower left",columnspacing=0.02,handletextpad=0.1,fontsize=fs-4)

s2.set_ylim([0.10,1.5])
s2.text(0.7,1.20,r"$\frac{\beta_x}{\beta_y}$",fontsize=fs+1)
s2.text(1.1,1.20,r"$\frac{410}{472}_{[\rm m]}$",fontsize=fs)
s2.text(6.7,1.08,r"$\frac{1055}{142}_{[\rm m]}$",fontsize=fs)
s2.text(0.7,1.42,r"$\epsilon_{Nx} = 0.352$ [mm*mrad]",fontsize=fs-4)
s2.text(0.7,1.36,r"$\epsilon_{Ny} = 0.365$ [mm*mrad]",fontsize=fs-4)

s2.set_ylabel(r"% Outside 3$\sigma$ at Target",fontsize=fs)
s2.set_xlabel(r"Beam Size Ratio at PBW, $\beta_X / \beta_Y$",fontsize=fs)
s2.set_title("Target Halo Growth \nwith Beam Roundness at PBW",fontsize = fs+2)
s2.legend(ncol=2,loc="lower right",columnspacing=0.02,handletextpad=0.1,fontsize=fs-4)

savename = "../../../Pictures/BeamSize,HaloDependence_"+datetime.now().strftime("%H-%M-%S")+".pdf"
print(savename)
#plt.show()
fig.savefig(savename,bbox_inches="tight")
plt.close()

