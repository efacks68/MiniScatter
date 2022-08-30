#read in
import numpy as np
import matplotlib.pyplot as plt
import csv
from datetime import datetime

fs=14
n=16
HistSigmaX = np.zeros(n)
HistSigmaY = np.zeros(n)
TargSigmaX = np.zeros(n)
TargSigmaY = np.zeros(n)
BeamSRatio = np.zeros(n)
OrigX = np.zeros(n)
OrigY = np.zeros(n)

with open('growthOutput.csv') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  line_count = 0
  z = 1
  for row in csv_reader:
    if line_count == 0:
      #print(f'Column names are {", ".join(row)}')
      line_count += 1
    elif line_count == 1 or line_count == 2 or line_count == 3 or line_count == 4:
      line_count += 1
      z +=1
    else:
      #print(row[0],row[6])
      #print(f"betaX: {row[0]}, betaY: {row[5]}, nemtY: {row[6]}")
      BeamSRatio[line_count-z] = float(row[0])
      HistSigmaX[line_count-z] = float(row[1])
      TargSigmaX[line_count-z] = float(row[2])
      OrigX[line_count-z] = float(row[3])
      HistSigmaY[line_count-z] = float(row[4])
      TargSigmaY[line_count-z] = float(row[5])
      OrigY[line_count-z] = float(row[6])
      line_count +=1
  #print("Processed ",line_count,"lines")
#print(datetime.now().strftime("%H-%M-%S"),"\n")

#for i in range(len(BeamSRatio)):
#  print(BeamSRatio[i],HistSigmaX[i],TargSigmaX[i],OrigX[i],HistSigmaY[i],TargSigmaY[i],OrigY[i])

nonzero = np.greater(BeamSRatio,0)
BeamSRatio = BeamSRatio[nonzero]
HistSigmaX = HistSigmaX[nonzero]
HistSigmaY = HistSigmaY[nonzero]
TargSigmaX = TargSigmaX[nonzero]
TargSigmaY = TargSigmaY[nonzero]
OrigX = OrigX[nonzero]
OrigY = OrigY[nonzero]

fig = plt.figure()
plt.tight_layout()
plt.scatter(BeamSRatio,HistSigmaX,c="aqua",label=r"GEANT4 $\sigma_x$")
plt.scatter(BeamSRatio,HistSigmaY,c="gold",label=r"GEANT4 $\sigma_y$")
plt.scatter(BeamSRatio,TargSigmaX,c="green",label=r"Müller $\sigma_x$")
plt.scatter(BeamSRatio,TargSigmaY,c="red",label=r"Müller $\sigma_y$")
plt.scatter(BeamSRatio,OrigX,c="steelblue",label=r"No PBW $\sigma_x$")
plt.scatter(BeamSRatio,OrigY,c="darkorange",label=r"No PBW $\sigma_y$")

plt.ylim([3,25.5])
plt.text(0.7,8.8,  r"$\frac{\beta_x}{\beta_y}$",fontsize=fs+1)
plt.text(1.0,8.8,  r"$\frac{410}{472}_{[\rm m]}$",fontsize=fs)
plt.text(6.7,5.6,  r"$\frac{1055}{142}_{[\rm m]}$",fontsize=fs)
plt.text(0.75,23.5,r"$\epsilon_{Nx} = 0.352$ [mm*mrad]",fontsize=fs-4)
plt.text(0.75,22  ,r"$\epsilon_{Ny} = 0.365$ [mm*mrad]",fontsize=fs-4)

plt.xlabel(r"Beam Size Ratio at PBW, $\beta_X / \beta_Y$",fontsize=fs)
plt.ylabel(r"$\sigma$ at Target [mm]",fontsize=fs)
plt.title(r"Target $\sigma$ Growth with PBW Beam Roundness",fontsize = fs+2)

plt.legend(ncol=3,loc="lower left",columnspacing=0.02,handletextpad=0.1,fontsize=fs-4)
savename = "../../../Pictures/BeamSizeDependence_"+datetime.now().strftime("%H-%M-%S")+".pdf"
print(savename)
#plt.show()
fig.savefig(savename)
plt.close()

