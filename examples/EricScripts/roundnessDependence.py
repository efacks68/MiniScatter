from datetime import datetime
origin = datetime.now().strftime("%H-%M-%S")
print(origin,"\n")
import numpy as np
import matplotlib.pyplot as plt
import os,sys,re
import csv

energy = 570
partZ = 1 #[C]
partA = 938.27209 #[MeV/c2]
um = 1e-6 #[m]
#energy = float(input("What is the beam Energy? "))
gamma_rel = 1 + energy/partA #from PrintTwissParameters
beta_rel = np.sqrt(gamma_rel*gamma_rel -1 )/gamma_rel
N=1e5
fs = 14 #set the axis label font size early
i=0

#first, make test arrays of Twiss to send
#Swinging back, read in CSV file Twiss and make everything a loop
fig = plt.figure()
Twiss = np.zeros((17,6))
#Orig = np.zeros((17,2))
with open('beam_sizes18Aug.csv') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  line_count = 0
  z=1
  for row in csv_reader:
    if line_count == 0:
      #print(f'Column names are {", ".join(row)}')
      line_count += 1
    #elif line_count ==1 or line_count == 2 or line_count == 3 or line_count == 4: #
    #  line_count += 1
    #  z += 1
    else:
      #print(row[0],row[6])
      #print(f"betaX: {row[0]}, betaY: {row[5]}, nemtY: {row[6]}")
      Twiss[line_count-z] = [row[0],row[1],row[2],row[4],row[5],row[6]]
      #Orig[line_count-z] = [row[3],row[7]]
      #print(Twiss[line_count-z])
      line_count += 1
  print("Processed ",line_count,"lines")
csv_file.close()

print(datetime.now().strftime("%H-%M-%S"),"\n")

n=line_count-z
HistSigmaX = np.zeros(n)
HistSigmaY = np.zeros(n)
EqSigmaX = np.zeros(n)
EqSigmaY = np.zeros(n)
BeamSRatio = np.zeros(n)
OrigX = np.zeros(n)
OrigY = np.zeros(n)
PBWpOut3sigX = np.zeros(n)
PBWpOut3sigY = np.zeros(n)
VacpOut3sigX = np.zeros(n)
VacpOut3sigY = np.zeros(n)
#print(Twiss[2])

#Twiss = [betaX,alphaX,emtX,betaY,alphaY,emtY]
#Twiss = np.array([144,-8,0.351,88,-1,0.365])
#Twiss = np.array([941.25,-58.81,0.11315,120.03,-7.46,0.12155])
for i in range(line_count-1): #6):#
  if Twiss[i][0] == 0.0:
    continue
  print("line",i)
  # pass values to scripts
  cmd = "python3 SimpleDemoPBWscript.py"
  cmd = cmd + " " + str(N) #necessary!
  #print(len(Twiss[i]),Twiss[i])
  for j in range(len(Twiss[i])):
    cmd = cmd + " " + str(Twiss[i][j])
  cmd += ' 0' #for thickness, the case for MagnetPBW
  print(cmd)

  #Get scripts to run
  #os.system(cmd)
  stream = os.popen(cmd)
  output = stream.read()
  print(output)

  #make accessible the values I want (from MiniScatter getData function)
  if re.search("PBW_570MeV",output): #the case of MagnetPBW
    #print(re.search("(?<=(Histogram Real PBW, Eq 8 sigma X: )).+(?=(mm))",output))
    HistSigmaX[i] = float(re.search("(?<=(Histogram Real PBW, Eq 8 sigma X: )).+(?=(mm))",output).group(0))
    HistSigmaY[i] = float(re.search("(?<=(Histogram Real PBW, Eq 8 sigma Y: )).+(?=(mm))",output).group(0))
    EqSigmaX[i] = float(re.search("(?<=(Fit       Real PBW, Eq 8 sigma X: )).+(?=(mm))",output).group(0))
    EqSigmaY[i] = float(re.search("(?<=(Fit       Real PBW, Eq 8 sigma Y: )).+(?=(mm))",output).group(0))
    PBWpOut3sigX[i] = float(re.search("(?<=(Real PBW, Eq 8 )).+(?=( % outside 3 sigma X))",output).group(0))
    PBWpOut3sigY[i] = float(re.search("(?<=(Real PBW, Eq 8 )).+(?=( % outside 3 sigma Y))",output).group(0))
    OrigX[i] =      float(re.search("(?<=(Fit       No PBW sigma X: )).+(?=(mm))",output).group(0)) #this is init Twiss drifted to Target
    OrigY[i] =      float(re.search("(?<=(Fit       No PBW sigma Y: )).+(?=(mm))",output).group(0))

    #print(Orig[i])

  # do stuff with values
  #make BeamSRatio array
  BeamSRatio[i] = Twiss[i][0]/Twiss[i][3]
  print("BetaX / BetaY = ",BeamSRatio[i])
  print("GEANT4:",HistSigmaX[i],HistSigmaY[i],"\nEq8:",EqSigmaX[i],EqSigmaY[i])
  print("NoPBW: ",OrigX[i],OrigY[i],"\n% Outside 3sigma:",PBWpOut3sigX[i],PBWpOut3sigY[i])

  #print(datetime.now().strftime("%H-%M-%S"),"\n")
  print("Now no PBW")
  if Twiss[i][0] == 0.0:
    continue
  print("line",i)
  # pass values to scripts
  cmd = "python3 SimpleDemoPBWscript.py"
  cmd = cmd + " " + str(N) #necessary!
  #print(len(Twiss[i]),Twiss[i])
  for j in range(len(Twiss[i])):
    cmd = cmd + " " + str(Twiss[i][j])
  cmd += ' 0.1' #for thickness, the case for no PBW
  print(cmd)

  #Get scripts to run
  #os.system(cmd)
  stream = os.popen(cmd)
  output = stream.read()
  print(output)
  if re.search("simplePBW_0.1mmVacuum",output): #the case of Vacuum
    #add in 2nd y axis to show % outside 3 sigma on the right! (figure out how to do with the change?)
    VacpOut3sigX[i] = float(re.search("(?<=(No PBW, Eq 8 )).+(?=( % outside 3 sigma X))",output).group(0))
    VacpOut3sigY[i] = float(re.search("(?<=(No PBW, Eq 8 )).+(?=( % outside 3 sigma Y))",output).group(0))
  
  #print(BeamSRatio[i],HistSigmaX[i],EqSigmaX[i],OrigX[i],HistSigmaY[i],EqSigmaY[i],OrigY[i],PBWpOut3sigX[i],PBWpOut3sigY[i],VacpOut3sigX[i],VacpOut3sigY[i],"\n\n")

noPBW = np.array([0.8])
VacpOut3sigX = np.array([0.269])
VacpOut3sigY = np.array([0.25])

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

savename = "../../../Pictures/BeamSizeDependence_"+datetime.now().strftime("%H-%M-%S")+".svg"
print(savename)
fig.savefig(savename,bbox_inches="tight")
plt.close()

#change this to write to CSV!
with open('roundnessDependence30Aug.csv',mode = 'w') as round_file:
  round_writer = csv.writer(round_file,delimiter = ',')
  round_writer.writerow(["BeamSRatio","HistSigmaX","EqSigmaX","OrigX","PBWpOut3sigX","HistSigmaY","EqSigmaY","OrigY","PBWpOut3sigY"])
  round_writer.writerow([noPBW[0],VacpOut3sigX[0],VacpOut3sigY[0]])
  for i in range(len(BeamSRatio)):
    round_writer.writerow([BeamSRatio[i],HistSigmaX[i],EqSigmaX[i],OrigX[i],PBWpOut3sigX[i],HistSigmaY[i],EqSigmaY[i],OrigY[i],PBWpOut3sigY[i]])
round_file.close()

print(origin)
print("now",datetime.now().strftime("%H-%M-%S"))