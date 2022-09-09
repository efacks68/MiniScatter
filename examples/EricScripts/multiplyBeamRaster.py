import csv, sys, math
import numpy as np
from decimal import *
print(getcontext())

print(sys.argv[1])
picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
beamFile = picPWD+sys.argv[1] #"PBW_570MeV_eX113um,eY122um_bX941m,bY120m_aX-59,aY-7_N1e+03"

#find number of particles
with open(beamFile+".csv") as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  N = 0
  for row in csv_reader:
    N += 1
csv_file.close()
print("There are ",N,"particles")

particle = []
x = np.zeros(N)
px = np.zeros(N)
y = np.zeros(N)
py = np.zeros(N)
z = np.zeros(N)
E = np.zeros(N)

#read in particles
with open(beamFile+".csv") as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  N = 0
  for row in csv_reader:
    particle.append(row[0])
    x[N] = Decimal(row[1]) #give up at this point, since it doesn't read in the correct values!
    px[N] = Decimal(row[2])
    y[N] = Decimal(row[3])
    py[N] = Decimal(row[4])
    z[N] = Decimal(row[5])
    E[N] = Decimal(row[6])

print("{:.17f}".format(x[1]),E[1])
#beam sizes in [mm]
sigx = 3.47 #[mm]
sigy = 3.404 #[mm]

a=4 #number of beams to make
evod = a % 2
if evod == 0: #even
  xShift = sigx * 0.5 #[mm]
  yShift = sigy * 0.5 #[mm]

  outFile = beamFile + "_mult" + str(a)
  print(len(E))
  with open(outFile+".csv",mode = 'w') as part_file:
    part_writer = csv.writer(part_file,delimiter = ',')
    print("+X,+Y")
    for i in range(len(E)):
      part_writer.writerow([particle[i], "{:.19f}".format(math.fsum([x[i],xShift*(evod + 1)])), "{:.19f}".format(px[i]), "{:.19f}".format(math.fsum([y[i],yShift*(evod + 1)])), "{:.19f}".format(py[i]), "{:.3f}".format(z[i]), "{:.3f}".format(E[i])])
      if i == 1:
        print(particle[i], "{:.19f}".format(math.fsum([x[i],xShift*(evod)])), "{:.19f}".format(px[i]), "{:.19f}".format(math.fsum([y[i],yShift*(evod)])), "{:.19f}".format(py[i]), "{:.3f}".format(z[i]), "{:.3f}".format(E[i]))
        print(particle[i], "{:.19f}".format(math.fsum([x[i],xShift*(evod + 1)])), "{:.19f}".format(px[i]), "{:.19f}".format(math.fsum([y[i],yShift*(evod + 1)])), "{:.19f}".format(py[i]), "{:.3f}".format(z[i]), "{:.3f}".format(E[i]))
    print("+X,-Y")
    for i in range(len(E)):
      part_writer.writerow([particle[i], "{:.19f}".format(math.fsum([x[i],xShift*(evod + 1)])), "{:.19f}".format(px[i]), "{:.19f}".format(math.fsum([y[i],-yShift*(evod + 1)])), "{:.19f}".format(py[i]), "{:.3f}".format(z[i]), "{:.3f}".format(E[i])])
    print("-X,+Y")
    for i in range(len(E)):
      part_writer.writerow([particle[i], "{:.19f}".format(math.fsum([x[i],-xShift*(evod + 1)])), "{:.19f}".format(px[i]), "{:.19f}".format(math.fsum([y[i],yShift*(evod + 1)])), "{:.19f}".format(py[i]), "{:.3f}".format(z[i]), "{:.3f}".format(E[i])])
    print("-X,-Y")
    for i in range(len(E)):
      part_writer.writerow([particle[i], "{:.19f}".format(math.fsum([x[i],-xShift*(evod + 1)])), "{:.19f}".format(px[i]), "{:.19f}".format(math.fsum([y[i],-yShift*(evod + 1)])), "{:.19f}".format(py[i]), "{:.3f}".format(z[i]), "{:.3f}".format(E[i])])
  part_file.close()

print(outFile)

