import numpy as np
import csv
#eps = .12 [mm-mrad]
#Beta : 100 to 2000, steps of 100
#Alpha: 0 to -100, steps of -5
#Produces sigmas: 1.1-24.5mm

#betas = [10,25,50,75,100, 150,200,250,500,1e3, 2e3,5e3] #12 -144
betas = [10,50,100,200,500,1e3,5e3] #7 -49
#betas = [10,50,100,200,500] #5 -25  
#alphas = np.arange(0,-100,-10) #10 - 100
alphas = np.arange(-5,-100,-20) #2 - 4
eps = 0.12 #[mm-mrad]
filename = "TwissRange0,12mm-mrad_"+str(len(betas))+"Bx-"+str(len(betas))+"By-"+str(len(alphas))+"Ax-"+str(len(alphas))+"Ay"
twissPwd = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/OpenXAL/OXALNotebooks/failureTwiss/"
print(filename,betas,alphas)

i = 0
Twiss = np.zeros([len(betas)**2*len(alphas)**2,6])
for bx in betas:
    for by in betas:
        for ax in alphas:
            for ay in alphas:
                Twiss[i] = [eps,bx,ax,eps,by,ay]
                i+=1
#print(Twiss)
with open(twissPwd+filename+".csv",mode = 'w',newline=None) as csv_file:
    csv_writer = csv.writer(csv_file,delimiter = ',')
    for i in range(len(Twiss[:,])):
        csv_writer.writerow([i,Twiss[i,0],Twiss[i,1],Twiss[i,2],Twiss[i,3],Twiss[i,4],Twiss[i,5]])
csv_file.close()
print(i)
