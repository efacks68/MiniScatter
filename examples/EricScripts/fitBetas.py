import csv,sys,re
from plotFit import findFit,gaussian
from numpy import zeros,mean,std,greater
import matplotlib.pyplot as plt

name=sys.argv[1]
from plotFit import numLines
n = numLines("../../../Pictures/"+name)
#print(n)
betaX=zeros(n)
betaY=zeros(n)
betaX.fill(-750)
betaY.fill(-750)

with open("../../../Pictures/"+name+".csv",mode='r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    next(csv_reader) #skips header
    rowNum = 0
    for row in csv_reader:
        #print(sample,rowNum,row[0],row[1],row[2])
        
        betaX[rowNum] = float(row[1])
        betaY[rowNum] = float(row[2])
        #print(rowNum,"Read in Betas:",betaX[rowNum],betaY[rowNum])
        rowNum+=1
#print(re.search("(([0-9]*)+(?=Pct))",name)[1])
if   int(re.search("(([0-9]*)+(?=Pct))",name)[1]) == 1:
    nBins = 8 #30 for 10pct, 8 for 1pct
elif int(re.search("(([0-9]*)+(?=Pct))",name)[1]) ==10:
    nBins = 30

nonzero = greater(betaX,-750)
betaX = betaX[nonzero]
betaY = betaY[nonzero]
#print(betaX.shape)

muX, sigmaX, amplX,intervalX = findFit(betaX,[len(betaX)/4,mean(betaX),std(betaX)],(0,[5e7,5e7,5e5]),nBins)
muY, sigmaY, amplY,intervalY = findFit(betaY,[len(betaY)/4,mean(betaY),std(betaY)],(0,[5e7,5e7,5e5]),nBins)

plt.clf()
fig = plt.figure(figsize=(15,5))
plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
s1 = fig.add_subplot(1,2,1)
s2 = fig.add_subplot(1,2,2)


fs=14
bgdbox=dict(pad=1,fc='w',ec='none')
propsR = dict(horizontalalignment="right",verticalalignment="top", backgroundcolor = 'w',bbox=bgdbox,fontsize=fs-2)

s1.hist(betaX,bins=nBins)
s1.plot(intervalX, gaussian(intervalX,amplX,muX,sigmaX), "r--", linewidth=2)
s1.set_xlim([750,1350])
if s1.get_xlim()[1] > 1700:
    #print()
    deltaX=1e-2
elif s1.get_xlim()[1] > 1200:
    deltaX=5e-3
else:
    deltaX=8e-4
s1.text(s1.get_xlim()[1]*(1-deltaX), s1.get_ylim()[1]*0.95, r"Nominal $\beta_{x}$= "+re.search("(([0-9]*)+(?=,))",name)[1]+" [m]", propsR)
    
s1.text(s1.get_xlim()[1]*(1-deltaX), s1.get_ylim()[1]*0.85,r"$\mu$="+"{:.3f} [m]".format(muX)+"\n"+r"$\sigma$="+"{:.3f}%".format(sigmaX/muX*100), propsR)
title1=r"Distribution of $\beta_x$ Values for "+re.search("(([0-9]*)+(?=Pct))",name)[1]+"% Variation"
#print(title1)
plt.setp(s1,title=title1,xlabel=r"$\beta$ [m]",ylabel="Counts [a.u.]")

s2.hist(betaY,bins=nBins)
s2.plot(intervalY, gaussian(intervalY,amplY,muY,sigmaY), "r--", linewidth=2)
s2.set_xlim([95,165])
if s2.get_xlim()[1] > 200:
    #print()
    deltaY=1e-2
elif s2.get_xlim()[1] > 150:
    deltaY=5e-3
else:
    deltaY=8e-4
s2.text(s2.get_xlim()[1]*(1-deltaY), s2.get_ylim()[1]*0.95, r"Nominal $\beta_{y}$= "+re.search("(([0-9]*)+(?=m_))",name)[1]+" [m]", propsR)
s2.text(s2.get_xlim()[1]*(1-deltaY), s2.get_ylim()[1]*0.85,r"$\mu$="+"{:.3f} [m]".format(muY)+"\n"+r"$\sigma$="+"{:.3f}%".format(sigmaY/muY*100), propsR)
title2=r"Distribution of $\beta_y$ Values for "+re.search("(([0-9]*)+(?=Pct))",name)[1]+"% Variation"
plt.setp(s2,title=title2,xlabel=r"$\beta$ [m]",ylabel="Counts [a.u.]")
plt.tight_layout()
plt.savefig("../../../Pictures/"+name+"Spread.png",bbox_inches='tight')
print("../../../Pictures/"+name+"Spread.png")

