#fitBetas.py
#python3 fitBetas.py HEBT-A2T_100pctField_1.0e-04Jitter_250x
#
import csv,sys,re
from plotFit import findFit,gaussian
from numpy import zeros,mean,std,greater
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

name=sys.argv[1]
from plotFit import numLines

if re.search("beta",name): #from plotFit getTwiss
    path = "../../../Pictures/"
    n = numLines(path+name)
    bX = 1
    bY = 2
elif re.search("Jitter",name): #from OpenXAL
    path = "../../../../OpenXAL/OXALNotebooks/failureTwiss/"
    n = numLines(path+name)
    bX = 2
    bY = 5

betaX=zeros(n)
betaY=zeros(n)

with open(path+name+".csv",mode='r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    next(csv_reader) #skips header
    rowNum = 0
    for row in csv_reader:
        #print(sample,rowNum,row[0],row[1],row[2])
        betaX[rowNum] = float(row[bX])
        betaY[rowNum] = float(row[bY])
        #print(rowNum,"Read in Betas:",betaX[rowNum],betaY[rowNum])
        rowNum+=1
        if rowNum == 200: break
#print(mean(betaX))
#print(betaX,"\n\n",betaY)
#Set Bins by %
if re.search("Jitter",name):
    if int(re.search("(([0-9]*)+(?=Jitter))",name)[1]) == 4:
        pct = 1
        nBins = 11
    elif int(re.search("(([0-9]*)+(?=Jitter))",name)[1]) == 3:
        pct = 10
        nBins = 12

if re.search("beta",name):
    if int(re.search("(([0-9]*)+(?=Pct))",name)[1]) == 1:
        pct = 1
        nBins = 10
    elif int(re.search("(([0-9]*)+(?=Pct))",name)[1]) == 10:
        pct = 10
        nBins = 30
nomBX = 1085.63
nomBY = 136.06

nonzero = greater(betaX,0)
betaX = betaX[nonzero]
betaY = betaY[nonzero]
#print(betaX.shape)

muX, sigmaX, amplX,intervalX = findFit(betaX,[len(betaX)/4,mean(betaX),std(betaX)],(0,[3.7e2,5e4,5e4]),nBins)
muY, sigmaY, amplY,intervalY = findFit(betaY,[len(betaY)/4,mean(betaY),std(betaY)],(0,[3.7e2,5e4,5e4]),nBins)

plt.clf()
fig = plt.figure(figsize=(15,5))
plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
s1 = fig.add_subplot(1,2,1)
s2 = fig.add_subplot(1,2,2)

fs=16
bgdbox=dict(pad=1,fc='w',ec='none')
propsR = dict(horizontalalignment="right",verticalalignment="top", backgroundcolor = 'w',bbox=bgdbox,fontsize=fs-3)

s1.hist(betaX,bins=nBins)
s1.plot(intervalX, gaussian(intervalX,amplX,muX,sigmaX), "r--", linewidth=2)
s1.set_xlim([980,1230])#750,1350
if s1.get_xlim()[1] > 1700:
    #print()
    deltaX=1e-2
elif s1.get_xlim()[1] > 1200:
    deltaX=2e-3
else:
    deltaX=8e-4
s1.text(s1.get_xlim()[1]*(1-deltaX), s1.get_ylim()[1]*0.98, r"Nominal $\beta_{x}$= "+str(nomBX)+" [m]", propsR)
#if pct == 1:
s1.xaxis.set_major_locator(ticker.MultipleLocator(20)) #to adjust the # of xticks to every 20
s1.text(s1.get_xlim()[1]*(1-deltaX), s1.get_ylim()[1]*0.85,r"$\mu$="+"{:.3f} [m]".format(muX)+"\n"+r"$\sigma$="+"{:.3f}[m]\n({:.3f}%)".format(sigmaX,sigmaX/muX*100), propsR)
title1=r"Distribution of $\beta_x$"+"\n for {:.0f}% QP Errors Around Nominal".format(pct)
#print(title1)
s1.set_title(title1,fontsize=fs+2)
s1.set_xlabel(r"$\beta$ [m]",fontsize=fs)
s1.set_ylabel("Counts [a.u.]",fontsize=fs)

s2.hist(betaY,bins=nBins)
s2.plot(intervalY, gaussian(intervalY,amplY,muY,sigmaY), "r--", linewidth=2)
s2.set_xlim([120,150])#95,165
if s2.get_xlim()[1] > 200:
    deltaY=1e-2
elif s2.get_xlim()[1] > 150:
    deltaY=5e-2
else:
    deltaY=2e-3
s2.text(s2.get_xlim()[1]*(1-deltaY), s2.get_ylim()[1]*0.98, r"Nominal $\beta_{y}$= "+str(nomBY)+" [m]", propsR)
#if pct == 1:
s2.xaxis.set_major_locator(ticker.MultipleLocator(2)) #to adjust the # of xticks to every 2
s2.text(s2.get_xlim()[1]*(1-deltaY), s2.get_ylim()[1]*0.85,r"$\mu$="+"{:.3f} [m]".format(muY)+"\n"+r"$\sigma$="+"{:.3f}[m]\n({:.3f}%)".format(sigmaY,sigmaY/muY*100), propsR)
title2=r"Distribution of $\beta_y$"+"\n for {:.0f}% QP Errors Around Nominal".format(pct)
s2.set_title(title2,fontsize=fs+2)
s2.set_xlabel(r"$\beta$ [m]",fontsize=fs)
s2.set_ylabel("Counts [a.u.]",fontsize=fs)
#plt.setp(s2.get_xticklabels(),rotation=45,ha="right")

#plt.tight_layout()
plt.savefig("../../../Pictures/"+name+"Spread.png",bbox_inches='tight',dpi=500)
print("../../../Pictures/"+name+"Spread.png")

