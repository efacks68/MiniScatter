#poissonEstimator.py
#from Haavard 27.3
#python3 poissonEstimator.py 2 #or 8

import numpy,os,sys
import matplotlib.pyplot as plt

if os.uname()[1] in {"tensor.uio.no", "heplab01.uio.no", "heplab04.uio.no","heplab03.uio.no"}:
    pwd = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"

n_images = int(1e4)
j_maxes = numpy.zeros(n_images)
j_maxes2 = numpy.zeros(n_images)
means = []
sigmas = []
means2 = []
sigmas2 = []

if sys.argv[1] == "2":
    end = 5e3    #5e3,5e4
    step = 2e2   #2e2,2e3
    binx = 400   #400,100
    biny = 200    #200,50
    binS = 2     #2, 8
    locx = 2.4e3   #2.4e3, 3e4
    locy = 1.11 #1.11, 1.025
    deltax = 0.05 #0.05, 0.02
    deltay = 0.05

    counts = [250,1200,2400,4800]
    vals = [1.30,1.1275,1.0896,1.06323]
elif sys.argv[1] == "8":
    end = 7e4; step = 3e3; binx = 100; biny = 50
    binS = 8; locx = 3e4; locy = 1.025; deltax = 0.05; deltay = 0.01
    counts = [2.4e3,1.5e4,3e4,6.5e4]
    vals = [1.07,1.0301,1.02128,1.0144]
lrange = range(int(step), int(end), int(step))  # from, to, in steps of
#lrange2 = range(int(step), int(end*a*a), int(step*a*a))  # from, to, in steps of

name="n{:.0e}_lambda{:.0e}x{:.0e}_bins{:.0f}x{:.0f}".format(n_images,end,step,binx,biny)
print(name)

#for finding and reading file if there or doing calculation if not:
###
import csv
if os.path.isfile(pwd+name+".csv"):
    print("Found data! Reading in!",pwd+name+".csv")
    from plotFit import numLines
    nLines = numLines(pwd+name)
    #print(nLines)
    lrange = numpy.zeros(nLines-1)
    means = numpy.zeros(nLines-1)
    sigmas = numpy.zeros(nLines-1)
    with open(pwd+name+".csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        i = 0
        z=1
        next(csv_reader) #skips header
        for row in csv_reader:
            lrange[i] = row[0]
            means[i] = row[1]
            sigmas[i] = row[2]
            #print(lrange[i],means[i],sigmas[i])
            i+=1
        csv_file.close()
else:
    print(lrange)
    for lamb in lrange:
        for i in range(n_images):
            img = numpy.random.poisson(lam=lamb, size=(binx, biny)) # size of flat top in pixels
            j_maxes[i] = img.max()
        means.append((j_maxes.mean()/lamb))
        sigmas.append(j_maxes.std()/lamb)
        print("lamb {}".format(lamb),(j_maxes.mean()/lamb))

    print("Writing CSV")
    with open(pwd+name+".csv",mode = 'w') as csv_file:
        csv_writer = csv.writer(csv_file,delimiter = ',')
        csv_writer.writerow(["lambda","mu","sigma"])
        for i in range(len(lrange)):
                csv_writer.writerow([lrange[i],means[i],sigmas[i]])#,VacpOutsideBoxes[i]])
        csv_file.close()
    print("CSV written",name)
###

#for lamb in lrange2:
#    for i in range(n_images):
#        img2 = numpy.random.poisson(lam=lamb, size = (int(binx/a),int(biny/a)))
#        j_maxes2[i] = img2.max()
#    means2.append((j_maxes.mean()/lamb))
#    sigmas2.append(j_maxes.std()/lamb)
#    print("lamb {}".format(lamb),(j_maxes.mean()/lamb))

fs=16
fig,ax = plt.subplots()
ax.plot(lrange, means,label=r"$\mu$"+" Bins {:.0f}x{:.0f}".format(binS,binS)+r" [mm$^2$]",c='b')
#ax.plot(lrange2, means2,label=r"$\mu$"+" {:.0f}x{:.0f} Bins".format(binx/a,biny/a),c='g')
#ax.set_yscale("log")
ax2=ax.twinx()
#plt.savefig("mu_lambda{:.0e}x{:.0e}_bins{:.0f}x{:.0f}".format(end,step,binx,biny))
ax2.plot(lrange, sigmas,label=r"$\sigma$"+" Bins {:.0f}x{:.0f}".format(binS,binS)+r" [mm$^2$]",c="orange")
#ax2.plot(lrange2, sigmas2,label=r"$\sigma$"+" {:.0f}x{:.0f} Bins".format(binx/a,biny/a),c="r")
#ax2.set_yscale("log")
ax.set_title("Poisson Noise Contribution to Peak Value",fontsize=fs+2)
ax.set_xlabel(r"$\lambda$",fontsize=fs+2)
ax.set_ylabel(r"Peak / $\lambda$",fontsize=fs+2)
ax2.set_ylabel(r"$\sigma$ / $\lambda$",fontsize=fs+2)

ax.scatter(counts,vals,color='r',s=40,label=r"Simulation $\lambda$")

ax.vlines(locx,ax.get_ylim()[0],locy+deltay/2,alpha=0.5,color="g")
ax.text(locx*(1+deltax),locy,"Flat-Top Counts / Bin",ha="left",fontsize=fs-2)
plt.text(.35, .98, "{:.0e} Samples".format(n_images), ha='center', va='top', transform=ax.transAxes,fontsize=fs-2) #best thing ever.
#ax2.annotate('Flat Top Count',
#        xy=(locx,locy), xycoords='data',
#        xytext=(locx,locy), textcoords='data',
#        arrowprops=dict(arrowstyle="->"))
order = [0,2,1]
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1+lines2
labels = labels1 + labels2
#print(lines)
ax2.legend([lines[idx] for idx in order],[labels[idx] for idx in order], loc="upper right",fontsize=fs-4)
#ax2.legend(lines1 + lines2, labels + labels2, loc="upper right",fontsize=fs-4)
#plt.show()
plt.savefig(pwd+name+".png",bbox_inches='tight',dpi=500)
print(pwd+name+".png")