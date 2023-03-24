#Script for searching Raster Magnet Amplitude space
# making 2D map of % outside box, with input from Carl Lindstrom.

def mapRADependence(args,Twiss,sample,paths,origBx,origBY):
    from datetime import datetime
    origin = datetime.now()
    print(origin,"\n")
    import numpy as np
    import matplotlib.pyplot as plt
    import os,csv#,argparse

    #Constants for running scripts
    physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ"
    zoff = "*-10" #[mm] with preappended * to keep covar defined at z=0
    options     = {'physList':args.physList, 'dependence':"RA", 'zoff':zoff, 'MiniRoot':True  }

    #Important things
    if args.t == 0:
        materials = ["G4_Al"] #overwrites potential user input. Needs work
    elif args.t == 0.1:
        materials = ["G4_Galactic"]
    boxes = [0]#,-.25,-.375,-.45,-.50]#,0.125,0.25,0.375] #make an args for 24.11.22

    #For the input amplitude range selection, 'short' or 'long'
    amplRatio = (args.aX + 0.001) / (args.aY + 0.001) #so no /zero
    if args.ampl == 's':
        #rXAmps = np.array([args.aX,49,50,51,52,53,54])
        rXAmps = np.array([args.aX,49,51,53,57])
        #Make Y amplitude values based on ratio (necessary?)
        rYAmps = np.zeros(len(rXAmps))
        for i in range(len(rXAmps)):
            rYAmps[i] = rXAmps[i]/amplRatio
        rXRange = rXAmps.max() - rXAmps.min()
        print(rXRange)
        legloc = "center right"
        twTloc = 0.55
        pLloc = 0.95
    elif args.ampl == 'l':
        #rXAmps = np.array([args.aX,15,25,30,35,40,45,50])
        rXAmps = np.array([args.aX,0.001,15,25,35])
        #Make Y amplitude values based on ratio (necessary?)
        rYAmps = np.zeros(len(rXAmps))
        for i in range(len(rXAmps)):
            rYAmps[i] = rXAmps[i]/amplRatio
        rXRange = rXAmps.max() - rXAmps.min()
        print(rXRange)
        legloc = "center left"
        twTloc = 0.85
        pLloc = 0.85
    elif args.ampl == 'srcy': #Short Range Constant Y amplitude
        rXAmps = np.array([args.aX,49,50,51,52,53,54])
        rYAmps = args.aY*0.9 * np.ones(len(rXAmps)) #90% amplitude constant
        rXRange = rXAmps.max() - rXAmps.min()
        print(rXRange)
        legloc = "center right"
        twTloc = 0.55
    elif args.ampl =="map":
        rXAmps = np.linspace(args.startX,args.eX,args.Nstep)#X
        rYAmps = np.linspace(args.startY,args.eY,args.Nstep)#Y
        rXRange = args.eX-args.startX
        print(args.startX,args.eX,args.Nstep,"\n",rXAmps,"\n",rYAmps)
        print("there are ",args.Nstep*args.Nstep,"points to plot. Expect that number/6 hours.")
        #print("there are ",args.NstepX*args.NstepY,"points to plot. Expect that number of minutes.")

    print(args.aX,args.aY)
    pOutsideBoxes = np.zeros([len(rYAmps),len(rXAmps)]) #figure out how to run N times/have spread to avg numbers
    coreJMeans = np.zeros([len(rYAmps),len(rXAmps)])
    #coreMeans = np.zeros([len(rYAmps),len(rXAmps)])
    #VacpOutsideBoxes = np.zeros(len(rXAmps))

    #Check if there is a CSV with the data already present. Speeds up plot modifications
    mapCsvPWD = paths['statsPWD']+"rAmplDependence/2DMap/"

    #now name has beta instead of emit!
    nPs = round((2*10+round(0.04*1/14*1e6)/1e3)*args.nP)*args.Nb
    name = "RasterAmplDependence_POutBox,coreJMean_nP{:.1e}_bX{:.0f},{:.0f}m_R{:.1f},{:.1f}mm_{:.0f}steps".format(nPs,Twiss[1],Twiss[4],rXRange,args.eY-args.startY,args.Nstep)
    if os.path.isfile(mapCsvPWD+name+".csv"):
        print("Found data! Reading in!",mapCsvPWD+name)
        #from plotFit import numLines
        #nLines = numLines(paths['csvPWD']+name)
        #print(nLines)
        with open(mapCsvPWD+name+".csv") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            i = 0
            j = 0
            k = 0
            z=1
            for row in csv_reader:
                if i == 0 and j ==0: #Header line
                    i += 1
                    j += 1
                else:
                    rXAmps[i-z] = row[0]
                    rYAmps[j-z] = row[1]
                    pOutsideBoxes[j-z][i-z] = row[2]
                    coreJMeans[j-z][i-z] = row[3]
                    #coreMeans[j-z][i-z] = row[4]
                    #VacpOutsideBoxes[line_count-z] = row[3]
                    #print(i-z,j-z,rXAmps[i-z],rYAmps[j-z],pOutsideBoxes[i-z][j-z])
                    j += 1
                    if j == len(rXAmps)+1:
                        i += 1
                        j = 0 + z
            csv_file.close()
    else:
        #print("Run simulations for RA ratio",amplRatio,"\t",name)
        for i in range(len(rXAmps)):
            for j in range(len(rYAmps)):
                print("\nline [",i,",",j,"]",rXAmps[i],rYAmps[j],"\t",datetime.now())
                #Create Rastered Beam file, runARasterMaker checks if the CSV is already present
                args.aX = rXAmps[i]
                args.aY = rYAmps[j]
                from runARasterMaker import runARasterMaker
                rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,paths['csvPWD'],options,sample,origBx,origBY)

                from simulation import setup
                #def                        setup(args,material,           beamFile,Twiss,options,paths):
                savename,simSetup_simple1 = setup(args,args.material,rasterBeamFile,Twiss,options,paths)

                import miniScatterDriver
                from plotFit import rasterImage

                TRYLOAD = True  #Try to load already existing data instead of recomputing, only if using getData_TryLoad function.
                (twiss_PBW, numPart_PBW, objects_PBW) = miniScatterDriver.getData_tryLoad(simSetup_simple1, tryload=TRYLOAD,getObjects=["tracker_cutoff_xy_PDG2212","init_xy"])
                [Jmax,pOutsideBoxes[j][i],beamArea,coreJMeans[j][i],centX,centY,rValue,chi2] = rasterImage(savename,"Target",objects_PBW["tracker_cutoff_xy_PDG2212"],simSetup_simple1["N"],args,Twiss,options,boxes,paths,sample)

                #rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,paths['csvPWD'],options,sample)
                ##Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
                #_,_,_,Jmaxes[j][i], pOutsideBoxes[j][i],beamArea,coreJMean,centX,centY,rValue,rDiff  = simulation(args,args.material,beamXAngle,beamYAngle,rasterBeamFile,Twiss,options,boxes,paths)
                ##Store % value for plotting
                ##Run with 0.1 thickness, which is the value that triggers the Vacuum rectangular 'PBW' for reference
                ##noPBWPOutBox = runPBW(energy,rasterBeamFile,beamType,0.1,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss,rXAmps[i],rYAmps[i],dependence)
                ##VacpOutsideBoxes[i] = noPBWPOutBox
                print(datetime.now().strftime("%H-%M-%S"))
                from plotFit import saveStats
                saveStats(paths['statsPWD'],Twiss,rasterBeamFile,Jmax,pOutsideBoxes[j][i],beamArea,coreJMeans[j][i],centX,centY,rValue,chi2,"MapRMA",args.reBin,args)
            print("time elapsed",datetime.now()-origin)

        #Save values for future quick use
        #if args.saveCsv:
        print("Writing CSV")
        with open(mapCsvPWD+name+".csv",mode = 'w') as csv_file:
            csv_writer = csv.writer(csv_file,delimiter = ',')
            csv_writer.writerow(["Raster X Amplitude","Raster Y Amplitude","POutBox PBW","Core J Mean [uA/mm2]"])
            for i in range(len(rXAmps)):
                for j in range(len(rYAmps)):
                    csv_writer.writerow([rXAmps[i],rYAmps[j],pOutsideBoxes[j][i],coreJMeans[j][i]])#,VacpOutsideBoxes[i]])
            csv_file.close()
        print("CSV written",mapCsvPWD+name)
        print("time elapsed",datetime.now()-origin)

    for i in range(len(rXAmps)):
        for j in range(len(rYAmps)):
            print(rXAmps[i],rYAmps[j],pOutsideBoxes[j][i],coreJMeans[j][i])

    #Plot for parameter search analysis
    fs=14
    minim = pOutsideBoxes.min()*0.9999
    maxim = pOutsideBoxes.max()*1.05
    plotRange = maxim-minim
    print(minim,maxim,plotRange)
    plt.close()
    plt.clf()
    X,Y = np.meshgrid(rXAmps,rYAmps)
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(15,6))
    plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
    from matplotlib.colors import LogNorm
    c = ax1.pcolor(X,Y,pOutsideBoxes,norm=LogNorm(vmin=minim, vmax=maxim),shading='auto',cmap='viridis')#,label=r"PBW % Outisde Box")
    ax1.set_xlabel("Horizontal Raster Amplitude [mm]",fontsize=fs)
    ax1.set_ylabel("Vertical Raster Amplitude [mm]",fontsize=fs)
    ax1.set_title("Rastered Beam % Outside Box on Target\nwith Raster Amplitude at PBW",fontsize = fs+2)
    cbarVals  = [minim,minim+plotRange*0.25,minim+plotRange*0.5,minim+plotRange*0.75,maxim] #make array for color bar values
    cbarLabels = ["{:.2f}".format(cbarVals[0]),"{:.2f}".format(cbarVals[1]),"{:.2f}".format(cbarVals[2]),
                "{:.2f}".format(cbarVals[3]),"{:.2f}".format(cbarVals[4])]#,"{:.1f}".format(cbarVals[5])] #make labels of Value
    cbarLabel = "% Outside Box"
    cbar = fig.colorbar(c, ax=ax1,pad=0.01,ticks=cbarVals)
    cbar.set_label(cbarLabel,labelpad=2,fontsize=fs-2)
    #cbar.set_ticks(cbarVals)
    cbar.set_ticklabels(cbarLabels)
    print("CbarVals",cbarVals)

    ##Max I Plot
    minimMax = coreJMeans.min()*0.9999
    maximMax = coreJMeans.max()*1.05
    plotRangeMax = maximMax-minimMax
    print(coreJMeans.min(),coreJMeans.max())
    d = ax2.pcolor(X,Y,coreJMeans,norm=LogNorm(vmin=minimMax, vmax=maximMax),shading='auto',cmap='viridis')
    ax2.set_xlabel("Horizontal Raster Amplitude [mm]",fontsize=fs)
    ax2.set_ylabel("Vertical Raster Amplitude [mm]",fontsize=fs)
    ax2.set_title("Mean Core Current Density on Target\nwith Raster Amplitude at PBW",fontsize = fs+2)
    cbarVals2  = [minimMax,minimMax+plotRangeMax*0.25,minimMax+plotRangeMax*0.5,minimMax+plotRangeMax*0.75,maximMax] #make array for color bar values
    cbarLabels2 = ["{:.1f}".format(cbarVals2[0]),"{:.1f}".format(cbarVals2[1]),"{:.1f}".format(cbarVals2[2]),
                "{:.1f}".format(cbarVals2[3]),"{:.1f}".format(cbarVals2[4])]
    cbarLabel2 = r"Current Density [$\mu$A/cm$^2$]"
    cbar2 = fig.colorbar(d, ax=ax2,pad=0.01,ticks=cbarVals2)
    cbar2.set_label(cbarLabel2,labelpad=2,fontsize=fs-2)
    #cbar2.set_ticks(cbarVals2)
    cbar2.set_ticklabels(cbarLabels2)
    from math import floor,ceil
    xlim1 = ax1.get_xlim()
    ylim1 = ax1.get_ylim()
    #22Jan23-you are adjusting the nominal markers, then runnning several different maps try to find an optimal ampl.
    #also, you have to adjust the original values in the Gaussian fitting as it keeps throwing 
    #  the "lower/upper bounds outside current parameter value." error which is about the initial value being outside the SetParamterLimit bounds

    #ax1.hlines(args.aY,floor(args.aX),ceil(args.aX),color='orange',lw=2)
    #ax1.vlines(args.aX,args.aY-(ylim1[1]-ylim1[0])*0.02,args.aY+(ylim1[1]-ylim1[0])*0.02,color='orange',lw=2)
    ax1.text(args.aX-0.3,args.aY,"Nominal",color="w",ha="right",va="center",fontsize=fs-2)
    #ax2.hlines(args.aY,floor(args.aX),ceil(args.aX),color='orange',lw=2)
    #ax2.vlines(args.aX,args.aY-(ylim1[1]-ylim1[0])*0.02,args.aY+(ylim1[1]-ylim1[0])*0.02,color='orange',lw=2)
    ax2.text(args.aX-0.3,args.aY,"Nominal",color="w",ha="right",va="center",fontsize=fs-2)

    ax1.plot(args.aX,args.aY,color="orange",marker="P",markersize=15)
    ax2.plot(args.aX,args.aY,color="orange",marker="P",markersize=15)
    for i in range(len(rXAmps)):
        for j in range(len(rYAmps)):
            ax2.text(rXAmps[i],rYAmps[j],"{:.0f}".format(coreJMeans[j][i]),color="w",ha="center",va="center",fontsize=fs+4)#,fontweight="bold")
            ax1.text(rXAmps[i],rYAmps[j],"{:.1f}".format(pOutsideBoxes[j][i]),color="w",ha="center",va="center",fontsize=fs+4)#,fontweight="bold")
    ##Set up texts to include with relevant info
    #xlim = plt.xlim()
    #ylim = plt.ylim()
    #plt.ylim([-0.7,ylim[1]])
    #ylim = plt.ylim()
    #plt.text(xlim[1]*0.98,ylim[0]+0.1,physList,fontsize = fs-5,color="k",horizontalalignment="right",verticalalignment="bottom")
    #plt.text(xlim[0]+0.2,ylim[1]*twTloc,"Beam Twiss at PBW:",fontsize=fs-4) #position Twiss print out depending on plot range
    #plt.text(xlim[0]+0.2,ylim[1]*(twTloc-0.07*1),r"$\epsilon_{Nx,Ny}$ = "+"{:.3f}, {:.3f}".format(Twiss[2],Twiss[5])+r"$_{[mm \cdot mrad]}$",fontsize=fs-4)
    #plt.text(xlim[0]+0.2,ylim[1]*(twTloc-0.07*2),r"$\beta_{x,y}$ = "+"{:.0f}, {:.0f}".format(Twiss[0], Twiss[3])+r"$_{[m]}$",fontsize=fs-4)
    #plt.text(xlim[0]+0.2,ylim[1]*(twTloc-0.07*3),r"$\alpha_{x,y}$ = "+"{:.1f}, {:.1f}".format(Twiss[1],Twiss[4]),fontsize=fs-4)
    #if args.ampl == 's' or args.ampl == 'l':
    #  plt.text(2*(xlim[1]-xlim[0])/7+xlim[0],ylim[1]*0.93,"Raster Amplitude H:V = 2.975",fontsize=fs-2)
    #plt.legend(loc=legloc)
    plt.tight_layout()
    dt = datetime.now()
    plt.savefig(mapCsvPWD+name+".png")#+dt.strftime("%H-%M-%S")
    print("saved",mapCsvPWD+name+".png")#+dt.strftime("%H-%M-%S")

    print(datetime.now()-origin)