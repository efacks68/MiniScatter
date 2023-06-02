#Script for searching Raster Magnet Amplitude space
# making 2D map of % outside box, with input from Carl Lindstrom.

def thicknessDependence(args,Twiss,sample,paths):
    from datetime import datetime
    origin = datetime.now()
    print(origin,"\n")
    import numpy as np
    import matplotlib.pyplot as plt
    import os,csv#,argparse

    #Constants for running scripts
    physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ"
    zoff        = "*-2" #[mm] with preappended * to keep covar defined at z=0
    options     = {'physList':physList, 'dependence':"Twiss", 'zoff':zoff, 'initTree':False,
                    'exitTree':False, 'targetTree':True, 'MCS':True, 'engPlot':False,
                    'mat3Plot':False, 'TwissFits':False, 'PBWT': 2.25 ,'MiniRoot':False}

    #Important things
    #ind = 1 #0 is for position distribution sigma, 1 is for angular distribution sigma
    if args.t == 0:
        args.materials = ["G4_Al"] #overwrites potential user input. Needs work
    elif args.t == 0.1:
        args.materials = ["G4_Galactic"]
    args.beamFile = ""
    args.compTargs = True
    if args.compTargs:
        options['targetTree'] = True
    #    options['exitTree'] = True
        options['MCS'] = True
    boxes = [0]#,-.25,-.375,-.45,-.50]#,0.125,0.25,0.375] #make an args for 24.11.22

    #Check if there is a CSV with the data already present. Speeds up plot modifications
    thickCsvPWD = paths['statsPWD']+"thicknessDependence/"
    
    #Define the array of different thicknesses
    PBWT = np.arange(args.minThick,args.maxThick+0.1,args.stepThick)
    options['MCS'] = True

    pOutsideBoxes = np.zeros(len(PBWT)) ; pOutsideBoxes.fill(-750) #figure out how to run N times/have spread to avg numbers
    Jmaxes = np.zeros(len(PBWT))        ; Jmaxes.fill(-750)
    e8TargxReal = np.zeros(len(PBWT))   ; e8TargxReal.fill(-750)
    e8TargyReal = np.zeros(len(PBWT))   ; e8TargyReal.fill(-750)
    targxTwissf = np.zeros(len(PBWT))   ; targxTwissf.fill(-750)
    targyTwissf = np.zeros(len(PBWT))   ; targyTwissf.fill(-750)

    #now name has beta instead of emit!
    from datetime import datetime
    if args.thickInd == 0:
        name = "tDependencePos8_beta{:.2f},{:.2f}m_R{:.1f},{:.1f}x{:.2f}mm_N{:.0e}".format(Twiss[1],Twiss[4],args.minThick,args.maxThick,args.stepThick,args.Nbeamlet)#+datetime.now().strftime("%H-%M-%S")
        a=1
    else: 
        name = "tDependenceAng8_beta{:.2f},{:.2f}m_R{:.1f},{:.1f}x{:.2f}mm_N{:.0e}".format(Twiss[1],Twiss[4],args.minThick,args.maxThick,args.stepThick,args.Nbeamlet)#+datetime.now().strftime("%H-%M-%S")
        a=1e3
    if os.path.isfile(thickCsvPWD+name+".csv"):
        print("Found data! Reading in!",thickCsvPWD,name)
        #from plotFit import numLines
        #nLines = numLines(paths['csvPWD']+name)
        #print(nLines)
        with open(thickCsvPWD+name+".csv") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            i = 0
            z=1
            next(csv_reader) #skips header
            for row in csv_reader:
                PBWT[i] = row[0]
                pOutsideBoxes[i] = row[1]
                e8TargxReal[i] = row[2]
                e8TargyReal[i] = row[3]
                targxTwissf[i] = row[4]
                targyTwissf[i] = row[5]
                i+=1
            csv_file.close()
    else:
        for i in range(len(PBWT)):
            print("\nline",i,PBWT[i])
            #Create Rastered Beam file, runARasterMaker checks if the CSV is already present

            options['PBWT'] = PBWT[i]
            from runPBW import runPBW
            beamFile = paths['csvPWD']+"tDependence_beta{:.2f},{:.2f}m_t{:.2f}mm_N{:.0e}".format(Twiss[1],Twiss[4],options['PBWT'],args.Nbeamlet)
            #args.t = PBWT[i] #if want simple, flat target of i thickness Al.
            [Jmaxes[i],pOutsideBoxes[i],beamArea,coreJMean,centX,centY,rValue,rDiff],e8SigsX,e8SigsY,targSigsX,targSigsY = runPBW(args,beamFile,Twiss,options,boxes,paths)
            #bc passing list for making smaller # of arguments, extract variables we want here
            e8TargxReal[i] = e8SigsX[args.thickInd]*a;e8TargyReal[i] = e8SigsY[args.thickInd]*a;targxTwissf[i]=targSigsX[args.thickInd]*a;targyTwissf[i]=targSigsY[args.thickInd]*a

            #rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,paths['csvPWD'],options,sample)
            ##Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
            #_,_,_,Jmaxes[i], pOutsideBoxes[i],beamArea,coreJMean,centX,centY,rValue,rDiff  = simulation(args,args.material,beamXAngle,beamYAngle,rasterBeamFile,Twiss,options,boxes,paths)
            ##Store % value for plotting
            ##Run with 0.1 thickness, which is the value that triggers the Vacuum rectangular 'PBW' for reference
            ##noPBWPOutBox = runPBW(energy,rasterBeamFile,beamType,0.1,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss,rXAmps[i],rYAmps[i],dependence)
            ##VacpOutsideBoxes[i] = noPBWPOutBox
            print(datetime.now().strftime("%H-%M-%S"))
            from plotFit import saveStats
            statsName = "thicknessStats"
            saveStats(paths['statsPWD'],Twiss,beamFile,Jmaxes[i],pOutsideBoxes[i],beamArea,coreJMean,centX,centY,rValue,rDiff,statsName,args.reBin,args)
        print("time elapsed",datetime.now()-origin)

        #Save values for future quick use
        #if args.saveCsv:
        print("Writing CSV")
        with open(thickCsvPWD+name+".csv",mode = 'w') as csv_file:
            csv_writer = csv.writer(csv_file,delimiter = ',')
            csv_writer.writerow(["thickness","POutBox PBW","Jmaxes [uA/mm2]"])
            for i in range(len(PBWT)):
                    csv_writer.writerow([PBWT[i],pOutsideBoxes[i],e8TargxReal[i],e8TargyReal[i],targxTwissf[i],targyTwissf[i]])#,VacpOutsideBoxes[i]])
            csv_file.close()
        print("CSV written",thickCsvPWD+name)
        print("time elapsed",datetime.now()-origin)

    from plotFit import getMoments
    #particle characteristic values
    if args.particle == "proton":
        partA = 938.27209 #[MeV/c2]
        partZ = 1
    elif args.particle == "electron":
        partA = 0.511 #[MeV/c2]
        partZ = 1
    gamma_rel = 1 + args.energy/partA #from PrintTwissParameters
    beta_rel = np.sqrt(gamma_rel*gamma_rel -1 )/gamma_rel
    #Get rid of Normalized Emittance!
    Igemtx = Twiss[0]/(beta_rel*gamma_rel)
    Igemty = Twiss[3]/(beta_rel*gamma_rel)
    #[beta,alpha,gemt,gamma]
    [sigXPBW,_] = getMoments([Twiss[1],Twiss[2],Igemtx,0])
    [sigYPBW,_] = getMoments([Twiss[4],Twiss[5],Igemty,0])

    nonEmpty = np.greater(pOutsideBoxes,-750) #remove unused elements
    pOutsideBoxes = pOutsideBoxes[nonEmpty]
    Jmaxes = Jmaxes[nonEmpty]
    e8TargxReal = e8TargxReal[nonEmpty]
    e8TargyReal = e8TargyReal[nonEmpty]
    targxTwissf = targxTwissf[nonEmpty]
    targyTwissf = targyTwissf[nonEmpty]

    print("For a beamlet",Twiss,"\nThick, % Out Box,\t Müller X, \tMüller Y,\t   Distrib X, \t    Distrib Y")
    for i in range(len(PBWT)):
        print(PBWT[i],",    {:.3f}".format(pOutsideBoxes[i]),",  ",e8TargxReal[i],e8TargyReal[i],targxTwissf[i],targyTwissf[i])

    #Plot for parameter search analysis
    fs=18
    fig= plt.figure(figsize=(17, 6))
    plt.subplots_adjust(wspace=0.25)
    ax = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    ax.scatter(PBWT,e8TargxReal,edgecolor='orange',label=r"Müller $\sigma_x$",facecolors="None",s=100,zorder=30)
    ax2.scatter(PBWT,e8TargyReal,color='orange',label=r"Müller $\sigma_y$",s=100,zorder=10)

    ax.scatter(PBWT,targxTwissf,edgecolor='b',label=r"GEANT4 $\sigma_x$",marker="s",facecolors="None",zorder=30)
    ax2.scatter(PBWT,targyTwissf,color='b',label=r"GEANT4 $\sigma_y$",marker="s",zorder=10)
    a = 0.35
    b = 0.7
    ind=0
    if args.thickInd == 0:
        #ax.set_xlim((args.minThick-0.2,args.maxThick+0.2))
        #ax.set_ylim((np.min(e8TargxReal)-1,np.max(targxTwissf)+1))
        #ax2.set_xlim((args.minThick-0.2,args.maxThick+0.2))
        #ax2.set_ylim((np.min(e8TargyReal)-1,np.max(e8TargyReal)+1))
        for idx, val in np.ndenumerate(PBWT):
            #print(idx,val)
            if round(val*1e2) / 1e2 == 2.25:
                ind = idx
                #print("ind=",idx)
        ax.vlines(2.25,targxTwissf[ind]-a-0.15,targxTwissf[ind]+a,color="m",zorder=-30)
        ax.text(1.8,13.5,"     Actual PBW Al\nThickness: 2.25mm\nGEANT4/Müller={:.3f}%".format((targxTwissf[ind]/e8TargxReal[ind]-1)*100),fontsize=fs-3)
        ax2.vlines(2.25,targyTwissf[ind]-b,targyTwissf[ind]+b,color="m",zorder=0)
        ax2.text(1.8,7,"     Actual PBW Al\nThickness: 2.25mm\nGEANT4/Müller={:.3f}%".format((targyTwissf[ind]/e8TargyReal[ind]-1)*100),fontsize=fs-3)
        """plt.annotate('Actual PBW Al', #up
        xy=(2.25, 11.6), xycoords='data',
        xytext=(-35, -25), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))
        plt.annotate('Thickness: 2.25mm',
        xy=(2.25, 9), xycoords='data', #down
        xytext=(-49, 20), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))"""
        
        ax.set_ylabel(r"$\sigma_x$ at Target [mm]",fontsize=fs,labelpad=25)
        ax.set_title(r"Target Position $\sigma_x$ Growth"+"\nwith PBW Thickness",fontsize=fs+2)
        ax2.set_ylabel(r"$\sigma_y$ at Target [mm]",fontsize=fs)
        ax2.set_title(r"Target Position $\sigma_y$ Growth"+"\nwith PBW Thickness",fontsize=fs+2)

        #xlim = plt.xlim()
        #ylim = plt.ylim()
        ax.text(0.01,0.95,"Beam Twiss at PBW:",fontsize=fs-5, transform=ax.transAxes)
        ax.text(0.01,0.90,r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[\mu m]}$",fontsize=fs-5, transform=ax.transAxes)
        ax.text(0.01,0.84,r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$",fontsize=fs-5, transform=ax.transAxes)
        ax.text(0.01,0.78,r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]),fontsize=fs-5, transform=ax.transAxes)
        ax2.text(0.01,0.95,"Beam Twiss at PBW:",fontsize=fs-5, transform=ax2.transAxes)
        ax2.text(0.01,0.90,r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[\mu m]}$",fontsize=fs-5, transform=ax2.transAxes)
        ax2.text(0.01,0.84,r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$",fontsize=fs-5, transform=ax2.transAxes)
        ax2.text(0.01,0.78,r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]),fontsize=fs-5, transform=ax2.transAxes)
    
    else:
        ax.set_ylabel(r"Angular $\sigma_x$ at Target [$\mu$rad]",fontsize=fs)
        ax.set_title(r"Target Anglular $\sigma_x$ Growth"+"\nwith PBW Thickness",fontsize=fs+2)
        ax.vlines(2.25,targxTwissf[ind]-a-0.2,targxTwissf[ind]+a,color="m")
        ax.text(1.8,ax.get_ylim()[1]*0.65,"     Actual PBW Al\nThickness: 2.25mm\nMüller/GEANT4={:.3f}%".format((e8TargxReal[ind]/targxTwissf[ind]-1)*100),fontsize=fs-3)
        #plt.annotate('   Actual PBW Al\nThickness: 2.25mm', #up
        #xy=(2.25, .00213), xycoords='data',
        #xytext=(-49, -60), textcoords='offset points',
        #arrowprops=dict(arrowstyle="->"))
        ax2.set_ylabel(r"Angular $\sigma_y$ at Target [$\mu$rad]",fontsize=fs)
        ax2.set_title(r"Target Anglular $\sigma_y$ Growth"+"\nwith PBW Thickness",fontsize=fs+2)
        ax2.vlines(2.25,targyTwissf[ind]-b,targyTwissf[ind]+b,color="m")
        ax2.text(1.8,ax2.get_ylim()[1]*0.65,"     Actual PBW Al\nThickness: 2.25mm\nMüller/GEANT4={:.3f}%".format((targyTwissf[ind]/e8TargyReal[ind]-1)*100),fontsize=fs-3)

        #xlim = plt.xlim()
        #ylim = plt.ylim()
        ax.text(0.01,0.95,"Beam Twiss at PBW:",fontsize=fs-5, transform=ax.transAxes)
        ax.text(0.01,0.90,r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[\mu m]}$",fontsize=fs-5, transform=ax.transAxes)
        ax.text(0.01,0.84,r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$",fontsize=fs-5, transform=ax.transAxes)
        ax.text(0.01,0.78,r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]),fontsize=fs-5, transform=ax.transAxes)
        ax2.text(0.01,0.95,"Beam Twiss at PBW:",fontsize=fs-5, transform=ax2.transAxes)
        ax2.text(0.01,0.90,r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[\mu m]}$",fontsize=fs-5, transform=ax2.transAxes)
        ax2.text(0.01,0.84,r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$",fontsize=fs-5, transform=ax2.transAxes)
        ax2.text(0.01,0.78,r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]),fontsize=fs-5, transform=ax2.transAxes)
    

    ax.set_xlabel("PBW Aluminium Thickness [mm]",fontsize=fs)
    ax2.set_xlabel("PBW Aluminium Thickness [mm]",fontsize=fs)

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    #([handles1[idx] for idx in order1],[labels1[idx] for idx in order1]
    #ax.legend(ncol=2,loc='lower right',columnspacing=0.02,handletextpad=0.1,fontsize=fs-7)
    
    ax.legend(lines, labels, ncol=2,loc='lower right',columnspacing=0.02,handletextpad=0.1,fontsize=fs-7)
    ax2.legend(lines2, labels2, ncol=2,loc='lower right',columnspacing=0.02,handletextpad=0.1,fontsize=fs-7)

    plt.tight_layout()
    dt = datetime.now()
    plt.savefig(thickCsvPWD+name+".png",dpi=1000)#+dt.strftime("%H-%M-%S")
    print("saved",thickCsvPWD+name+".png")#+dt.strftime("%H-%M-%S")

    print(datetime.now()-origin)