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
    zoff        = "*-1" #[mm] with preappended * to keep covar defined at z=0
    options     = {'physList':physList, 'dependence':"Twiss", 'zoff':zoff, 'initTree':False,
                    'exitTree':False, 'targetTree':True, 'MCS':True, 'engPlot':False,
                    'mat3Plot':False, 'TwissFits':False, 'PBWT': 2.25 ,'MiniRoot':False}

    #Important things
    #ind = 1 #0 is for position distribution sigma, 1 is for angular distribution sigma
    if args.t == 0:
        args.materials = ["G4_Al"] #overwrites potential user input. Needs work
    elif args.t == 0.1:
        args.materials = ["G4_Galactic"]
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
    if args.thickInd == 0:
        name = "tDependencePos_beta{:.2f},{:.2f}m_R{:.1f},{:.1f}x{:.2f}mm".format(Twiss[1],Twiss[4],args.minThick,args.maxThick,args.stepThick)
    else: 
        name = "tDependenceAng_beta{:.2f},{:.2f}m_R{:.1f},{:.1f}x{:.2f}mm".format(Twiss[1],Twiss[4],args.minThick,args.maxThick,args.stepThick)
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
            beamFile = paths['csvPWD']+"tDependence_beta{:.2f},{:.2f}m_t{:.1f}mm".format(Twiss[1],Twiss[4],options['PBWT'])
            [Jmaxes[i],pOutsideBoxes[i],beamArea,coreJMean,centX,centY,rValue,rDiff],e8SigsX,e8SigsY,targSigsX,targSigsY = runPBW(args,beamFile,Twiss,options,boxes,paths)
            #bc passing list for making smaller # of arguments, extract variables we want here
            e8TargxReal[i] = e8SigsX[args.thickInd];e8TargyReal[i] = e8SigsY[args.thickInd];targxTwissf[i]=targSigsX[args.thickInd];targyTwissf[i]=targSigsY[args.thickInd]

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

    plt.scatter(PBWT,e8TargxReal,c='green',label=r"Müller $\sigma_x$",s=65)
    plt.scatter(PBWT,e8TargyReal,c='red',label=r"Müller $\sigma_y$",s=65)
    plt.scatter(PBWT,targxTwissf,c='aqua',label=r"GEANT4 $\sigma_x$")
    plt.scatter(PBWT,targyTwissf,c='gold',label=r"GEANT4 $\sigma_y$")

    if args.thickInd == 0:
        plt.xlim((args.minThick-0.2,args.maxThick+0.2))
        plt.ylim((np.min(e8TargyReal)-5,np.max(e8TargxReal)+4))
        plt.annotate('Actual PBW Al', #up
        xy=(2.25, 14), xycoords='data',
        xytext=(-35, -20), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))
        plt.annotate('Thickness: 2.25mm',
        xy=(2.25, 10.6), xycoords='data', #down
        xytext=(-49, 15), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))
        
        plt.ylabel(r"$\sigma$ at Target [mm]",fontsize=fs)
        plt.title(r"Target Position $\sigma$ Growth"+"\nwith PBW Thickness",fontsize=fs+2)

        xlim = plt.xlim()
        ylim = plt.ylim()
        plt.text(xlim[0]+0.05,ylim[1]*0.95,"Beam Twiss at PBW:",fontsize=fs-5)
        plt.text(xlim[0]+0.05,ylim[1]*0.90,r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[\mu m]}$",fontsize=fs-5)
        plt.text(xlim[0]+0.05,ylim[1]*0.84,r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$",fontsize=fs-5)
        plt.text(xlim[0]+0.05,ylim[1]*0.78,r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]),fontsize=fs-5)
    
    else:
        plt.ylabel(r"$\sigma$ at Target [mrad]",fontsize=fs)
        plt.title(r"Target Anglular $\sigma$ Growth"+"\nwith PBW Thickness",fontsize=fs+2)
        #plt.vlines(2.25,plt.ylim()[0],plt.ylim()[1],color="m")
        #plt.text(2.3,plt.ylim()[0]+0.1,"Actual PBW Al Thickness: 2.25mm",ha="right")
        plt.annotate('   Actual PBW Al\nThickness: 2.25mm', #up
        xy=(2.25, .00213), xycoords='data',
        xytext=(-49, -60), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))

        xlim = plt.xlim()
        ylim = plt.ylim()
        plt.text(xlim[0]+0.05,ylim[1]*0.97,"Beam Twiss at PBW:",fontsize=fs-5)
        plt.text(xlim[0]+0.05,ylim[1]*0.94,r"$\epsilon_{Nx,Ny}$="+"{:.3f}, {:.3f}".format(Twiss[0],Twiss[3])+r"$_{[\mu m]}$",fontsize=fs-5)
        plt.text(xlim[0]+0.05,ylim[1]*0.91,r"$\beta_{x,y}$="+"{:.0f}, {:.0f}".format(Twiss[1], Twiss[4])+r"$_{[m]}$",fontsize=fs-5)
        plt.text(xlim[0]+0.05,ylim[1]*0.88,r"$\alpha_{x,y}$="+"{:.1f}, {:.1f}".format(Twiss[2],Twiss[5]),fontsize=fs-5)
    

    plt.xlabel("PBW Aluminium Thickness [mm]",fontsize=fs)

    #handles1, labels1 = plt.gca().get_legend_handles_labels()
    #print(labels1)
    #order1=[0,1,2,3]
    #([handles1[idx] for idx in order1],[labels1[idx] for idx in order1]
    plt.legend(ncol=2,loc='lower right',columnspacing=0.02,handletextpad=0.1,fontsize=fs-7)

    plt.tight_layout()
    dt = datetime.now()
    plt.savefig(thickCsvPWD+name+".png")#+dt.strftime("%H-%M-%S")
    print("saved",thickCsvPWD+name+".png")#+dt.strftime("%H-%M-%S")

    print(datetime.now()-origin)