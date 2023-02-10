#Script for searching Raster Magnet Amplitude space
# making 2D map of % outside box, with input from Carl Lindstrom.



#Command Line arguments for save control
#parser = argparse.ArgumentParser()
#parser.add_argument("--beamClass", type=str,    default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
#parser.add_argument("--l",         type=str,    help="Load Particles or not",   default="PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03")
#parser.add_argument("--twiss",     type=float,  nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
#parser.add_argument("--t",         type=float,  default=0,     help="PBW Thickness [mm], 0=>MagnetPBW, 0.1 = Vacuum, >0.1 => solid Al Xmm thick. Default=0")
#parser.add_argument("--energy",    type=float,  default=570,   help="Beam Energy [MeV]. Default=570")
#parser.add_argument("--Nb",        type=int,    default=10,    help="Number of macroparticles per beamlet. Default=10")
#parser.add_argument("--nP",        type=float,  default=1e2,   help="Numper of beamlets in pulse. Default=1e3")
#parser.add_argument("--rX",        type=float,  default=0,     help="X distance from beam axis [mm]. Default=0")
#parser.add_argument("--rY",        type=float,  default=0,     help="Y distance from beam axis [mm]. Default=0")
#parser.add_argument("--aX",        type=float,  default=54.65, help="RM X Amplitude [mm]. Default=54.65")
#parser.add_argument("--aY",        type=float,  default=18.37, help="RM Y Amplitude [mm]. Default=18.37")
#parser.add_argument("--failure",   type=float,  default=0,     choices = range(0,5),  help="Which RM Failure case, 0-4. Default=0")
#parser.add_argument("--magFails",  type=int,    default=2,     choices = range(0,5),  help="Number of Raster Magnets that fail, 1-4. Default=2")
#parser.add_argument("--xlim",      type=float,  default=150,   help="+/- value for horizontal axis of output rastered image [mm]. Default=150")
#parser.add_argument("--ylim",      type=float,  default=150,   help="+/- value for vertical axis of output rastered image [mm]. Default=150")
#parser.add_argument("--maxim",     type=float,  default=0  ,   help="Maximum current density value for output rastered imagem[uA/cm^2]. Default=0")
#parser.add_argument("--edgeRaster",action="store_true",  help="Only populate edges of raster. Default=False")
#parser.add_argument("--PBIP",      action="store_true",  default=False,   help="Is PBIP present? Default=False")
#parser.add_argument("--noText",    action="store_true",  default=False,   help="Turns off printed text when called. Default=False")
#parser.add_argument("--noBox",     action="store_true",  default=False,   help="Turns off printed box when called. Default=False")
#parser.add_argument("--savePics",  action="store_true",  default=False,   help="Saves Rastered Image. Default=False")
#parser.add_argument("--saveGrads", action="store_true",  default=False,   help="Plots gradients of beam at Target. Default=False")
#parser.add_argument("--saveEdges", action="store_true",  default=False,   help="Plots edges on Raster Image. Default=False")
#parser.add_argument("--gaussFit",  action="store_true",  default=False,   help="Computes sum of Gaussian fits for central axis projection. Default=False")
#parser.add_argument("--saveFits",  action="store_true",  default=False,   help="Saves plots of Gaussian Fitting. Default=False")
#parser.add_argument("--saveHist",  action="store_true",  default=False,   help="Saves Histogram of proton density at target. Default=False")
#parser.add_argument("--saveRaster",action="store_true",  default=False,   help="Saves plot of rastered beam. Default=False")

#parser.add_argument("--ampl",   type=str,     default='map', help="Range of amplitudes to use: short(nominal-10% less) or large(nominal-70% less)")
#parser.add_argument("--eX",     type=int,     default=55,    help="End ampl X")
#parser.add_argument("--eY",     type=int,     default=20,    help="End ampl Y")
#parser.add_argument("--startX", type=int,     default=30,    help="Start ampl for X")
#parser.add_argument("--startY", type=int,     default=10,    help="Start ampl for Y")
#parser.add_argument("--NstepX", type=int,     default=6,     help="N steps for X")
#parser.add_argument("--NstepY", type=int,     default=6,     help="N steps for Y")
#args = parser.parse_args()
#print(args)

def mapRADependence(args,Twiss):
    from datetime import datetime
    origin = datetime.now()
    print(origin,"\n")
    import numpy as np
    import matplotlib.pyplot as plt
    import os,csv#,argparse
    from runARasterMaker import runARasterMaker
    from simulation import simulation

    #Constants for running scripts
    physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ"
    zoff = "*-10" #[mm] with preappended * to keep covar defined at z=0
    options     = {'physList':physList, 'dependence':"Twiss", 'zoff':zoff }

    #options     = {'noText':args.noText, 'noBox':args.noBox, 'wide':True, 'physList':physList, 'dependence':"RA",
    #                            'xlim':args.xlim, 'ylim':args.ylim, 'maxim':args.maxim, 'saveHist':args.saveHist,
    #                            'PBIP':args.PBIP, 'beamClass':args.beamClass, 'Nb':args.Nb, 'failure':args.failure,
    #                            'magFails':args.magFails, 'saveRaster':args.saveRaster, 'saveFits':args.saveFits,
    #                            'saveGrads':args.saveGrads, 'saveEdges':args.saveEdges, 'gaussFit':args.gaussFit }

    #Important things
    if args.t == 0:
        materials = ["G4_Al"] #overwrites potential user input. Needs work
    elif args.t == 0.1:
        materials = ["G4_Galactic"]
    boxes = [0]#,-.25,-.375,-.45,-.50]#,0.125,0.25,0.375] #make an args for 24.11.22

    if os.uname()[1] == "tensor.uio.no":
        csvPWD = "/scratch2/ericdf/PBWScatter/CSVs/"
    elif os.uname()[1] == "mbef-xps-13-9300":
        csvPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/"
    else: print("Help! Unknown build directory!, scatterPBW.py l 61")

    #For the input amplitude range selection, 'short' or 'long'
    amplRatio = (args.aX + 0.001) / (args.aY + 0.001) #so no /zero
    defaultRMAmplX = 54.65
    defaultRMAmplY = 18.37
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
        rXAmps = np.linspace(args.startX,args.eX,args.NstepX)
        rYAmps = np.linspace(args.startY,args.eY,args.NstepY)
        rXRange = args.eX-args.startX
        #print(start,end,step,"\n",rXAmps,"\n",rYAmps)
        print("there are ",args.NstepX*args.NstepY,"points to plot. Expect that number of minutes.")

    pOutsideBoxes = np.zeros([len(rYAmps),len(rXAmps)])
    Jmaxes = np.zeros([len(rYAmps),len(rXAmps)])
    #coreMeans = np.zeros([len(rYAmps),len(rXAmps)])
    #VacpOutsideBoxes = np.zeros(len(rXAmps))

    #Check if there is a CSV with the data already present. Speeds up plot modifications
    mapCsvPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/rAmplDependence/2DMap/"
    statsPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"

    #now name has beta instead of emit!
    name = "RasterAmplDependence_POutBox,Imax_bX{:.0f}m_R{:.1f},{:.1f}mm".format(Twiss[1],rXRange,args.eY-args.startY)
    if os.path.isfile(mapCsvPWD+name+".csv"):
        print("Found data! Reading in!",name)
        #from plotFit import numLines
        #nLines = numLines(csvPWD+name)
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
                    Jmaxes[j-z][i-z] = row[3]
                    #coreMeans[j-z][i-z] = row[4]
                    #VacpOutsideBoxes[line_count-z] = row[3]
                    #print(i-z,j-z,rXAmps[i-z],rYAmps[j-z],pOutsideBoxes[i-z][j-z])
                    j += 1
                    if j == len(rXAmps)+1:
                        i += 1
                        j = 0 + z
            csv_file.close()
    else:
        print("Run simulations for RA ratio",amplRatio,"\t",name)
        for i in range(len(rXAmps)):
            for j in range(len(rYAmps)):
                print("\nline [",i,",",j,"]",rXAmps[i],rYAmps[j])
                #Create Rastered Beam file, runARasterMaker checks if the CSV is already present
                rasterBeamFile, beamXAngle, beamYAngle = runARasterMaker(args,Twiss,csvPWD,options)
                #Send raster beam file to runPBW which simulates with MiniScatter or opens already run data. Full PBW model
                _,_,_,Jmaxes[j][i], pOutsideBoxes[j][i], dispY,dispX,rValue = simulation(args,args.material,False,beamXAngle,beamYAngle,rasterBeamFile,Twiss,options,boxes)
                #Store % value for plotting
                #Run with 0.1 thickness, which is the value that triggers the Vacuum rectangular 'PBW' for reference
                #noPBWPOutBox = runPBW(energy,rasterBeamFile,beamType,0.1,beamXAngle,beamYAngle,PBIP,args.savePics,physList,Twiss,rXAmps[i],rYAmps[i],dependence)
                #VacpOutsideBoxes[i] = noPBWPOutBox
                print(datetime.now().strftime("%H-%M-%S"))
                from plotFit import saveStats
                saveStats(statsPWD,rasterBeamFile,Jmaxes[j][i],pOutsideBoxes[j][i],dispY,dispX,rValue)
            print("time elapsed",datetime.now()-origin)

        #Save values for future quick use
        #if args.saveCsv:
        print("Writing CSV")
        with open(mapCsvPWD+name+".csv",mode = 'w') as csv_file:
            csv_writer = csv.writer(csv_file,delimiter = ',')
            csv_writer.writerow(["Raster X Amplitude","Raster Y Amplitude","POutBox PBW","Jmaxes [uA/mm2]"])
            for i in range(len(rXAmps)):
                for j in range(len(rYAmps)):
                    csv_writer.writerow([rXAmps[i],rYAmps[j],pOutsideBoxes[j][i],Jmaxes[j][i]])#,VacpOutsideBoxes[i]])
            csv_file.close()
        print("CSV written",mapCsvPWD+name)
        print("time elapsed",datetime.now()-origin)

    for i in range(len(rXAmps)):
        for j in range(len(rYAmps)):
            print(rXAmps[i],rYAmps[j],pOutsideBoxes[j][i],Jmaxes[j][i])

    #Plot for parameter search analysis
    fs=14
    minim = pOutsideBoxes.min()+0.01
    maxim = pOutsideBoxes.max()*1.1
    plotRange = maxim-minim
    print(minim,maxim,plotRange)
    picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/rAmplDependence/2DMap/"
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
    cbarVals  = [minim,minim+plotRange*0.1,minim+plotRange*0.4,minim+plotRange*0.7,maxim] #make array for color bar values
    cbarLabels = ["{:.1f}".format(cbarVals[0]),"{:.1f}".format(cbarVals[1]),"{:.1f}".format(cbarVals[2]),
                "{:.1f}".format(cbarVals[3]),"{:.1f}".format(cbarVals[4])]#,"{:.1f}".format(cbarVals[5])] #make labels of Value
    cbarLabel = "% Outside Box"
    cbar = fig.colorbar(c, ax=ax1,pad=0.01,ticks=cbarVals)
    cbar.set_label(cbarLabel,labelpad=2,fontsize=fs-2)
    #cbar.set_ticks(cbarVals)
    cbar.set_ticklabels(cbarLabels)
    print("CbarVals",cbarVals)

    ##Max I Plot
    minimMax = Jmaxes.min()*0.9999
    maximMax = Jmaxes.max()*1.0001
    plotRangeMax = maximMax-minimMax
    print(Jmaxes.min(),Jmaxes.max())
    d = ax2.pcolor(X,Y,Jmaxes,norm=LogNorm(vmin=minimMax, vmax=maximMax),shading='auto',cmap='viridis')
    ax2.set_xlabel("Horizontal Raster Amplitude [mm]",fontsize=fs)
    ax2.set_ylabel("Vertical Raster Amplitude [mm]",fontsize=fs)
    ax2.set_title("Peak Current Density on Target\nwith Raster Amplitude at PBW",fontsize = fs+2)
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
    #22Jan-you are adjusting the nominal markers, then runnning several different maps
    #try to find an optimal ampl.
    #also, you have to adjust the original values in the Gaussian fitting as it keeps throwing 
    #  the "lower/upper bounds outside current parameter value." error which is about the initial
    #  value being outside the SetParamterLimit bounds

    ax1.hlines(defaultRMAmplY,floor(defaultRMAmplX),ceil(defaultRMAmplX),color='o')
    ax1.vlines(defaultRMAmplX,defaultRMAmplY-(ylim1[1]-ylim1[0])*0.02,defaultRMAmplY+(ylim1[1]-ylim1[0])*0.02,color='o')
    ax1.text(defaultRMAmplX+1,defaultRMAmplY+1,"Nominal Y, Nominal X",dict(horizontalalignment="right",verticalalignment="bottom"))
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
    plt.savefig(picPWD+name+".png")#+dt.strftime("%H-%M-%S")
    print("saved",picPWD+name+".png")#+dt.strftime("%H-%M-%S")

    print(datetime.now()-origin)