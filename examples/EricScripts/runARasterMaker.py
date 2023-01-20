
def runARasterMaker(energy,NperBunch,nPulses,envXatBPM94,envYatBPM94,edges,Twiss,rasterXAmplitude0,rasterYAmplitude0,csvPWD,options):
    import numpy as np
    from math import pi, asin, sin
    from datetime import datetime
    start = datetime.now()
    #distance 1e0 unit in this script is mm, so 1e3=1meter
    um = 1e-6
    mm = 1e-3

    #Twiss 
    nemtX = Twiss[0] #[mm-mrad]
    betaX = Twiss[1] #[m]
    alphX = Twiss[2]
    nemtY = Twiss[3] #[mm-mrad]
    betaY = Twiss[4] #[m]
    alphY = Twiss[5]
    #Calculate Geometric Emittance for Covariance Matrix
    partA = 938.27209 #[MeV/c2]
    partZ = 1
    gamma_rel = 1 + energy / partA #from PrintTwissParameters
    beta_rel = np.sqrt(gamma_rel * gamma_rel - 1 ) / gamma_rel
    gemtX = nemtX*um / (beta_rel * gamma_rel) #[m]
    gemtY = nemtY*um / (beta_rel * gamma_rel) #[m]

    t_pulse = round(0.04 * 1/14 * 1e6) # mus
    pulse_start = 10 #why have a delay?
    t_end = (2* pulse_start + t_pulse) /1000# - 2#3 ms
    time_length = round(t_end * nPulses ) #number of pulses, nominal = 2.86e3

    t = np.linspace(0,t_end,time_length) #array of steps of length time_length
    N_t = len(t) # number of time samples
    n_tii  = 10 #number of positions per us
    totX = np.zeros([N_t*n_tii*NperBunch,2])
    totY = np.zeros([N_t*n_tii*NperBunch,2])
    #print("t_pulse {:.2f}, t_end {:.2f}, time_length {:.2f}, N_t {:.2f}".format(t_pulse,t_end,time_length,N_t))

    ##Raster Constants
    Fx = 39.953*1e3 #[kHz]
    Fy = 28.7051*1e3 #[kHz]
    pb = 1
    dt =  np.mean(np.diff(t)) #[s]
    delta_t = np.linspace(0,dt,n_tii) #[s]
    z = -10 #[mm] for z location to generate protons at in MiniScatter

    #Assumes Protons
    covX = gemtX/mm * np.asarray([[betaX/mm,-alphX],[-alphX,(1+alphX**2)/(betaX/mm)]]) #[mm]
    covY = gemtY/mm * np.asarray([[betaY/mm,-alphY],[-alphY,(1+alphY**2)/(betaY/mm)]]) #[mm]

    #Calculate Envelope Center Angle
    dBPM93to94 = 3031 #[mm] from OpenXAL(?)
    dBPM93toPBW = 16822 #[mm] PBW 4400mm upstream of Target: https://gitlab.esss.lu.se/ess-bp/ess-lattice/-/blob/HEBT_RASTER_V29/9.0_HEBT/Beam_Physics/lattice.dat
    dBPM93toTarg = 21222 #[mm] from Synoptic Viewer https://confluence.esss.lu.se/pages/viewpage.action?pageId=222397499
    dPBWtoTarg = 4400 #[mm] from lattice and Synoptic
    envXAngle = envXatBPM94 / dBPM93to94 #x' = distance from beamline axis at BPM94, assume Cross Over at BPM93 / distance BPM 93 to 94
    envYAngle = envYatBPM94 / dBPM93to94 #y' not radians, as per Kyrre 2.11.22
    beamletXAngle = 0 #default
    beamletYAngle = 0 #default

    #Envelope Center Offset, i.e. raster centre on BEW
    envXCenterOffset = envXAngle * dBPM93toPBW #[mm] Takes into account angular drift since Cross Over
    envYCenterOffset =  envYAngle * dBPM93toPBW #[mm]
    envXCenOff = envXCenterOffset * np.ones(N_t) #[mm]
    envYCenOff = envYCenterOffset * np.ones(N_t) #[mm]

    #Since generating particles just before PBW, must scale a0 by dPBWtoTarg * envAngle = (1- dPBWtoTarg / dBPM93toTarg)
    #amplScale = 1 - dPBWtoTarg / dBPM93toTarg #double check you account for Z before PBW in MiniScatter! beam production plane in GEANT, not exact PBW center!
    amplScale = 1 #Cyrille said the RM Amplitude is already scaled
    sRasterXAmpl = rasterXAmplitude0 * amplScale * np.ones(N_t) #[mm]
    sRasterYAmpl = rasterYAmplitude0 * amplScale * np.ones(N_t) #[mm]

    #For weighting edges case (--edges argument)
    Right = 50
    Top = 17
    i=0 #for Bunch number iterator

    #Pick name based on beam; default: "PBW_570MeV_beta1007,130m_RMamp55,18mm_N2.9e+05_NpB10_NPls1e+03"
    if options['beamClass'] == "ESS" or options['beamClass'] == "Yngve":
        name = "PBW_{:.0f}MeV_beta{:.0f},{:.0f}m_RMamp{:.0f},{:.0f}mm_N{:.1e}_NpB{:.0f}_NPls{:.0e}".format(energy,betaX,betaY,rasterXAmplitude0,rasterYAmplitude0,len(totX[:,0]),NperBunch,nPulses)
    elif options['beamClass'] == "pencil":
        name = "PBW_{:.0f}MeV_pencilBeam_RMamp{:.0f},{:.0f}mm_N{:.1e}_NpB{:.0f}_NPls{:.1e}".format(energy,rasterXAmplitude0,rasterYAmplitude0,len(totX[:,0]),NperBunch,nPulses)
    elif options['beamClass'] == "Twiss":
        name = "PBW_{:.0f}MeV_eX{:.2f},eY{:.2f}um_bX{:.0f},bY{:.0f}m_aX{:.0f},aY{:.0f}_RMamp{:.0f},{:.0f}mm_N{:.1e}_NpB{:.0f}_NPls{:.0e}".format(energy,
                    nemtX,nemtY,betaX,betaY,alphX,alphY,rasterXAmplitude0,rasterYAmplitude0,len(totX[:,0]),NperBunch,nPulses)
    if envXatBPM94 != 0:
        name = name + "_X{:.0f}mrad".format(envXAngle*1e3)
    if envYatBPM94 != 0:
        name = name + "_Y{:.0f}mrad".format(envYAngle*1e3)

    #Raster Magnet Failure Options
    idcy = round(N_t * (1 - options['magFails'] / 4 )) #produces which quarter to end is 0
    if options['failure'] == 1: #Horizontal Fail
        sRasterXAmpl[idcy:] = 0 # no RM amplitude for magFails-th quarter of pulse
        name = name + "_failure1-" + str(options['magFails'])+"f"
    elif options['failure'] == 2: # Vertical Fail
        sRasterYAmpl[idcy:] = 0
        name = name + "_failure2-" + str(options['magFails'])+"f"
    elif options['failure'] == 3: # H & V Fails
        sRasterXAmpl[idcy:] = 0
        sRasterYAmpl[idcy:] = 0
        name = name + "_failure3-" + str(options['magFails'])+"f"
    elif options['failure'] == 4: # Correlated Motion
        Fy = Fx
        name = name + "_failure4-" + str(options['magFails'])+"f"
    if edges:
        name = name + "_edges"

    #Calculate periods  
    periodX = pb/Fx * np.ones(N_t) #[s] used for beamlet center calculation
    periodY = pb/Fy * np.ones(N_t) #[s]

    #If file found, don't make again!
    outname = csvPWD + name
    from os.path import isfile
    if isfile(outname+".csv"):
        print("Found: ",outname,".csv",sep="")
    else:
        print("CSV not found. Making: ",outname,".csv",sep="")
        centroids = np.zeros([N_t*n_tii,2])
        for jj in range(N_t):
            for ii in range(n_tii):
                tjj_ii = t[jj] + delta_t[ii]

                #Calculate the Raster Magnet contribution to Beamlet Center Location relative to beamline axis, as per Cyrille Thomas
                beamletX = envXCenOff[jj] + 2 * sRasterXAmpl[jj] / pi * asin(sin(2 * pi / periodX[jj] * tjj_ii )) #[mm]
                beamletY = envYCenOff[jj] + 2 * sRasterYAmpl[jj] / pi * asin(sin(2 * pi / periodY[jj] * tjj_ii )) #[mm]

                #Calculate the total Beamlet Angle = Envelope Center Angle + the angle given to each beamlet by the Raster Magnets
                beamletXAngle = beamletX / dBPM93toTarg + envXAngle #[mm]
                beamletYAngle = beamletY / dBPM93toTarg + envYAngle #[mm]

                #save total beamlet position
                centroids[i,0] = beamletX
                centroids[i,1] = beamletY
                NperBunch = NperBunch

                #In case of weighting edges of raster 
                if edges:
                    if beamletX > -Right and beamletX < Right and beamletY < Top and beamletY > -Top: #set weight depending on position
                        NperBunch = 5 #decrease center NperBunch
                    else:
                        NperBunch = NperBunch #Edges get full NperBunch

                    #Generate beamlet distributions
                    rng = np.random.default_rng()
                    ptsX = rng.multivariate_normal([beamletX,beamletXAngle],covX,size = NperBunch) #mean is [pos,ang]!
                    ptsY = rng.multivariate_normal([beamletY,beamletYAngle],covY,size = NperBunch)
                    
                    for k in range(NperBunch): #put this beamlet into total. Could just be written, figure that out later.
                        totX[NperBunch*i+k,0] = ptsX[k,0]
                        totX[NperBunch*i+k,1] = ptsX[k,1]
                        totY[NperBunch*i+k,0] = ptsY[k,0]
                        totY[NperBunch*i+k,1] = ptsY[k,1]
                    i+=1

        #Check on output parameters
        print("Centroid X max: {:.2f}mm; Particle X max: {:.2f}mm".format(np.max(centroids[:,0]),np.max(totX[:,0])),"; Shape:",np.shape(totX))
        #Remove 0,0 particles, should be none except for weighted edges case
        nonzero = np.not_equal(totX[:,0],0) #be careful!
        totX = totX[nonzero]
        totY = totY[nonzero]
        #print("After masking out 0s: ",np.shape(totX))

        import csv
        with open(outname+".csv",mode = 'w',newline=None) as part_file:
            part_writer = csv.writer(part_file,delimiter = ',')
            for i in range(len(totX)):
                part_writer.writerow(["proton", totX[i,0], totX[i,1], totY[i,0], totY[i,1], z, energy])
        part_file.close()

        finish = datetime.now()
        print(name,"took: ",finish-start,"s")

        if options['saveRaster']:
            #found the below method: https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
            import matplotlib.pyplot as plt
            import mpl_scatter_density
            plt.close()
            fig = plt.figure()
            s1 = fig.add_subplot(1,1,1,projection="scatter_density")
            x=totX[:,0]
            y=totY[:,0]
            density = s1.scatter_density(x,y,cmap='jet')
            fig.colorbar(density,label=r"Protons/mm^2")
            s1.set_xlabel("X [mm]")
            s1.set_ylabel("Y [mm]")
            s1.set_xlim([-100,100])
            s1.set_ylim([-50,50])
            s1.set_title("Rastered Beam Number Density\n{:.1e} protons {:.2f}ms".format(len(totX),time_length*1e-3))
            import os
            if os.uname()[1] == "tensor.uio.no":
                picPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
            plt.savefig(picPWD+name+".png")
            print(picPWD+name+".png")
            plt.close()

    return outname, envXAngle,envYAngle
  