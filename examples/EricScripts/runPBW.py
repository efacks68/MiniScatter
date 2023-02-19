#runPBW.py #need better name
#Eric Fackelman
#5 Nov 2022
#Uses the MiniScatter Python interface to model the ESS beam interaction with the PBW

def runPBW(args,beamFile,Twiss,options,boxes):
    #Code setup
    import numpy as np
    from math import floor,log10,ceil
    from plotFit import plotFit,findFit,numLines
    from simulation import simulation
    from datetime import datetime

    #Define something preliminarily
    materials = ["G4_Al"] #need to move this
    N = args.Nbeamlet
    #print("You've entered: {:f}mm thick".format(thick),materials,", {:.0e} protons, ".format(N),betax,alphx,nemtx,betay,alphy,nemty)
    mag=floor(log10(N)) #for dynamically scaling the halo plots

    #Compute Angle as in runARasterMaker bc used for MCS calculation
    dBPM93to94 = 3031 #[mm] from OpenXAL(?)
    beamXAngle = args.rX / dBPM93to94
    beamYAngle = args.rY / dBPM93to94

    #Opening figure only if doing material plotting
    if not args.matPlots:
        savename,xtarg,ytarg,Jmax,pOutsideBoxes,dispY,dispX,rValue = simulation(args,args.material,beamXAngle,beamYAngle,beamFile,Twiss,options,boxes)
    elif args.matPlots:
        materials = ["G4_Galactic","G4_Al","G4_Au"] #,"G4_AIR"full list as of 1.4.22
        import matplotlib.pyplot as plt
        #Create the fig before the loop with 2 plots side by side
        fig = plt.figure(figsize=(15,8))
        plt.subplots_adjust(wspace=0.25) #increase width space to not overlap
        s1 = fig.add_subplot(1,2,1)
        s2 = fig.add_subplot(1,2,2)
        fs = 14 #setting the axis label font size early

        #Dictionaries for sigma and % outside of 3 Sigma for each material
        sigmax = {}
        sigmay = {}
        pOut3sigx = {}
        pOut3sigy = {}

        #Start the text boxes for displaying info in the graphs
        percenttextx = r"% outside 3$\sigma_x$:"
        percenttexty = r"% outside 3$\sigma_y$:"
        sigmatextx = r"$\sigma_x$:"
        sigmatexty = r"$\sigma_y$:"
    
        for material in materials:
            #function for preparing the run and running miniScatterDriver functions
            #savename,xtarg,ytarg,targPOutBox,targImax, targCoreMeanI =  simulation(args,mat,beamXAngle,beamYAngle,beamFile,Twiss,options,boxes)
            savename,xtarg,ytarg,Jmax,pOutsideBoxes,dispY,dispX,rValue = simulation(args,mat,   args.rX,   args.rY,beamFile,Twiss,options,boxes)
            #Now plot the distributions with various views depending on the material
            if material == "G4_Galactic" or material == "G4_AIR":
                if options['mat3Plot']:
                    print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xtarg),np.max(ytarg)))
                    xmax = ceil(np.max(xtarg)/100)*100
                    #For plotting x,y plots for core view and halo zoom-ins
                    #plotFit(xs,    ys, savename,xlim,   ylim,material)
                    #xlim<=10 => [-xlim*sigma,xlim*sigma]; xlim>10 => [-xlim,xlim]
                    #ylim==0 => with plot; ylim!=0 => [0,ylim] (5 particles /(mag of N) is for halo)
                    #plotFit(xtarg,ytarg,savename,  3,      0,material,thick) #3 sigma core
                    plotFit(xtarg,ytarg,savename, xmax,5/(10**(mag+0)),material,args.t) #full range halo, with dynamic halo zoom
            elif material == "G4_Al":
                if options['mat3Plot']:
                    print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xtarg),np.max(ytarg)))
                    xmax = ceil(np.max(xtarg)/100)*100
                    print(xmax)
                    #plotFit(xs,    ys, savename,xlim,   ylim,material)
                    #plotFit(xtarg,ytarg,savename,  3,      0,material,thick) #3 sigma core
                    plotFit(xtarg,ytarg,savename, xmax,5/(10**(mag+0)),material,args.t) #full range halo
                    #plotFit(xtarg,ytarg,savename,  10,10/(10**(mag+0)),material,thick) #10 sigma range halo
            elif material == "G4_Au":
                if options['mat3Plot']:
                    print("Max in x: {:.3f} and in y: {:.3f}".format(np.max(xtarg),np.max(ytarg)))
                    xmax = ceil(np.max(xtarg)/100)*100
                    #plotFit(xs,    ys, savename,xlim,   ylim,material) 
                    #plotFit(xtarg,ytarg,savename,  3,      0,material,thick) #3 sigma core
                    plotFit(xtarg,ytarg,savename, xmax,15/(10**(mag+0)),material,args.t) #full range halo
                    #plotFit(xtarg,ytarg,savename,  10,15/(10**(mag+0)),material,thick) #10 sigma range halo
        
            #Multi-Material Plot section, continues the material loop to plot data.
            #For creating the PDF of particle distribution from findFit function optimized mu and sigma.
            from scipy.stats import norm 
            
            #Define plotting characteristics depending on material
            if material == "G4_Galactic":
                mat = "Vacuum"
                color = "magenta"
                dash = "m--"
                tsp = 0.2 #transparency
            elif material == "G4_Al":
                mat ="Al"
                color = "blue"
                dash = "b--"
                tsp = 0.5
            elif material == "G4_AIR":
                mat ="Air"
                color = "cyan"
                tsp = 0.2
                dash = "c--"
            elif material == "G4_Au":
                mat = "Au"
                color = "gold"
                tsp = 0.5
                dash = "r--"

            #Use Scipy.optimize.curve_fit in my findFit function to get mus and sigmas:
            mux, sigmax[material], xinterval = findFit(xtarg) #dynamically gets parameters AND histogram intervals!
            muy, sigmay[material], yinterval = findFit(ytarg)

            #Find range of particles that are within 3 sigma for each material
            sigx=np.abs(xtarg)>3*sigmax[material]# and np.abs(xtarg)<10*sigma)
            sigy=np.abs(ytarg)>3*sigmay[material]

            #Utilize dictionary to find % of particles within 3 sigma 
            pOut3sigx[material] = len(xtarg[sigx])/len(xtarg)*100 
            pOut3sigy[material] = len(ytarg[sigy])/len(ytarg)*100
            #print(material," gives ",pOut3sigx[material],"% outisde 3 sigma in X")
            #print(material," gives ",pOut3sigy[material],"% outisde 3 sigma in Y")

            #Update the texts to include this materials sigma and % outside 3 sigma
            percenttextx +="\n"+mat+" = "+"{:.2f}".format(pOut3sigx[material])+"%"
            percenttexty +="\n"+mat+" = "+"{:.2f}".format(pOut3sigy[material])+"%"
            sigmatextx +="\n"+mat+" = "+"{:.2f}".format(sigmax[material])+"mm"
            sigmatexty +="\n"+mat+" = "+"{:.2f}".format(sigmay[material])+"mm"

            #Make the histogram of the full energy distrubtion for X
            nx, binsx, patchesx = s1.hist(xtarg, xinterval, density=True, facecolor=color, alpha=tsp,label=mat)
            
            #Add the "best fit" line using the earlier norm.fit mus and sigmas for X
            y1 = norm.pdf(binsx, mux, sigmax[material])
            l1 = s1.plot(binsx, y1, dash, linewidth=1,label=mat) #important to keep it as l# or it won't work

            #Make the histogram of the full energy distrubtion for Y
            ny, binsy, patchesy = s2.hist(ytarg, yinterval, density=True, facecolor=color, alpha=tsp,label=mat)
            
            #Add the "best fit" line using the earlier norm.fit mus and sigmas for Y
            y2 = norm.pdf(binsy, muy, sigmay[material])
            l2 = s2.plot(binsy, y2, dash, linewidth=1,label=mat) #labeled by material short-name

            #Set Plot characterists
            s=10 # number of sigma width to plot
            sigvals = sigmax.values() #get values of dictionary
            xlim1 = s*max(sigvals) #use s*max as the xlim
            s1.set_xlim([-xlim1,xlim1]) #seems to show all distributions, 30.3.22
            s1.set_title("Fitted X Distributions",fontsize=fs)
            s1.set_xlabel("X Position [mm]",fontsize=fs)
            s1.set_ylabel("Probability Density",fontsize=fs)
            s1.legend()
            ylim1=s1.get_ylim() #dynamically get the ylimits
            s1.text(-xlim1*0.95,ylim1[1]*0.75,sigmatextx) #set the texts at 3/4 and 1/2 of ylim
            s1.text(-xlim1*0.95,ylim1[1]*0.5,percenttextx) #xlim is fixed

            sigvals = sigmay.values() #get values of dictionary
            xlim2 = s*max(sigvals) #use s*max as the xlim
            s2.set_xlim([-xlim2,xlim2])
            s2.set_title("Fitted Y Distributions",fontsize=fs)
            s2.set_xlabel("Y Position [mm]",fontsize=fs)
            s2.set_ylabel("Probability Density",fontsize=fs)
            s2.legend()
            ylim2=s2.get_ylim()
            s2.text(-xlim2*0.95,ylim2[1]*0.75,sigmatexty) #set the texts at 3/4 and 1/2 of ylim
            s2.text(-xlim2*0.95,ylim2[1]*0.5,percenttexty) #x position is fixed
            fig.suptitle(rf"Distributions at ESS Target of 10$^{{:d}}$ Protons".format(mag)+" Through Various Material PBWs\n"+
                    rf"Initial beam of $\epsilon_x={{:.1e}}, \beta_x={{:.0e}}, \alpha_x={{:.1f}}$, ".format(Twiss[0],Twiss[1],Twiss[2]) +
                rf" $\epsilon_y={{:.1e}}, \beta_y={{:.0e}}, \alpha_y={{:.1f}}$ ".format(Twiss[3],Twiss[4],Twiss[5]),fontsize=18,y=0.99) #fontweight="bold",
            #suptitle can have values inside math if use 2 sets of {{}} - fix from "linuxtut.com" blog post

            dt = datetime.now()
            name = savename+"_"+dt.strftime("%H-%M-%S")+"_multi."+args.picFormat #update the savename to not overwrite others
            #print(name)
            if args.savePics:
                fig.savefig(name,bbox_inches="tight")
            plt.close() #be sure to close the plot
        #print(datetime.now().strftime("%H-%M-%S"),"\t",datetime.now()-start)

    #return Jmax,pOutsideBoxes,dispY,dispX,rValue
