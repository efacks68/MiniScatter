#betaPOutFit.py
#python3 spreadTwiss.py --source twiss --twissFile HEBT-A2T_100pctField_1.0e-03Jitter_200x --beamClass jitter --samples 200 --saveFits
#with Confidence Ellipse fitting from https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html

def betaPOutFit(args,Twiss,paths,origBX,origBY,beamFile,axis):
    import csv,re,matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    from betaPOutFit import confidence_ellipse
#    from matplotlib.pyplot import scatter,savefig,close,plot,xlim,ylim,text,title,xlabel,ylabel,tight_layout,xticks
    from numpy import greater,zeros,mean,std
    from scipy.stats import linregress
    #from plotFit import findFit, gaussian
    from math import floor,log10

    betas = zeros(args.samples)
    jMaxes = zeros(args.samples)
    pOuts = zeros(args.samples)
    betas.fill(-750) #allow better filter than nonzero?
    jMaxes.fill(-750)
    pOuts.fill(-750)
    #print(Twiss,beamFile)
    i=0
    #axis="x"#x"
    if axis in {"Y","y"}:
        ind = 5  #beta-5, emitt-4, alpha-6
    elif axis in {"X","x"}:
        ind = 2 #beta-2, emitt-2,alpha-3

    if re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(Jitter)))",args.twissFile)[1] == "-03":
        deltaX = .999
        deltaY = 1.002
        pct=10
    elif re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(Jitter)))",args.twissFile)[1] == "-04":
        deltaX = .999
        deltaY = 1.0003
        pct=1
    #How to select which entries? Make keys to search for

    pbwKey = "PBW_{:.0f}MeV".format(args.energy)
    nBkey = "{:.0f}_NpB{:.0f}".format(floor(log10(2.88e4*args.Nb)),args.Nb)
    betaKey = "beta{:.2f},{:.2f}m".format(Twiss[1],Twiss[4])
    origKey = "OrigbX{:.2f},bY{:.2f}".format(origBX,origBY)
    failKey = "failure{:.0f}-{:.0f}f".format(args.failure,args.magFails)
    qpKey = "QP"+args.qpNum
    print(qpKey)
    #print(nBkey,pctKey,origKey,betaKey)
    #print(re.search(nBkey,beamFile) , re.search(pctKey,beamFile) , re.search(origKey,beamFile), re.search(betaKey,beamFile) ,"\n\n")
    with open(paths['statsPWD']+args.statsFile+".csv",mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        i = 0
        #requires differentiating between type of run we're looking at. Can't have both in one line
        if args.betaSpread != 0: #sampleIn files
            for row in csv_reader:
                #if i == 6: #for finding indices out of 7sigma
                #    print("index",i,"failed file:",row[0])
                #print(re.search(nBkey,row[0]) , re.search(pctKey,row[0]) , (re.search(origKey,row[0])))
                if re.search(nBkey,row[0]) and re.search(pctKey,row[0]) and re.search(origKey,row[0]):
                    betas[i] = float(row[ind])
                    jMaxes[i] = float(row[7])
                    pOuts[i] = float(row[8])
                    if i == args.samples-1:
                        lenbreak = True
                        #print("lenbreak 1")
                        break
                    i+=1
        elif args.failure != 0: #RM fails, though only PBW... files at this time
            for row in csv_reader:
                if re.search(nBkey,row[0]) and re.search(failKey,row[0]) and re.search(pbwKey,row[0]) and re.search(betaKey,row[0]) :
                    betas[i] = float(row[ind])
                    jMaxes[i] = float(row[7])
                    pOuts[i] = float(row[8])
                    if i == args.samples-1:
                        lenbreak = True
                        #print("lenbreak 1")
                        break
                    i+=1
        elif args.beamClass == "jitter": #for jitter studies from OXAL
            pctKey = re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(pct)))",args.twissFile)[1]+"pctField"
            jitterKey = re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(Jitter)))",args.twissFile)[1]+"Jitter"#_{:.0f}x".format(args.samples)
            print(pctKey,nBkey,jitterKey,origKey)
            for row in csv_reader:
                #if i == 6: #for finding indices out of 7sigma
                #    print("index",i,"failed file:",row[0])
                #print(re.search(nBkey,row[0]) , re.search(jitterKey,row[0]) , (re.search(origKey,row[0])))
                if re.search(nBkey,row[0]) and re.search(args.twissFile,row[0]) and re.search(pctKey,row[0]):
                    betas[i] = float(row[ind])
                    jMaxes[i] = float(row[7])
                    pOuts[i] = float(row[8])
                    if i == args.samples-1:
                        #print("lenbreak 1")
                        lenbreak = True
                        break
                    i+=1
        elif args.beamClass == "qpFail": #for individual QP spread given in qpNum from OXAL twissFile 
            pctKey = re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(pct)))",args.twissFile)[1]+"pctField"
            jitterKey = re.search("(([-+]?[0-9]*\.?[0-9]*)(?=(Jitter)))",args.twissFile)[1]+"JitterJL"
            print("qpFail",pctKey,nBkey,origKey,qpKey,jitterKey)
            for row in csv_reader:
                #if i == 6: #for finding indices out of 7sigma
                #    print("index",i,"failed file:",row[0])
                #print(re.search(nBkey,row[0]) , re.search(jitterKey,row[0]) , (re.search(origKey,row[0])))
                if re.search(nBkey,row[0]) and re.search(origKey,row[0]) and re.search(qpKey,row[0]) and re.search(pctKey,row[0]):# and re.search(jitterKey,row[0]):
                    betas[i] = float(row[ind])
                    jMaxes[i] = float(row[7])
                    pOuts[i] = float(row[8])
                    if i == args.samples-1:
                        lenbreak = True
                        print(row[0])
                        #print("lenbreak 1")
                        break
                    i+=1
        else: #nominal
            i = 0
            for row in csv_reader:
                #if i == 44: #for finding indices out of 7sigma
                #    print("index",i,"failed file:",row[0])
                #print("nominal",re.search(failKey,row[0]))#re.search(nBkey,row[0]) , re.search(pbwKey,row[0]) , (re.search(betaKey,row[0])))
                if re.search(nBkey,row[0]) and re.search(pbwKey,row[0]) and re.search(betaKey,row[0]) and not re.search("failure",row[0]):
                    betas[i] = float(row[ind])
                    jMaxes[i] = float(row[7])
                    pOuts[i] = float(row[8])
                    if i == args.samples-1:
                        #print("lenbreak 1")
                        lenbreak = True
                        break
                    i+=1
    csv_file.close()

    nonEmpty = greater(betas,-750) #remove unused elements
    betas = betas[nonEmpty]
    pOuts =pOuts[nonEmpty]

    fig, ax_nstd = plt.subplots(figsize=(6, 6))

    mu = mean(betas), mean(pOuts)
    sigma = std(betas), std(pOuts)

    print(mu[0],ax_nstd.get_ylim()[0],ax_nstd.get_ylim()[1])
    print(mu[1],ax_nstd.get_xlim()[0],ax_nstd.get_xlim()[1])
    ax_nstd.axvline(mu[0],c='grey', lw=1,alpha=0.4)
    ax_nstd.axhline(mu[1],c='grey', lw=1,alpha=0.4)

    ax_nstd.scatter(betas, pOuts, s=1)

    confidence_ellipse(betas, pOuts, ax_nstd, n_std=1,
                    label=r'$1\sigma$', edgecolor='firebrick')
    confidence_ellipse(betas, pOuts, ax_nstd, n_std=2,
                    label=r'$2\sigma$', edgecolor='fuchsia', linestyle='--')
    confidence_ellipse(betas, pOuts, ax_nstd, n_std=3,
                    label=r'$3\sigma$', edgecolor='blue', linestyle=':')

    ax_nstd.scatter(mu[0], mu[1], c='red', s=3)
    #ax_nstd.set_title('Different standard deviations')

    fs=16
    slope, intercept, r, p, se = linregress(betas, pOuts)
    plt.plot(betas,slope*betas+intercept,c='g',alpha=0.5,label="Fit")
    if axis in {"Y","y"}:
        plt.title(r"$\beta_y$ vs. % Outside Target Area"+"\nfor {:.0f}% ".format(pct)+r"$\beta$"+" Variations Around Nominal",fontsize=fs+2)
        plt.xlabel(r"$\beta_y$ [m]",fontsize=fs)
    elif axis in {"X","x"}:
        plt.title(r"$\beta_x$ vs. % Outside Target Area"+"\nfor {:.0f}% ".format(pct)+r"$\beta$"+" Variations Around Nominal",fontsize=fs+2)
        plt.xlabel(r"$\beta_x$ [m]",fontsize=fs)
    plt.ylabel("% Outside Target Area",fontsize=fs)

    plt.text(0.99,0.01,r"R$^2$"+" = {:.4f}\n% = {:.3e}".format(r**2,slope)+r"$\beta$"+" + {:.3f}".format(intercept),ha="right",va="bottom",fontsize=fs-2,transform=ax_nstd.transAxes)

    ax_nstd.legend(loc="upper left",fontsize=fs-2)
    plt.tight_layout()
    print(paths['statsPWD']+args.twissFile+"_{:.0f}pBetaVpOut".format(pct)+axis+".png")
    plt.savefig(paths['statsPWD']+args.twissFile+"_{:.0f}pBetaVpOut".format(pct)+axis+".png",bbox_inches='tight',dpi=args.dpi)

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)