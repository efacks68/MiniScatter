#plot spread of same Twiss

def spreadHist(args,Twiss,n):
    from os import uname
    import csv
    from matplotlib.pyplot import hist,savefig,close
    from numpy import greater,zeros

    if uname()[1] == "tensor.uio.no":
        statsPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    elif uname()[1] == "mbef-xps-13-9300":
        statsPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"

    readPOut = zeros(n)
    #print(readPOut.shape)
    i=0
    lenbreak=False
    with open(statsPWD+"EvalStats.csv",mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            if float(row[1]) == Twiss[0] and float(row[2]) == Twiss[1] and float(row[3]) == Twiss[2] and float(row[4]) == Twiss[3] and float(row[5]) == Twiss[4] and float(row[6]) == Twiss[5]:
                #print("found Twiss")
                readPOut[i] = float(row[8]) #[mm-mrad]
                #print(i,readPOut[i])
                if i+1 == n:
                    lenbreak = True
                    break
                i+=1
                #print(lenbreak)
            if lenbreak:
                break
    csv_file.close()

    nonzero = greater(readPOut,0)
    readPOut = readPOut[nonzero]
    print(readPOut.shape)
    hist(readPOut)
    name=statsPWD+"beamPOutHist_"+str(args.pOff)+"pOffNom.png"
    savefig(name)
    close()
    print(name)

#from argparse import ArgumentParser

#parser = ArgumentParser()
#parser.add_argument("--beamClass", type=str,   default="ESS", help="Determines beam Twiss: 'ESS', 'Yngve', or 'pencil. If other, just do --twiss. Default='ESS'")
#parser.add_argument("--twiss",     type=float, nargs=6,       help="Twiss parameters in form: NemtX[mm*mrad],BetaX[m],AlphX,NemtY[mm*mrad],BetaY[m],AlphY")
#parser.add_argument("--pOff", type=int,default=0)
#args = parser.parse_args()

#if args.beamClass == 'Yngve': #smallest from Yngve
#    Twiss = [0.3519001,144.15027172522036,-8.184063058768368,0.3651098,88.04934327630778,-1.0382192928960423]
#elif args.beamClass == 'ESS': #from my OpenXAL calculation
#    Twiss = [0.11315,1006.80,-60.44,0.12155,129.72,-7.72]
    #pOff = 10 #negative to have less than
#    Twiss[1] = Twiss[1] * (1 + args.pOff/100)
#    Twiss[4] = Twiss[4] * (1 + args.pOff/100)
#elif args.beamClass == 'pencil': #"pencil" beam of ~0 emittance
#    Twiss = [0.0001,0.15,0,0.0001,0.15,0]
#if args.twiss:
#        Twiss = [args.twiss[0],args.twiss[1],args.twiss[2],args.twiss[3],args.twiss[4],args.twiss[5]]
#        args.beamClass = "Twiss"

#n=20

#run(args,Twiss,n)