#beamletScatter.py
#run various methods on beamlet

def beamletScatter(args,Twiss):
    from datetime import datetime
    origin = datetime.now()
    print(origin)
    import numpy as np
    #import matplotlib.pyplot as plt
    from os import uname
    #from sys import path as sysPath
    #from argparse import ArgumentParser
    from runARasterMaker import runARasterMaker
    from runPBW import runPBW

    #Constants for running scripts
    physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ" or "QGSP_BERT__SS"
    zoff = "*-10" #[mm] with preappended * to keep covar defined at z=0

    options     = {'physList':physList, 'dependence':"Twiss", 'zoff':zoff }

    #Important things
    if args.t == 0:
        args.material = "Al" #overwrites potential user input. Needs work
    elif args.t == 0.1:
        args.material = "Vac"
    boxes = [0]#,-.25,-.375,-.45,-.50]#,0.125,0.25,0.375] #make an args for 24.11.22

    if uname()[1] == "tensor.uio.no":
        csvPWD = "/scratch2/ericdf/PBWScatter/CSVs/"
        statsPWD = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    elif uname()[1] == "mbef-xps-13-9300":
        csvPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/scatterPBWFiles/"
        statsPWD = "/home/efackelman/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
    else:
        csvPWD = input("Path from home to direction you like to save root files to: ")
        statsPWD = "."

    #Send to runPBW to simulate with MiniScatter or opens already run data
    runPBW(args,args.beamFile,Twiss,options,boxes)



