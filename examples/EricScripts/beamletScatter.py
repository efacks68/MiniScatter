#beamletScatter.py
#run running analysis on beamlet

def beamletScatter(args,Twiss,i,paths):
    from datetime import datetime
    origin = datetime.now()
    print(origin)
    import numpy as np
    from os import uname
    from runARasterMaker import runARasterMaker
    from runPBW import runPBW

    #print(args,Twiss)
    #Constants for running scripts
    physList    = "QGSP_BERT_EMZ" # "QGSP_BERT_EMZ" or "FTFP_BERT_EMZ" or "QGSP_BERT__SS"
    zoff        = "*-1" #[mm] with preappended * to keep covar defined at z=0
    options     = {'physList':args.physList, 'dependence':"Twiss", 'zoff':zoff, 'initTree':False,
                    'exitTree':False, 'targetTree':False, 'MCS':False, 'engPlot':False,
                    'mat3Plot':False, 'TwissFits':False,'MiniRoot':False }

    #Important things
    if args.t == 0:
        args.material = "Al" #overwrites potential user input. Needs work
    elif args.t == 0.1:
        args.material = "Vac"

    args.beamFile = ""#/scratch2/ericdf/PBWScatter/ESS/tDependence_beta1085.63,136.06m_t2.25mm_N1e+05_runW_QBZ"

    if args.gaussFit and not args.compTargs:
        options['MiniRoot'] = True
    if args.compTargs:
        options['targetTree'] = True
    #    options['exitTree'] = True
        options['MCS'] = True
    if options['TwissFits']:
        options['initTree'] = True
        options['exitTree'] = True
        options['targetTree'] = True
    if args.saveParts:
        options['MiniRoot'] = False
        options['initTree'] = True
    boxes = [0]#,-.25,-.375,-.45,-.50]#,0.125,0.25,0.375] #make an args for 24.11.22

    #Send to runPBW to simulate with MiniScatter or opens already run data
    PMASreturn,e8SigsX,e8SigsY,targSigsX,targSigsY = runPBW(args,args.beamFile,Twiss,options,boxes,paths)

    print("Simulation took ",datetime.now()-origin,"s long",sep="")