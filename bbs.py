#!/usr/bin/env python
import os, argparse, time, skymodel, pyrap.tables as pt
from colors import *

def skyf(ms,n=5,gap=1,st=[1,0,0,0],type='POINT'):
    model = skymodel.write(ms = ms, n=n, gap=gap, st=st, type=type)
    return model

def timeavg(msin, msout, tstep):
    t = open('timeavg.parset', 'w')
    t.write("""msin = %s
    msout = %s
    steps = [avg]
    avg.type = average
    avg.timestep = %s""" % (msin, msout, tstep) )
    t.close()
    os.system('DPPP timeavg.parset')

def newcol(ms, icol, ocol):
    t = open('newcol.parset', 'w')
    t.write("""msin = %s
    msin.datacolumn = %s
    msout = .
    msout.datacolumn = %s
    steps = []"""%(ms, icol, ocol))
    t.close()
    os.system('DPPP newcol.parset')
    
def predict(ms, sky='', col='MODEL_DATA', force='-f', parmdb='instrument', \
            E='F', G='F', mode='DEFAULT', baselines='CS* &', sources='', chunksize='195'):
    t1=time.time()
    try: pt.taql('update %s set %s=0+0i'%(ms,col))
    except: None
    pt.taql('UPDATE %s set FLAG=False' % ms)
    pt.taql('UPDATE %s set FLAG_ROW=False' % ms)
    parset = 'predict.parset'
    t = open(parset, 'w')
    t.write("""Strategy.InputColumn = DATA
Strategy.TimeRange = []
Strategy.ChunkSize = %s
Strategy.Baselines = %s
Strategy.Steps = [predict]

Step.predict.Operation = PREDICT
Step.predict.Model.Gain.Enable = %s
Step.predict.Model.Beam.Enable = %s
Step.predict.Model.Beam.Mode = %s
Step.predict.Model.Sources = [%s]
Step.predict.Output.Column = %s""" % (chunksize, baselines, G, E, mode, sources, col) )
    t.close()
    log = ms+'.'+col+'.log'
    if sky=='': sky = skyf(ms,n,gap,st,type)
    command = "calibrate-stand-alone -t 16 -v %s --parmdb-name %s %s %s %s &> bbs_log.txt" \
              % (force, parmdb, ms, parset, sky)
    print '--> %s' % command
    os.system(command)
    t2=time.time()
    print colr.red + '--> time taken: %s min' % str((t2-t1)/60.)[0:4] + colr.end

def correct(ms, sky, icol, ocol, force='-f', parmdb='instrument', E='T', G='F', mode='DEFAULT', \
            baselines='[CR]S*&', sources='', chunksize='195'):
    parset = 'correct.parset'
    try: pt.taql('update %s set %s=0+0i'%(ms,ocol))
    except: None
    t = open(parset, 'w')
    t.write("""Strategy.InputColumn = %s
Strategy.TimeRange = []
Strategy.ChunkSize = %s
Strategy.Baselines = %s
Strategy.UseSolver = F
Strategy.Correlations = []
Strategy.Steps = [correct]

Step.correct.Operation = CORRECT
Step.correct.Model.Cache.Enable = T
Step.correct.Model.Gain.Enable = %s
Step.correct.Model.Beam.Enable = %s
Step.predict.Model.Beam.Mode = %s
Step.correct.Model.Sources = [%s]
Step.correct.Output.Column = %s""" % (icol, chunksize, baselines, G, E, mode, sources, ocol) )
    t.close()
    log = ms+'.'+ocol+'.log'
    command = "calibrate-stand-alone -t 16 -v %s --parmdb-name %s %s %s %s &> bbs_log.txt" \
              % (force, parmdb, ms, parset, sky)
    print command
    os.system(command)

def calibrate(ms, sky, icol, ocol, force='-f', parmdb='instrument', E='T', G='F', \
              baselines='[CR]S*&', ssources='', csources='', chunksize='300', nmiter='100'):
    t1=time.time()
    try: pt.taql('update %s set %s=0+0i'%(ms,ocol))
    except: None
    parset = 'calibrate.parset'
    t = open(parset, 'w')
    t.write("""Strategy.InputColumn = %s
Strategy.ChunkSize = %s
Strategy.Baselines = %s
Strategy.UseSolver = F
Strategy.Correlations = []
Strategy.Steps = [solve,correct]

Step.solve.Operation = SOLVE
Step.solve.Model.Cache.Enable = T
Step.solve.Model.Gain.Enable = %s
Step.solve.Model.Beam.Enable = %s
Step.solve.Model.Ionosphere.Enable = F
Step.solve.Model.Sources = [%s]
Step.solve.Solve.Parms = ["Gain:*","Gain:*"]
Step.solve.Solve.ExclParms = []
Step.solve.Solve.CalibrationGroups = []
Step.solve.Solve.CellSize.Freq = 0
Step.solve.Solve.CellSize.Time = 1
Step.solve.Solve.CellChunkSize = 75
Step.solve.Solve.UVRange = []
Step.solve.Solve.PropagateSolutions = F
Step.solve.Solve.Options.MaxIter = %s
Step.solve.Solve.Options.EpsValue = 1e-9
Step.solve.Solve.Options.EpsDerivative = 1e-9
Step.solve.Solve.Options.ColFactor = 1e-9
Step.solve.Solve.Options.LMFactor = 1.0
Step.solve.Solve.Options.BalancedEqs = F
Step.solve.Solve.Options.UseSVD = T

Step.correct.Operation = CORRECT
Step.correct.Model.Cache.Enable = T
Step.correct.Model.Gain.Enable = %s
Step.correct.Model.Beam.Enable = %s
Step.correct.Model.Sources = [%s]
Step.correct.Output.Column = %s""" % (icol, chunksize, baselines, G, E, ssources, nmiter, \
                                       G, E, csources, ocol) )
    t.close()
    log = ms+'.'+ocol+'.log'
    command = "calibrate-stand-alone -t 16 -v %s --parmdb-name %s %s %s %s | tee %s" \
              % (force, parmdb, ms, parset, sky, log)
    print '--> %s'%command
    os.system(command)
    t2=time.time()
    print colr.red + '--> time taken: %s min' % str((t2-t1)/60.)[0:4] + colr.end

def main():
    parser=argparse.ArgumentParser(description='Run BBS predict, solve and/or correct')
    parser.add_argument('-t', help='task: pred, cal, cor', required=True)
    parser.add_argument('-i', help='MS name', required=True)
    parser.add_argument('-s', help='Sky model name', required=True)
    parser.add_argument('-c', help='Column name, or input column', \
                        default='MODEL_DATA', required=False)
    parser.add_argument('-o', help='Output column', default='', required=False)
    parser.add_argument('-f', help="force -f/''", default='-f', required=False)
    parser.add_argument('-p', help='Parameter database name', default='instrument', \
                        required=False)
    parser.add_argument('-e', help='Beam enabled: T/F', default='F', required=False)
    parser.add_argument('-g', help='Gain enabled: T/F', default='F', required=False)
    parser.add_argument('-b', help='Baselines: casa selection syntax', \
                        default='[CR]S*&', required=False)
    parser.add_argument('-so', help='Names of sources; sources to solve for in case \
    of calibration', default='', required=False)
    parser.add_argument('-sc', help='Names of sources to correct for after calibration', \
                        default='', required=False)
    parser.add_argument('-ch', help='Chunksize, amount of data loaded at once; \
    in units of #timeslots', default='100', required=False)
    args=parser.parse_args()
    if args.t == 'pred':
        predict(ms=args.i, sky=args.s, col=args.c, force=args.f, parmdb=args.p, \
                E=args.e, G=args.g, baselines=args.b, sources=args.so, chunksize=args.ch)
    elif args.t == 'cor':
        correct(ms=args.i, sky=args.s, icol=args.c, ocol=args.o, force=args.f, \
                parmdb=args.p, E=args.e, G=args.g, baselines=args.b, sources=args.so, \
                chunksize=args.ch)
    elif args.t == 'cal':
        calibrate(ms=args.i, sky=args.s, icol=args.c, ocol=args.o, force=args.f, \
                  parmdb=args.p, E=args.e, G=args.g, baselines=args.b, ssources=args.so, \
                  csources=args.sc, chunksize=args.ch)

if __name__=='__main__':
    main()
