#! /usr/bin/env python
import os, re, argparse, pyrap.tables as pt, pyrap.images as pim, time, numpy as np
from colors import *

def Image(ms='',table='',col='',force='T',img='',fits='',selectdata='T',uvrange='.03~3klambda',\
	   niter='0',psfmode='clark',imsize='1024',imagermode='', cell='0.5arcmin',\
	   weight='uniform',gridmode='widefield',wprojplanes='64',st='IQUV',noise=None):
    t1 = time.time()
    # prepare columns
    print ms, img
    print uvrange, imsize, cell, weight
    if force=='T':
        if col!='CORRECTED_DATA':
            table = pt.table(ms, readonly=False, ack=False)
            d = table.getcol(col)
            if noise!=None:
                #rms = noise * np.sqrt(len(d))
                rms = noise * 1425. # for upto 800 lambda
                noise = np.random.normal(0,rms,d.shape) + 1j * np.random.normal(0,rms,d.shape)
                d = d + noise
                print '--> Added rms noise of %s Jy' % str(rms)
            print "Copying data from %s into CORRECTED_DATA for imaging"% col
            table.putcol('CORRECTED_DATA', d)
            table.close()
    elif force == 'F':
        table = pt.table(ms, readonly=False, ack=False)
	print "Backuping CORRECTED_DATA in DATA"
	table.putcol('DATA', table.getcol('CORRECTED_DATA'))
	print "Copying data from %s into CORRECTED_DATA for imaging"% col
	table.putcol('CORRECTED_DATA', table.getcol(col))
    else: None
    if img=='': img='img/'+ms[3:-3]+'.'+col+'.'+str(imsize)+'.'+cell+'.'+weight

    t=open('casathon.py','w')
    t.write("clean(vis='"+ms+"', imagename='"+img+"', selectdata="+selectdata+", \
    uvrange='"+uvrange+"',\n mode='mfs', niter="+niter+", psfmode='"+psfmode+"', \
    imagermode='"+imagermode+"', imsize="+imsize+", cell='"+cell+"',\n\
    stokes='"+st+"', weighting='"+weight+"', robust=-1, gridmode='"+gridmode+"', \
    wprojplanes="+wprojplanes+")\n")
    t.close()
    print "CASA parset file casathon.py created"
    if img!='': os.system('rm -r %s.*' % img)
    os.system('casapy --nologger -c casathon.py')
    if fits=='': fits=img+'.fits'
    pim.image(img+'.image').tofits(fits, overwrite=True)

    if force=='F':
        print "Restoring CORRECTED_DATA from DATA"
	table.putcol('CORRECTED_DATA', table.getcol('DATA'))
	table.close()
    t2 = time.time()
    print colr.red + '--> time taken:\t%s min' % str((t2-t1)/60)[0:5] + colr.end
    return fits

def noise_test():
    ms = 'ms/L80273_SAP000_SB179_uv.MS.dppp.1ch.dppp'
    t = pt.table(ms, readonly=False)
    d = t.getcol('CORRECTED_DATA')
    noise = np.random.normal(0,1,d.shape) + 1j * np.random.normal(0,1,d.shape)
    t.putcol('CORRECTED_DATA', noise)
    t.close()
    Image(ms=ms, col='CORRECTED_DATA', img='tests/noise_test', imsize='480', cell='0.5arcmin', \
               uvrange='0.03~0.8klambda', st='IQUV', niter='0', weight='uniform')

def main():
     parser=argparse.ArgumentParser(description='Imaging by CASA from a \
     specific column of a MS')
     parser.add_argument('-i',help='MS name',required=True)
     parser.add_argument('-o',help='IMG name',required=False)
     parser.add_argument('-c',help='Column name',required=True)
     parser.add_argument('-f',help='force T/F, use T only if CORRECTED_DATA \
     column is not needed',required=False)
     parser.add_argument('-im',default='400', help='image size in pixels',required=False)
     parser.add_argument('-cs',default='1.5arcmin', help='cell size in arcmin',\
                         required=False)
     parser.add_argument('-n',default='0', help='number of iterations',required=False)
     args=parser.parse_args()
     Image(ms=args.i, col=args.c, force=args.f, imsize=args.im, \
           cell=args.cs, img=args.o, niter=args.n)

if __name__=='__main__':
     main()
     #noise_test()
