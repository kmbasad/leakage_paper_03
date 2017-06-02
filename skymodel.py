#!/usr/bin/env python
import numpy as np, pyrap.tables as pt
from sys import argv

def hms2deg(inp):
    inp = str(inp).split(' ')
    if len(inp) == 3:
        return str(((float(inp[2])/60. + float(inp[1]))/60. + float(inp[0]))*15.)
    else:
        x = float(inp[0])/15.
	h = int(x)
	m = int((x-h)*60.)
	s = (((x-h)*60.)-m)*60.
	return str(h)+':'+str(m)+':'+str(s)

def dms2deg(inp):
    inp = str(inp).split(' ')
    if len(inp) == 3:
        return inp[0][0:1]+str(((float(inp[2])/60. + float(inp[1]))/60. + float(inp[0])))
    else:
        if inp<0: sign = '-'
        else: sign = '+'
        x = float(inp[0])
	d = int(x)
	m = int((x-d)*60.)
	s = (((x-d)*60.)-m)*60.
	return sign+str(d)+'.'+str(m)+'.'+str(s)	

def write(ms, n=6, gap=0.3, st=[10,0,0,0], type='POINT', name=''):
    t = pt.table(ms, readonly=True, ack=False)
    field = pt.table(t.getkeyword('FIELD'), readonly=True, ack=False)
    refdir = field.getcell('REFERENCE_DIR', 0)[0]
    cra = refdir[0]*(180/np.pi)
    cdec = refdir[1]*(180/np.pi)

    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    freq = tsp.getcell('REF_FREQUENCY',0)
    freq = str(round(freq/1e6,1))+'e6'

    if cra==0.0: x = np.linspace(0., 360.-(360./n), n)
    else: x = np.linspace(cra-n/2*gap, cra+n/2*gap, n)
    if cdec==90.0: y = np.linspace(cdec-n*gap, cdec, n)
    else: y = np.linspace(cdec-n/2*gap, cdec+n/2*gap, n)
    mesh = np.meshgrid(x,y)

    if name=='': name = 'sky/gridded_sky_'+str(n)+'-'+str(gap)
    t=open(name,'w')
    t.write("# (Name, Type, RA, DEC, I, Q, U, V, MajorAxis, MinorAxis, Orientation, ReferenceFrequency='150e6', SpectralIndex='[-0.75]') = format\n")

    c=0
    st=map(str,st)
    for i in range(n):
        for j in range(n):
            c+=1
            ra = hms2deg(mesh[0][i][j])
            dec = dms2deg(mesh[1][i][j])
            t.write('S'+str(c)+', '+type+', '+ra+', '+dec+', '+st[0]+', '+st[1]+', '+\
                    st[2]+', '+st[3]+', '+'0, 0, 0, '+freq+', -0.75\n')
    t.close()
    print 'saved as %s'%name
    return name

if __name__=="__main__":
    write(argv[1])
