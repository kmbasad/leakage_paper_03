#!/usr/bin/env python
import numpy as np, pyfits as pf, scipy.constants as sc, scipy.integrate as si, plots, sys, os
import matplotlib.pyplot as pl, scipy.interpolate as sint
from matplotlib.colors import LogNorm
from matplotlib.mlab import griddata

def radial_profile(data, c=None):
    y, x = np.indices(data.shape)
    if not c: c = map(int, [x.max()/2., y.max()/2.])
    r = np.hypot(x - c[0], y - c[1])
    r = r.astype(np.int)
    ind = np.where(r<=data.shape[0]/2.)
    tbin = np.bincount(r[ind], data[ind])
    nr = np.bincount(r[ind])
    rp = tbin / nr
    return rp

def binning(d,n,c=None):
    x = np.indices(d.shape)
    if not c: c = int(x.max()/2.)
    r = abs(x-c)
    s,e = np.histogram(r,weights=d.reshape(1,len(d)),bins=n)
    c,e = np.histogram(r,bins=n)
    return s/c

def spherical_profile(data, nbins, c=None, lcut=0, ucut=0):
    z, y, x = np.indices(data.shape)
    if not c: c = map(int, [x.max()/2., y.max()/2., z.max()/2.])
    r = np.hypot(x-c[0], y-c[1], z-c[2])
    sum, e = np.histogram(r, weights=data, bins=nbins)
    count, edges = np.histogram(r, bins=nbins)
    PS = sum/count
    l = len(e[e <= lcut])
    if ucut==0: ucut=data.shape[0]/2.
    u = len(e[e <= ucut])
    print len(PS), len(PS[l:u])
    return PS[l:u], e[l:u]

def comoving(z,dNu=0.2e6):
    H_0 = 100.*1e3 # Hubble constant at z=0 in units of h where h = H_0/100 km/s/Mpc
    H_0 = 67.8*1e3
    omega_m = 0.3 # matter density
    omega_l = 0.7 # cosmological constant
    DH = sc.c/H_0 # Hubble distance in Mpc
    fHI = sc.c/0.21
    Ezinv = lambda rs: 1./np.sqrt(omega_m*(1+rs)**3+omega_l) # inverse dimensionless Hubble parameter function
    Dc = DH * si.quad(Ezinv, 0, z)[0]

    H_z = H_0 * np.sqrt(omega_m*(1+z)**3 + omega_l) # Hubble constant at redshift z
    kzE = 2*np.pi*H_z*fHI/(sc.c*(1+z)**2)

    return Dc, kzE

def comoving_grid(umin, umax, xbin, zbin, dNu='', nu0='', Nz='', h='', B='', nu_avg=None, fov=3.78/2):
    try:
        Nx, Ny, Nz = h['NAXIS1'], h['NAXIS2'], h['NAXIS3']
        dNu = h['CDELT3'] # in Hz
        nu0 = h['CRVAL3'] # starting frequency
    except: None
    B = Nz*dNu # bandwidth
    fHI = sc.c/0.21
    nu = nu0+B/2. # central frequency
    z = fHI/nu - 1. # central redshift

    # comoving Mpc and k_parallel/Eta at redshift=z; Eta=small_bandwidth
    DM, kzE = comoving(z)
    
    # k_perp kr[0,1] and k_para kr[2,3] ranges
    kr = np.zeros((4), dtype='float32') # k-ranges will be stored in this list
    kr[0], kr[1] = 2*np.pi*umin/DM, 2*np.pi*umax/DM
    kr[2], kr[3] = kzE/B, kzE/dNu
    
    # calculate x and z grid
    k_perp = np.linspace(kr[0], kr[1], xbin)
    k_para = np.linspace(kr[2], kr[3], zbin)

    # the wedge
    k_wedge = wedge(k_perp, z, fov*(np.pi/180.) )

    return k_perp, k_para, k_wedge

def wedge(ku,z,fov):
    H_0 = 67.8*1e3
    omega_m = 0.3 # matter density
    omega_l = 0.7 # cosmological constant
    DH = sc.c/H_0 # Hubble distance in Mpc
    Ez = np.sqrt(omega_m*(1+z)**3+omega_l)
    DM, kzE = comoving(z)
    wedge = np.sin(fov) * DM*Ez/(DH*(1+z)) * ku
    return wedge

def threeD_wedge():
    p = np.load('3C196_diffuse_P.npy')
    h = pf.getheader('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.fits')
    zbin, xbin = p.shape[0], p.shape[1]
    x, y, w = comoving_grid(h, umin=30., umax=800., xbin=xbin, zbin=zbin)
    print w.min(), w.max()
    #plots.pcolr(p, x=x, y=y, title='', name='dum.pdf', sh=True, w=w)

def main(ft_factor=1, pixel=None, umin=30., umax=800., xbin=50, zbin=20, kbin=20, title='', \
           name='', fits='', shw=False, cube=''):
    cube = pf.getdata(fits)
    h = pf.getheader(fits)
    nu0, dNu = h['CRVAL3'], h['CDELT3']
    print nu0, dNu
    Nz = cube.shape[0]
    ft_shape = (cube.shape[0], cube.shape[1]*ft_factor, cube.shape[2]*ft_factor)
    print '--> Performing n-D Fourier transform with shape %s' % str(ft_shape)
    cube_ft = np.fft.fftn(cube, ft_shape) / np.sqrt(cube.size) # nD FFT
    cube_ft = np.fft.fftshift(cube_ft) # shift to center and normalize

    # check if the FT went okay
    # http://stackoverflow.com/questions/19444373/normalization-of-2d-fft
    pp_ratio = np.sum(np.abs(cube)**2) / np.sum(np.abs(cube_ft)**2)
    vm_ratio = np.std(cube)**2 / np.mean(np.abs(cube_ft)**2)
    print '--> Power-power ratio = %f'% pp_ratio
    print '--> Variance-integral ratio = %f'% vm_ratio

    # cut upto the longest baseline/shortest scale
    N = cube_ft.shape
    if pixel==None:
        fmin = 1/(np.abs(h['CDELT1'])*(np.pi/180.)*N[1]) # spatial resolution corresponding to the pixel size
    else: fmin = 1/(np.abs(pixel/60.)*(np.pi/180.)*N[1])
    px = int(umax/fmin)
    print 'cut: %s'%px
    c = N[1]/2
    print umin, umax, c

    cube_ft = cube_ft[:, c-px:c+px, c-px:c+px]
    print '--> Cut everything above umax, shape: %s'% str(cube_ft.shape)

    # calculate 3D power spectrum
    cube_ft = np.abs(cube_ft)**2

    # radial averaging in k_perp plane
    lcut = umin/fmin # lower uv cut
    dum, r = radial_profile(cube_ft[0,:,:], nbins=xbin, lcut=lcut) # find the actual bin size in xy-plane
    print '--> Radial binning in k_perp for every z using %s bins'% len(dum)
    ps_ubin = np.zeros((cube_ft.shape[0], len(dum)), dtype='float32')
    for i in range(cube_ft.shape[0]):
        ps_ubin[i,:], r = radial_profile(cube_ft[i,:,:], nbins=xbin, lcut=lcut)
    xbin = len(dum)

    # binning in k_para direction
    print '--> Binning in k_para using %s bins'% zbin
    ps_uzbin = np.zeros((zbin, xbin), dtype='float32')
    for i in range(xbin):
        ps_uzbin[:,i], r = binning(ps_ubin[:,i], zbin)
    np.save(name[:-4]+'.npy', ps_uzbin)

    # Calculate comoving grid and plot dimensional power spectrum
    print '--> Calculating power spectrum'
    #x, y, w = ps.comoving_grid(umin=umin, umax=umax, xbin=xbin, zbin=zbin, h=h)
    x, y, w = comoving_grid(umin=umin, umax=umax, xbin=xbin, zbin=zbin, dNu=dNu, nu0=nu0, Nz=Nz)
    PS = ps_uzbin

    if name!='': plots.pcolr(PS, x=x, y=y, title=title, name=name, sh=shw)
    
    return PS, x, y

def per_baseline(ms):
    t = pt.table(ms, readonly=True, ack=False)
    
    ts = t.query(sortlist='TIME',columns='TIME')
    firstTime = ts.getcell("TIME", 0)
    lastTime  = ts.getcell("TIME", t.nrows()-1)
    ts.close()
    intTime = t.getcell("INTERVAL", 0)
    print 'Integration time:\t%f sec' % (intTime)
    nTimeslots = (lastTime - firstTime) / intTime
    print 'Number of timeslots:\t%d' % (nTimeslots)
    timeslots = [0,nTimeslots]
    
    tant = pt.table(t.getkeyword('ANTENNA'), readonly=True, ack=False)
    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    numChannels = len(tsp.getcell('CHAN_FREQ',0))
    print 'Number of channels:\t%d' % (numChannels)
    print 'Reference frequency:\t%5.2f MHz' % (tsp.getcell('REF_FREQUENCY',0)/1.e6)

    antList = tant.getcol('NAME')
    antToPlot = range(2)

    tsel = t.query('TIME >= %f AND TIME <= %f AND ANTENNA1 IN %s AND ANTENNA2 IN %s' % \
                   (firstTime+timeslots[0]*intTime, firstTime+timeslots[1]*intTime, \
                    str(antToPlot), str(antToPlot) ) )
    x = np.linspace(0,(nTimeslots*intTime)/3600., 390)
    
    for tpart in tsel.iter(["ANTENNA1","ANTENNA2"]):
        ant1 = tpart.getcell("ANTENNA1", 0)
        ant2 = tpart.getcell("ANTENNA2", 0)
        if ant1 != ant2:
            print antList[ant1], antList[ant2]
            mc = tpart.getcol('MODEL_DATA_C')
            pl.plot(x, mc[:,0,0], 'r-', markersize=2, label='XX')
            pl.plot(x, mc[:,0,1], 'g-', markersize=2, label='XY')
            pl.plot(x, mc[:,0,2], 'b-', markersize=2, label='YX')
            pl.plot(x, mc[:,0,3], 'y-', markersize=2, label='YY')
            pl.xlim([0,x.max()])
            pl.legend()
            pl.show()
            pl.close()

def vis2ps(DIR, col, umin, umax, ubin, zbin, name, sTime=None, iTime=None):
    MSs = np.sort(os.listdir(DIR))
    sb = len(MSs)
    t = pt.table(DIR+MSs[sb/2], ack=False)
    
    # Time selection (optional)
    fT = t.getcell("TIME", 0)
    lT  = t.getcell("TIME", t.nrows()-1)
    iT = t.getcell("INTERVAL", 0)
    mT = (lT - fT)/2.
    T = t.getcol('TIME')
    if sTime=='first': ind = np.where(T<=(fT+iTime*3600.))[0]
    elif sTime=='last': ind = np.where(T>=(lT-iTime*3600.))[0]
    elif sTime=='mid':
        m0, m1 = fT+mT-((iTime/2.)*3600.), fT+mT+((iTime/2.)*3600.)
        ind = np.where(np.logical_and(T>=m0, T<=m1))[0]
    uvw = t.getcol('UVW')
    b = np.hypot(uvw[:,0], uvw[:,1])
    if sTime in ['first','mid','last']:
        b = b[ind]
        print 'Selected %s %s hours, %s columns'%(sTime,str(iTime),str(len(b)))

    V = np.zeros((sb,ubin,4), dtype=np.cfloat) # V array
    t = pt.table(DIR+MSs[0], ack=False)
    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    nu0 = tsp.getcell('REF_FREQUENCY',0)
    t.close()
    t = pt.table(DIR+MSs[1], ack=False)
    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    nu1 = tsp.getcell('REF_FREQUENCY',0)
    t.close()
    dNu = nu1 - nu0
    for i in range(sb):
        print MSs[i]
        t = pt.table(DIR+MSs[i], ack=False)
        tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
        f = tsp.getcell('REF_FREQUENCY',0)
        fact = sc.c/f
        Umin, Umax = umin*fact, umax*fact
        on = np.linspace(Umin, Umax, ubin+1)
        d = t.getcol(col)
        for j in range(ubin):
            ind = np.where(np.logical_and(b >= on[j], b <= on[j+1]))[0]
            for k in range(4):
                V[i,j,k] = np.mean(np.abs(d[ind,0,k]))
        t.close()
    #np.save('NCP_GRF_visibilities.npy', V)
    S = np.zeros((sb,ubin,5), dtype=np.cfloat) # S array
    S[:,:,0] = (V[:,:,0] + V[:,:,3])/2.
    S[:,:,1] = (V[:,:,0] - V[:,:,3])/2.
    S[:,:,2] = (V[:,:,1] + V[:,:,2])/2.
    S[:,:,3] = (V[:,:,1] - V[:,:,2])/2.*1j
    S[:,:,4] = np.abs(S[:,:,1] + (1j * S[:,:,2]))

    S_ft = np.fft.fft(S, axis=0) / sb
    S_ft = np.fft.fftshift(S_ft, axes=0)

    PS = np.abs(S_ft[sb/2:,:,:]) ** 2
    PS_full = np.abs(S_ft) ** 2

    np.save(name+'_PS.npy', PS)

    x, y, w = comoving_grid(umin, umax, xbin=ubin, zbin=sb/2, dNu=dNu, nu0=nu0, Nz=50)
    PSDL = np.zeros(PS.shape, dtype='float32')
    for k in range(5):
        for i in range(len(y)):
            for j in range(len(x)):
                PSDL[i,j,k] = (PS[i,j,k] * x[j]**2 * y[i]) / (2*np.pi)**2
    np.save(name+'_PSDL.npy', PSDL)

    x, y, w = comoving_grid(umin, umax, xbin=ubin, zbin=sb/2, dNu=dNu, nu0=nu0, Nz=50)
    PSDL_full = np.zeros(PS_full.shape, dtype='float32')
    for k in range(5):
        for i in range(len(y)):
            for j in range(len(x)):
                PSDL_full[i,j,k] = (PS_full[i,j,k] * x[j]**2 * y[i]) / (2*np.pi)**2
    np.save(name+'_PSDL_full.npy', PSDL_full)
    print "--> Saved as%s"%name+'_PSDL.npy'


def vis2ps_v2(DIR, col, umin, umax, ubin, zbin, name, sTime=None, iTime=None):
    MSs = np.sort(os.listdir(DIR))
    sb = len(MSs)
    t = pt.table(DIR+MSs[sb/2], ack=False)
    
    # Time selection (optional)
    fT = t.getcell("TIME", 0)
    lT  = t.getcell("TIME", t.nrows()-1)
    iT = t.getcell("INTERVAL", 0)
    mT = (lT - fT)/2.
    T = t.getcol('TIME')
    if sTime=='first': ind = np.where(T<=(fT+iTime*3600.))[0]
    elif sTime=='last': ind = np.where(T>=(lT-iTime*3600.))[0]
    elif sTime=='mid':
        m0, m1 = fT+mT-((iTime/2.)*3600.), fT+mT+((iTime/2.)*3600.)
        ind = np.where(np.logical_and(T>=m0, T<=m1))[0]
    uvw = t.getcol('UVW')
    b = np.hypot(uvw[:,0], uvw[:,1])
    if sTime in ['first','mid','last']:
        b = b[ind]
        print 'Selected %s %s hours, %s columns'%(sTime,str(iTime),str(len(b)))

    V = np.zeros((sb,ubin,4), dtype=np.cfloat) # V array
    t = pt.table(DIR+MSs[0], ack=False)
    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    nu0 = tsp.getcell('REF_FREQUENCY',0)
    t.close()
    t = pt.table(DIR+MSs[1], ack=False)
    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    nu1 = tsp.getcell('REF_FREQUENCY',0)
    t.close()
    dNu = nu1 - nu0
    S = np.zeros((d.shape[0],4), dtype=np.cfloat) # S array
    for i in range(sb):
        print MSs[i]
        t = pt.table(DIR+MSs[i], ack=False)
        tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
        f = tsp.getcell('REF_FREQUENCY',0)
        fact = sc.c/f
        Umin, Umax = umin*fact, umax*fact
        on = np.linspace(Umin, Umax, ubin+1)
        d = t.getcol(col)
        S[:,0] = (d[:,0,0] + d[:,0,3])/2.
        S[:,1] = (d[:,0,0] - d[:,0,3])/2.
        S[:,2] = (d[:,0,1] + d[:,0,2])/2.
        S[:,3] = S[:,1] + (1j * S[:,2])
        for j in range(ubin):
            ind = np.where(np.logical_and(b >= on[j], b <= on[j+1]))[0]
            for k in range(4):
                V[i,j,k] = np.mean(S[ind,k])
        t.close()
    #np.save('NCP_GRF_visibilities.npy', V)

    S_ft = np.fft.fft(V, axis=0) / sb
    S_ft = np.fft.fftshift(S_ft, axes=0)

    PS = np.abs(S_ft[sb/2:,:,:]) ** 2
    PS_full = np.abs(S_ft) ** 2

    np.save(name+'_PS.npy', PS)

    x, y, w = comoving_grid(umin, umax, xbin=ubin, zbin=sb/2, dNu=dNu, nu0=nu0, Nz=50)
    PSDL = np.zeros(PS.shape, dtype='float32')
    for k in range(4):
        for i in range(len(y)):
            for j in range(len(x)):
                PSDL[i,j,k] = (PS[i,j,k] * x[j]**2 * y[i]) / (2*np.pi)**2
    np.save(name+'_PSDL.npy', PSDL)

    x, y, w = comoving_grid(umin, umax, xbin=ubin, zbin=sb, dNu=dNu, nu0=nu0, Nz=50)
    PSDL_full = np.zeros(PS_full.shape, dtype='float32')
    for k in range(4):
        for i in range(len(y)):
            for j in range(len(x)):
                PSDL_full[i,j,k] = (PS_full[i,j,k] * x[j]**2 * y[i]) / (2*np.pi)**2
    np.save(name+'_PSDL_full.npy', PSDL_full)
    print "--> Saved as%s"%name+'_PSDL.npy'

def vis2ps_gridded(DIR, col, umin, umax, ubin, zbin, name, sTime=None, iTime=None):
    MSs = np.sort(os.listdir(DIR))
    sb = len(MSs)
    t = pt.table(DIR+MSs[sb/2], ack=False)
    
    # Time selection (optional)
    fT = t.getcell("TIME", 0)
    lT  = t.getcell("TIME", t.nrows()-1)
    iT = t.getcell("INTERVAL", 0)
    mT = (lT - fT)/2.
    T = t.getcol('TIME')
    if sTime=='first': ind = np.where(T<=(fT+iTime*3600.))[0]
    elif sTime=='last': ind = np.where(T>=(lT-iTime*3600.))[0]
    elif sTime=='mid':
        m0, m1 = fT+mT-((iTime/2.)*3600.), fT+mT+((iTime/2.)*3600.)
        ind = np.where(np.logical_and(T>=m0, T<=m1))[0]
    uvw = t.getcol('UVW')
    b = np.hypot(uvw[:,0], uvw[:,1])
    if sTime in ['first','mid','last']:
        b = b[ind]
        print 'Selected %s %s hours, %s columns'%(sTime,str(iTime),str(len(b)))

    V = np.zeros((sb,ubin,4), dtype=np.cfloat) # V array
    t = pt.table(DIR+MSs[0], ack=False)
    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    nu0 = tsp.getcell('REF_FREQUENCY',0)
    t.close()
    t = pt.table(DIR+MSs[1], ack=False)
    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    nu1 = tsp.getcell('REF_FREQUENCY',0)
    t.close()
    dNu = nu1 - nu0
    for i in range(sb):
        print MSs[i]
        t = pt.table(DIR+MSs[i], ack=False)
        tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
        f = tsp.getcell('REF_FREQUENCY',0)
        fact = sc.c/f
        Umin, Umax = umin*fact, umax*fact
        on = np.linspace(Umin, Umax, ubin+1)
        d = t.getcol(col)
        for j in range(ubin):
            ind = np.where(np.logical_and(b >= on[j], b <= on[j+1]))[0]
            for k in range(4):
                V[i,j,k] = np.mean(d[ind,0,k])
        t.close()
    #np.save('NCP_GRF_visibilities.npy', V)
    S = np.zeros((sb,ubin,5), dtype=np.cfloat) # S array
    S[:,:,0] = (V[:,:,0] + V[:,:,3])/2.
    S[:,:,1] = (V[:,:,0] - V[:,:,3])/2.
    S[:,:,2] = (V[:,:,1] + V[:,:,2])/2.
    S[:,:,3] = (V[:,:,1] - V[:,:,2])/2.*1j
    S[:,:,4] = S[:,:,1] + (1j * S[:,:,2])

    S_ft = np.fft.fft(S, axis=0) / sb
    S_ft = np.fft.fftshift(S_ft, axes=0)

    PS = np.abs(S_ft[sb/2:,:,:]) ** 2
    PS_full = np.abs(S_ft) ** 2

    np.save(name+'_PS.npy', PS)

    x, y, w = comoving_grid(umin, umax, xbin=ubin, zbin=sb/2, dNu=dNu, nu0=nu0, Nz=50)
    PSDL = np.zeros(PS.shape, dtype='float32')
    for k in range(5):
        for i in range(len(y)):
            for j in range(len(x)):
                PSDL[i,j,k] = (PS[i,j,k] * x[j]**2 * y[i]) / (2*np.pi)**2
    np.save(name+'_PSDL.npy', PSDL)

    x, y, w = comoving_grid(umin, umax, xbin=ubin, zbin=sb/2, dNu=dNu, nu0=nu0, Nz=50)
    PSDL_full = np.zeros(PS_full.shape, dtype='float32')
    for k in range(5):
        for i in range(len(y)):
            for j in range(len(x)):
                PSDL_full[i,j,k] = (PS_full[i,j,k] * x[j]**2 * y[i]) / (2*np.pi)**2
    np.save(name+'_PSDL_full.npy', PSDL_full)
    print "--> Saved as%s"%name+'_PSDL.npy'

#vis2ps_gridded('MS2/', 'MODEL_DATA_15d', 30., 175., 30, 25, 'vis2ps_dum', sTime=sTime,iTime=.5)

def vis2ps_grid(dr, col, umin, umax, grid, zbin=None, name=''):
    MSs = np.sort(os.listdir(dr))
    sb = len(MSs)
    stokes = 5

    # calculate nu0, dNu
    t = pt.table(dr+MSs[0], ack=False)
    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    nu0 = tsp.getcell('REF_FREQUENCY',0)
    t.close()
    t = pt.table(dr+MSs[1], ack=False)
    tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
    nu1 = tsp.getcell('REF_FREQUENCY',0)
    dNu = nu1 - nu0

    size = ((umax-grid/2)*2.)/grid+1
    ug = vg = np.linspace(-umax, umax, size)
    Sg = np.zeros((sb,size,size,5), dtype=np.cfloat) # S array
    Sga = np.zeros((sb,size,size,5), dtype='float32') # S array

    for i in range(sb):
        print MSs[i]
        t = pt.table(dr+MSs[i], ack=False)

        # frequency
        tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
        f = tsp.getcell('REF_FREQUENCY',0)
        fact = sc.c/f

        # Stokes visibilities
        V = t.getcol(col)[:,0,:]
        S = np.zeros((V.shape[0],5), dtype=np.cfloat) # S array
        S[:,0] = (V[:,0]+V[:,3])/2.
        S[:,1] = (V[:,0]-V[:,3])/2.
        S[:,2] = (V[:,1]+V[:,2])/2.
        S[:,3] = (V[:,1]-V[:,2])/2.*1j
        S[:,4] = S[:,1] + 1j*S[:,2]

        # UV cut
        Umin, Umax = umin*fact, umax*fact
        uvw = t.getcol('UVW')
        u, v = uvw[:,0], uvw[:,1]
        b = np.hypot(u,v)
        ind = np.where(np.logical_and(b>=Umin, b<=Umax))[0]
        S = S[ind,:]
        u, v = u[ind], v[ind]

        # Gridding
        ug, vg = ug*fact, vg*fact
        Sg_r = sint.griddata((u, v), S.real, (ug[None,:], vg[:,None]), method='nearest')
        Sg_i = sint.griddata((u, v), S.imag, (ug[None,:], vg[:,None]), method='nearest')
        #Sg_a = sint.griddata((u, v), np.abs(S), (ug[None,:], vg[:,None]), method='nearest')
        
        # Save in the full-SB array
        Sg[i,:,:,:] = Sg_r + 1j*Sg_i
        #Sga[i,:,:,:] = Sg_a
        t.close()

    S_ft = np.fft.fft(Sg, axis=0) / sb
    S_ft = np.fft.fftshift(S_ft, axes=0)

    # Create the 3D PS
    PS = np.abs(S_ft) ** 2
    np.save(name+'_PS_3D.npy', PS)
    print '--> 3D PS saved with shape:'
    print PS.shape

    # Create 2D PS
    # while radially averaging, lcut=2 to remove the inner 2*grid lambdas
    lcut = 2 # the central 2*2 pixels are < 32-lambda, thus rejected
    ubin = len(radial_profile(PS[0,:,:,0])[lcut:])
    PS_ubin = np.zeros((sb,ubin,stokes), dtype='float32')
    for s in range(stokes):
        for f in range(sb):
            PS_ubin[f,:,s] = radial_profile(PS[f,:,:,s])[lcut:]

    if zbin==None: zbin = len(PS_ubin[sb/2:,0,0])
    PS_uzbin = np.zeros((zbin,ubin,stokes))
    for s in range(stokes):
        for u in range(ubin):
            PS_uzbin[:,u,s] = binning(PS_ubin[:,u,s], zbin)
            #PS_uzbin[:,u,s] = PS_ubin[sb/2:,u,s]

    np.save(name+'_PS_2D_dimensional.npy', PS_uzbin)

    PS_uzbin_dl = np.zeros(PS_uzbin.shape, dtype='float32')
    x, y, w = comoving_grid(umin, umax, xbin=ubin, zbin=zbin, dNu=dNu, nu0=nu0, Nz=sb)
    for s in range(stokes):
        for f in range(zbin):
            for u in range(ubin):
                PS_uzbin_dl[f,u,s] = (PS_uzbin[f,u,s] * x[u]**2 * y[f]) / (2.*np.pi)**2

    np.save(name+'_PS_2D_dimensionless.npy', PS_uzbin_dl)
    print '--> 2D PS saved with shape:'
    print PS_uzbin_dl.shape

def grid_excon(dr, col, umin=30, umax=180, name=''):
    SBs = range(179,229)
    d = np.zeros((len(SBs),27,27,4), dtype=np.cfloat)
    for i in range(len(SBs)):
        MS = dr+str(SBs[i])+'_uv_002.MS.dppp'
        os.system('excon -m %s -c %s -w 32 -l %s -u %s -x 4 -p 34 -d 300'\
            % (MS, col, str(umin), str(umax)) )
        I = pf.getdata(MS+'_GR.fits')+1j*pf.getdata(MS+'_GI.fits')
        d[i,:,:,0] = I[0,0,137:164,137:164]
        Q = pf.getdata(MS+'_GRQ.fits')+1j*pf.getdata(MS+'_GIU.fits')
        d[i,:,:,1] = Q[0,0,137:164,137:164]
        U = pf.getdata(MS+'_GRU.fits')+1j*pf.getdata(MS+'_GIU.fits')
        d[i,:,:,2] = U[0,0,137:164,137:164]
        V = pf.getdata(MS+'_GRV.fits')+1j*pf.getdata(MS+'_GIV.fits')
        d[i,:,:,3] = V[0,0,137:164,137:164]
    np.save('Gridded_vis_%s_%s.npy'%(name,col), d)


def vis2ps_grid_excon(file, umin=30, umax=175, zbin=12, name=''):
    d = np.load(file)
    d[:,:,:,3] = d[:,:,:,1] + 1j*d[:,:,:,2]

    sb = d.shape[0]
    stokes = d.shape[3]
    S_ft = np.fft.fft(d, axis=0) / sb
    S_ft = np.fft.fftshift(S_ft, axes=0)

    # Create the 3D PS
    PS = np.abs(S_ft) ** 2
    np.save(name+'_PS_3D.npy', PS)
    print '--> 3D PS saved with shape:'
    print PS.shape

    # Create 2D PS
    # while radially averaging, lcut=2 to remove the inner 2*grid lambdas
    lcut = 2 # the central 2*2 pixels are < 32-lambda, thus rejected
    ucut = 11
    ubin = len(radial_profile(PS[0,:,:,0])[lcut:ucut])
    PS_ubin = np.zeros((sb,ubin,stokes), dtype='float32')
    for s in range(stokes):
        for f in range(sb):
            PS_ubin[f,:,s] = radial_profile(PS[f,:,:,s])[lcut:ucut]

    if zbin==None: zbin = len(PS_ubin[sb/2:,0,0])
    PS_uzbin = np.zeros((zbin,ubin,stokes))
    for s in range(stokes):
        for u in range(ubin):
            PS_uzbin[:,u,s] = binning(PS_ubin[:,u,s], zbin)
            #PS_uzbin[:,u,s] = PS_ubin[sb/2:,u,s]

    np.save(name+'_PS_2D_dimensional.npy', PS_uzbin)

    PS_uzbin_dl = np.zeros(PS_uzbin.shape, dtype='float32')
    x, y, w = comoving_grid(umin, umax, xbin=ubin, zbin=zbin, dNu=0.2e6, nu0=150e6, Nz=sb)
    for s in range(stokes):
        for f in range(zbin):
            for u in range(ubin):
                PS_uzbin_dl[f,u,s] = (PS_uzbin[f,u,s] * x[u]**2 * y[f]) / (2.*np.pi)**2

    np.save(name+'_PS_2D_dimensionless.npy', PS_uzbin_dl)
    print '--> 2D PS saved with shape:'
    print PS_uzbin_dl.shape

    return
    PS = PS_uzbin_dl
    L = np.sqrt(PS[:,:,0]/PS[:,:,3])*1e2
    print L.max(), np.median(L)
    plots.pcolr(L, x=x, y=y, title='', name='', sh=True, vmin=1e-1,vmax=0.7)


#vis2ps_grid_excon(file='Gridded_vis_3C196_MODEL_DATA_15d_C.npy')
