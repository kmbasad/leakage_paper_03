#!/usr/bin/env python
import pyfits as pf, radial_profile as rp, scipy.integrate as si, sys, \
       scipy.constants as sc, numpy as np, plots, sys, matplotlib.pyplot as pl
from sys import argv
from numpy import *
from pylab import *
sys.path.append('../')

def radial_profile(data, nbins, c=None, lcut=0, ucut=0):
    y, x = np.indices(data.shape)
    if not c: c = map(int, [x.max()/2., y.max()/2.])
    r = np.hypot(x - c[0], y - c[1])
    sum, e = np.histogram(r, weights=data, bins=nbins)
    count, edges = np.histogram(r, bins=nbins)
    PS = sum/count
    l = len(e[e <= lcut])
    if ucut==0: ucut=data.shape[0]/2.
    u = len(e[e <= ucut])
    return PS[l:u], e[l:u]

def binning(d,n,c=None):
    x = np.indices(d.shape)
    if not c: c = int(x.max()/2.)
    r = abs(x-c)
    s,e = histogram(r,weights=d.reshape(1,len(d)),bins=n)
    c,e = histogram(r,bins=n)
    return s/c, e

def spherical_profile(data, nbins, c=None, lcut=0, ucut=0):
    z, y, x = np.indices(data.shape)
    if not c: c = map(int, [x.max()/2., y.max()/2., z.max()/2.])
    r = np.sqrt((x-c[0])**2+(y-c[1])**2+(z-c[2])**2)
    sum, e = np.histogram(r, weights=data, bins=nbins)
    count, edges = np.histogram(r, bins=nbins)
    PS = sum/count
    l = len(e[e <= lcut])
    if ucut==0: ucut=data.shape[0]
    u = len(e[e <= ucut])
    print len(PS), len(PS[l:u])
    return PS[l:u], e[l:u]

def twoD(inp,h=None, nbin=50,n=0,umin=30.,umax=800.,sh=False,square=True, name=None):
    try:
        im = pf.getdata(inp)[0,n,:,:]
        h = pf.getheader(inp)
    except: im, h = inp, h
    im_ft = np.fft.fftn(im)
    im_ft = np.fft.fftshift(im_ft) / sqrt(im.size)
    #im_ft = np.fft.fftshift(im_ft)
    # cut upto the longest baseline/shortest scale
    N = im.shape[0]
    px = abs(h['CDELT1'])*(pi/180.)
    fmin = 1/(px*N) # minimum spatial frequency
    fmax = 1/(px*2)
    c = N/2
    ucut = umax/fmin
    lcut = umin/fmin
    #lcut=0
    if square==True: ps2d = np.abs(im_ft)**2
    else: ps2d = np.abs(im_ft)
    ps1d = radial_profile(ps2d, nbin, lcut=lcut, ucut=ucut)[0]

    # calculate comoving grid
    nu = h['RESTFRQ']
    kmin = 2*pi*umin/comoving(1420e6/nu-1)[0]
    kmax = 2*pi*umax/comoving(1420e6/nu-1)[0]
    k = np.linspace(kmin, kmax, len(ps1d))
    figure(figsize=(5,4))
    plot(k, ps1d,'o-', color='black')
    xscale('log')
    yscale('log')
    xlim(min(k),max(k))
    xlabel('$k$')
    ylabel('$P(k)$ [mK]$^2$')
    grid()
    if sh==True: show()
    else: savefig('/home/users/khan/plots/'+name, bbox_inches='tight', dpi=40)
    close()
    return ps1d, k

def threeD(cube, header='', name='', ft_factor=1, title='', xbin=30, zbin=20, kbin=30, umin=30., \
           umax=800., st='', shw=False, pixel=1.17, B='', dNu='', nu0='', vmin='', vmax=''):
    # open file and perform 3D FFT, shift large scales to the center and normalize
    try:
        fits = pf.open(cube)
        cube, h = fits[0].data, fits[0].header
        try: cd1 = h['CDELT1']
        except: h = pf.getheader('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.fits')
    except: cube, h = cube, header
    cube = cube*1e3 # convert from K to mK
    ft_shape = (cube.shape[0], cube.shape[1]*ft_factor, cube.shape[2]*ft_factor)
    print '--> Performing n-D Fourier transform with shape %s' % str(ft_shape)
    cube_ft = np.fft.fftn(cube, ft_shape) / sqrt(cube.size) # nD FFT
    cube_ft = np.fft.fftshift(cube_ft) # shift to center and normalize
    #np.save('cube_ft.npy', cube_ft)
    #cube_ft = np.load('cube_ft.npy')

    # check if the FT went okay
    # http://stackoverflow.com/questions/19444373/normalization-of-2d-fft
    #pp_ratio = np.sum(np.abs(cube)**2) / np.sum(np.abs(cube_ft)**2)
    vm_ratio = std(cube)**2 / mean(abs(cube_ft)**2)
    #print '--> Power-power ratio = %f'% pp_ratio
    print '--> Variance-integral ratio = %f'% vm_ratio

    #exit()

    # cut upto the longest baseline/shortest scale
    N = cube_ft.shape
    try: fmin = 1/(abs(h['CDELT1'])*(pi/180.)*N[1]) # spatial resolution corresponding to the pixel size
    except: fmin = 1/(abs(pixel/60.)*(pi/180.)*N[1])
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
    #np.save(name+'PS.npy', ps_uzbin)

    # Calculate comoving grid and plot dimensional cyl PS
    print '--> Calculating power spectrum'
    x, y, w1 = comoving_grid(umin=umin, umax=umax, xbin=xbin, zbin=zbin, h=h, fov=3.5)
    x, y, w2 = comoving_grid(umin=umin, umax=umax, xbin=xbin, zbin=zbin, h=h, fov=45.)
    w=[]
    w.append(w1)
    w.append(w2)
    #except: x, y, w = comoving_grid(h, umin=umin, umax=umax, xbin=xbin, zbin=zbin, B=B, dNu=dNu, nu0=nu0)
    PS = ps_uzbin
    #plots.pcolr(PS, x=x, y=y, title=title+': [mK$^2$Mpc$^{-3}$]', name=name+'_CPS.pdf', sh=shw, \
    #            vmin=vmin, vmax=vmax)

    # Calculate and plot dimensionless cyl PS
    
    print '--> Calculating dimensionless power spectrum x^2*y/(2*pi)^2'
    PSDL = np.zeros(PS.shape, dtype='float32')
    for i in range(len(x)):
        for j in range(len(y)):
            #PSDL[j,i] = PS[j,i] * 2*pi * abs(x[1]-x[0]) * x[i] * abs(y[1]-y[0])
            PSDL[j,i] = (PS[j,i] * x[i]**2 * y[j]) / (2*pi)**2
            #k = hypot(x[i], y[j])
            #PSDL[j,i] = PS[j,i] * (k**3/(2*pi**2))
    np.save(name+'.npy', PSDL)
    np.save(name+'_x.npy', x)
    np.save(name+'_y.npy', y)
    np.save(name+'_w.npy', w)
    print PSDL.shape, x.shape, y.shape
    #return PSDL, x, y

    #plots.pcolr(PSDL, x=x, y=y, title=title+': $\\Delta^2$ [mK]$^2$', name=name+'_CPS_dl.pdf', sh=shw, \
    #            vmin=vmin, vmax=vmax,w=w)
    # integral-variance ratio
    vm_ratio = np.std(cube)**2 / np.mean(PSDL)
    print '--> CylPS Variance-mean ratio = %f'% vm_ratio
    #return

    # spherical PS
    
    print '--> Calculating spherical PS using %s bins' % zbin
    ps_sph = spherical_profile(cube_ft, nbins=zbin)[0]
    ps_sph_dl = np.zeros((len(ps_sph)), dtype='float32')
    x, y, w = comoving_grid(umin=umin, umax=umax, xbin=zbin, zbin=zbin, h=h)
    #except: x, y, w = comoving_grid(h, umin=umin, umax=umax, xbin=kbin, zbin=kbin, B=B, dNu=dNu, nu0=nu0)
    K = []
    for i in range(len(ps_sph)):
        k = hypot(x[i], y[i])
        K.append(k)
        ps_sph_dl[i] = ps_sph[i] * k**3/(2*pi**2)
    np.save(name+'_spherical.npy', ps_sph_dl)
    np.save(name+'_spherical_K.npy', K)
    
    vm_ratio = np.std(cube)**2 / np.mean(ps_sph_dl)
    print '--> SphPS Variance-mean ratio = %f'% vm_ratio

    return K, ps_sph_dl

def threeD_wedge():
    p = np.load('3C196_diffuse_P.npy')
    h = pf.getheader('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.fits')
    zbin, xbin = p.shape[0], p.shape[1]
    x, y, w = comoving_grid(h, umin=30., umax=800., xbin=xbin, zbin=zbin)
    print w.min(), w.max()
    #plots.pcolr(p, x=x, y=y, title='', name='dum.pdf', sh=True, w=w)

def comoving(z,dNu=0.2e6):
    H_0 = 100.*1e3 # Hubble constant at z=0 in units of h where h = H_0/100 km/s/Mpc
    H_0 = 67.8*1e3
    omega_m = 0.3 # matter density
    omega_l = 0.7 # cosmological constant
    DH = sc.c/H_0 # Hubble distance in Mpc
    fHI = sc.c/0.21
    Ezinv = lambda rs: 1./sqrt(omega_m*(1+rs)**3+omega_l) # inverse dimensionless Hubble parameter function
    Dc = DH * si.quad(Ezinv, 0, z)[0]

    H_z = H_0 * sqrt(omega_m*(1+z)**3 + omega_l) # Hubble constant at redshift z
    kzE = 2*pi*H_z*fHI/(sc.c*(1+z)**2)

    return Dc, kzE

def wedge(ku,z,fov):
    H_0 = 67.8*1e3
    omega_m = 0.3 # matter density
    omega_l = 0.7 # cosmological constant
    DH = sc.c/H_0 # Hubble distance in Mpc
    Ez = sqrt(omega_m*(1+z)**3+omega_l)
    DM, kzE = comoving(z)
    wedge = fov * ((DM*Ez)/(DH*(1+z))) * ku
    return wedge

def comoving_grid(umin, umax, xbin, zbin, h='', B='', dNu='', nu0='', Nz='', nu_avg=None, fov=4.):
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
    kr[0], kr[1] = 2*pi*umin/DM, 2*pi*umax/DM
    kr[2], kr[3] = kzE/B, kzE/dNu
    
    # calculate x and z grid
    k_perp = linspace(kr[0], kr[1], xbin)
    k_para = linspace(kr[2], kr[3], zbin)

    # the wedge
    k_wedge = wedge(k_perp, z, fov*(pi/180.) )

    #plot(k_perp, k_wedge)
    #xscale('log')
    #yscale('log')
    #show()

    return k_perp, k_para, k_wedge

def ps3d_3C196():
    I = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.fits')
    In = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.noisy.fits')
    Q = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Qcube.K.fits')
    Qn = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Qcube.K.noisy.fits')
    U = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Ucube.K.fits')
    Un = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Ucube.K.noisy.fits')
    V = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Vcube.K.fits')
    P = Q + 1j*U
    Pn = Qn + 1j*Un
    h = pf.getheader('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.fits')
    #comoving_grid(h, 30., 800., 50, 20)
    threeD(P, h, name='3C196_diffuse_P')
    #threeD(Pn, h, name='3C196_diffuse_P_noisy')
    #threeD(I, h, name='3C196_diffuse_I')
    #threeD(In, h, name='3C196_diffuse_I_noisy')
    #threeD(V, h, name='3C196_diffuse_V')

def ps3d_3C196_3deg():
    I = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.fits')[:,60:420,60:420]
    In = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.noisy.fits')[:,60:420,60:420]
    Q = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Qcube.K.fits')[:,60:420,60:420]
    Qn = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Qcube.K.noisy.fits')[:,60:420,60:420]
    U = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Ucube.K.fits')[:,60:420,60:420]
    Un = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Ucube.K.noisy.fits')[:,60:420,60:420]
    V = pf.getdata('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Vcube.K.fits')[:,60:420,60:420]
    P = Q + 1j*U
    Pn = Qn + 1j*Un
    h = pf.getheader('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.fits')
    #comoving_grid(h, 30., 800., 50, 20)
    threeD(P, h, name='3C196_3deg_diffuse_P')
    #threeD(Pn, h, name='3C196_3deg_diffuse_P_noisy')
    #threeD(I, h, name='3C196_3deg_diffuse_I')
    #threeD(In, h, name='3C196_3deg_diffuse_I_noisy')
    #threeD(V, h, name='3C196_3deg_diffuse_V')

def ps3d_NCP():
    Q = pf.getdata('L86762_natural10-800-sc2cl_Qcube.fits')
    h = pf.getheader('L86762_natural10-800-sc2cl_Qcube.fits')
    threeD(Q, h, name='NCP_PS.pdf')

def gmca():
    I = pf.getdata('gmca_res_I_noiseless.fits')
    In = pf.getdata('gmca_res_I_noisy.fits')
    h = pf.getheader('L80273_SAP000_SB099_uv.MS.dppp.1ch.dppp.real.Icube.K.fits')
    threeD(I, h, name='gmca_3C196', title='Post-GMCA')
    threeD(In, h, name='gmca_3C196_noisy', title='Post-GMCA noise')

if __name__=='__main__':
    #threeD(cube='rest_GMCA_Res_04sources.fits', pixel='11.94', umin=30, umax=250, xbin=10, zbin=10, kbin=20)
    #sys.exit()
    if argv[1]=='3c': ps3d_3C196()
    elif argv[1]=='3c3': ps3d_3C196_3deg()
    elif argv[1]=='2d': twoD(inp=argv[2],n=int(argv[3]))
    elif argv[1]=='g': gmca()
    elif argv[1]=='w': threeD_wedge()
    elif argv[1]=='sta': threeD(argv[2], shw=True)
    elif argv[1]=='ncp': ps3d_NCP()
    else: threeD(cube=argv[1], umin=float(argv[2]), umax=float(argv[3]), pixel=float(argv[4]), \
                 xbin=int(argv[5]), zbin=int(argv[6]), kbin=int(argv[7]), name='PS', B=10e6, dNu=0.5e6, \
                 nu0=170e6)
