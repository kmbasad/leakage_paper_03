#!/usr/bin/env python
import sys
sys.path.append('../')
import os, pyfits as pf, numpy as np, plots
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt, cm, patches
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_grid
import power_spectrum as ps

pre = '/Users/kmbasad/Data/paper_03/'
dp = '../paper_03_latex/'

def power_spectra(filename):
    Q4=pre+filename+'Q_4d.fits'
    U4=pre+filename+'U_4d.fits'
    Q9=pre+filename+'Q_9d.fits'
    U9=pre+filename+'U_9d.fits'
    
    #q4, x, y = ps.threeD(cube=Q4, name=Q4[:-5])
    #u4, x, y = ps.threeD(cube=U4, name=U4[:-5])
    #q9, x, y = ps.threeD(cube=Q9, name=Q9[:-5])
    #u9, x, y = ps.threeD(cube=U9, name=U9[:-5])

    q4 = np.load(Q4[:-5]+'.npy')
    u4 = np.load(U4[:-5]+'.npy')
    q9 = np.load(Q9[:-5]+'.npy')
    u9 = np.load(U9[:-5]+'.npy')
    x, y, w = np.load(Q4[:-5]+'_x.npy'), np.load(Q4[:-5]+'_y.npy'), np.load(Q4[:-5]+'_w.npy')

    labels = ['$Q,4^\circ$', '$U,4^\circ$', '$Q,9^\circ$', '$U,9^\circ$']
    if 'NCP' in filename: name='ncp'
    elif '3C196' in filename: name='3c196'

    norm = LogNorm(vmin=5e2, vmax=5e5)
    cmap = cm.cubehelix_r
    fig = plt.figure(1, figsize=(11,4))
    ax = axes_grid.ImageGrid(fig, 111, nrows_ncols=(1,4), axes_pad=0.05, add_all=True,\
                            share_all=True, aspect=False, label_mode='L', \
                            cbar_mode='single',cbar_size='8%',cbar_pad=0.05, cbar_location='right')

    X, Y = np.meshgrid(x,y)
    im = ax[0].pcolormesh(X,Y,q4, cmap=cmap, norm=norm)
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[1].pcolormesh(X,Y,u4, cmap=cmap, norm=norm)
    ax[2].pcolormesh(X,Y,q9, cmap=cmap, norm=norm)
    ax[3].pcolormesh(X,Y,u9, cmap=cmap, norm=norm)
    print q4.min(), q4.max()

    # Colorbar
    fs = 9 # fontsize
    cb = ax.cbar_axes[0].colorbar(im, format='%.2f')
    cb.ax.minorticks_on()
    cb.ax.set_yscale('log')
    cb.ax.set_ylabel('Power $\\Delta^2$ [mK]$^2$')
    cb.ax.set_yticks([1e3,1e4,1e5])
    cb.ax.set_yticklabels(['$10^3$', '$10^4$','$10^5$'], fontsize=fs)

    for i in range(4):
        ax[i].plot(x,w[0], color='k', linestyle='solid', linewidth=.6)
        ax[i].plot(x,w[1], color='k', linestyle='dashed', linewidth=.6)

        ax[i].set_xlabel('$k_\perp$  [Mpc$^{-1}$]')
        ax[i].text(.023,1.5,'%s: %s'%(filename[4:-8], labels[i]))
        ax[i].set_ylabel('$k_\parallel$ [Mpc$^{-1}$]')
        ax[i].set_yticks([.1,1])
        ax[i].set_yticklabels([.1,1], fontsize=fs)
        ax[i].set_xticks([.03,.1,.5])
        ax[i].set_xticklabels([.03,.1,.5], fontsize=fs)
        ax[i].set_xlim([x.min(), x.max()])
        ax[i].set_ylim([y.min(), y.max()])
        if i!=0: ax[i].get_yaxis().set_visible(False)
    #plt.show()
    plt.savefig("%sobserved_ps_%s_2d.pdf"%(dp,name), bbox_inches='tight')
    plt.close()

    return

    # ================
    # Power spectra averaged over some k_perp values
    # ================
    fig = plt.figure(figsize=(5,4))

    q4m = np.mean(q4[:,0:4],axis=1)
    u4m = np.mean(u4[:,0:4],axis=1)
    q9m = np.mean(q9[:,0:4],axis=1)
    u9m = np.mean(u9[:,0:4],axis=1)

    plt.plot(y,q4m, 'g-', label=labels[0], markersize=5)
    plt.plot(y,u4m, 'g--', label=labels[1], markersize=5)
    plt.plot(y,q9m, 'b-', label=labels[2], markersize=5)
    plt.plot(y,u9m, 'b--', label=labels[3], markersize=5)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e3,2e4])
    plt.legend()
    plt.xlabel('$k_\parallel$ Mpc$^{-1}$')
    plt.ylabel('Power $\\Delta^2$ [mK]$^2$')
    plt.xticks([.1,1], ['0.1', '1.0'])
    plt.grid()
    #plt.savefig(dp+'observed_ps_%s_1d.pdf'%name, bbox_inches='tight')
    #plt.show()
    plt.close()

#power_spectra(filename='obs/NCP_L86762_')
#power_spectra(filename='obs/3C196_L80508_')

def calculate_ps(Q,U='',P='',name=''):
    """
    The old version that calculates both cyl and sphe
    deprecated on 7 Apr 2017
    """
    kq, q = ps.threeD(cube=Q, title='$Q$', name=name+'_Q',vmin=.3e3,vmax=2e5)
    ku, u = ps.threeD(cube=U, title='$U$', name=name+'_U',vmin=.3e3,vmax=2e5)
    kp, p = ps.threeD(cube=P, title='$|Q+iU|$', name=name+'_P',vmin=.3e3,vmax=2e5)
    #np.save(name+'_kp.npy', kp)
    #return
    pl.plot(np.sort(kq), q, 'ro-', label='$Q$')
    pl.plot(np.sort(ku), u, 'bo-', label='$U$')
    pl.plot(np.sort(kp), p, 'go-', label='$Q+iU$')
    pl.xscale('log')
    pl.yscale('log')
    #xlim(min(k), max(K))
    pl.xlabel('$k$ [Mpc$^{-1}$]', fontsize=16)
    pl.ylabel('$\\Delta^2(k)$ [mK]$^2$', fontsize=16)
    pl.ylim([1e2,1e7])
    pl.xlim([4e-2,1.5e0])
    pl.grid()
    pl.legend(loc='upper left')
    pl.savefig(name+'_SPS_dl.pdf', dpi=80, bbox_inches="tight")
    #pl.show()
    pl.close()
    
    #PS2, x, y = PS.threeD(fits='Pcube_10MHz_zeroF.fits')

#calculate_ps(Q=pre+'cubes/L86762_natural10-800-sc2cl_Qcube_10MHz_4d_K.fits', U=pre+'cubes/L86762_natural10-800-sc2cl_Ucube_10MHz_4d_K.fits', P=pre+'cubes/L86762_natural10-800-sc2cl_Pcube_10MHz_4d_K.fits', name=pre+'plots/'+'NCP_4d')
#calculate_ps(d2+'cubes/3C196_L80508_Qcube_10MHz_4d_K.fits', d2+'cubes/3C196_L80508_Ucube_10MHz_4d_K.fits', d2+'cubes/3C196_L80508_Pcube_10MHz_4d_K.fits',name='3C196_4d')
#calculate_ps(Q=d1+'cubes/L86762_natural10-800-sc2cl_Qcube_10MHz_9d_K.fits', U=d1+'cubes/L86762_natural10-800-sc2cl_Ucube_10MHz_9d_K.fits', P=d1+'cubes/L86762_natural10-800-sc2cl_Pcube_10MHz_9d_K.fits',name='NCP_9d')
#calculate_ps(d2+'cubes/3C196_L80508_Qcube_10MHz_9d_K.fits', d2+'cubes/3C196_L80508_Ucube_10MHz_9d_K.fits', d2+'cubes/3C196_L80508_Pcube_10MHz_9d_K.fits',name='3C196_9d')

def ps_ratio(path):
    N4=path+'NCP_L86762_P_4d.fits'
    N9=path+'NCP_L86762_P_9d.fits'
    T4=path+'3C196_L80508_P_4d.fits'
    T9=path+'3C196_L80508_P_9d.fits'
    
    #n4, k = ps.threeD(cube=N4, name=N4[:-5])
    #n9, k = ps.threeD(cube=N9, name=N9[:-5])
    #t4, k = ps.threeD(cube=T4, name=T4[:-5])
    #t9, k = ps.threeD(cube=T9, name=T9[:-5])
    #sys.exit()

    n4 = np.load(N4[:-5]+'_spherical.npy')
    n9 = np.load(N9[:-5]+'_spherical.npy')
    t4 = np.load(T4[:-5]+'_spherical.npy')
    t9 = np.load(T9[:-5]+'_spherical.npy')

    x, y = np.load(N4[:-5]+'_x.npy'), np.load(N4[:-5]+'_y.npy')
    k4 = np.load(N4[:-5]+'_spherical_K.npy')
    k9 = np.load(N9[:-5]+'_spherical_K.npy')
    k34 = np.load(T4[:-5]+'_spherical_K.npy')
    k39 = np.load(T9[:-5]+'_spherical_K.npy')
    lim = len(x[x<0.1])
    print x[:lim]

    n4m = (np.sqrt(n4) * 0.27e-2)**2
    n9m = (np.sqrt(n9) * 0.27e-2)**2
    t4m = (np.sqrt(t4) * 0.35e-2)**2
    t9m = (np.sqrt(t9) * 0.35e-2)**2
    #t4m = (np.sqrt( np.mean(t4[:,0:lim],axis=1) ) * 0.34e-2)**2
    #t9m = (np.sqrt( np.mean(t9[:,0:lim],axis=1) ) * 0.36e-2)**2

    labels = ['NCP, $4^\circ$', 'NCP, $9^\circ$', '3C196, $4^\circ$', '3C196, $9^\circ$']
    plt.plot(k4,n4m, 'og-', label=labels[0], markersize=5)
    plt.plot(k9,n9m, 'og--', label=labels[1], markersize=5)
    plt.plot(k34,t4m, 'ob-', label=labels[2], markersize=5)
    plt.plot(k39,t9m, 'ob--', label=labels[3], markersize=5)
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim([1e3,2e4])
    plt.legend(ncol=2)
    plt.xlabel('$k$ [Mpc$^{-1}$]')
    plt.ylabel('Power $\\Delta^2$ [mK]$^2$')
    #plt.xticks([.1,1], ['0.1', '1.0'])
    #plt.yticks([1e-2,5e-2], ['$10^{-2}$',5e-2])
    plt.grid()
    plt.savefig(dp+'predicted_leakage_ps.pdf', bbox_inches='tight')
    #plt.show()
    plt.close()

ps_ratio(pre+'obs/')

def copy_MS():
    d = '/data1/users/lofareor/khan/P03/'
    l = open(d+'L86762_avg1ch10s_BBS.ref').readlines()[5:]
    x = [i.strip().split() for i in l]
    for i in range(179,180):
        print i
        os.system('cp -r /net/'+x[i][3]+'/'+x[i][0]+' '+d+'NCP/MS/')

def cube_10MHz(Q,U):
    # 10 deg, 1200 pix, 0.5 arcmin, 302 SBs
    q = pf.getdata(Q)
    u = pf.getdata(U)
    h = pf.getheader(Q)
    fov = h['NAXIS1'] * h['CDELT2']
    nu0, dnu = h['CRVAL3'], h['CDELT3']
    #return
    # take 191 to 241 SBs, and inner 600 pixels for 5 deg
    nu0 = nu0 + dnu * 191
    h['CRVAL3'] = nu0
    q9, u9 = q[190:240,60:1140,60:1140], u[190:240,60:1140,60:1140]
    q4, u4 = q[190:240,360:840,360:840], u[190:240,360:840,360:840]
    pf.writeto(Q[:-5]+'_10MHz_9d.fits', q9, h, clobber=True)
    pf.writeto(U[:-5]+'_10MHz_9d.fits', u9, h, clobber=True)
    pf.writeto(Q[:-5]+'_10MHz_4d.fits', q4, h, clobber=True)
    pf.writeto(U[:-5]+'_10MHz_4d.fits', u4, h, clobber=True)
    return
    t = open('freq.txt', 'w')
    for i in range(50):
        t.write(str(nu0 + dnu * i)+'\n')
    t.close()

#cube_10MHz(d1+'cubes/L86762_natural10-800-sc2cl_Qcube.fits', d1+'cubes/L86762_natural10-800-sc2cl_Ucube.fits')
#cube_10MHz(d2+'cubes/3C196_L80508_Qcube.fits', d2+'cubes/3C196_L80508_Ucube.fits')


def beam():
    M = np.load('NCP_L86762_SB189.45deg.image_M.npy')
    I = M[79:129,:,:,0,0]
    L = np.sqrt(M[79:129,:,:,0,1]**2 + M[79:129,:,:,0,2]**2)
    PSi, x, y = cyl_ps(cube=I, pixel=21.)
    PSl, x, y = cyl_ps(cube=L, pixel=21.)
    plots.pcolr((PSl/PSi)*1e2, x=x, y=y, title='', name='beam_leakage_PS.pdf', sh=False)

def desystematize(fits,fro,to):
    P0 = pf.getdata(fits)
    h = pf.getheader(fits)
    P_ft = np.fft.fft(P0, axis=0)
    P_fts = np.fft.fftshift(P_ft, axes=0)
    P_fts[fro:to+1,:,:] = P_fts[0,:,:]
    P = np.fft.ifft(P_fts, axis=0)
    pf.writeto(fits[:-5]+'_0f.fits', P.real, h, clobber=True)
    #plots.imshw(np.mean(P0, axis=0)*1e3, sh=True, name=fits[:-5]+'.pdf', vmin=-14,vmax=4, \
    #            title='Polarized flux [mJy]')
    #plots.imshw(np.mean(P.real, axis=0)*1e3, sh=True, name=fits[:-5]+'_FT.pdf', vmin=-0.8, vmax=0.8, \
    #            title='Polarized flux [mJy]')

#desystematize('NCP_Q_5deg_10MHz.fits',25,26)
#desystematize('NCP_U_5deg_10MHz.fits',25,26)
#desystematize('L86762_natural10-800-sc2cl_Qcube.fits',149,153)
#desystematize('3C196_L80508_Qcube_10MHz.fits',25,26)
#desystematize('3C196_L80508_Ucube_10MHz.fits',25,26)
