#!/usr/bin/env python
import numpy as np, pyfits as pf, time, sys
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid import axes_grid
import scipy.stats as ss
from sys import argv
import power_spectrum_cyl as psc, plots, glob
#import GRF_sim_plot as gsp

# make 64^3 GRF using Vibor's code in my laptop, copy the cube into target
# make sure you have 50 MS of NCP: SB179 to SB228
# average the MSs down to 2 minutes interval for faster calculation
# add CORRECTED_DATA columns to the 2min MSs, CASA might not work
# make an image of 18.75arcmin resolution, 50 pix, 15 deg
# copy the GRF into the image to make a model for prediction
# make a BBS model using the CASA .image file
# predict using the specified MS and SKY. remove flags before predicting

d1 = '/data1/users/lofareor/khan/P03/NCP/'
d2 = '/data1/users/lofareor/khan/P03/3C196/'
d3 = '/data1/users/lofareor/khan/P03/3C295/'
dp = '/home/users/khan/plots/'
dp = '../paper03_latex/'
path = '/Users/kmbasad/Data/paper_03/sim/'

def MS_timeavg_newcol(inp, out):
    MSs = np.sort(os.listdir(inp))
    for i in range(1):
        msin = inp+MSs[i]
        print msin
        msout = out+MSs[i]
        bbs.timeavg(msin, msout, '12')

#MS_timeavg_newcol(inp='/data3/users/lofareor/khan/3C196/ms/50SBs/', out='/data1/users/lofareor/khan/P03/3C196/MS2/')
#MS_timeavg_newcol(inp='/data3/users/lofareor/khan/3C295/L104068_BSAGE/50SBs/', out='/data1/users/lofareor/khan/3C295/MS2/')

def img_for_bbs_model(MS, IMG, imsize, cell):
    #bbs.newcol(MS, 'MODEL_DATA', 'CORRECTED_DATA')
    #pt.taql('update %s set FLAG=False'%MS)
    #pt.taql('update %s set FLAG_ROW=False'%MS)
    casa.Image(ms=MS, img=IMG, col='CORRECTED_DATA', uvrange='30~1000lambda', imsize=imsize,\
               cell=cell, st='IQUV')

#img_for_bbs_model(MS='MS2/L80273_SAP000_SB179_uv.MS.dppp.1ch.dppp', IMG='3C196_image_15d', imsize='50', cell='18.75arcmin')
#img_for_bbs_model(MS='MS2/L80273_SAP000_SB179_uv.MS.dppp.1ch.dppp', IMG='3C196_image_9d', imsize='28', cell='18.75arcmin')
#img_for_bbs_model(MS='MS2/L80273_SAP000_SB179_uv.MS.dppp.1ch.dppp', IMG='3C196_image_4d', imsize='16', cell='15arcmin')
#img_for_bbs_model(MS='MS2/L86762_SAP000_SB185_uv_002.MS.dppp', IMG='NCP_image_15d', imsize='50', cell='18.75arcmin')
#img_for_bbs_model(MS='MS2/L86762_SAP000_SB185_uv_002.MS.dppp', IMG='NCP_image_9d', imsize='28', cell='18.75arcmin')
#img_for_bbs_model(MS='MS2/L86762_SAP000_SB185_uv_002.MS.dppp', IMG='NCP_image_4d', imsize='16', cell='15arcmin')

def make_model(Q, U):
    q, u = pf.getdata(Q), pf.getdata(U)
    pf.writeto(Q[:-7]+'28.fits', q[:,11:39,11:39], clobber=True)
    pf.writeto(U[:-7]+'28.fits', u[:,11:39,11:39], clobber=True)
    pf.writeto(Q[:-7]+'16.fits', q[:,17:33,17:33], clobber=True)
    pf.writeto(U[:-7]+'16.fits', u[:,17:33,17:33], clobber=True)

#make_model('fsim_modelB_Q_50.fits', 'fsim_modelB_U_50.fits')

def bbs_models(im, Q, U, sk):
    global d
    n = pf.getdata(im+'.fits').shape[2]
    im = pim.image(im+'.image')
    Q = pf.getdata(Q)
    U = pf.getdata(U)
    # n=50 for 15 deg, and 28 for 9 deg, and 16 for 4 deg (15arcmin)
    nu0, dNu = 150000000, 195312
    for i in range(50):
        print i+1
        sb = str(179+i)
        t = open(sk+sb+'.txt', 'w')
        nu = nu0 + i*dNu
        t.write("# (Name, Type, RA, DEC, I, Q, U, V, MajorAxis, MinorAxis, Orientation, ReferenceFrequency='"+str(nu)+"', SpectralIndex='[0.]') = format\n")
        c = 0
        for j in range(n):
            for k in range(n):
                c += 1
                f,ch,dec,ra = im.toworld((0,0,j,k))
                ra = sky.hms2deg(ra*(180./np.pi))
                dec = sky.dms2deg(dec*(180./np.pi))
                t.write('S'+str(c)+', POINT, '+str(ra)+', '+str(dec)+', 0, '+str(Q[i,j,k])+', '+\
                    str(U[i,j,k])+', 0, 0, 0, 0, '+str(nu)+'., 0.\n')
        t.close()

#bbs_models('model/3C196_image_4d', 'model/fsim_modelB_Q_16.fits', 'model/fsim_modelB_U_16.fits', 'sky/4deg/SB')
#bbs_models('model/3C196_image_9d', 'model/fsim_modelB_Q_28.fits', 'model/fsim_modelB_U_28.fits', 'sky/9deg/SB')
#bbs_models('model/3C196_image_15d', 'model/fsim_modelB_Q_50.fits', 'model/fsim_modelB_U_50.fits', 'sky/15deg/SB')

#bbs_models('model/NCP_image_4d', 'model/fsim_modelB_Q_16.fits', 'model/fsim_modelB_U_16.fits', 'sky/4deg/SB')
#bbs_models('model/NCP_image_9d', 'model/fsim_modelB_Q_28.fits', 'model/fsim_modelB_U_28.fits', 'sky/9deg/SB')
#bbs_models('model/NCP_image_15d', 'model/fsim_modelB_Q_50.fits', 'model/fsim_modelB_U_50.fits', 'sky/15deg/SB')

def bbs_predict_correct(msd, skd, col):
    MSs = np.sort(os.listdir(msd))
    Skys = np.sort(os.listdir(skd))
    t1 = time.time()
    for i in range(len(MSs)):
        MS = msd+MSs[i]
        #bbs.newcol(MS, 'DATA', 'CORRECTED_DATA')
        sky = skd+Skys[i]
        print MS, sky
        bbs.predict(ms=MS, col=col, sky=sky, G='F', E='T', baselines='CS* &')
        bbs.correct(ms=MS, icol=col, ocol=col+'_C', sky=sky, G='F', E='T', \
                    baselines='CS* &', chunksize='195')
    t2 = time.time()
    print '---> Time taken: %s' % str((t2-t1)/60.)[0:4]

#bbs_predict_correct(d2+'MS2/', d2+'sky/4deg/', 'MODEL_DATA_4d')
#bbs_predict_correct(d2+'MS2/', d2+'sky/9deg/', 'MODEL_DATA_9d')
#bbs_predict_correct(d2+'MS2/', d2+'sky/15deg/', 'MODEL_DATA_15d')

#bbs_predict_correct(d1+'MS2/', d1+'sky/4deg/', 'MODEL_DATAN_4d')
#bbs_predict_correct(d1+'MS2/', d1+'sky/9deg/', 'MODEL_DATAN_9d')
#bbs_predict_correct(d1+'MS2/', d1+'sky/15deg/', 'MODEL_DATAN_15d')

def casa_image():
    global d
    MSs = np.sort(os.listdir(d+'NCP/MS2/'))
    t1 = time.time()
    for i in range(len(MSs)):
        print MSs[i]
        MS = d+'NCP/MS2/'+MSs[i]
        img = d+'GRF/10Jy_9deg/'+MSs[i]
        casa.Image(ms=MS, img=img, col='MODEL_DATA_10Jy_9d_C', uvrange='30~1000lambda', imsize='28',\
                   cell='1.875arcmin', st='IQUV', weight='uniform')
    t2 = time.time()
    print '---> Time taken: %s' % str((t2-t1)/60.)[0:4]
#casa_image()

def model_ps():
    q = pf.getdata('fsim_modelB_Q_50.fits')
    u = pf.getdata('fsim_modelB_U_50.fits')
    h = pf.getheader('3C196_image_15d.fits')
    p = np.abs(q+1j*u)
    PS.twoD(p[24,:,:]*1e3, h=h, nbin=30, umin=30., umax=200., sh=False, name='FG_3C196_model_1dps.pdf')
    #plots.imshw(q[24,:,:]*1e3, sh=False, name='FG_3C196_model_Q.pdf', title='$Q$ [mK]')
    #plots.imshw(u[24,:,:]*1e3, sh=False, name='FG_3C196_model_U.pdf', title='$U$ [mK]')
    return
    f = np.linspace(150,160,p.shape[0])
    pl.figure(figsize=(5,4))
    pl.plot(f, q[:,24,24]*1e3, label='$Q$')
    pl.plot(f, u[:,24,24]*1e3, label='$U$')
    pl.legend()
    pl.xlabel('Frequency [MHz]')
    pl.ylabel('Polarized intensity [mK]')
    pl.ylim([-62,62])
    pl.savefig(dp+'FG_3C196_model_fspectra.pdf', bbox_inches='tight', dpi=30)

#model_ps()


# -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- #


def grid(dr, col, name):
    degs = ['_15d', '_9d', '_4d']
    for i in range(3):
        coln = col+degs[i]+'_C'
        psc.grid_excon(dr, coln, name=name)

#grid('3C196/MS2/L80273_SAP000_SB', 'MODEL_DATA', '3C196')
#grid('NCP/MS2/L86762_SAP000_SB', 'MODEL_DATAN', 'NCP')

def power_spectra_vis(file, col, name, umin=30.,umax=175.,grid=16.,zbin=10):
    degs = ['_15d', '_9d', '_4d']
    for i in range(3):
        coln = col+degs[i]+'_C'
        namen = name+degs[i]
        psc.vis2ps_grid_excon(file=file+'%s.npy'%coln, name=namen)
        #psc.vis2ps_grid(msd, coln, umin,umax,grid,zbin, name=namen)
    return

#power_spectra_vis('Gridded_vis_3C196_', 'MODEL_DATA', d2+'PS/FG_3C196_excon')
#power_spectra_vis('Gridded_vis_NCP_', 'MODEL_DATAN', d1+'PS/FG_NCP_excon')
#power_spectra_vis(d2+'MS2/', 'MODEL_DATA', d2+'PS/FG_3C196_zbinned')
#power_spectra_vis(d1+'MS2/', 'MODEL_DATAN', d1+'PS/FG_NCP_zbinned')


# -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- #


def power_spectra_vis_pcolr_grid(name, cbar=True,xtl=True,ytl=True,sh=False,dr=''):
    fig = pl.figure(1, figsize=(10,10))
    ax = axes_grid.ImageGrid(fig, 111, nrows_ncols=(3,3), axes_pad=0.05, add_all=True,\
                            share_all=True, aspect=False, label_mode='L', \
                            cbar_mode='none')
    degs = ['15d', '9d', '4d']
    fs = 12
    c = 0
    for i in range(len(degs)):
        PS = np.load(dr+name+degs[i]+'_PS_2D_dimensionless.npy')
        I, P = PS[:,:,0], PS[:,:,3]
        L = np.sqrt(I/P) *100.
        print PS.shape
        ubin, zbin = PS.shape[1], PS.shape[0]
        nu0, dNu = 150.e6, 0.2e6
        x, y, w = psc.comoving_grid(30.,175., ubin, zbin, dNu, nu0, 50)
        X, Y = np.meshgrid(x,y)
        unit = '[mJy]$^2$'
        #d = np.zeros(P.shape+(3,), dtype='float')
        cmap = pl.cm.cubehelix_r
        imP = ax[3*i+0].pcolormesh(X,Y,P*1e6,cmap=cmap, norm=LogNorm(vmin=.5e-3,vmax=5e1))
        imI = ax[3*i+1].pcolormesh(X,Y,I*1e6,cmap=cmap, norm=LogNorm(vmin=.5e-8,vmax=5e-4))
        imL = ax[3*i+2].pcolormesh(X,Y,L,cmap=cmap, norm=LogNorm(vmin=.5e-1,vmax=5e0))
        
        fs = 9
        if i==0:
            cax = plots.make_colorbar(ax[0], fs)
            cb = pl.colorbar(imP, cax=cax, orientation='horizontal')
            #cb.set_ticks([1e-9,1e-7,1e-5])
            cb.set_label('$P_P$: $\Delta^2$ [mK]$^2$')
            cb.ax.xaxis.set_label_position('top')
            cb.ax.xaxis.set_ticks_position('top')
            
            cax = plots.make_colorbar(ax[1], fs)
            cb = pl.colorbar(imI, cax=cax, orientation='horizontal')
            #cb.set_ticks([1e-14,1e-12,1e-10])
            cb.set_label('$P_I$ $\Delta^2$ [mK]$^2$')
            cb.ax.xaxis.set_label_position('top')
            cb.ax.xaxis.set_ticks_position('top')
            
            cax = plots.make_colorbar(ax[2], fs)
            cb = pl.colorbar(imL, cax=cax, orientation='horizontal')
            cb.set_label('$L_I$ [$\%$]')
            cb.ax.xaxis.set_label_position('top')
            cb.ax.xaxis.set_ticks_position('top')
            cb.ax.minorticks_off()

        for j in range(3):
            ax[c].set_yticks([.2,.6,1,1.4,1.8])
            ax[c].set_xlabel('$k_\perp$ [Mpc$^{-1}$]')
            ax[c].text(0.026,1.65,degs[i][:-1]+'$^\circ$, '+name[3:-7], color='black')
            ax[c].set_ylabel('$k_\parallel$ [Mpc$^{-1}$]')
            c+=1

    pl.savefig(dr+name+'PSDLs.pdf', bbox_inches='tight')
    #pl.show()
    pl.close()

#power_spectra_vis_pcolr_grid('FG_NCP_excon_', dr=path)
#power_spectra_vis_pcolr_grid('FG_3C196_excon_', dr=path)
#power_spectra_vis_pcolr_grid('FG_3C196_zbinned_', dr=d2)

# -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- #



def power_spectra_vis_hist_grid(name,dr,sh=True):
    fig = pl.figure(1, figsize=(13,5))
    g = axes_grid.ImageGrid(fig, 111, nrows_ncols=(1,3), axes_pad=0.05, add_all=True,\
                            share_all=True, aspect=False, label_mode='L', cbar_mode='none')

    degs = ['15d', '9d', '4d']
    dlab = ['$15^{\circ}$', '$9^{\circ}$', '$4^{\circ}$']
    fs = 12
    stats = np.zeros((4,6))
    for i in range(3):
        PS = np.load(dr+name+degs[i]+'_PS_3D.npy')
        I, P = PS[:,:,:,0], PS[:,:,:,3]
        #print np.sqrt(np.mean(I)/np.mean(P))*100.
        #I, P = PS[:,:,0], PS[:,:,4]
        L = np.sqrt(I/P)*100.
        L = L[np.logical_not(np.isnan(L))]
        if '3C196' in name: Ln = L; nbins = 4000
        elif 'NCP' in name: Ln = L[L<10.]; nbins = 1000
        n, bins, patches = g[i].hist(Ln, alpha=0.5, bins=nbins, color='steelblue', normed=True, histtype='stepfilled')
        mx = bins[np.where(n==n.max())][0]
        stats[:,i] = (np.mean(L), np.median(L), np.std(L), mx)
        #g[i].axvline(mx, color='r', linestyle='dashed', linewidth=1, label=dlab[i]+', Peak: $%.2f$'%mx)
        g[i].axvline(np.median(L), color='g', linestyle='dotted', linewidth=1, label=dlab[i]+', Median = $%.2f$'%np.median(L))
        #g[i].axvline(L.mean(), color='b', linestyle='dashed', linewidth=1, label='Mean = $%.2f$'%np.mean(L))
    
        # Fit Rice distribution to the RMS histogram
        #L1 = Ln[Ln<0.7]
        #b, l, s = ss.rice.fit(L1)
        #pdf = ss.rice.pdf(bins, b, loc=l, scale=s)
        #m, v = ss.rice.stats(b, loc=l, scale=s, moments='mv')
        #med = ss.rice.median(b, loc=l, scale=s)
        #mx = bins[np.where(pdf==pdf.max())][0]
        #stats[:,3+i] = (m, med, np.sqrt(v), mx)
        #mad = np.median(np.abs(med-L1))
        
        #med, mad = str(round(med,2)), str(round(mad,2))
        #pdf = pdf*(n.max()/pdf.max())
        #g[i].plot(bins, pdf, color='k', label='Rice: ($%s, \ %s$)'%(med,mad))

        # Fit Pareto distribution with shape parameter 1
        b, l, s = ss.pareto.fit(Ln, 1)
        pdf = ss.pareto.pdf(bins, b, l, s)
        median, mean = ss.pareto.median(b,l,s), ss.pareto.mean(b,l,s)
        g[i].plot(bins, pdf, color='black', linestyle='dashed', label='Fitted Pareto distribution')
        g[i].axvline(median, color='r', linestyle='dotted', linewidth=1, label='Pareto median = $%.2f$'%np.median(L))
        print median, mean
        
        g[i].legend(prop={'size':12})
        g[i].set_xlabel('$L_I$ [$\%$]', fontsize=fs+1)
        g[i].set_ylabel('Frequency of occurence', fontsize=fs)
        g[i].set_xticks([1,2,3,4,5])
        #g[i].set_yticks([.2,.6,1,1.4])
        g[i].set_xlim(-0.01,5.2)
        #g[i].set_ylim(0,1.7)

        # Find Median absolute deviation
        
    if sh==False: pl.savefig(dp+name+'L_gridded_hist.pdf', bbox_inches='tight')
    else: pl.show()
    pl.close()
    return
    stats = np.around(stats,2)
    t = open(name+'stats.txt', 'w')
    labels = ['Mean', 'Median', 'Standard deviation', 'Peak']
    for i in range(4):
        t.write('\multicolumn{1}{|l|}{%s}'%labels[i])
        for j in range(6): t.write(' & '+str(stats[i,j]))
        t.write(' \\\ \n')
    t.close()

    return

power_spectra_vis_hist_grid('FG_NCP_excon_', dr=path, sh=False)
power_spectra_vis_hist_grid('FG_3C196_excon_', dr=path, sh=False)

#power_spectra_vis_hist_grid('FG_3C196_linear_', dr=d2, sh=True)


# -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- # -- #



def awi_predict_image():
    MSs = np.sort(os.listdir('NCP/MS/'))
    t1 = time.time()
    for i in [10]:
        MS = 'NCP/MS/'+MSs[i]
        img = 'GRF/awmodel/'+MSs[i]
        make_model(MS, img, i)
        awi_predict(MS, img)
        store_data(MS, 'MODEL_DATA', 'MODEL_DATA_GRF', 'mdata1')
        img = 'GRF/awimg/'+MSs[i]
        casa.Image(ms=MS, img=img, col='MODEL_DATA_GRF', uvrange='30~800lambda', imsize='480',\
                   cell='2.5arcmin', st='IQUV')
    t2=time.time()

def awi_predict(ms, img):
    pt.taql('UPDATE %s set MODEL_DATA=0+0i' % ms)
    awi.Image(ms=ms, col='MODEL_DATA', img=img, npix='480', cell='2.5arcmin', uvdist='0.03~3klambda', \
              twe='0.2',tw='300.0', st='IQUV', niter='0', op='predict', element='1', model=img, pbc='1e-6')

def make_model(ms, img, SB=191):
    awi.Image(ms=ms, col='CORRECTED_DATA', img=img, npix='480', cell='2.5arcmin', st='IQUV', \
              op='image', niter='0', wp='2', element='0')
    Q = pf.getdata('GRF/GRF_480_pi-2.7_Q.fits')[SB,:,:]
    U = pf.getdata('GRF/GRF_480_pi-2.7_U.fits')[SB,:,:]

    mod = pim.image(img)
    mod_data = mod.getdata()
    data = np.zeros(mod_data.shape)
    data[0,1,:,:] = Q
    data[0,2,:,:] = U
    mod.putdata(data)

def store_data(ms, icol, ocol, name):
    t = pt.table(ms, readonly=False)
    try: t.removecols(ocol)
    except: None
    m = t.getcol(icol)
    coldmi = t.getdminfo(icol)
    coldmi["NAME"] = name
    t.addcols (pt.maketabdesc(pt.makearrcoldesc(ocol, 0.,  valuetype='complex', shape=[1,4])),coldmi)
    print '--> Created the column %s' % ocol
    t.putcol(ocol, m)
    print '--> Copied data from %s column to %s column' % (icol, ocol)
    t.close()

def create_column(ms, col):
    t = pt.table('../L86762_SAP000_SB180_uv_002.MS.dppp')
    coldmi = t.getdminfo('DATA')
    t.close()
    t = pt.table(ms, readonly=False)
    try: t.removecols(col)
    except: None
    coldmi["NAME"] = 'cdata1'
    t.addcols (pt.maketabdesc(pt.makearrcoldesc(col, 0.,  valuetype='complex', shape=[1,4])),coldmi)
    t.close()





def power_spectra_vis_pcolr(name, cbar=True,xtl=True,ytl=True,sh=False,dr=''):
    fig = pl.figure(1, figsize=(10,13))
    g = axes_grid.ImageGrid(fig, 111, nrows_ncols=(3,3), axes_pad=0.1, add_all=True,\
                            share_all=True, aspect=False, label_mode='L', \
                            cbar_mode='none')
    degs = ['15d', '9d', '4d']
    fs = 12
    c = 0
    for i in range(len(degs)):
        PS = np.load(dr+'PS/'+name+degs[i]+'_PS_2D_dimensionless.npy')
        I, P = PS[:,:,0], PS[:,:,4]
        L = np.sqrt(I/P) *100.
        print PS.shape
        ubin, zbin = PS.shape[1], PS.shape[0]
        nu0, dNu = 150.e6, 190000.
        x, y, w = psc.comoving_grid(30.,175., ubin, zbin, dNu, nu0, 50)
        X, Y = np.meshgrid(x,y)
        unit = '[mJy]$^2$'
        #d = np.zeros(P.shape+(3,), dtype='float')
        cmap = pl.cm.rainbow
        imP = g[3*i+0].pcolormesh(X,Y,P*1e3,cmap=cmap, norm=LogNorm(vmin=1e-10,vmax=3e-4))
        imI = g[3*i+1].pcolormesh(X,Y,I*1e3,cmap=cmap, norm=LogNorm(vmin=1e-15,vmax=5e-9))
        imL = g[3*i+2].pcolormesh(X,Y,L,cmap=cmap, norm=LogNorm(vmin=.2e-1,vmax=1e1))
        if i==0:
            cax = plots.make_colorbar(g[0], fs)
            cb = pl.colorbar(imP, cax=cax, orientation='horizontal')
            cb.set_ticks([1e-9,1e-7,1e-5])
            cb.set_label('$P_P$: $\Delta^2$ [mJy]$^2$')
            cb.ax.xaxis.set_label_position('top')
            cax = plots.make_colorbar(g[1], fs)
            cb = pl.colorbar(imI, cax=cax, orientation='horizontal')
            cb.set_ticks([1e-14,1e-12,1e-10])
            cb.set_label('$P_I$ $\Delta^2$ [mJy]$^2$')
            cb.ax.xaxis.set_label_position('top')
            cax = plots.make_colorbar(g[2], fs)
            cb = pl.colorbar(imL, cax=cax, orientation='horizontal')
            cb.set_label('$L_I$ [$\%$]')
            cb.ax.xaxis.set_label_position('top')
        for j in range(3):
            g[c].set_xscale('log')
            g[c].set_yscale('log')
            g[c].set_xlim([x.min(), x.max()])
            g[c].set_ylim([y.min(), y.max()])
            g[c].set_xlabel('$k_\perp$ [Mpc$^{-1}$]')
            c+=1
        kp = 'deg$^2$) $k_\parallel$ [Mpc$^{-1}$]'
        g[0].set_ylabel('($15^2$%s'%kp)
        g[3].set_ylabel('($9^2$%s'%kp)
        g[6].set_ylabel('($4^2$%s'%kp)
    #pl.savefig(dp+name+'PSDLs.pdf', dpi=80, bbox_inches='tight')
    pl.show()
    pl.close()

#power_spectra_vis_pcolr('FG_NCP_', dr=d1)
#power_spectra_vis_pcolr('FG_3C196_', dr=d2)


def power_spectra_vis_hist(name,dr):
    fig = pl.figure(1, figsize=(13,5))
    g = axes_grid.ImageGrid(fig, 111, nrows_ncols=(1,3), axes_pad=0.05, add_all=True,\
                            share_all=True, aspect=False, label_mode='L', cbar_mode='none')

    degs = ['15d', '9d', '4d']
    dlab = ['$15^{2\circ}$', '$9^{2\circ}$', '$4^{2\circ}$']
    fs = 12
    for i in range(3):
        print degs[i]
        PS = np.load(dr+'PS/'+name+degs[i]+'_PS_3D.npy')
        I, Q, U, P = PS[:,:,0], PS[:,:,1], PS[:,:,2], PS[:,:,4]
        q = np.sqrt(I/Q)*100.
        L = np.sqrt(I/P)*100.
        Ln = L[L<3.5]
        n, bins, patches = g[i].hist(Ln, alpha=0.5, bins=100, color='grey', normed=True, histtype='stepfilled')
        mx = bins[np.where(n==n.max())][0]
        print np.mean(L), np.median(L), np.std(L), mx
        g[i].axvline(mx, color='r', linestyle='dashed', linewidth='1', label=dlab[i]+', Peak: $%s$'%str(mx)[:4])
        g[i].axvline(np.median(L), color='g', linestyle='dashed', linewidth='1', label='Median = $%s$'%str(np.median(L))[:4])
        g[i].axvline(L.mean(), color='b', linestyle='dashed', linewidth='1', label='Mean = $%s$'%str(np.mean(L))[:4])
    
        # Fit Rice distribution to the RMS histogram
        L1 = Ln[Ln<1.]
        b, l, s = ss.rice.fit(L1)
        pdf = ss.rice.pdf(bins, b, loc=l, scale=s)
        m, v = ss.rice.stats(b, loc=l, scale=s, moments='mv')
        med = ss.rice.median(b, loc=l, scale=s)
        mx = bins[np.where(pdf==pdf.max())][0]
        print m, med, np.sqrt(v), mx
        
        m, sd = str(round(m,2)), str(round(np.sqrt(v),2))
        pdf = pdf*(n.max()/pdf.max())
        g[i].plot(bins, pdf, color='k', label='Rice: ($%s, \ %s$)'%(m,sd))
        g[i].legend(prop={'size':12})
        g[i].set_xlabel('$L_I$ [$\%$]', fontsize=fs+1)
        g[i].set_ylabel('Frequency of occurence', fontsize=fs)
        g[i].set_xticks([0,1,2,3])
        #g[i].set_yticks([.2,.6,1,1.4])
        g[i].set_xlim(-.1,3.5)
        #g[i].set_ylim(0,1.7)
    #pl.savefig(dp+name+'L_hist.pdf', bbox_inches='tight', dpi=80)
    pl.show()
    pl.close()
    return

#power_spectra_vis_hist('FG_3C196_', dr=d2)
#power_spectra_vis_hist('FG_NCP_', dr=d1)
#power_spectra_vis_hist('FG_3C196_last_', dr=d2)
