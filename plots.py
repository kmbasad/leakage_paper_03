import numpy as np, os, sys, matplotlib
from mpl_toolkits.axes_grid import axes_grid
from mpl_toolkits.axes_grid import make_axes_locatable
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib.colors import LogNorm
from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as pl

def make_colorbar(g, fs, bba=(0, 0, 1, 1.07)):
    ax = g
    cax = inset_axes(ax,
                     width="100%", # width = 10% of parent_bbox width
                     height="5%", # height : 50%
                     loc=1,
                     bbox_to_anchor=bba,
                     bbox_transform=ax.transAxes,
                     borderpad=0
                     )
    mpl.artist.setp(getp(gca(),'xticklabels'), 'color', 'black', fontsize=fs)
    for i in ['top','bottom','left','right']: cax.spines[i].set_linewidth(0.01)
    return cax

def pcolr(d,x=[], y=[], ticks=None, fs=18, name='', title='', sh=True, cbmm=False,\
          vmin=None,vmax=None, cticks=None, cbar=True,xtl=True,ytl=True, w=None):
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    fig = figure(1)
    g = axes_grid.ImageGrid(fig, 111, nrows_ncols=(1,1), axes_pad=0.01, add_all=True,\
                            share_all=False, aspect=False, label_mode='1', cbar_mode='none')
    cmap = cm.rainbow
    X, Y = meshgrid(x,y)
    #print X.shape, Y.shape, d.shape
    im = g[0].pcolormesh(X,Y,d, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))

    if cbar==True:
        cax = make_colorbar(g[0], fs)
        cb = colorbar(im, cax=cax, orientation='horizontal')
        cb.set_label(title, fontsize=fs)
        if ticks!=None: cb.set_ticks(ticks)
        cb.ax.xaxis.set_label_position('top')

    g[0].set_xscale('log')
    g[0].set_yscale('log')
    g[0].set_xlim([x.min(), x.max()])
    g[0].set_ylim([y.min(), y.max()])

    #majorFormatter = FormatStrFormatter('%1.1g')
    #g[0].xaxis.set_major_formatter(majorFormatter)
    #g[0].yaxis.set_major_formatter(majorFormatter)
    #cb.ax.xaxis.set_major_formatter(majorFormatter)

    if xtl==True: g[0].set_xlabel('$k_\perp$ [Mpc$^{-1}$]', fontsize=fs+2)
    if ytl==True: g[0].set_ylabel('$k_\parallel$ [Mpc$^{-1}$]', fontsize=fs+2)
    g[0].set_xticks([1e-1], minor=False)
    g[0].set_xticklabels(['0.1'], minor=False, fontsize=fs+1)
    g[0].set_yticks([1e-1,1e0], minor=False)
    g[0].set_yticklabels(['0.1', '1.0'], minor=False, fontsize=fs+1)

    g[0].set_aspect(0.9)

    if w!=None:
        for i in range(len(w)): g[0].plot(x,w[i], linestyle='dashed')

    if sh==True: show()
    else:
        savefig(name, dpi=60, bbox_inches="tight")
        print '--> Saved as %s' % name
    close()

def image(d, ticks=[], fs=14, name='', title='', x=[], y=[], sh=True):
    fig = figure(1)
    g = axes_grid.ImageGrid(fig, 111, nrows_ncols=(1,1), axes_pad=0.01, add_all=True,\
                            share_all=False, aspect=False, label_mode='1', cbar_mode='none')
    im = g[0].imshow(d, origin='lower', cmap=cm.rainbow, norm=LogNorm())

    # colorbar
    cax = make_colorbar(g[0], fs-1)
    cb = colorbar(im, cax=cax, orientation='horizontal')
    cb.set_label(title, fontsize=fs-1)
    cb.ax.xaxis.set_label_position('top')

    # ticks
    xt = range(5,d.shape[1],6)
    g[0].set_xticks(xt)
    g[0].set_xticklabels(np.round(x[xt[:]],1), fontsize=fs-2)
    yt = range(5,d.shape[0],6)
    g[0].set_yticks(yt)
    g[0].set_yticklabels(np.round(y[yt[:]],1), fontsize=fs-2)

    # labels
    g[0].set_xlabel('$k_\perp$ (Mpc$^{-1}$)', fontsize=fs)
    g[0].set_ylabel('$k_\parallel$ (Mpc$^{-1}$)', fontsize=fs)

    if sh==True: show()
    else: savefig('/home/users/khan/plots/'+name, dpi=128, bbox_inches="tight")
    close()

def imshw(d, ticks=None, fs=14, name='', title='', x=[], y=[], sh=True, vmin=None, vmax=None):
    fig = figure(1)
    g = axes_grid.ImageGrid(fig, 111, nrows_ncols=(1,1), axes_pad=0.01, add_all=True,\
                            share_all=False, aspect=False, label_mode='1', cbar_mode='none')
    im = g[0].imshow(d, origin='lower', cmap=cm.rainbow, vmin=vmin, vmax=vmax)

    # colorbar
    cax = make_colorbar(g[0], fs-1)
    cb = colorbar(im, cax=cax, orientation='horizontal')
    cb.set_label(title, fontsize=fs-1)
    cb.ax.xaxis.set_label_position('top')
    if ticks!=None: cb.set_ticks(ticks)

    # ticks
    tks = [9,19,29,39]
    tls = [3,6,9,12]
    g[0].set_xticks(tks)
    g[0].set_xticklabels(tls)
    g[0].set_yticks(tks)
    g[0].set_yticklabels(tls)
    g[0].set_xlabel('$\\phi$ [degrees]')
    g[0].set_ylabel('$\\theta$ [degrees]')

    if sh==True: show()
    else: savefig('/home/users/khan/plots/'+name, dpi=80, bbox_inches="tight")
    close()

def imshw_multi(d, nrc=None, fs=14, name='', sh=True, vmin=None, vmax=None, ticks=None):
    if nrc==None: nrc = (d.shape[0]/2, d.shape[0]/2)
    ha = linspace(-d.shape[0]/2., d.shape[0]/2., d.shape[0])
    fig = figure(1, figsize=(10,4))
    g = axes_grid.ImageGrid(fig, 111, nrows_ncols=nrc, axes_pad=0., add_all=True,\
                            share_all=True, aspect=True, label_mode='L', cbar_mode='single', \
                            cbar_size='5%')
    if nrc==(1,1): it=1
    else: it = d.shape[0]
    for i in range(it):
        if nrc==(1,1): dn = d
        else: dn = d[i,:,:]
        im = g[i].imshow(dn, origin='lower', cmap=cm.rainbow,vmin=vmin,vmax=vmax)
        g[i].set_xticks([])
        g[i].set_yticks([])
        g[i].text(5,40,str(round(ha[i],1)), fontsize=13)
    cb = g[0].cax.colorbar(im)
    if ticks!=None: cb.ax.set_yticks(ticks)

    if sh==True: show()
    else: savefig('/home/users/khan/plots/'+name, dpi=80, bbox_inches="tight")
    close()

def line(x,y,styles,labels,xlbl,ylbl,fname,widths=[],xlimit=[],leg=['best',10,1,'',0.5,0.5],\
         scalex='log',scaley='log',ylimit=[],fs=11,lw=0.5,sh=False,\
         tickformat=[True,True], stformat='%1.2g'):
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    figure(num=None, figsize=(5,4))
    if widths==[]: widths = ones((len(y)))
    for i in range(len(y)):
        plot(x, y[i,:], styles[i], linewidth=widths[i], label=labels[i])
    yscale(scaley)
    xscale(scalex)
    if xlimit==[]: xlim([min(x), max(x)])
    else: xlim(xlimit)
    if ylimit!=[]: ylim(ylimit)
    lg = legend(loc=leg[0], prop={'size':leg[1]}, ncol=leg[2], title=leg[3], borderpad=leg[4],\
                labelspacing=leg[5])
    xlabel(xlbl, fontsize=fs)
    ylabel(ylbl, fontsize=fs)
    for j in ['top','bottom','left','right']:
        gca().spines[j].set_linewidth(lw)
    setp(getp(gca(),'yticklabels'), 'color', 'black', fontsize=fs)
    setp(getp(gca(),'xticklabels'), 'color', 'black', fontsize=fs)
    setp(lg.get_title(), fontsize=fs-1)
    lg.get_frame().set_linewidth(lw)
    majorFormatter = FormatStrFormatter(stformat)
    ax = gca()
    if tickformat[0]==True: ax.xaxis.set_major_formatter(majorFormatter)
    if tickformat[1]==True: ax.yaxis.set_major_formatter(majorFormatter)
    if sh==True: show()
    else: savefig('/home/users/khan/plots/'+fname, dpi=80, bbox_inches="tight")
    close()


"""
if i in [0,1,2]:
            cax = plots.make_colorbar(g[0], fs)
            cb = pl.colorbar(imP, cax=cax, orientation='horizontal')
            #cb.set_ticks([1e-9,1e-7,1e-5])
            cb.set_label('$P_P$: $\Delta^2$ [mK]$^2$')
            cb.ax.xaxis.set_ticks_position('top')
            cb.ax.xaxis.set_label_position('top')
            cax = plots.make_colorbar(g[1], fs)
            cb = pl.colorbar(imI, cax=cax, orientation='horizontal')
            #cb.set_ticks([1e-14,1e-12,1e-10])
            cb.set_label('$P_I$ $\Delta^2$ [mK]$^2$')
            cb.ax.xaxis.set_label_position('top')
            cax = plots.make_colorbar(g[2], fs)
            cb = pl.colorbar(imL, cax=cax, orientation='horizontal')
            cb.set_label('$L_I$ [$\%$]')
            cb.ax.xaxis.set_label_position('top')

"""
