'''
Run a slew of various plots to assess the quality of the resulting SDSS 
morphology catalogs I create with galaxy.py
'''

import numpy as np
from astropy.table import Table, vstack, join
from collections import OrderedDict, defaultdict
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pdb 
import pandas as pd

import matplotlib as mpl
mpl.rcParams.update({'font.size': 24, 
                    'font.family': 'STIXGeneral', 
                    'mathtext.fontset': 'stix',
                    'xtick.labelsize':16,
                    'ytick.labelsize':16,
                    'xtick.major.width':2,
                    'ytick.major.width':2,
                    'axes.linewidth':2,
                    'lines.linewidth':3,
                    'legend.fontsize':18})


def lotz_mergers(m20):
    return -0.14*m20+0.33

def lotz_separation(m20):
    return  0.14*m20+0.8

def compare_parameters(dat, fout='mycat'):
    cm = plt.cm.get_cmap('viridis')
    #fsmooth = "t01_smooth_or_features_a01_smooth_weighted_fraction"
    fsmooth = "t01_smooth_or_features_a01_smooth_debiased"
    colors = dat[fsmooth]

    fig = plt.figure(figsize=(8,8))
    gs = plt.GridSpec(4,4)
    gs.update(wspace=0.005, hspace=0.005)

    #'''
    ax10 = fig.add_subplot(gs[3,0])
    ax10.scatter(dat['A'], dat['C'], c=colors, marker='.', alpha=0.25, cmap=cm)
    ax10.set_xlabel(r'$A$')
    ax10.set_ylabel(r'$C$')
    ax10.set_ylim(0., 6.)
    ax10.set_xlim(0., .8)
    ax10.set_xticks([0.1, 0.3, 0.5, 0.7])
    ax10.set_yticks([1,2,3,4,5])
    #ax10.tick_params(axis='x', pad=8)
    #plt.text(0.5, 0.5, '10', fontsize=20, color='red', transform=ax10.transAxes)
    
    ################
    
    ax8 = fig.add_subplot(gs[3,1], sharey=ax10)
    plt.setp(ax8.get_yticklabels(), visible=False)
    ax8.scatter(dat['G'], dat['C'], c=colors, marker='.', alpha=0.25, cmap=cm)
    ax8.set_xlabel(r'$G$')
    #ax8.set_ylabel('C')
    ax8.set_xlim(0.2, 0.8)
    ax8.set_ylim(0., 6.)
    ax8.set_xticks([0.3, 0.5, 0.7])
    #plt.text(0.5, 0.5, '8', fontsize=20, color='red', transform=ax8.transAxes)
    
    ax9 = fig.add_subplot(gs[2,1], sharex=ax8)
    plt.setp(ax9.get_xticklabels(), visible=False)
    ax9.scatter(dat['G'], dat['A'], c=colors, marker='.', alpha=0.25, cmap=cm)
    #ax9.set_xlabel('G')
    ax9.set_ylabel(r'$A$')
    ax9.set_xlim(0.2, 0.8)
    ax9.set_ylim(0., .8)
    ax9.set_yticks([0.1, 0.3, 0.5, 0.7])
    #plt.text(0.5, 0.5, '9', fontsize=20, color='red', transform=ax9.transAxes)
    
    ################
    
    ax5 = fig.add_subplot(gs[3,2], sharey=ax10)
    plt.setp(ax5.get_yticklabels(), visible=False)
    ax5.scatter(dat['M20'], dat['C'], c=colors, marker='.', alpha=0.25, cmap=cm)
    ax5.set_xlabel(r'M$_{20}$')
    #ax5.set_ylabel('C')
    ax5.set_xlim(0., -3.5)
    ax5.set_ylim(0., 6.)
    ax5.set_xticks([-2.5, -1.5, -0.5])
    #plt.text(0.5, 0.5, '5', fontsize=20, color='red', transform=ax5.transAxes)
    
    ax6 = fig.add_subplot(gs[2,2], sharex=ax5, sharey=ax9)
    plt.setp((ax6.get_xticklabels(), ax6.get_yticklabels()), visible=False)
    ax6.scatter(dat['M20'], dat['A'], c=colors, marker='.', alpha=0.25, cmap=cm)
    #ax6.set_xlabel('M20')
    #ax6.set_ylabel('A')
    ax6.set_xlim(0., -3.5)
    ax6.set_ylim(0., .7)
    #plt.text(0.5, 0.5, '6', fontsize=20, color='red', transform=ax6.transAxes)
    
    ax7 = fig.add_subplot(gs[1,2], sharex=ax5)
    plt.setp(ax7.get_xticklabels(), visible=False)
    ax7.scatter(dat['M20'], dat['G'], c=colors, marker='.', alpha=0.25, cmap=cm)
    #ax7.set_xlabel('M20')
    ax7.set_ylabel(r'$G$')
    ax7.set_xlim(0., -3.5)
    ax7.set_ylim(0.2, 0.8)
    ax7.set_yticks([0.3, 0.5, 0.7])
    #plt.text(0.5, 0.5, '7', fontsize=20, color='red', transform=ax7.transAxes)
    
    ################
    
    ax1 = fig.add_subplot(gs[3,3], sharey=ax10)
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax1.scatter(dat['elipt'], dat['C'], c=colors, marker='.', alpha=0.25, cmap=cm)
    ax1.set_xlabel(r'$1-b/a$')
    #ax1.set_ylabel('C')
    ax1.set_ylim(0., 6.)
    ax1.set_xlim(0., .9)
    ax1.set_xticks([0.2, 0.4, 0.6, 0.8])
    #plt.text(0.5, 0.5, '1', fontsize=20, color='red', transform=ax1.transAxes)
    
    ax2 = fig.add_subplot(gs[2,3], sharex=ax1, sharey=ax9)
    plt.setp((ax2.get_xticklabels(), ax2.get_yticklabels()), visible=False)
    ax2.scatter(dat['elipt'], dat['A'], c=colors, marker='.', alpha=0.25, cmap=cm)
    #ax2.set_xlabel('elipt')
    #ax2.set_ylabel('A')
    ax2.set_ylim(0., .8)
    ax2.set_xlim(0., .9)
    #plt.text(0.5, 0.5, '2', fontsize=20, color='red', transform=ax2.transAxes)
    
    ax3 = fig.add_subplot(gs[1,3], sharex=ax1, sharey=ax7)
    plt.setp((ax3.get_xticklabels(),ax3.get_yticklabels()), visible=False)
    ax3.scatter(dat['elipt'], dat['G'], c=colors, marker='.', alpha=0.25, cmap=cm)
    #ax3.set_xlabel('elipt')
    #ax3.set_ylabel('G')
    ax3.set_xlim(0., .9)
    ax3.set_ylim(0.2, 0.8)
    #plt.text(0.5, 0.5, '3', fontsize=20, color='red', transform=ax3.transAxes)
    
    ax4 = fig.add_subplot(gs[0,3])
    plt.setp(ax4.get_xticklabels(), visible=False)
    im = ax4.scatter(dat['elipt'], dat['M20'], c=colors, marker='.', alpha=0.25, cmap=cm)
    #ax4.set_xlabel('elipt')
    ax4.set_ylabel(r'M$_{20}$')
    ax4.set_ylim(0., -3.5)
    ax4.set_xlim(0., .9)
    ax4.set_yticks([-2.5, -1.5, -0.5])
    #plt.text(0.5, 0.5, '4', fontsize=20, color='red', transform=ax4.transAxes)

    # add colorbar
    cbaxes = fig.add_axes([0.13, 0.87, 0.45, 0.02])
    cbax = fig.colorbar(mappable=im, cax=cbaxes, orientation='horizontal')
    cbax.ax.set_xlabel(r'$f_{\mathrm{smooth}}$')

    #gs.tight_layout(fig)
    plt.savefig('morph_params_'+fout+'.png', bbox_inches='tight')
    plt.show()
    plt.close()

def compare_parameters_density(dat, fout='mycat'):
    
    colors = 'black'

    fig = plt.figure(figsize=(8,8))
    gs = plt.GridSpec(4,4)
    gs.update(wspace=0.005, hspace=0.005)

    #'''
    ax10 = fig.add_subplot(gs[3,0])
    ax10.set_adjustable('box-forced')
    lims = [0, .8, 0, 6.]
    grid, extent, aspect = density(dat['A'], dat['C'], lims, 25)
    test=grid.T[grid.T !=0.]
    print np.min(test), np.max(test)
    #pdb.set_trace()
    im1 = ax10.imshow(grid.T, origin='low', extent=lims, 
                     interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                     aspect=aspect, cmap=plt.get_cmap('Blues'))
    ax10.set_xlabel('A', fontsize=20)
    ax10.set_ylabel('C', fontsize=20)
    ax10.set_ylim(0., 6.)
    ax10.set_xlim(0., .8)
    ax10.set_xticks([0.2, 0.4, 0.6])
    ax10.set_yticks([1,2,3,4])

    ################
    
    ax8 = fig.add_subplot(gs[3,1], sharey=ax10)
    ax8.set_adjustable('box-forced')
    plt.setp(ax8.get_yticklabels(), visible=False)
    lims = [0.2, .8, 0., 6.]
    grid, extent, aspect = density(dat['G'], dat['C'], lims, 25)
    im = ax8.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.), 
                    aspect=aspect, cmap=plt.get_cmap('Blues'))
    ax8.set_xlabel('G', fontsize=20)
    ax8.set_xlim(0.2, 0.8)
    ax8.set_ylim(0., 6.)
    ax8.set_xticks([0.4, 0.5, 0.6])
    
    ax9 = fig.add_subplot(gs[2,1], sharex=ax8)
    ax9.set_adjustable('box-forced')
    plt.setp(ax9.get_xticklabels(), visible=False)
    lims = [0.2, .8, 0., .8]
    grid, extent, aspect = density(dat['G'], dat['A'], lims, 25)
    im = ax9.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('Blues'))
    ax9.set_ylabel('A', fontsize=20)
    ax9.set_xlim(0.2, 0.8)
    ax9.set_ylim(0., .8)
    ax9.set_yticks([0.2, 0.4, 0.6])
    
    ################
    
    ax5 = fig.add_subplot(gs[3,2], sharey=ax10)
    ax5.set_adjustable('box-forced')
    plt.setp(ax5.get_yticklabels(), visible=False)
    lims = [-3.5, 0., 0., 6.]
    grid, extent, aspect = density(dat['M20'], dat['C'], lims, 25)
    im = ax5.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('Blues'))
    ax5.set_xlabel('M20', fontsize=20)
    ax5.set_xlim(0., -3.5)
    ax5.set_ylim(0., 6.)
    ax5.set_xticks([-2.5, -1.5, -0.5])
    
    ax6 = fig.add_subplot(gs[2,2], sharex=ax5, sharey=ax9)
    ax6.set_adjustable('box-forced')
    plt.setp((ax6.get_xticklabels(), ax6.get_yticklabels()), visible=False)
    lims = [-3.5, 0., 0., .8]
    grid, extent, aspect = density(dat['M20'], dat['A'], lims, 25)
    im = ax6.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=1./aspect, cmap=plt.get_cmap('Blues'))
    ax6.set_xlim(0., -3.5)
    ax6.set_ylim(0., .8)
    
    ax7 = fig.add_subplot(gs[1,2], sharex=ax5)
    ax7.set_adjustable('box-forced')
    plt.setp(ax7.get_xticklabels(), visible=False)
    lims = [-3.5, 0., 0.2, .8]
    grid, extent, aspect = density(dat['M20'], dat['G'], lims, 25)
    im = ax7.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=1./aspect, cmap=plt.get_cmap('Blues'))
    ax7.set_ylabel('G', fontsize=20)
    ax7.set_xlim(0., -3.5)
    ax7.set_ylim(0.2, 0.8)
    ax7.set_yticks([0.4, 0.5, 0.6])
    
    ################
    
    ax1 = fig.add_subplot(gs[3,3], sharey=ax10)
    ax1.set_adjustable('box-forced')
    plt.setp(ax1.get_yticklabels(), visible=False)
    lims = [0., .9, 0., 6.]
    grid, extent, aspect = density(dat['elipt'], dat['C'], lims, 25)
    im = ax1.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('Blues'))
    ax1.set_xlabel('elipt', fontsize=20)
    ax1.set_xlim(0., .9)
    ax1.set_ylim(0., 6.)
    ax1.set_xticks([0.2, 0.4, 0.6, 0.8])
    
    ax2 = fig.add_subplot(gs[2,3], sharex=ax1, sharey=ax9)
    ax2.set_adjustable('box-forced')
    plt.setp((ax2.get_xticklabels(), ax2.get_yticklabels()), visible=False)
    lims = [0., .9, 0., .8]
    grid, extent, aspect = density(dat['elipt'], dat['A'], lims, 25)
    im = ax2.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=1./aspect, cmap=plt.get_cmap('Blues'))
    ax2.set_xlim(0., .9)
    ax2.set_ylim(0., .8)
    
    ax3 = fig.add_subplot(gs[1,3], sharex=ax1, sharey=ax7)
    ax3.set_adjustable('box-forced')
    plt.setp((ax3.get_xticklabels(),ax3.get_yticklabels()), visible=False)
    lims = [0., .9, 0.2, .8]
    grid, extent, aspect = density(dat['elipt'], dat['G'], lims, 25)
    im = ax3.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=1./aspect, cmap=plt.get_cmap('Blues'))
    ax3.set_xlim(0., .9)
    ax3.set_ylim(0.2, 0.8)
    
    ax4 = fig.add_subplot(gs[0,3])
    ax4.set_adjustable('box-forced')
    plt.setp(ax4.get_xticklabels(), visible=False)
    lims = [0., .9, -3.5, 0.]
    grid, extent, aspect = density(dat['elipt'], dat['M20'], lims, 25)
    im = ax4.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('Blues'))
    ax4.set_ylabel('M20', fontsize=20)
    ax4.set_xlim(0., .9)
    ax4.set_ylim(0., -3.5)
    ax4.set_yticks([-2.5, -1.5, -0.5])

    # add colorbar
    cbaxes = fig.add_axes([0.13, 0.87, 0.45, 0.02])
    cbax = fig.colorbar(mappable=im1, cax=cbaxes, orientation='horizontal')
	#ticks=[-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1],
    cbax.ax.set_xlabel('number of galaxies', fontsize=16)

    plt.savefig('morph_params_'+fout+'.png')
    plt.show()
    plt.close()


def density(x, y, limits, numbins):

    notnans = np.where((~np.isnan(x) & ~np.isnan(y)))
    x, y = x[notnans], y[notnans]

    xsize = np.abs(limits[1]-limits[0])
    ysize = np.abs(limits[3]-limits[2])

    if xsize > ysize:
        aspect = ysize/xsize
    else:
        aspect = xsize/ysize

    gridx = np.linspace(limits[0], limits[1], numbins)
    gridy = np.linspace(limits[2], limits[3], numbins)

    grid, xedges, yedges = np.histogram2d(x, y, bins=[gridx, gridy])
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    #plt.figure()
    #plt.imshow(grid.T, origin='low', extent=extent, interpolation='nearest',
    #           aspect=aspect)
    #plt.show()
    return grid, extent, aspect


def gini_m20_compare(dat, fout):

    fig = plt.figure(figsize=(20,8))

    x1 = np.arange(-3.5, 1., .1)
    y1 = lotz_mergers(x1)
    x2 = np.arange(-3.5, -1.6, .1)
    y2 = lotz_separation(x2)

    ax1 = fig.add_subplot(121)
    lims = [-3.5, 1., .2, .8]
    grid, extent, aspect = density(dat['M20'], dat['G'], lims, 25)
    im = ax1.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=1./aspect, cmap=plt.get_cmap('Reds'))
    cbar = fig.colorbar(im)
    ax1.plot(x1, y1,'k',lw=2, label='Merger Line (Lotz 08)')
    ax1.plot(x2, y2,'k--',lw=2, label='Morph Line (Lotz 08)')
    ax1.set_title('Elliptical Apertures')
    ax1.set_ylabel('Gini', fontsize=14)
    ax1.set_xlabel('M20', fontsize=14)
    ax1.set_xlim(1., -3.5)
    ax1.set_ylim(.2, .8)

    ax2 = fig.add_subplot(122)
    lims = [-3.5, 1., .2, .8]
    grid, extent, aspect = density(dat['M20_c'], dat['G_c'], lims, 25)
    im = ax2.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=1./aspect, cmap=plt.get_cmap('Reds'))
    cbar = fig.colorbar(im)
    ax2.plot(x1, y1,'k',lw=2, label='Merger Line (Lotz 08)')
    ax2.plot(x2, y2,'k--', lw=2, label='Morph Line (Lotz 08)')
    ax2.set_title('Circular Apertures')
    ax2.set_xlabel('M20', fontsize=14)
    ax2.set_xlim(1., -3.5)
    ax2.set_ylim(.2, .8)
    
    #plt.tight_layout()
    
    plt.savefig('G-M20_'+fout+'.png')

    plt.show()
    plt.close()

def Rp_compare(dat, fout):
    fig = plt.figure(figsize=(8,5))

    ax1 = fig.add_subplot(111)
    lims = [0., 80., 0., 75.]
    grid, extent, aspect = density(dat['Rp_corr']*.396, dat['petrorad_i'], lims, 25)
    im = ax1.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('Reds'))

    #divider = make_axes_locatable(ax1)
    #cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(im, fraction=0.035, pad=0.04)
    cbar.set_label("Galaxy number density", fontsize=18)

    ax1.plot([0,100], [0,100], 'k--',lw=1.5)
    ax1.set_xlabel(r'R$_p$ [arcsec]')
    ax1.set_ylabel(r'SDSS R$_p$ [arcsec]')
    ax1.set_xlim(0,80.)
    ax1.set_ylim(0,75.)

    """
    lowrp = dat[dat['Rp_corr']*0.396<=20]
    coeff = np.polyfit(lowrp['Rp_corr']*.396, lowrp['petrorad_i'], deg=1)
    x = np.linspace(0, 20, 100)
    y = coeff[0]*x + coeff[1]
    print y
    pdb.set_trace()
    """

    """
    ax2 = fig.add_subplot(212)
    lims = [0., 100., 0., 75.]
    grid, extent, aspect = density(dat['Rp_old']*.396, dat['petrorad_i'], lims, 25)
    im = ax2.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('Reds'))
    cbar = fig.colorbar(im)
    ax2.plot([0,100], [0,100], 'k--',lw=1.5, label='1-to-1')
    #ax2.plot(x, y, ls='-.', color='r', label='linear fit')
    #ax2.text(67, 10, r'$y={0:.2f}x+{1:.2f}$'.format(coeff[0], coeff[1]), color='r')

    ax2.set_xlabel(r'My R$_p$ (elliptical: old) [arcsec]')
    ax2.set_ylabel(r'SDSS R$_p$ (circular) [arcsec]')
    ax2.set_xlim(0,100.)
    ax2.set_ylim(0,75.)
    ax2.legend(loc='upper left', frameon=False)    
	"""

    plt.savefig('PetroRad_compare_'+fout+'.png', bbox_inches='tight')
    plt.show()
    plt.close()

def elliptical_circular_compare_SDSS(dat, fout):

    fig = plt.figure(figsize=(10,10))

    ax1 = fig.add_subplot(221)
    lims = [0., 6., 0., 6.]
    grid, extent, aspect = density(dat['C'], dat['C_c'], lims, 25)
    im = ax1.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('Oranges'))
    cbar = fig.colorbar(im)
    #ax1.plot(dat['C'], dat['C_c'], 'y.')
    ax1.plot([0., 6.], [0., 6.], 'k--', lw=2)
    ax1.set_xlim(0., 6.)
    ax1.set_ylim(0., 6.)
    #ax1.set_xlabel('Elliptical')
    ax1.set_ylabel('Circular')
    ax1.set_title('Concentration')

    ax2 = fig.add_subplot(222)
    #ax2.plot(dat['A'], dat['A_c'], 'b.')
    lims = [-0.5, 1.5, -0.5, 1.5]
    grid, extent, aspect = density(dat['A'], dat['A_c'], lims, 25)
    im = ax2.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('BuPu'))
    cbar = fig.colorbar(im)
    ax2.plot([-0.5, 1.5], [-0.5, 1.5], 'k--', lw=2)
    ax2.set_xlim(-0.5, 1.5)
    ax2.set_ylim(-0.5, 1.5)
    ax2.set_title('Asymmetry')

    ax3 = fig.add_subplot(223)
    #ax3.plot(dat['G'], dat['G_c'], 'r.')
    lims = [0., 1.0, 0., 1.0]
    grid, extent, aspect = density(dat['G'], dat['G_c'], lims, 25)
    im = ax3.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('Reds'))
    cbar = fig.colorbar(im)
    ax3.plot([0., 1.], [0., 1.], 'k--', lw=2)
    ax3.set_xlabel('Elliptical')
    ax3.set_ylabel('Circular')
    ax3.set_title('Gini')

    ax4 = fig.add_subplot(224)
    #ax4.plot(dat['M20'], dat['M20_c'], 'g.')
    lims = [-3.5, 0., -3.5, 0.]
    grid, extent, aspect = density(dat['M20'], dat['M20_c'], lims, 25)
    im = ax4.imshow(grid.T, origin='low', extent=lims, 
                    interpolation='nearest', norm=LogNorm(vmin=1.,vmax=15000.),
                    aspect=aspect, cmap=plt.get_cmap('Greens'))
    cbar = fig.colorbar(im)
    ax4.plot([-3.5, 0.], [-3.5, 0.], 'k--', lw=2)
    ax4.set_xlabel('Elliptical')
    ax4.set_title('M20')

    plt.tight_layout()
    plt.savefig('ell_circ_'+fout+'.png')
    plt.show()
    plt.close()

def distributions(dat, fout):
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax1.hist(dat['REDSHIFT'], bins=50, color='red', alpha=.5, 
             normed=True, range=(0,.5))
    ax1.set_xlabel('z')
    
    ax2 = fig.add_subplot(223)
    ax2.hist(dat['MEDIAN'], bins=50, color='red', alpha=.5, 
             normed=True, range=(5, 15))
    ax2.set_xlabel('log[Stellar Mass]')
    
    ax3 = fig.add_subplot(222)
    ax3.hist(dat['petrorad_i'], bins=50, color='red', alpha=.5, 
             normed=True,range=(0,75))
    ax3.set_xlabel('SDSS i band Petro Rad [arcseconds]')
    
    ax4 = fig.add_subplot(224)
    ax4.hist(dat['Rp_c']*.396, bins=50, color='red', alpha=.5, 
             normed=True, range=(0,75))
    ax4.set_xlabel('My Petro Rad [arcseconds]')
    
    plt.tight_layout()
    plt.savefig('z_Rp_mass_dist_'+fout+'.png')
    plt.show()
    plt.close()

def z_mag_rp(dat, fout):
    fig = plt.figure(figsize=(10,8))
    ax1 = fig.add_subplot(211)
    ax1.plot(dat['REDSHIFT'], dat['PETROMAG_I'], 'k.')
    ax1.set_ylabel('i-band Mag')

    ax2 = fig.add_subplot(212)
    ax2.plot(dat['REDSHIFT'], dat['Rp']*.396, 'k.')
    ax2.set_ylabel('Petrosian Rad [arcsec]')
    ax2.set_xlabel('z')

    plt.savefig('z_mag_rp_'+fout+'.png')
    plt.show()
    plt.close()


def remove_flags(dat):
    colnames = dat.colnames
    flags = [name for name in colnames if 'flag' in name and name[0]!='t']

    newdat = dat.copy()
    for f in flags:
        newdat = newdat[np.where(newdat[f] == 0)]

    return newdat
    
def main():

    #dat = Table.read("GZ2All_ancillary_morphology_urls.fits")
    #dat = dat[dat.colnames[:100]]
    #dat['objid'] = dat['dr7objid']

    dat = Table.read("GZmain_petrorads_highzgal.csv")
    dat['objid'] = dat['objID']

    #morphold = Table.read("../SDSSmorphology_catalogs/92415/SDSSmorphology_full_catalog_92415.fits")
    #morphold = morphold[['objid', 'Rp', 'Rp_c']]
    #morphold['Rp_old'] = morphold['Rp']

    morphnew = Table.read("../SDSSmorphology_full_catalog_110817.csv")
    #morphnew['Rp_new'] = morphnew['Rp']


    clean_mask = ((morphnew['bflag']==0) & (morphnew['oflag']==0) & \
                 (morphnew['uflag']==0) & (morphnew['Rpflag']==0))    
    morphclean = morphnew[clean_mask]
    #newdat = Table.read('../measure_morph_bin4_goodapertures_newRpdef.csv')

    fulldat = join(morphclean, dat, keys='objid')

    idx = np.random.choice(range(len(fulldat)), size=5000, replace=False)
    sample = fulldat[idx]

    fout = 'random_sample_fsmooth_colored_clean'
    #compare_parameters_density(dat, fout='fullGZ2_ell_noflags_density')

    #compare_parameters(sample, fout=fout)

    #elliptical_circular_compare_SDSS(dat, fout)
    
    Rp_compare(fulldat, fout='cleanSample_forThesis')
    
    #distributions(dat, fout)
    #gini_m20_compare(dat, fout)
    #z_mag_rp(dat, fout)

    pdb.set_trace()


if __name__ == '__main__':
    main()
