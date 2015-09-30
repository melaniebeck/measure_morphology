'''
Run a slew of various plots to assess the quality of the resulting SDSS 
morphology catalogs I create with galaxy.py
'''

import os
import string
import argparse
import numpy as np
from astropy.table import Table
from scipy.stats import nanmedian 
from collections import OrderedDict, defaultdict
import matplotlib.pyplot as plt
from math import pi
import pdb 

def lotz_mergers(m20):
    return -0.14*m20+0.33

def lotz_separation(m20):
    return  0.14*m20+0.8

def compare_parameters(dat, fout='mycat'):
    
    colors = 'black'

    fig = plt.figure(figsize=(10,8))
    gs = plt.GridSpec(4,4)
    gs.update(wspace=0.005, hspace=0.005)

    #'''
    ax10 = fig.add_subplot(gs[3,0])
    ax10.scatter(dat['A'], dat['C'], color=colors, marker='.', alpha=0.25)
    ax10.set_xlabel('A', fontsize=20)
    ax10.set_ylabel('C', fontsize=20)
    ax10.set_ylim(0., 6.)
    ax10.set_xlim(0., .8)
    ax10.set_xticks([0.2, 0.4, 0.6])
    ax10.set_yticks([1,2,3,4])
    #ax10.tick_params(axis='x', pad=8)
    #plt.text(0.5, 0.5, '10', fontsize=20, color='red', transform=ax10.transAxes)
    '''
    ax10.text(0,19, 'Ellipticals', fontsize=16, color='red')
    ax10.text(0,18, 'Disks', fontsize=16, color='purple')
    ax10.text(0,17, 'Edge-on Disks', fontsize=16, color='green')
    ax10.text(0,16, 'Mergers', fontsize=16, color='yellow')
    '''
    
    ################
    
    ax8 = fig.add_subplot(gs[3,1], sharey=ax10)
    plt.setp(ax8.get_yticklabels(), visible=False)
    ax8.scatter(dat['G'], dat['C'], color=colors, marker='.', alpha=0.25)
    ax8.set_xlabel('G', fontsize=20)
    #ax8.set_ylabel('C')
    ax8.set_xlim(0.2, 0.8)
    ax8.set_ylim(0., 6.)
    ax8.set_xticks([0.4, 0.5, 0.6])
    #plt.text(0.5, 0.5, '8', fontsize=20, color='red', transform=ax8.transAxes)
    
    ax9 = fig.add_subplot(gs[2,1], sharex=ax8)
    plt.setp(ax9.get_xticklabels(), visible=False)
    ax9.scatter(dat['G'], dat['A'], color=colors, marker='.', alpha=0.25)
    #ax9.set_xlabel('G')
    ax9.set_ylabel('A', fontsize=20)
    ax9.set_xlim(0.2, 0.8)
    ax9.set_ylim(0., .8)
    ax9.set_yticks([0.2, 0.4, 0.6])
    #plt.text(0.5, 0.5, '9', fontsize=20, color='red', transform=ax9.transAxes)
    
    ################
    
    ax5 = fig.add_subplot(gs[3,2], sharey=ax10)
    plt.setp(ax5.get_yticklabels(), visible=False)
    ax5.scatter(dat['M20'], dat['C'], color=colors, marker='.', alpha=0.25)
    ax5.set_xlabel('M20', fontsize=20)
    #ax5.set_ylabel('C')
    ax5.set_xlim(0., -3.5)
    ax5.set_ylim(0., 6.)
    ax5.set_xticks([-2.5, -1.5, -0.5])
    #plt.text(0.5, 0.5, '5', fontsize=20, color='red', transform=ax5.transAxes)
    
    ax6 = fig.add_subplot(gs[2,2], sharex=ax5, sharey=ax9)
    plt.setp((ax6.get_xticklabels(), ax6.get_yticklabels()), visible=False)
    ax6.scatter(dat['M20'], dat['A'], color=colors, marker='.', alpha=0.25)
    #ax6.set_xlabel('M20')
    #ax6.set_ylabel('A')
    ax6.set_xlim(0., -3.5)
    ax6.set_ylim(0., .8)
    #plt.text(0.5, 0.5, '6', fontsize=20, color='red', transform=ax6.transAxes)
    
    ax7 = fig.add_subplot(gs[1,2], sharex=ax5)
    plt.setp(ax7.get_xticklabels(), visible=False)
    ax7.scatter(dat['M20'], dat['G'], color=colors, marker='.', alpha=0.25)
    #ax7.set_xlabel('M20')
    ax7.set_ylabel('G', fontsize=20)
    ax7.set_xlim(0., -3.5)
    ax7.set_ylim(0.2, 0.8)
    ax7.set_yticks([0.4, 0.5, 0.6])
    #plt.text(0.5, 0.5, '7', fontsize=20, color='red', transform=ax7.transAxes)
    
    ################
    
    ax1 = fig.add_subplot(gs[3,3], sharey=ax10)
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax1.scatter(dat['elipt'], dat['C'], color=colors, marker='.', alpha=0.25)
    ax1.set_xlabel('elipt', fontsize=20)
    #ax1.set_ylabel('C')
    ax1.set_ylim(0., 6.)
    ax1.set_xlim(0., .9)
    ax1.set_xticks([0.2, 0.4, 0.6, 0.8])
    #plt.text(0.5, 0.5, '1', fontsize=20, color='red', transform=ax1.transAxes)
    
    ax2 = fig.add_subplot(gs[2,3], sharex=ax1, sharey=ax9)
    plt.setp((ax2.get_xticklabels(), ax2.get_yticklabels()), visible=False)
    ax2.scatter(dat['elipt'], dat['A'], color=colors, marker='.', alpha=0.25)
    #ax2.set_xlabel('elipt')
    #ax2.set_ylabel('A')
    ax2.set_ylim(0., .8)
    ax2.set_xlim(0., .9)
    #plt.text(0.5, 0.5, '2', fontsize=20, color='red', transform=ax2.transAxes)
    
    ax3 = fig.add_subplot(gs[1,3], sharex=ax1, sharey=ax7)
    plt.setp((ax3.get_xticklabels(),ax3.get_yticklabels()), visible=False)
    ax3.scatter(dat['elipt'], dat['G'], color=colors, marker='.', alpha=0.25)
    #ax3.set_xlabel('elipt')
    #ax3.set_ylabel('G')
    ax3.set_xlim(0., .9)
    ax3.set_ylim(0.2, 0.8)
    #plt.text(0.5, 0.5, '3', fontsize=20, color='red', transform=ax3.transAxes)
    
    ax4 = fig.add_subplot(gs[0,3])
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax4.scatter(dat['elipt'], dat['M20'], color=colors, marker='.', alpha=0.25)
    #ax4.set_xlabel('elipt')
    ax4.set_ylabel('M20', fontsize=20)
    ax4.set_ylim(0., -3.5)
    ax4.set_xlim(0., .9)
    ax4.set_yticks([-2.5, -1.5, -0.5])
    #plt.text(0.5, 0.5, '4', fontsize=20, color='red', transform=ax4.transAxes)

    #gs.tight_layout(fig)
    plt.savefig('compare_parameters_'+fout+'.png')
    plt.show()
    plt.close()


def bin_contour():

    #'''
    x, y = dat['A'], dat['C']
    notnans = np.where((~np.isnan(x) & ~np.isnan(y)))
    x, y = x[notnans], y[notnans]

    gridx = np.linspace(0,.8, 10)
    gridy = np.linspace(2,6,20)

    grid, xedges, yedges = np.histogram2d(x, y, bins=[gridx, gridy]) 
    #normed=True)

    plt.figure()
    myextent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.imshow(grid.T, origin='low', extent=myextent, interpolation='nearest',
               aspect=.3)
    #plt.pcolormesh(gridx,gridy,grid)
    #plt.plot(x,y,'bo', alpha=.2)
    plt.colorbar()
    plt.show()

    pdb.set_trace()

   # '''
def gini_m20_compare(cat1, cat2, titles, outfile):

    fig = plt.figure(figsize=(20,8))

    x1 = np.arange(-3.5, 1., .1)
    y1 = lotz_mergers(x1)
    x2 = np.arange(-3.5, -1.6, .1)
    y2 = lotz_separation(x2)

    ax1 = fig.add_subplot(131)

    ax1.plot(cat1['M20_2'], cat1['G_2'], 'r.', alpha=.5) 
    ax1.plot(x1, y1,'k',lw=2, label='Merger Line (Lotz 08)')
    ax1.plot(x2, y2,'k--',lw=2, label='Morph Line (Lotz 08)')
    #ax.plot([0,-3.], [.4, .4], 'k--')   

    ax1.set_ylabel('Gini', fontsize=14)
    ax1.set_xlabel('M20', fontsize=14)
    ax1.set_title('Ellipt Apers, M20/G bad [Original]')
    ax1.set_xlim(1., -3.5)
    ax1.set_ylim(.2, .8)

    ax2 = fig.add_subplot(132)

    ax2.plot(cat1['M20'], cat1['G'], 'r.',  alpha=.5) 
    ax2.plot(x1, y1,'k',lw=2, label='Merger Line (Lotz 08)')
    ax2.plot(x2, y2,'k--', lw=2, label='Morph Line (Lotz 08)')
    #ax2.plot([0,-3.], [.4, .4], 'k--')   

    #ax2.set_ylabel('Gini', fontsize=14)
    ax2.set_xlabel('M20', fontsize=14)
    ax2.set_title('Ellipt Apers, M20/G fixed')
    ax2.set_xlim(1., -3.5)
    ax2.set_ylim(.2, .8)
    
    ax3 = fig.add_subplot(133)
    ax3.plot(cat2['M20'], cat2['G'], 'r.',  alpha=.5) 
    ax3.plot(x1, y1,'k',lw=2, label='Merger Line (Lotz 08)')
    ax3.plot(x2, y2,'k--', lw=2, label='Morph Line (Lotz 08)')
    #ax3.plot([0,-3.], [.4, .4], 'k--')   

    #ax3.set_ylabel('Gini', fontsize=14)
    ax3.set_xlabel('M20', fontsize=14)
    ax3.set_title('Circular Apers, M20/G fixed')
    ax3.set_xlim(1., -3.5)
    ax3.set_ylim(.2, .8)

    plt.tight_layout()
    
    plt.savefig('G-M20_SDSS.png')

    plt.show()
    plt.close()

def Rp_compare(cat1, cat2, xylabels, histlabels, outname):
    try:
        cat1['Rp'] = cat1['Rp_1']
        cat2['Rp'] = cat2['Rp_1']
    except:
       pass

    fig = plt.figure(figsize=(8,8))

    ax1 = fig.add_subplot(111)
    ax1.plot(cat1['Rp_2']*.396, cat2['petrorad_i'], 'r.', 
             label='Ell Aper; Original cleaning')
    ax1.plot(cat1['Rp']*.396, cat2['petrorad_i'], 'b.', 
             label='Ell Aper; better cleaning')
    ax1.plot(cat2['Rp']*.396, cat2['petrorad_i'], 'y.', 
             label='Circ Aper; better cleaning')
    ax1.plot([0,100], [0,100], 'k--',lw=1.5, 
             label='1-to-1')
    ax1.set_xlabel('My Petrosian Radii [arcsec]')
    ax1.set_ylabel('SDSS Petrosian Radius [arcsec]')
    #ax1.set_title('Petrosian Radius [pixels]'
    ax1.set_ylim(0,55)
    ax1.legend(loc='best')
    
    '''
    ax2 = fig.add_subplot(122)
    ax2.hist(cat1['Rp'], bins=100, normed=True, color='blue', alpha=.5,
             label=histlabels[0])
    ax2.hist(cat2['Rp'], bins=100, normed=True, color='red', alpha=.5, 
             label=histlabels[1])
    ax2.set_xlabel('Rp [pixels]')
    ax2.legend(loc='best')
    #'''

    plt.savefig('Rp_SDSS.png')
    plt.show()
    plt.close()


def elliptical_circular_compare_SDSS():
    cat = Table.read('test_galchanges.fits')

    fig = plt.figure(figsize=(16,16))
    ax1 = fig.add_subplot(221)

    #'''
    ax1.plot(cat['C'], cat['C_c'], 'ys')
    ax1.plot([2., 5.], [2., 5.], 'k--', lw=1.5)
    #ax1.set_xlabel('Elliptical')
    ax1.set_ylabel('Circular')
    ax1.set_title('Concentration')
    '''
    ax1.plot(cat1['Mx'], cat2['Mx'], 'ys')
    ax1.plot([2., 5.], [2., 5.], 'k--', lw=1.5)
    #ax1.set_xlabel('Elliptical')
    ax1.set_ylabel('Circular')
    ax1.set_title('Concentration')
    '''
    ax2 = fig.add_subplot(222)
    ax2.plot(cat['A'], cat['A_c'], 'bs')
    ax2.plot([-.1, .5], [-.1, .5], 'k--', lw=1.5)
    #ax2.set_xlabel('Elliptical')
    #ax2.set_ylabel('Circular')
    ax2.set_title('Asymmetry')

    ax3 = fig.add_subplot(223)
    ax3.plot(cat['G2'], cat['G2_c'], 'rs')
    ax3.plot([.25, .75], [.25, .75], 'k--', lw=1.5)
    ax3.set_xlabel('Elliptical')
    ax3.set_ylabel('Circular')
    ax3.set_title('Gini')

    ax4 = fig.add_subplot(224)
    ax4.plot(cat['M20'], cat['M20_c'], 'gs')
    ax4.plot([-3., -1.], [-3., -1.], 'k--', lw=1.5)
    ax4.set_xlabel('Elliptical')
    #ax4.set_ylabel('Circular')
    ax4.set_title('M20')

    plt.tight_layout()
    plt.savefig('morphparams_ell_circ_SDSS2.png')
    plt.show()
