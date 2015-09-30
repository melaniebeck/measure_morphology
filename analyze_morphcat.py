'''
Run a slew of various plots to assess the quality of the resulting SDSS 
morphology catalogs I create with galaxy.py
'''

import numpy as np
from astropy.table import Table, vstack
from collections import OrderedDict, defaultdict
import matplotlib.pyplot as plt
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
    plt.savefig('morph_params_'+fout+'.pdf')
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

def gini_m20_compare(dat, fout):

    fig = plt.figure(figsize=(20,8))

    x1 = np.arange(-3.5, 1., .1)
    y1 = lotz_mergers(x1)
    x2 = np.arange(-3.5, -1.6, .1)
    y2 = lotz_separation(x2)

    ax1 = fig.add_subplot(121)

    ax1.plot(dat['M20'], dat['G'], 'r.', alpha=.5) 
    ax1.plot(x1, y1,'k',lw=2, label='Merger Line (Lotz 08)')
    ax1.plot(x2, y2,'k--',lw=2, label='Morph Line (Lotz 08)')
    
    ax1.set_title('Elliptical Apertures')
    ax1.set_ylabel('Gini', fontsize=14)
    ax1.set_xlabel('M20', fontsize=14)
    ax1.set_xlim(1., -3.5)
    ax1.set_ylim(.2, .8)

    ax2 = fig.add_subplot(122)

    ax2.plot(dat['M20_c'], dat['G_c'], 'r.',  alpha=.5) 
    ax2.plot(x1, y1,'k',lw=2, label='Merger Line (Lotz 08)')
    ax2.plot(x2, y2,'k--', lw=2, label='Morph Line (Lotz 08)')
    
    ax2.set_title('Circular Apertures')
    ax2.set_xlabel('M20', fontsize=14)
    ax2.set_xlim(1., -3.5)
    ax2.set_ylim(.2, .8)
    
    plt.tight_layout()
    
    plt.savefig('G-M20_'+fout+'.pdf')

    plt.show()
    plt.close()

def Rp_compare(dat, fout):
    fig = plt.figure(figsize=(8,8))

    ax1 = fig.add_subplot(111)
    ax1.plot(dat['Rp']*.396, dat['petrorad_i'], 'r.', 
             label='Elliptical')

    ax1.plot(dat['Rp_c']*.396, dat['petrorad_i'], 'b.', 
             label='Circular')

    ax1.plot([0,100], [0,100], 'k--',lw=1.5, 
             label='1-to-1')

    ax1.set_xlabel('My Petrosian Radii [arcsec]')
    ax1.set_ylabel('SDSS Petrosian Radius [arcsec]')

    ax1.set_ylim(0,55)
    ax1.legend(loc='best')
    
    plt.savefig('PetrosianRad_compare_'+fout+'.pdf')
    plt.show()
    plt.close()

def elliptical_circular_compare_SDSS(dat, fout):

    fig = plt.figure(figsize=(16,16))

    ax1 = fig.add_subplot(221)
    ax1.plot(dat['C'], dat['C_c'], 'y.')
    ax1.plot([1., 7.], [1., 7.], 'k--', lw=2)
    #ax1.set_xlabel('Elliptical')
    ax1.set_ylabel('Circular')
    ax1.set_title('Concentration')

    ax2 = fig.add_subplot(222)
    ax2.plot(dat['A'], dat['A_c'], 'b.')
    ax2.plot([-.1, .7], [-.1, .7], 'k--', lw=2)
    #ax2.set_xlabel('Elliptical')
    #ax2.set_ylabel('Circular')
    ax2.set_title('Asymmetry')

    ax3 = fig.add_subplot(223)
    ax3.plot(dat['G'], dat['G_c'], 'r.')
    ax3.plot([.2, .8], [.2, .8], 'k--', lw=2)
    ax3.set_xlabel('Elliptical')
    ax3.set_ylabel('Circular')
    ax3.set_title('Gini')

    ax4 = fig.add_subplot(224)
    ax4.plot(dat['M20'], dat['M20_c'], 'g.')
    ax4.plot([-3.5, 0.], [-3.5, 0.], 'k--', lw=2)
    ax4.set_xlabel('Elliptical')
    ax4.set_title('M20')

    plt.tight_layout()
    plt.savefig('ell_circ_'+fout+'.pdf')
    plt.show()

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
    plt.savefig('z_Rp_mass_dist_'+fout+'.pdf')
    plt.show()


def main():

    #data1 = Table.read('GZ2Photoz_ancillary_morphology_masses.fits')
    data2 = Table.read('GZ2Specz_ancillary_morphology_masses.fits')

    #dat = vstack((data1, data2))


    fout = 'fullGZ2'
    #compare_parameters(data2, fout='fullGZ2_ell')
    #elliptical_circular_compare_SDSS(data2, fout)
    #Rp_compare(data2, fout)
    #distributions(data2, fout)

    gini_m20_compare(data2, fout)

if __name__ == '__main__':
    main()
