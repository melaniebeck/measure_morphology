# A random script for doing random things that need doing
# Happy Birthday Jaybirdy. I miss you more than I can bear. bare? 

import os
import string
import argparse
import pyfits as fits
import numpy as np
from astropy.table import Table
from scipy.spatial import cKDTree
from scipy.interpolate import interp1d
from scipy.stats import nanmedian 
from collections import OrderedDict, defaultdict
from random import gauss
import matplotlib.pyplot as plt
from math import pi
import pdb #"""for doing an IDL-like stop"""
from utils import find_closest
import run_sextractor


def comparecatalogs():
    #'''
    # read in sample of galaxies -- I need to compare some parameters
    #testgals = Table.read('bright_gals.fits')#,format='ascii.ipac', guess=False)
    cassata = Table.read('Cassata_catalog.tbl',format='ascii.ipac')
    zest_sub = Table.read('bright_gals_zest.fits')
    
    testcoords = zip(zest_sub['ra'], zest_sub['dec'])
    
    coords = defaultdict(list)         
    for testcoord in testcoords:    
        (ra, dec), index = find_closest(testcoord, zip(cassata['ra'], cassata['dec']))
        coords['ra'].append(ra)
        coords['dec'].append(dec)
        coords['idx'].append(index)
    #'''
    
    casssub = cassata[coords['idx']]
    casssub.write('bright_gals_cassata.fits', clobber=True)
    

def plotzestcassatamine():
    cassdat = Table.read('bright_gals_cassata.fits')
    zestdat = Table.read('bright_gals_zest.fits')
    mydat = Table.read('bright_gals_newpetro.fits')
    
    zestdat['rpet'][np.where(zestdat['rpet'] > 1000.)] = 0.

    #xaxis = np.linspace(0., 250., len(mydat['PETRO_RAD']))

    plt.plot(cassdat['r_petro'], zestdat['rpet'], 'go', label='Cassata vs. ZEST')
    plt.plot(mydat['PETRO_RAD'], zestdat['rpet'], 'bo', label='Mine vs. ZEST')
    plt.plot([0,250], [0,250], 'k--')
    plt.legend()
    plt.xlabel('ZEST')
    plt.ylabel('Mine/Cassata')
    plt.title('Petrosian Radius Comparion - 3 catalogs')
    plt.savefig('rpet_comparison_3cats.png')
    plt.show()
    
def selectnewsample():
    #bright = fits.open('bright_gals_zest.fits')
    zest = Table.read('ZEST_catalog.fits')
    racut = np.where((zest['ra'] >= 149.95) & (zest['ra'] <= 150.065))
    cut = np.where( (zest['dec'][racut] >= 2.6) & (zest['dec'][racut] <= 2.7))
    
    sample = zest[cut]
    sample.write('bigsample_zest.fits', overwrite=True)
    pdb.set_trace()
    #indexes = []
    #for coord in zip(bright['ra'], bright['dec']):
    #    (blah1, blah2), index = find_closest(coord, zip(zest['ra'][cut],zest['dec'][cut]))
    #    indexes.append(index)
    
    #print len(cut[0])
    #print max(zest['acs_mag_auto'][cut]), min(zest['acs_mag_auto'][cut])
    #plt.hist(zest['acs_mu_class'][cut])
    #plt.show()
    if len(cut[0]) < 1000.:
        coords = Table([zest['ra'][cut], zest['dec'][cut]], names=('ra', 'dec'))
        coords.write('bigsample.dat', format='ascii.csv')
        #pdb.set_trace()
    else: 
        print 'make it smaller'

    
def plotpetrorad_largesample():
    dat = Table.read('bigsample_mycat_v10.fits')
    zest = Table.read('bigsample_zest.fits')

    zrad = zest['rpet']
    mrad = dat['rpet']

    mu = np.nanmean(mrad/zrad)
    sig = np.nanstd(mrad/zrad)
    
    # Histogram of magnitude distribution
    plt.figure()
    plt.hist(zest['acs_mag_auto'], bins=40, normed=1)
    ax = plt.gca()
    ax.invert_xaxis()
    plt.title('Magnitude distribution for 971 galaxies')
    plt.xlabel('Magnitude F814')
    #plt.savefig('magdist_bigsample.png')
    #plt.show()
    
    #exit()
    #pdb.set_trace()
    
    # Histogram of MRAD/ZRAD ratio with mean/sigma
    plt.figure()
    plt.hist(mrad/zrad, bins=100, range=[0,3], color='blue')
    plt.plot([mu,mu], [0,130], lw=4)
    plt.axvline(x=mu, linewidth=3, color='k')
    plt.axvline(x=mu+sig, linewidth=3, color='k', linestyle='--')
    plt.axvline(x=mu-sig, linewidth=3, color='k', linestyle='--')
    #plt.plot([mu+sig, mu+sig], [0,130], 'k--', linewidth=4)
    #plt.plot([mu-sig, mu-sig], [0,130], 'k--', linewidth=4)
    #plt.hist(zrad, bins=20, range=[0,200], color='blue', alpha=0.5, label='Claudias')
    plt.legend(loc='upper right')
    plt.xlabel('Radius [pixels]')
    
    #plt.savefig('petro_ratio_bigsample_v3.png')
    #plt.show()
    
    #exit()

    # My radius vs. ZEST radius and ratio
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(zrad, mrad, 'ro', [0,200], [0,200], 'k--')
    plt.title('Petrosian Radius')
    plt.xlabel('ZEST [pixels]')
    plt.ylabel('Mine [pixels]')
    
    plt.subplot(2,1,2)
    plt.plot(zrad, mrad/zrad, 'ro')
    plt.plot([0,200], [1,1], 'k--')
    plt.yscale('log')
    plt.xlabel('Radius [pixels]')
    plt.ylabel('Mine/ZEST')
    
    '''
    subset = np.where(mrad/zrad < .5)
    i = 0
    for xy in zip(zrad[subset],mrad[subset]/zrad[subset]):
        plt.annotate(str(subset[0][i]), xy=xy)
        i += 1
    #'''
    plt.savefig('rpet_bigsample.png')
    
    plt.show()
    plt.close()
    
    #exit()
    '''
    # create a table with pertinent infos for the subset
    mr = mrad[subset[0]]
    ra = dat['ra'][subset[0]]
    dec = dat['dec'][subset[0]]
    
    zr = zrad[subset[0]]
    flag = zest['flagrpet'][subset[0]]
    bigtable = Table([subset[0], mr, zr, ra, dec, flag], 
                names=('imgnum', 'myrad', 'zestrad', 'ra', 'dec', 'flag'))
    bigtable.write('toolow_bigsample_v3.txt', format='ascii.fixed_width')
    
    exit()
    '''

def plotasymmetry(cat1, cat2, outname):
    a1 = cat1['aa']
    a2 = cat2['A']    

    plt.plot(a1, a2, 'ro', [-.2,.5], [-.2,.5], 'k')
    plt.axvline(x=0., linestyle='--', color='k')
    plt.axhline(y=0., linestyle='--', color='k')
    plt.ylim(-0.3, .8)
    plt.xlim(-.2, .5)
    plt.title('Asymmetry')
    plt.xlabel('ZEST')
    plt.ylabel('Mine')
    plt.savefig(outname)
    plt.show()

    #pdb.set_trace()


def plotconcentration():
    dat1 = Table.read('bigsample_mycat_v10.fits') # circulur
    dat2 = Table.read('bigsample_mycat_v9.fits') # elliptical
    zest = Table.read('bigsample_zest.fits')
    
    myc1 = dat1['C'] # circular
    myc2 = dat2['C'] # elliptical
    zc = zest['cc']
    '''
    mavg = np.nanmean(myc)
    mstd = np.nanstd(myc)
    zavg = np.nanmean(zc)
    zstd = np.nanstd(zc)
    
    print "my avg:", mavg, "my std:", mstd
    print "ZEST avg:", zavg, "ZEST std:", zstd
    '''
    ratio1 = myc1/zc
    ratio2 = myc2/zc
    #ravg = np.nanmean(ratio)
    #rstd = np.nanstd(ratio)
    '''
    plt.figure()
    #plt.hist(myc, bins=20, range=(1,5), alpha=0.5)
    #plt.hist(zc, bins=20, range=(1,5), color='r', alpha=0.5)
    plt.hist(myc/zc, bins=25, range=(0.5,1.5))
    plt.axvline(x=ravg, color='k', linewidth=3)
    plt.axvline(x=ravg+rstd, color='k', linestyle='--')
    plt.axvline(x=ravg-rstd, color='k', linestyle='--')
    plt.title('My C / ZEST C')
    plt.savefig('Cratio_hist_circ.png')
    plt.close()
    '''
    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(zc, myc1, 'ro', [1,5], [1,5], 'k--')
    plt.ylim(1., 6.)
    plt.xlabel('ZEST C')
    plt.ylabel('My C')
    plt.title('Concentration using Circular Annuli')

    plt.subplot(2,2,2)
    plt.plot(zc, myc2, 'ro', [1,5], [1,5], 'k--')
    plt.ylim(1.0, 6.)
    plt.xlabel('ZEST C')
    plt.ylabel('My C')
    plt.title('Concentration using Elliptical Annuli')

    plt.subplot(2,2,3)
    plt.plot(zc, ratio1, 'bo')
    plt.axhline(y=1.0, linestyle='--', color='k')
    plt.ylim(.4, 2.)
    plt.xlabel('ZEST')
    plt.ylabel('My C / ZEST C')

    plt.subplot(2,2,4)
    plt.plot(zc, ratio2, 'bo')
    plt.axhline(y=1.0, linestyle='--', color='k')
    plt.ylim(.4, 2.)
    plt.xlabel('ZEST C')
    plt.ylabel('My C / ZEST C')    
    '''
    toohigh = np.where( ratio > ravg+2*rstd )
    i = 0
    for xy in zip(zc[toohigh], ratio[toohigh]):
        #print xy
        plt.annotate( '%i' %(toohigh[0][i]+1), xy=xy)
        i += 1
    toolow = np.where( ratio < ravg-rstd )
    i = 0
    for xy in zip(zc[toolow], ratio[toolow]):
        #print xy
        plt.annotate( '%i' %(toolow[0][i]+1), xy=xy)
        i += 1
    #'''
    plt.savefig('C_compare_apshape_v2.png')
    plt.show()
    plt.close()

    exit()

    mdev = dat[toohigh]
    zdev = zest[toohigh]
    bigtable = Table([mdev['name'], mdev['C'], zdev['cc'], mdev['rpet'], \
                          zdev['rpet'], zdev['flagrpet']],\
                         names=('name', 'myC', 'zestC', 'myrpet', 'zestrpet', 'flag'))
    bigtable.write('concentration_toohigh.txt', format='ascii.fixed_width')

    f = open('concentration_toohigh.sh', 'w+')
    f.write('ds9 ')
    for name in mdev['name']: 
        f.write('output/%s.fits[1] '%name)
    f.write('-single -lock frame image -scale zscale -zoom to fit')
    f.write('-single -lock scale yes')
    f.close() 

    mdev = dat[toolow]
    zdev = zest[toolow]
    bigtable = Table([mdev['name'], mdev['C'], zdev['cc'], mdev['rpet'], \
                          zdev['rpet'], zdev['flagrpet']],\
                         names=('name', 'myC', 'zestC', 'myrpet', 'zestrpet', 'flag'))
    bigtable.write('concentration_toolow.txt', format='ascii.fixed_width')

    f = open('concentration_toolow.sh', 'w+')
    f.write('ds9 ')
    for name in mdev['name']: 
        f.write('output/%s.fits[1] '%name)
    f.write('-single -lock frame image -scale zscale -zoom to fit')
    f.write('-single -lock scale yes')
    f.close()

    #pdb.set_trace()

def testbackground():
    dat = Table.read('bigsample_mycat_v9.fits')
    
    #pdb.set_trace()
    plt.figure()
    plt.hist(dat['med'], bins=20, range=(-0.001, 0.001))
    plt.title('Background median')
    plt.savefig('bkgmed_hist.png')
    plt.figure()
    plt.hist(dat['rms'], bins=20, range=(0.001, 0.006))
    plt.title('Background RMS')
    plt.savefig('bkgrms_hist.png')
    #plt.show()

    rmscut= np.where(dat['rms'] < 0.005)
    ra, dec, rms = dat['ra'][rmscut], dat['dec'][rmscut], dat['rms'][rmscut]
    fig,ax = plt.subplots(1,1, figsize=(10,6))
    im = ax.scatter(ra, dec, s=50, c=rms, alpha=.4, cmap=plt.cm.RdYlBu)
    ax.set_xlabel('RA [degrees]')
    ax.set_ylabel('Dec [degrees]')
    ax.set_ylim(1.55, 2.9)
    ax.set_xlim(149.55,150.8 )
    plt.colorbar(im,ax=ax)
    #fig.show()
    fig.savefig('position_rms.png')
    pdb.set_trace()

def plotgini(sample1, sample2, outname):

    plt.figure()
    plt.plot(sample1, sample2, 'ro', [0,1], [0,1], 'k--', zgwut, wut, 'bo')
    plt.ylim(0., 6.)
    plt.xlim(.2, .8)
    plt.xlabel('Sample 1')
    plt.ylabel('Sample 2 (Mine)')
    plt.savefig(outname)

    '''
    plt.figure()
    plt.hist(mg/zg, bins=60, range=(0., 4.5))
    plt.title('My Gini / ZEST Gini')
    #plt.hist(mg, bins=40, range=(0,1), alpha=0.5)
    #plt.hist(zg, bins=20, range=(0,1), alpha=0.5)
    plt.savefig('histGini_bad.png')
    #'''

def plotM20(zest, mine, outname):
    zm = zest['m20']
    mm = mine['M1']

    plt.figure()
    plt.plot(zm, mm, 'ro', [-3,0], [-3,0])
    plt.title("M20 comparison between ZEST and mine")
    plt.xlabel('ZEST M20')
    plt.ylabel('My M20')
    plt.savefig(outname)

    '''
    plt.figure()
    plt.plot(zm, mm/zm, 'ro')
    plt.axhline(y=1., ls='--')
    plt.savefig(outname)
    plt.show()
    #'''

def analyse_cleaning():
    '''
    Flag Key: (old)
    0:  didn't pass through anything. Oops
    2:  DIST < 10 & Bdist < 10
    3:  DIST > 10 & Bdist < 10 & Barea > 0.75 Farea
    4:  DIST > 10 & Bdist < 10 
    5:  DIST > 10 & Bdist > 10 & Fdist < 10
    6:  DIST > 10 & Bdist > 10 & Fdist > 10
    '''
    flagtype = [1,2,3,4,5,6,7,8]
    dat = Table.read('data7.txt', format='ascii')
    
    ff = dat['Flag']

    plt.figure()
    weights= np.ones_like(ff, dtype='float')/len(ff)
    plt.hist(ff, weights=weights)
    plt.xlabel('Flag value')
    plt.ylabel('Proportion of Flag Value for Sample')
    plt.savefig('flag_freq_new.png')
    plt.close()
    #plt.show()


    datf = {}
    for f in flagtype:
        datf["flag"+str(f)] = dat['name', 'Fdist', 'Bdist', \
                                  'Fidx', 'Farea', 'Barea'][ff == f]

    colors=['red', 'blue', 'green','yellow', 
            'cyan', 'purple','black', 'magenta']
    '''
    plt.figure()
    for idx, f in enumerate(flagtype):
        datname = 'flag'+str(f)
        plt.plot(datf[datname]['Fdist'], datf[datname]['Bdist'], 
                 color=colors[idx],  marker='o', ls='None', label=datname)
    plt.legend()
    plt.xlabel('Fdist')
    plt.ylabel('Bdist')
    plt.legend()
    plt.savefig('FBdist_flags_data4.png')
    plt.close()
    #plt.show()
    
    for idx, f in enumerate(flagtype):
        plt.figure()
        datname = 'flag'+str(f)
        weights = np.ones_like(datf[datname]['Farea'])/len(datf[datname]['Farea'])
        try:
            bins = np.arange(min(datf[datname]['Farea']), 
                             max(datf[datname]['Farea'])+100, 100)
            plt.hist(datf[datname]['Farea'], weights=weights, bins=bins, 
                     alpha=0.5, color=colors[idx], label=datname)
        except ValueError:
            pass
        plt.legend()
        #plt.xlim(0,30000)
        plt.xlabel('Farea')
        plt.ylabel('Proportion of subclass: Flag')
        plt.close()
    #plt.show()

    #'''
    plt.figure()
    plt.hist(datf['flag1']['Fdist'], color='red', alpha=.5)
    plt.hist(datf['flag1']['Bdist'], color='blue', alpha=.5)
    plt.show()
    pdb.set_trace()

    # test SE 
    images = datf['flag1']['name']
    for idx, i in enumerate(images):
        if (idx % 50 == 0):
            f=open('category1_'+str(idx)+'.sh', 'wb')
            f.write('ds9s ')
            f.write('output/datacube7/f_'+i+'[1] ')
        else:
            f.write('output/datacube7/f_'+i+'[1] ')
            if (idx+1) % 50 == 0:
                f.close()
    f.close()

    pdb.set_trace()


    #f = open('testB2dist_med.sh', 'wb')
    
    #for n in med:
    #    f.write('ds9m output/datacube2/f_%s\n'%n['name'])
    #f.close()

    f = open('Fdist_large.txt', 'wb')
    fdist = dat[dat['Fdist']>20.]
    for n in fdist:
        f.write('output/datacube2/f_%s\n'%n['name'])
    f.close()

def compare_catalogs(z, m, fout='mycat'):

    fig = plt.figure(figsize=(10,8))
    gs = plt.GridSpec(2,4)

    cut = np.where(m['Rp'] != c)
    # plot Petrosian Radius
    ax1 = fig.add_subplot(gs[0,0])
    ax1.plot(z['rpet'][cut], m['Rp'][cut],'ro', [0,200],[0,200], 'k--')
    ax1.set_ylim(-5., 200.)
    ax1.set_xlim(-5., 200.)
    ax1.set_title('Petrosian Radius')

    # plot Concentration
    ax2 = fig.add_subplot(gs[0,1])
    ax2.plot(z['cc'][cut], m['C'][cut], 'ro', [1,5], [1,5], 'k--')
    ax2.set_ylim(1., 6.)
    ax2.set_title('Concentration')

    # plot Asymmetry
    ax3 = fig.add_subplot(gs[0,2])
    ax3.plot(z['aa'][cut], m['A'][cut], 'ro', [-.2,.5], [-.2,.5], 'k')
    ax3.vlines(0., -.3, .8, linestyle='--', color='k')
    ax3.hlines(0., -0.2, .5, linestyle='--', color='k')
    ax3.set_ylim(-0.3, .8)
    ax3.set_xlim(-.2, .5)
    ax3.set_title('Asymmetry')

    # plot Gini 
    ax4 = fig.add_subplot(gs[1,0])
    ax4.plot(z['gg'][cut], m['G1'][cut], 'ro',  [0,1], [0,1], 'k--')
    ax4.set_ylim(0., 1.)
    ax4.set_xlim(0., 1.)
    ax4.set_title('Gini')

    ax5 = fig.add_subplot(gs[1,1])
    ax5.plot(z['gg'][cut], m['G2'][cut], 'ro',  [0,1], [0,1], 'k--')
    ax5.set_ylim(0., 1.)
    ax5.set_xlim(0., 1.)
    ax5.set_title('Gini')

    # plot M201
    ax6 = fig.add_subplot(gs[1,2])
    ax6.plot(z['m20'][cut], m['M1'][cut], 'ro', [-3.,0.], [-3.,0.], 'k--')
    ax6.set_title('M20')

    # BONUS PLOT
    ax7 = fig.add_subplot(gs[1,3])
    ax7.plot(z['m20'][cut], m['M2'][cut], 'ro', [-3.,0.], [-3.,0.], 'k--')
    ax7.set_title('M20')

    gs.tight_layout(fig)
    plt.savefig('compare_catalogs_zest&'+fout+'.png')
    plt.show()
    plt.close()

def compare_parameters(dat, fout='mycat'):
    try:
        dat['A']= dat['A_1a']
    except:
        pass

    try:
        colors = dat['t_color']
        colors = [c.rstrip() for c in colors]
    except:
        colors = 'black'

    #pdb.set_trace()

    fig = plt.figure(figsize=(10,8))
    gs = plt.GridSpec(4,4)
    gs.update(wspace=0.005, hspace=0.005)

    #'''
    ax10 = fig.add_subplot(gs[3,0])
    ax10.scatter(dat['A'], dat['C'], color=colors, marker='.', alpha=0.25)
    ax10.set_xlabel('A', fontsize=20)
    ax10.set_ylabel('C', fontsize=20)
    ax10.set_ylim(0., 5.)
    ax10.set_xlim(0., .8)
    ax10.set_xticks([0.2, 0.4, 0.6])
    ax10.set_yticks([1,2,3,4])
    #ax10.tick_params(axis='x', pad=8)
    #plt.text(0.5, 0.5, '10', fontsize=20, color='red', transform=ax10.transAxes)

    ax10.text(0,19, 'Ellipticals', fontsize=16, color='red')
    ax10.text(0,18, 'Disks', fontsize=16, color='purple')
    ax10.text(0,17, 'Edge-on Disks', fontsize=16, color='green')
    ax10.text(0,16, 'Mergers', fontsize=16, color='yellow')

    
    ################
    
    ax8 = fig.add_subplot(gs[3,1], sharey=ax10)
    plt.setp(ax8.get_yticklabels(), visible=False)
    ax8.scatter(dat['G'], dat['C'], color=colors, marker='.', alpha=0.25)
    ax8.set_xlabel('G', fontsize=20)
    #ax8.set_ylabel('C')
    ax8.set_xlim(0.3, 0.7)
    ax8.set_ylim(0., 5.)
    ax8.set_xticks([0.4, 0.5, 0.6])
    #plt.text(0.5, 0.5, '8', fontsize=20, color='red', transform=ax8.transAxes)
    
    ax9 = fig.add_subplot(gs[2,1], sharex=ax8)
    plt.setp(ax9.get_xticklabels(), visible=False)
    ax9.scatter(dat['G'], dat['A'], color=colors, marker='.', alpha=0.25)
    #ax9.set_xlabel('G')
    ax9.set_ylabel('A', fontsize=20)
    ax9.set_xlim(0.3, 0.7)
    ax9.set_ylim(0., .8)
    ax9.set_yticks([0.2, 0.4, 0.6])
    #plt.text(0.5, 0.5, '9', fontsize=20, color='red', transform=ax9.transAxes)
    
    ################
    
    ax5 = fig.add_subplot(gs[3,2], sharey=ax10)
    plt.setp(ax5.get_yticklabels(), visible=False)
    ax5.scatter(dat['M20'], dat['C'], color=colors, marker='.', alpha=0.25)
    ax5.set_xlabel('M20', fontsize=20)
    #ax5.set_ylabel('C')
    ax5.set_xlim(0., -3.)
    ax5.set_ylim(0., 5.)
    ax5.set_xticks([-2.5, -1.5, -0.5])
    #plt.text(0.5, 0.5, '5', fontsize=20, color='red', transform=ax5.transAxes)
    
    ax6 = fig.add_subplot(gs[2,2], sharex=ax5, sharey=ax9)
    plt.setp((ax6.get_xticklabels(), ax6.get_yticklabels()), visible=False)
    ax6.scatter(dat['M20'], dat['A'], color=colors, marker='.', alpha=0.25)
    #ax6.set_xlabel('M20')
    #ax6.set_ylabel('A')
    ax6.set_xlim(0., -3.)
    ax6.set_ylim(0., .8)
    #plt.text(0.5, 0.5, '6', fontsize=20, color='red', transform=ax6.transAxes)
    
    ax7 = fig.add_subplot(gs[1,2], sharex=ax5)
    plt.setp(ax7.get_xticklabels(), visible=False)
    ax7.scatter(dat['M20'], dat['G'], color=colors, marker='.', alpha=0.25)
    #ax7.set_xlabel('M20')
    ax7.set_ylabel('G', fontsize=20)
    ax7.set_xlim(0., -3.)
    ax7.set_ylim(0.3, 0.7)
    ax7.set_yticks([0.4, 0.5, 0.6])
    #plt.text(0.5, 0.5, '7', fontsize=20, color='red', transform=ax7.transAxes)
    
    ################
    
    ax1 = fig.add_subplot(gs[3,3], sharey=ax10)
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax1.scatter(dat['elipt'], dat['C'], color=colors, marker='.', alpha=0.25)
    ax1.set_xlabel('elipt', fontsize=20)
    #ax1.set_ylabel('C')
    ax1.set_ylim(0., 5.)
    ax1.set_xticks([0.2, 0.4, 0.6, 0.8])
    #plt.text(0.5, 0.5, '1', fontsize=20, color='red', transform=ax1.transAxes)
    
    ax2 = fig.add_subplot(gs[2,3], sharex=ax1, sharey=ax9)
    plt.setp((ax2.get_xticklabels(), ax2.get_yticklabels()), visible=False)
    ax2.scatter(dat['elipt'], dat['A'], color=colors, marker='.', alpha=0.25)
    #ax2.set_xlabel('elipt')
    #ax2.set_ylabel('A')
    ax2.set_ylim(0., .8)
    #plt.text(0.5, 0.5, '2', fontsize=20, color='red', transform=ax2.transAxes)
    
    ax3 = fig.add_subplot(gs[1,3], sharex=ax1, sharey=ax7)
    plt.setp((ax3.get_xticklabels(),ax3.get_yticklabels()), visible=False)
    ax3.scatter(dat['elipt'], dat['G'], color=colors, marker='.', alpha=0.25)
    #ax3.set_xlabel('elipt')
    #ax3.set_ylabel('G')
    ax3.set_ylim(0.3, 0.7)
    #plt.text(0.5, 0.5, '3', fontsize=20, color='red', transform=ax3.transAxes)
    
    ax4 = fig.add_subplot(gs[0,3])
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax4.scatter(dat['elipt'], dat['M20'], color=colors, marker='.', alpha=0.25)
    #ax4.set_xlabel('elipt')
    ax4.set_ylabel('M20', fontsize=20)
    ax4.set_ylim(0., -3.)
    ax4.set_yticks([-2.5, -1.5, -0.5])
    #plt.text(0.5, 0.5, '4', fontsize=20, color='red', transform=ax4.transAxes)

    #gs.tight_layout(fig)
    plt.savefig('compare_parameters_'+fout+'.png')
    plt.show()
    plt.close()

def lotz_mergers(m20):
    return -0.14*m20+0.33

def lotz_separation(m20):
    return  0.14*m20+0.8

def investigate_classification_contamination(cat):
    cat['t_color'] = np.array([color.rstrip() for color in cat['t_color']])

    true_neg = cat[np.where(cat['correct'] & (cat['probability']<.5))]
    strongTypeI = cat[np.where((~cat['correct']) & (cat['probability']>.5))]
    titles1 = ['NE classified as E', 'NE classified as NE']

    true_pos = cat[np.where((cat['correct']) & (cat['probability']>.5))]
    strongTypeII = cat[np.where((~cat['correct']) & (cat['probability']<.5))]
    titles2 = ['E classified as NE', 'E classified as E']

    gini_m20_compare(strongTypeI, true_neg, titles1, filename='notEllipticals')
    gini_m20_compare(strongTypeII, true_pos, titles2, filename='Ellipticals')

def select_balanced_sample(cat):
    '''
    NOT FLESHED OUT -- probably won't work
    '''
    # I fucking hate topcat.
    # Select out the sample I want WITHOUT that piece of shit:
    color = cat['t_color']
    smooth = cat['t01_smooth_or_features_a01_smooth_debiased']

    balance_cut = np.where((color=="purple") | (color=="yellow") | 
                           (color=="green") | ((color=="red")&(smooth>=0.68)))
    balanced = cat[balance_cut]


    color = ['purple', 'red', 'green', 'yellow']
    for c in color:
        cut = np.where(balanced['t_color']==c)
        num = len(cat[cut])
        per = float(len(balanced[cut]))/len(balanced)
        print num, per

def select_galaxy_classes(cat):

    merger_cut = np.where(cat['G'] > -0.115*cat['M20']+0.384)
    merger = cat[merger_cut]

    merger_disk_cut = np.where(merger['t_color']=='purple')
    merger_disk = merger[merger_disk_cut]
    
    high_m20_cut = np.where( (merger_disk['M20'] < -1.) & 
                             (merger_disk['M20'] > -1.5) )
    high_m20 = merger_disk[high_m20_cut]


    elliptical_cut = np.where((cat['G'] < -0.115*cat['M20']+0.384) & 
                              (cat['G'] > 0.115*cat['M20']+0.769))

    disk_ell_cut = np.where(cat[elliptical_cut]['t_color']=='purple')
    disk_ell = cat[disk_ell_cut]

    irregular_cut = np.where((cat['G'] > 0.115*cat['M20']+0.697) &
                          (cat['G'] < 0.115*cat['M20']+0.769) &
                          (cat['G'] < -0.115*cat['M20']+0.384))

    spiral_cut = np.where((cat['G'] < -0.115*cat['M20']+0.384) &
                       (cat['G'] < 0.115*cat['M20']+0.697))
    ell_disk_cut = np.where(cat[spiral_cut]['t_color']=='red')
    ell_disk = cat[ell_disk_cut]
    
    print len(disk_ell), len(ell_disk)


def gini_m20_compare(cat1, cat2, titles, outfile):
    try:
        cat1['M20']=cat1['M20_1']
        cat1['G']=cat1['G_1']
        cat2['M20']=cat2['M20_1']
        cat2['G']=cat2['G_1']
    except:
        pass

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

def gini_m20_compare_morphtypes(cat1, cat2):
    
    colors = ['purple', 'red', 'green', 'yellow']
    labels = ['Disks', 'Ellipticals', 'Edge On', 'Mergers']
    alphas = [.5, .3, .4, .75]
    i = 1
    for c,l,a in zip(colors, labels, alphas):

        cut = np.where(cat1['t_color']==c)
        x1 = np.arange(-3., 0., .1)
        y1 = lotz_mergers(x1)
        x2 = np.arange(-3., -1.6, .1)
        y2 = lotz_separation(x2)

        # plot ZEST's G/M20 plane
        ax1 = fig.add_subplot(220+i)

        ax1.scatter(cat1['M20'][cut], cat1['G'][cut], 
                    color=cat1['t_color'][cut], 
                    marker='^',  alpha=a, label=l)

        ax1.plot(x1, y1,'k',lw=2, label='Merger Line (Lotz 08)')
        ax1.plot(x2, y2,'k--', lw=2, label='Morph Line (Lotz 08)')
        ax1.plot([0,-3.], [.4, .4], 'k--')

        ax1.set_ylabel('Gini', fontsize=14)
        ax1.set_xlabel('M20', fontsize=14)
        ax1.set_title(l, fontsize=16)

        ax1.set_xlim(0., -3.)
        ax1.set_ylim(.3, .75)

        i+=1

    fig.tight_layout()
    plt.savefig('G_M20_18K.png')
    #ax1.legend(loc='best')

    
def A_M20_distances(data):
    cuts = np.where((data['oflag']==0) & (data['uflag']==0) & 
                    (data['M2']!=np.inf))
    cat = data[cuts]
    dist = np.sqrt((cat['Mcx2']-cat['Acx'])**2+(cat['Mcy2']-cat['Acy'])**2)

    bleh = np.where(dist>=10)
    
    print len(cat[bleh])
    print cat['name'][bleh]
    #pdb.set_trace()
    plt.figure()
    binsize=1.
    bins = np.arange(min(dist), max(dist)+binsize, binsize)
    plt.hist(dist, bins=bins)
    plt.show()

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


def main():
    
    parser = argparse.ArgumentParser(description='Compare my catalog with ZEST\nPlots Concentration, Asymmetry, Gini, and M20 of the two catalogs')
    parser.add_argument('catalog1', type=str, 
        help='ZEST catalog name')
    parser.add_argument('catalog2', type=str,
        help='My catalog name')
    args = parser.parse_args()


    cat1 = Table.read(args.catalog1)
    cat2 = Table.read(args.catalog2)
    #cat = cat.filled()


    labels=['Elliptical', 'Circular']
    #labels=['Original', 'New (M20/G fixed)']
    outtag = 'circularRp'

    #gini_m20_compare(cat1, cat2, titles=labels, 
    #                 outfile='SDSS_compare_%s'%outtag)

    Rp_compare(cat1, cat2, xylabels=labels, histlabels=labels, 
               outname='SDSS_compare_%s'%outtag)
    exit()
    
    #compare_parameters(cat, fout='18Ksample')

    f = open('incorrect_gals_2classes.txt', 'w+')
    prefix = 'scp beck@ramon-1.spa.umn.edu:/data/extragal/beck/output/datacube/'
    badgals = cat[~cat['correct']]
    for gal in badgals:
        f.write(prefix+gal['name'].rstrip()+'.fits sdss_toclean/\n')
    f.close()

    exit()

    #bigz = Table.read('catalogs/ZEST_catalog.fits')
    #colorcode(bigz, outname='catalogs/ZEST_catalog_colors.fits')
    #plotM20(zest, orig, 'm20_comp_orig.png')
    #plotM20(zest, new, 'm20_comp_new.png')
    #plotasymmetry(zest, new, 'asym_comp_new.png')
    #plotasymmetry(zest, orig, 'asym_comp_orig.png')



if __name__ == '__main__':
    main()
