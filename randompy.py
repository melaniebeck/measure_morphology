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

def plotasymmetry():
    dat1 = Table.read('bigsample_mycatv7.fits') 
    dat2 = Table.read('bigsample_mycat_v10.fits') 
    zest = Table.read('bigsample_zest.fits')
    #pdb.set_trace()
    za = zest['aa']
    ma1 = dat1['asym']/2.    
    ma2 = dat2['A']/2. 
    '''
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(za, ma1, 'ro', [-.2,.5], [-.2,.5], 'k')
    plt.ylim(-0.3, .8)
    plt.xlim(-.2, .5)
    plt.axvline(x=0., linestyle='--', color='k')
    plt.axhline(y=0., linestyle='--', color='k')
    #plt.axvline(x=.2, linestyle='--', color='b')
    #plt.axhline(y=.2, linestyle='--', color='b')
    plt.title('Asym with 1.5*Rp')
    plt.xlabel('ZEST')
    plt.ylabel('Mine')
    '''
    #plt.subplot(1,2,2)
    plt.plot(za, ma2, 'ro', [-.2,.5], [-.2,.5], 'k')
    plt.axvline(x=0., linestyle='--', color='k')
    plt.axhline(y=0., linestyle='--', color='k')
    #plt.axvline(x=.2, linestyle='--', color='b')
    #plt.axhline(y=.2, linestyle='--', color='b')
    plt.ylim(-0.3, .8)
    plt.xlim(-.2, .5)
    plt.title('Asymmetry within 1 Petrosian radius')
    plt.xlabel('ZEST')
    plt.ylabel('Mine')
    '''
    subset = np.where((za >= 0.2) & (ma2 <= 0.2))
    
    i = 0
    for xy in zip(za[subset], ma2[subset]):
        #print xy
        plt.annotate( '%.4f' %dat2['med'][subset[0][i]], xy=xy)
        i += 1
    '''
    plt.show()
    plt.savefig('asym_bigsample.png')
    pdb.set_trace()

    '''
    delx, dely=[],[]
    for idx, c in  enumerate(dat1['center']):
        #pdb.set_trace()
        delx.append(dat1['x'][idx]-dat1['center'][0][0])
        dely.append(dat1['y'][idx]-dat1['center'][0][1])

    delx = np.array(delx)
    dely = np.array(dely)
    #pdb.set_trace()
    

    msub = dat1[subset]
    zsub = zest[subset]
    #delx = delx[subset]
    #dely = dely[subset]
    zasym = zrad[subset]
    zrpet = zest['rpet'][subset]
    flag = zest['flagrpet'][subset]
    bigtable = Table([msub['name'], mrad1[subset], zsub['aa'], msub['rpet'], zsub['rpet'], \
                          zsub['flagrpet'], msub['e'], msub['theta']*180./pi], 
                     names=('name', 'myasym', 'zestasym', 'myrpet', 'zestrpet', 'flag',  \
                                'myelongation', 'mytheta'))
    bigtable.write('asym_1rpet.txt', format='ascii.fixed_width')

    f = open('asym_1rpet.sh', 'w+')
    f.write('ds9 ')
    for name in msub['name']: 
        f.write('%s.fits[1] '%name)
    f.close()
    '''

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

def plotgini(zestdata, mydata):
    zest = fits.getdata(zestdata)
    dat = fits.getdata(mydata)

    zg = zest['gg']
    mg = dat['G']
    
    zgwut = zg[mg > 1.]
    wut = mg[mg > 1.]
    #pdb.set_trace()

    plt.figure()
    plt.plot(zg, mg, 'ro', [0,1], [0,1], 'k--', zgwut, wut, 'bo')
    plt.ylim(0., 6.)
    plt.xlim(.2, .8)
    plt.xlabel('ZEST Gini')
    plt.ylabel('My Gini')
    plt.savefig('compareGini_bad.png')

    plt.figure()
    plt.hist(mg/zg, bins=60, range=(0., 4.5))
    plt.title('My Gini / ZEST Gini')
    #plt.hist(mg, bins=40, range=(0,1), alpha=0.5)
    #plt.hist(zg, bins=20, range=(0,1), alpha=0.5)
    plt.savefig('histGini_bad.png')

    #plt.show()
    #pdb.set_trace()

    #'''
    ratio = mg/zg
    sample_b = dat[np.where(ratio > 1.15)]
    sample_s = dat[np.where(ratio < 0.9)]
    f = open('Gtest_b.sh', 'w+')
    f.write('ds9 ')
    for name in sample_b['name']: 
        f.write('output/datacube/%s.fits[1] '%name)
    f.write('-single -lock frame image -scale zscale -zoom to fit')
    f.close
    f = open('Gtest_s.sh', 'w+')
    f.write('ds9 ')
    for name in sample_s['name']: 
        f.write('output/datacube/%s.fits[1] '%name)
    f.write('-single -lock frame image -scale zscale -zoom to fit')
    f.close
    #'''

def plotM20(zestdata, dataset):
    zest = fits.getdata(zestdata)
    dat = fits.getdata(dataset)

    #pdb.set_trace()

    zm = zest['m20']
    mm = dat['M20']

    plt.figure()
    plt.plot(zm, mm, 'ro', [-3,0], [-3,0])
    plt.title("M20 comparison between ZEST and mine")
    plt.xlabel('ZEST M20')
    plt.ylabel('My M20')
    plt.savefig('M20_test.png')

    plt.figure()
    plt.plot(zm, mm/zm, 'ro')
    plt.axhline(y=1., ls='--')
    plt.savefig('M20_ratio.png')
    #plt.show()

    ratio = mm/zm

    #mtest = dat[mm > -0.5]
    mtest = dat[ratio < 0.5]
    f = open('M20_bad2.sh', 'w+')
    f.write('ds9 ')
    for name in mtest['name']:
        f.write('output/datacube/%s.fits[1] '%name)
    f.write('-single -lock frame image -scale zscale -zoom to fit')
    f.close()

    f = open('M20_bad_cont2.sh', 'w+')
    for name in mtest['name']:
        f.write('evince output/m20/%s_m20.pdf \n' %name)
    f.close()

    #pdb.set_trace()

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
    

    
def main():
    
    parser = argparse.ArgumentParser(description='Compare my catalog with ZEST\nPlots Concentration, Asymmetry, Gini, and M20 of the two catalogs')
    parser.add_argument('catalog1', type=str, 
        help='ZEST catalog name')
    parser.add_argument('catalog2', type=str,
        help='My catalog name')
    args = parser.parse_args()

    #comparecatalogs()
    #plotzestcassatamine()
    #selectnewsample()
    #plotpetrorad_largesample()
    #plotasymmetry()
    #plotconcentration()
    #testbackground()
    #plotgini(zestdata, dataset)
    #plotM20(zestdata, dataset)
    #analyse_cleaning()
    #pdb.set_trace()

    z = Table.read(args.catalog1, format='fits')
    m = Table.read(args.catalog2, format='ascii', delimiter=' ')

    names = ['cat', 'oflag', 'uflag', 'ra', 'dec', 'e', 'x', 'y', 'a', 
             'b', 'theta', 'elipt', 'kron', 'Rp', 'Rpflag', 'Rp_SB', 
             'r20', 'r80', 'C', 'A', 'G1', 'G2', 'M1', 'M2']
    names2 = ['Ac', 'Mc1', 'Mc2', 'med', 'rms' ]

    # rename the columns cuz I messed up saving the table
    for i,n in enumerate(names):
        m.rename_column('col'+str(i+1), n)

    fig = plt.figure(figsize=(10,8))
    gs = plt.GridSpec(2,3)

    cut = np.where(m['uflag'] != 0)

    pdb.set_trace()

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
    ax4.plot(z['gg'][cut], m['G2'][cut], 'ro',  [0,1], [0,1], 'k--')
    ax4.set_ylim(-0.25, 1.)
    ax4.set_xlim(0., 1.)
    ax4.set_title('Gini')

    # plot M201
    ax5 = fig.add_subplot(gs[1,1])
    ax5.plot(z['m20'][cut], m['M1'][cut], 'ro', [-3.,0.], [-3.,0.], 'k--')
    ax5.set_title('M20')

    # BONUS PLOT
    ax6 = fig.add_subplot(gs[1,2])
    ax6.plot(z['m20'][cut], m['M2'][cut], 'ro', [-3.,0.], [-3.,0.], 'k--')
    ax6.set_title('M20')


    gs.tight_layout(fig)
    plt.savefig('compare_catalogs.png')
    plt.show()
    plt.close()

    pdb.set_trace()

if __name__ == '__main__':
    main()


'''
test = np.genfromtxt(r'testgal17_mel_prof_IDL2.dat')
crpix = test[:,0]
csb = test[:,1]
cavgsb = test[:,2]


# now we need to find the intersection of u(R)/<u(R)> with 0.2
f = interp1d(sb/avgsb, apix[1:len(apix)-1], 
bounds_error=False, assume_sorted=False)
rpet = f(0.2)  

fc = interp1d(csb/cavgsb, crpix, 
bounds_error=False, assume_sorted=False) 
crpet = fc(0.2)

print rpet, crpet

plt.subplot(3,1,1)
plt.title('Surface Brightness')
plt.plot(apix[1:len(apix)-1], sb, 'ro', crpix, csb,'ko') # rpix, sb2, 'go',
#plt.axis([xlims[0], xlims[1], ylims[0], ylims[1]])
plt.xscale('log')

plt.subplot(3,1,2)
plt.title('Avg Surface Brightness')
plt.plot(apix[1:len(apix)-1], avgsb, 'r^', crpix, cavgsb, 'k^') #rpix, avgsb1, 'r^', 
#plt.axis([xlims[0], xlims[1], ylims[0], ylims[1]])        
plt.xscale('log')

plt.subplot(3,1,3)
plt.title('SB/Avg_SB')
plt.plot(apix[1:len(apix)-1], sb/avgsb, 'rs', 
crpix, csb/cavgsb, 'ks', [1,250],[.2, .2], 'k--', 
[crpet, crpet], [0,1.2], 'k--', [rpet, rpet], [0, 1.2], 'r--')
#plt.axis([xlims[0], xlims[1], 0., 1.1])        
plt.xscale('log')

plt.savefig('testgal17.png')
plt.show()

pdb.set_trace()
#'''
