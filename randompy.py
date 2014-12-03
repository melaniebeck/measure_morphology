# A random script for doing random things that need doing

import os
import string
import pyfits as fits
import numpy as np
from astropy.table import Table
from scipy.spatial import cKDTree
from scipy.interpolate import interp1d
from collections import OrderedDict, defaultdict
from random import gauss
import matplotlib.pyplot as plt
from math import pi
import pdb #"""for doing an IDL-like stop"""
from utils import find_closest



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
    dat = Table.read('bigsample_mycatv3.fits')
    zest = Table.read('bigsample_zest.fits')
    
    # a few of the large sample just couldn't be found by SExtractor in my code
    # mask them out of the ZEST sample for now
    '''
    zest['num']np.arange(1,len(zest)+1)

    flagged = np.array([19, 21, 62, 356, 357, 361, 478])
    idx = np.arange(len(zest))
    mask = np.ones(idx.shape, dtype=bool)
    mask[flagged] = False
    
    included = idx[mask]
    excluded = idx[~mask]  
   
    zest = zest[included]
    #'''

    zrad = zest['rpet']
    mrad = dat['rpet']

    mu = np.nanmean(mrad/zrad)
    sig = np.nanstd(mrad/zrad)
    
    # Histogram of magnitude distribution
    plt.figure()
    plt.hist(zest['acs_mag_auto'], bins=40, normed=1)
    ax = plt.gca()
    ax.invert_xaxis()
    plt.title('Magnitude distribution for 964 galaxies')
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
    
    plt.savefig('petro_ratio_bigsample_v3.png')
    #plt.show()
    
    #exit()

    # My radius vs. ZEST radius and ratio
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(zrad, mrad, 'ro', [0,200], [0,200], 'k--')
    plt.title('PetroRad comparison 964 galaxies')
    plt.xlabel('ZEST [pixels]')
    plt.ylabel('Mine [pixels]')
    
    plt.subplot(2,1,2)
    plt.plot(zrad, mrad/zrad, 'ro')
    plt.plot([0,200], [1,1], 'k--')
    plt.yscale('log')
    plt.xlabel('Radius [pixels]')
    plt.ylabel('Mine/ZEST')
    
       
    subset = np.where(mrad/zrad < .5)
    i = 0
    for xy in zip(zrad[subset],mrad[subset]/zrad[subset]):
        plt.annotate(str(subset[0][i]), xy=xy)
        i += 1
    #'''
    plt.savefig('petro_bigsample_label_v3.png')

    print len(subset[0])
    
    #plt.show()
    plt.close()
    #plt.clf()
    
    #exit()
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
    
    high = np.where(mrad/zrad >=1.3)
    low = np.where(mrad/zrad <=0.7)
    
    totes = len(high[0])+len(low[0])
    
    ratio = float(totes)/len(mrad)*100.   
    print ratio

    #exit()

def plotasymmetry():
    dat = Table.read('bigsample_mycatv5.fits')
    zest = Table.read('bigsample_zest.fits')

    zrad = zest['aa']
    mrad = dat['asym'][:,0]-dat['asym'][:,1]

    mu = np.nanmean(mrad/zrad)
    sig = np.nanstd(mrad/zrad)

    plt.figure()
    plt.hist(mrad, bins=20, range=[-.1,1.6])
    plt.title('My Asymmetry')
    plt.figure()
    plt.hist(zrad, bins=20, range=[-.1, 1.])
    plt.title('ZEST Asymmetry')
    plt.show()

    # Histogram of MRAD/ZRAD ratio with mean/sigma
    plt.figure()
    plt.hist(mrad/zrad, bins=100, range=[-10.,20.], color='blue')
    #plt.plot([mu,mu], [0,1], lw=4)
    #plt.axvline(x=mu, linewidth=3, color='k')
    #plt.axvline(x=mu+sig, linewidth=3, color='k', linestyle='--')
    #plt.axvline(x=mu-sig, linewidth=3, color='k', linestyle='--')
    plt.legend(loc='upper right')
    plt.xlabel('Asymmetry')
    plt.savefig('asymmetry_bigsample_v1.png')

    # My radius vs. ZEST radius and ratio
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(zrad, mrad, 'ro', [0,2], [0,2], 'k--')
    plt.title('Asymmetry comparison 971 galaxies')
    plt.xlabel('ZEST')
    plt.ylabel('Mine')
    
    plt.subplot(2,1,2)
    plt.plot(zrad, mrad/zrad, 'ro')
    plt.plot([0,2], [1,1], 'k--')
    plt.ylim([-.5, 200.])
    #plt.yscale('log')
    plt.xlabel('ZEST Asymmetry')
    plt.ylabel('Mine/ZEST')
    plt.savefig('asymmetry_ratio_bigsample_v1.png')
    
    #pdb.set_trace()
    plt.show()

    exit()

def main():

    #comparecatalogs()
    #plotzestcassatamine()
    
    #selectnewsample()
    #plotpetrorad_largesample()
    
    plotasymmetry()

    #pdb.set_trace()
    
    
    
if __name__ == '__main__':
    main()
    
    '''   
    # construct a new table? 
    c1 = fits.Column(name='RA', format='D10', array=infos['ra'])
    c2 = fits.Column(name='DEC', format='D10', array=infos['dec'])
    c3 = fits.Column(name='A_IMAGE', format='D10', array=infos['a'])
    c4 = fits.Column(name='B_IMAGE', format='D10', array=infos['b'])
    c5 = fits.Column(name='THETA', format='D10', array=infos['theta'])
    c6 = fits.Column(name='ELONGATION', format='D10', array=infos['e'])
    c7 = fits.Column(name='PETRO_RAD', format='D10', array=infos['rpet'])
       
    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7])
    tbhdu.writeto('.fits', clobber=True)                                

    #for gal in galaxies:
    #    rpet.append(gal.rpet) 

#    with open('prad2.txt','w') as f:
#        for thing in rpet:
#            f.write("%f\n" %thing)
    #'''
    
# construct a dictionary? 
#infos = dict( zip(['ra', 'dec', 'a', 'b', 'theta', 'elongation', 'petro_rad'],
#                 [g.ra, g.dec, g.a, g.b, g.theta, g.e, g.rpet]) for g in galaxies)  
    
# for now let's explore the catalog of ZEST classified galaxies
# to determine which ones we should make postage stamps of

#f = open('IpacTableFromSource.tbl',r)
#mytable = Table.read('IpacTableFromSource.tbl',format='ascii.ipac')

#smalltable = mytable[0:10]
#smalltable.write('test.tex',format='ascii.latex
#print mytable.colnames
'''
bright = mytable[mytable['acs_mag_auto'] < 20.]
bright = bright[bright['acs_mag_auto'] > 19.95]
bright = bright[bright['acs_clean'] == 1.]
bright = bright[bright['acs_mu_class'] == 1.]
'''
#print len(bright)
#print min(bright['aa']), max(bright['aa'])

#pdb.set_trace()

#bright.write('bright_gals.fits')#,format='ascii.ipac')

#rpet =  bright['rpet']

#print rpet

# I've created some cutouts to begin working on calculating their asymmetry
# These cutouts are galaxies that have previously been classified by ZEST
# and have all parameters measured for comparison. 

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
