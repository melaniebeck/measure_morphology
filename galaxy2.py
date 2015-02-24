

import glob
import argparse
import os
import string
import pdb
import pyfits as fits
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from scipy.interpolate import interp1d
from astropy.table import Table
from collections import defaultdict
from random import gauss
from photutils import aperture_photometry, EllipticalAnnulus, \
                              EllipticalAperture
import scipy.ndimage.interpolation as sp_interp
#import scipy.ndimage.interpolation.rotate as sp.rotate
import utils2

class Galaxy(object):

    def __init__(self, hdulist, filename): #, cdat
        # the HDUList should contain:
        #   original image      hdulist[0]
        #   cleaned image       hdulist[1]
        #   bright seg map      hdulist[2]
        #   faint segmap        hudlist[3]
        #   SE catalog          hdulist[4]

        clean_dat = hdulist[1].data
        segmap = hdulist[3].data
        cat = hdulist[4].data[hdulist[4].header['SECATIDX']]
        
        name = os.path.basename(filename)
        name = string.split(name,'.fits')
        self.name = name[0]

        # checks? are performed outside of the galaxy initialization in __init__
        # but additional ones should probably be performed here as well (someday...)
                
        # initialize SExtractor parameters        
        self.e = cat['ELONGATION']
        self.x, self.y = cat['X_IMAGE'], cat['Y_IMAGE']
        self.kron = cat['KRON_RADIUS']
        self.a, self.b = cat['A_IMAGE'], cat['B_IMAGE']
        self.theta = cat['THETA_IMAGE']*pi/180. # in radians?
        self.ra, self.dec = cat['ALPHA_J2000'], cat['DELTA_J2000']

        # initialize morphological parameters
        self.rpet, self.rpetflag = self.get_petro(clean_dat)
        imgcenter = [clean_dat.shape[0]/2, clean_dat.shape[1]/2]

        if not np.isnan(self.rpet):
            aperture = EllipticalAperture(imgcenter, 1.5*self.rpet, 
                                          1.5*self.rpet/self.e, self.theta)                            
            self.asym, self.center = self.get_asymmetry(clean_dat, aperture, segmap)
            #self.gini = self.get_gini(clean_dat, aperture)
            #self.conc = self.get_concentration(clean_dat, aperture)
            #self.m20 = self.get_m20(clean_dat, aperture)
        else:
            self.asym, self.center = np.nan, (np.nan, np.nan)
            #self.gini = np.nan
            #self.conc = np.nan
            #self.m20 = np.nan
        self.elipt = cat['ELLIPTICITY']            

        #dir(self) in cmd line to see all the hidden shits
        
    def get_petro(self, clean_dat):

        r_flag = 1
        
        imgsize = clean_dat.shape[0]/2.
        se_apsize = self.a*self.kron 
        positions = [self.x, self.y]
        
        # condition of np.log10(imgsize/constant) ensures that the maximum
        # radius will never exceed the size of the image
        apix = 10*np.logspace(-1.0, np.log10(imgsize/10.), num=20)

        bpix = apix/self.e

        flux1, flux2, ans1, ans2 = [], [], [], []
        for idx in range(1,len(apix)-1):
            # for SB, annuli from radius "before" to radius "after" the current
            # radius in apix
            an1 = EllipticalAnnulus(positions, apix[idx-1], apix[idx+1], 
                                bpix[idx+1], self.theta)
            ans1.append(an1)
            flux1.append(aperture_photometry(clean_dat, an1, method='exact'))
            
            # for Avg SB, annuli from radius "before" to current radius in apix
            an2 = EllipticalAnnulus(positions, apix[idx-1], apix[idx], 
                                bpix[idx], self.theta)
            ans2.append(an2) 
            flux2.append(aperture_photometry(clean_dat, an2, method='exact'))

        phot_table1 = np.hstack(flux1)
        counts1 = phot_table1['aperture_sum']
        annuli1 = np.hstack(ans1)  
              
        phot_table2 = np.hstack(flux2)
        counts2 = phot_table2['aperture_sum']
        annuli2 = np.hstack(ans2)  
              
        areas1, areas2 = [], []      
        for an1, an2 in zip(annuli1, annuli2):
            #an1.plot(color='blue', lw=1.5, alpha=0.5)
            areas1.append(an1.area()) 
            areas2.append(an2.area())                      
        
        # calculate the average surface brightness
        csum, asum = [], []    
        for idx in range(1,len(counts2)+1):
            csum.append(sum(counts2[0:idx]))
            asum.append(sum(areas2[0:idx]))
        csum = np.array(csum)
        asum = np.array(asum)
               
        # SB and avg SB as a function of radius
        sb = counts1/areas1      
        avgsb = csum/asum  
        
        # make some plots! (Eventually make this optional - user controlled)
        #self.petro_plot(clean_dat, annuli1, apix, sb, avgsb)
 
        # now we need to find the intersection of u(R)/<u(R)> with 0.2
        x = apix[1:-1]
        radii = np.linspace(np.min(x), np.max(x), num=1000)
        f = interp1d(x, sb/avgsb, kind='cubic')        
        ratios = f(radii)

        rpet = np.nan
        if not np.any(np.isnan(ratios)):
            dirty = ratios-0.2
            # if any value in dirty is non negative, shit be rockin ... 
            if (dirty < 0).any():
                for idx,thing in enumerate(dirty):
                    if thing < 0:
                        rpet = np.mean([radii[idx-1],radii[idx]])
                        #print rpet
                        return rpet, r_flag
            # if not, then shit be broked ... 
            else: 
                print "Ratio never crosses 0.2. Cannot calculate Petrosian radius."
                r_flag = 2 
        # if there are nans, shit be broked ... 
        else:
            print "Interpoation failed! Cannot calculate Petrosian radius."
            r_flag = 0
        
        return rpet, r_flag     
    
    def petro_plot(self, clean_dat, annuli, radii, sb, avgsb):

        # plot the annuli overlaid on cleaned image
        plt.figure()
        for an in annuli:
            an.plot(color='blue', lw=1.5, alpha=0.5)
        imgplot = plt.imshow(clean_dat, cmap='gray_r', origin='lower')
        imgplot.set_clim(-0.009, 0.022)
        plt.title(self.name)        
        plt.savefig('output/apertures/'+self.name+'_aper.png', bbox_inches='tight')
        plt.close()
        
        # plot the SB, AvgSB and Ratio as a fcn of radius
        rpix = radii[1:len(radii)-1]
        zeros = np.zeros(len(rpix))
        eta = []
        for x in range(0,len(rpix)):
            eta.append(0.2)   
    
        xlims = [rpix[0], rpix[len(rpix)-2]]
        ylims = [-.2*max(avgsb), max(avgsb) + 0.2*max(avgsb)] 
        
        plt.figure()
        plt.subplot(211)
        plt.title('Surface Brightness Profile')
        plt.plot(rpix, sb, 'ro', label='Surface Brightness')
        plt.plot(rpix, avgsb, 'g^', label='Average SB')
        plt.plot(rpix, zeros, 'k--')
        plt.xscale('log')
        plt.axis([xlims[0], xlims[1], ylims[0], ylims[1]])
        plt.legend(loc='upper right')
        
        plt.subplot(212)
        plt.title('u(R)/<u(R)>')
        plt.plot(rpix, sb/avgsb, 'bo', label='SB/Avg_SB')
        plt.plot(rpix, eta, 'k--', label='eta=0.2')  
        plt.xscale('log')
        plt.axis([xlims[0], xlims[1], 0., 1.1])
        plt.legend(loc='upper right')    
        plt.savefig('output/profiles/'+self.name+'_prof.png', bbox_inches='tight')
        
        plt.close() 
               
        return 0
        
        
    def get_asymmetry(self, clean_dat, ap, segmap):
        ''' 
        In this one, we're doing it my way
        #'''

        galcenter = np.array([self.x, self.y])
        imgcenter = np.array([clean_dat.shape[0]/2., clean_dat.shape[1]/2.])
        delta = imgcenter - galcenter

        # to look at how shift works or to re-verify that I did it right...
        #shift_test()

        scale = ap.area()/(clean_dat.shape[0]*clean_dat.shape[1]-ap.area())
        asyms = defaultdict(list)
        prior_points = []
        bkgmin = 100.

        while True:
            '''
            this loop creates two things:
            1. a dictionary that continuously grows in size where each key is a 
               point on the image grid  and an associated asymmetry value 
            2. a list of asymmetry values which will always contain only 9 values.
               these 9 values are the asymmetries for that particular run. the first 
               value in this list will always be the "central" point of the current
               3x3 grid being tested
            '''
            # deltas is a np.array of lists
            # points is a list of tuples
            
            ba = []
            ga = []
            deltas, points = utils2.generate_deltas(imgcenter, .3, delta)

            for d, p in zip(deltas, points): 
                # if the point already exists in the dictionary, don't run asym codes!
                if p not in asyms: 
                    # Doing it my way -------->
                    newdata = sp_interp.shift(clean_dat, d)
                    rotdata = sp_interp.rotate(newdata, 180.)
                    residual = np.abs(newdata-rotdata)
                    numerator = aperture_photometry(residual, ap)
                    denominator = aperture_photometry(np.abs(newdata), ap)
                    num, den = numerator['aperture_sum'], denominator['aperture_sum']

                    bkgasym = float(((np.sum(residual) - num)/den)*scale)
                    galasym = float(num/den)

                    asym = float(num/den - bkgasym)

                    # create an array of asyms ... 
                    ba.append(bkgasym)
                    ga.append(galasym)
                    # ... and a dictionary that maps each asym to a point on the image grid
                    asyms[p].append([galasym, bkgasym])

                # just take the value that's already in the dictionary for that point
                else:
                    ga.append(asyms[p][0][0])
                    ba.append(asyms[p][0][1])
            
            # want to find the min background asymmetry, regardless of which position it's at
            bb = ba[np.where(ba == np.min(ba))[0]]
            if bb < bkgmin:
                bkgmin = bb

            # if the asymmetry found at the original center (first delta in deltas) 
            # is the minimum, we're done
            if ga[0] == np.min(ga):
                center = imgcenter - deltas[0]
                self.asym_plot(newdata, residual, ap, imgcenter, center)
                #pdb.set_trace()
                return ga[0]-bkgmin, center
            else:
                minloc = np.where(ga == np.min(ga))[0]
                #pdb.set_trace()
                delta = deltas[minloc[0]]
                prior_points = list(points)


    def shift_test(self, clean_dat, imgcenter):
        # [row, col] --> [y, x]
        clean_dat[252,255] = .1
        # ----> [251., 251.] - [252., 255.] = [-1., -4.] => [dy, dx]
        # shift works by putting it in [dy,dx] order 
        newdata = sp.shift(clean_dat, [-1, -4])
        rotdata = sp.rotate(newdata, 180.)
        residual = np.abs(newdata - rotdata)

        plt.figure()
        plt.imshow(clean_dat)
        plt.title('original')
        plt.plot(imgcenter[0], imgcenter[1], 'k+', mew=2, ms=10)
        plt.plot(255, 252, 'r+', mew=2, ms=10)
        #plt.xlim(imgcenter[0]-2*self.rpet, imgcenter[0]+2*self.rpet)
        #plt.ylim(imgcenter[1]-2*self.rpet, imgcenter[1]+2*self.rpet)
        
        plt.figure()
        plt.imshow(newdata)
        plt.title('shifted')
        #plt.xlim(imgcenter[0]-2*self.rpet, imgcenter[0]+2*self.rpet)
        #plt.ylim(imgcenter[1]-2*self.rpet, imgcenter[1]+2*self.rpet)
        plt.plot(imgcenter[0], imgcenter[1], 'k+', mew=2, ms=10)
        plt.plot(255, 252, 'r+', mew=2, ms=10)

        plt.figure()
        plt.imshow(rotdata)
        plt.title('rotated')
        #plt.xlim(imgcenter[0]-2*self.rpet, imgcenter[0]+2*self.rpet)
        #plt.ylim(imgcenter[1]-2*self.rpet, imgcenter[1]+2*self.rpet)
        plt.plot(imgcenter[0], imgcenter[1], 'k+', mew=2, ms=10)
        plt.plot(255, 252, 'r+', mew=2, ms=10)
        
        plt.figure()
        plt.imshow(residual)
        plt.title('residual')
        #plt.xlim(imgcenter[0]-2*self.rpet, imgcenter[0]+2*self.rpet)
        #plt.ylim(imgcenter[1]-2*self.rpet, imgcenter[1]+2*self.rpet)
        plt.plot(imgcenter[0], imgcenter[1], 'k+', mew=2, ms=10)
        plt.plot(255, 252, 'r+', mew=2, ms=10)               
        
        plt.show()
        pdb.set_trace()
        

    #def asym_plot(self, clean_dat, shift, center, aperture):
    def asym_plot(self, shifted, residual, aperture, imgcenter, galcenter):
        '''
        plt.figure()
        plt.imshow(org)
        plt.title('original')
        plt.plot(imgcenter[0], imgcenter[1], 'k+', mew=2, ms=10)
        plt.plot(galcenter[1], galcenter[0], 'r+', mew=2, ms=10)
        
        plt.figure()
        plt.imshow(shift)
        plt.title('shifted')
        plt.plot(imgcenter[0], imgcenter[1], 'k+', mew=2, ms=10)
        plt.plot(galcenter[1], galcenter[0], 'r+', mew=2, ms=10)

        plt.figure()
        plt.imshow(rot)
        plt.title('rotated')
        plt.plot(imgcenter[0], imgcenter[1], 'k+', mew=2, ms=10)
        plt.plot(galcenter[1], galcenter[0], 'r+', mew=2, ms=10)
        
        plt.figure()
        plt.imshow(resid)
        plt.title('residual')
        plt.plot(imgcenter[0], imgcenter[1], 'k+', mew=2, ms=10)
        plt.plot(galcenter[1], galcenter[0], 'r+', mew=2, ms=10) 
        ap.plot()
        plt.show()
        '''
        #imgcenter = [clean_dat.shape[0]/2., clean_dat.shape[1]/2.]
        #newdata = utils2.shift_image(clean_dat, shift[1], shift[0])
        #residual = newdata-np.rot90(newdata,2)
        low = imgcenter[0]-2*self.rpet
        high = imgcenter[1]+2*self.rpet
        plt.figure()
        plt.imshow(shifted)
        plt.plot(galcenter[1], galcenter[0], 'k+', mew=2, ms=10)
        plt.plot(imgcenter[1], imgcenter[0], 'r+', mew=2, ms=10)
        aperture.plot()
        plt.xlim(low, high)
        plt.ylim(low, high)
        plt.title('Shifted Cleaned Image')
        plt.savefig('output/asyms/'+self.name+'_asym1.png')
        plt.close()

        plt.figure()
        plt.imshow(residual)
        plt.plot(galcenter[1], galcenter[0], 'k+', mew=2, ms=10, label="Asymmetry Center")
        plt.plot(imgcenter[1], imgcenter[0], 'r+', mew=2, ms=10, label="Image Center")
        aperture.plot()
        plt.xlim(low, high)
        plt.ylim(low, high)
        plt.title('Asymmetry Residuals (I - I180)')
        plt.legend()
        plt.savefig('output/asyms/'+self.name+'_asym2.png')
        plt.close()
        return 0.

    def get_gini(self, clean_dat, ap):
        return 0.
            
    def get_concentration(self, clean_dat, ap):
        return 0.
        
    def get_m20(self, clean_dat, ap):
        return 0.
        
    def show(self, clean_dat):
        return 0.


####################### main ############################

def main():
    
    parser = argparse.ArgumentParser(description='Perform LLE/PCA/whatevs')
    parser.add_argument('directory', type=str, 
        help='Directory of fits images on which to run LLE.')
    #parser.add_argument('SEconfig', type=str,
    #    help='Specify the config.sex file for sextractor configuration parameters')
    args = parser.parse_args()
    
    fitsfiles = sorted(glob.glob(args.directory+'*.fits'))
    #fitsfiles = sorted(glob.glob(args.directory))
        
    galaxies = []
    for idx, f in enumerate(fitsfiles): 
        filename = 'output/f_'+os.path.basename(f)
        if not os.path.isfile(filename):
            print "File not found! Running SExtractor before proceeding."
            print "Cleaning ", os.path.basename(f)
            utils2.clean_frame(f)

        #else:
        # run everything else
        print "Running", os.path.basename(f)
        hdulist = fits.open(filename, memmap=True)
        galaxies.append(Galaxy(hdulist,filename)) #
        hdulist.close()
        #pdb.set_trace()
    
    pdb.set_trace()
    info = Table(rows=[g.__dict__ for g in galaxies])
    info.write('bigsample_mycatv5.fits', overwrite=True)

    exit()  


if __name__ == '__main__':
    main()
    

    '''
    baa = []
    # minimize the background first
    for d, p in zip(deltas, points):
    if p not in asyms: 
    
    bkgdata = sp_interp.shift(bkg_img, d)
    bkgrot = sp_interp.rotate(bkgdata, 180.)
    b_resid = np.abs(bkgdata-bkgrot)
    b_num = aperture_photometry(b_resid, aperture)
    b_den = aperture_photometry(np.abs(bkgdata), aperture)
    bn, bd = b_num['aperture_sum'], b_den['aperture_sum']
    basym = bn/(2*bd)
    baa.append(basym)
    # ... and a dictionary that maps each asym to a point on the image grid
    basyms[p].append(basym)
    
    if baa[0] == np.min(baa):
    center = imgcenter - deltas[0]
    #self.asym_plot()
    return baa[0], center
    else:
    minloc = np.where(baa == np.min(baa))[0]
    delta = deltas[minloc[0]]
    prior_points = list(points)
    print delta
    
    print baa
    pdb.set_trace()
    
    # Doing it Claudia's way
    galdata = sp_interp.shift(smaller, d)
    galrot = sp_interp.rotate(galdata, 180.)
    g_resid = np.abs(galdata-galrot)
    
    
    g_num = aperture_photometry(g_resid, aperture)
    g_den = aperture_photometry(np.abs(galdata), aperture)
    b_num = aperture_photometry(b_resid, aperture)
    b_den = aperture_photometry(np.abs(bkgdata), aperture)
    
    gn, gd = g_num['aperture_sum'], g_den['aperture_sum']
    
    asym = gn/(2*gd) - bn/(2*bd)
    
    #self.asym_plot(smaller, galdata, galrot, g_resid, aperture, imgcenter_sm, p)
    self.asym_plot(bkg_img, bkgdata, bkgrot, b_resid, aperture, imgcenter_sm, p)
    pdb.set_trace()
    
    '''
    
