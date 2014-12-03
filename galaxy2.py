

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
from photutils import aperture_photometry, EllipticalAnnulus, \
                              EllipticalAperture
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
            self.asym, self.center = self.get_asymmetry(clean_dat, aperture)
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
        
        
    def get_asymmetry(self, clean_dat, ap):
        galcenter = np.array([self.x, self.y])
        imgcenter = np.array([clean_dat.shape[0]/2., clean_dat.shape[1]/2.])
        '''
        # create aperture at center of image
        ap1 = utils2.EllipticalAperture( imgcenter, 1.5*self.rpet, 
                                         1.5*self.rpet/self.e, self.theta, clean_dat)
        # create galaxy and background masks from aperture
        apmask = ap1.aper.astype('float')
        bkgmask = np.logical_not(apmask).astype('float')
        scale = utils2.scale(apmask, clean_dat.shape)
        #'''

        scale = ap.area()/(clean_dat.shape[0]*clean_dat.shape[1]-ap.area())
        delta = imgcenter - galcenter
 
        prior_points = []
        asym = np.zeros((9,2))

        counter = 0
        while True:
            # generate shifts from initial pixel (delta) to 8 surrounding "pixels"
            deltas, points = utils2.generate_deltas(imgcenter, .2, delta)
            #print "new points that *should* match with rotated asym:\n", np.array(points)
            if counter == 0:
                indexes = [0,1,2,3,4,5,6,7,8]
            else:
                # find which values of points and asym were calculated on the prior run
                #for p in points:
                for a,pp in prior_run.iteritems():
                    if pp not in points:
                        print pp, a
                        newdata = utils2.shift_image(clean_dat, 
                                                     deltas[idx][1], deltas[idx][0])
                        image = np.abs(newdata-np.rot90(newdata,2))
                        numerator = aperture_photometry(image, ap, method='exact')
                        denominator = aperture_photometry(newdata, ap, method='exact')
                        n, d = numerator['aperture_sum'], denominator['aperture_sum']
                        bkgasym = ((np.sum(image) - n)/d)*scale
                        asym = n/d-bkgasym, bkgasym
                        pdb.set_trace()
                #indexes = [i for i, item in enumerate(points) if item not in prior_points]

            print indexes
            pdb.set_trace()
            # find asymmetry in the new set of "pixels" -- need to optimize this!!!
            #for idx, d in enumerate(deltas):
            for idx in indexes: 
                # using my codez
                # measure the asymmetry for 9 locations at and around the original delta
                #asym[idx] = utils2.measure_asymmetry(clean_dat, apmask, bkgmask, d, scale)

                print deltas[idx]
                # one way of doing it using photutils
                #newdata = utils2.shift_image(clean_dat, d[1], d[0])
                newdata = utils2.shift_image(clean_dat, deltas[idx][1], deltas[idx][0])
                image = np.abs(newdata-np.rot90(newdata,2))
                numerator = aperture_photometry(image, ap, method='exact')
                denominator = aperture_photometry(newdata, ap, method='exact')
                n, d = numerator['aperture_sum'], denominator['aperture_sum']
                bkgasym = ((np.sum(image) - n)/d)*scale
               
                asym[idx] = n/d-bkgasym, bkgasym
 
            minloc = np.where(asym[:,0] == asym[:,0].min())[0]
            minasym = asym[minloc[0]]
            mindelta = deltas[minloc[0]]

            print minloc #, mindelta

            # if the asymmetry found at the original center is the minimum, we're done
            if asym[0,0] == asym[:,0].min():
                center = imgcenter - mindelta
                asymmetry = minasym
                self.asym_plot(clean_dat, mindelta, center, ap)
                pdb.set_trace()
                return asymmetry, center

            # if not, repeat the process until we find the asymmetry minimum
            else:
                delta = mindelta
                center = points[minloc[0]]
                prior_run = dict((a,p) for a, p in zip(asym[:,0], points))
                #pdb.set_trace()
                #asym = np.array(list(asym[minloc[0]:]) + list(asym[0:minloc[0]]))
                #pdb.set_trace()
                #print "asymmetry after:\n", np.array(asym)
                counter += 1
                prior_points = list(points)
                #print "now asym has been rotated for the new set of points:\n"
                #print np.array(asym)


        #return asymmetry, center

    def asym_plot(self, clean_dat, shift, center, aperture):
        imgcenter = [clean_dat.shape[0]/2., clean_dat.shape[1]/2.]
        newdata = utils2.shift_image(clean_dat, shift[1], shift[0])
        residual = newdata-np.rot90(newdata,2)
        plt.figure()
        plt.imshow(newdata)
        plt.plot(center[0], center[1], 'k+', mew=2, ms=10)
        plt.plot(imgcenter[0], imgcenter[1], 'r+', mew=2, ms=10)
        aperture.plot()
        plt.xlim(imgcenter[0]-3*self.rpet, imgcenter[0]+3*self.rpet)
        plt.ylim(imgcenter[1]-3*self.rpet, imgcenter[1]+3*self.rpet)
        plt.title('Shifted Cleaned Image')
        plt.savefig('output/asyms/'+self.name+'_asym1.png')

        plt.close()
        plt.figure()
        plt.imshow(residual)
        plt.plot(center[0], center[1], 'k+', mew=2, ms=10)
        plt.plot(imgcenter[0], imgcenter[1], 'r+', mew=2, ms=10)
        aperture.plot()
        plt.xlim(imgcenter[0]-3*self.rpet, imgcenter[0]+3*self.rpet)
        plt.ylim(imgcenter[1]-3*self.rpet, imgcenter[1]+3*self.rpet)
        plt.title('Asymmetry Residuals (I - I180)')
        plt.savefig('output/asyms/'+self.name+'_asym2.png')
        #exit()

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
    
    fitsfiles = sorted(glob.glob(args.directory+'01*.fits'))
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
    
    pdb.set_trace()
    info = Table(rows=[g.__dict__ for g in galaxies])
    info.write('bigsample_mycatv5.fits', overwrite=True)

    exit()  


if __name__ == '__main__':
    main()
    
'''
        shape = clean_dat.shape
        galcenter = np.array([self.x, self.y])
        imgcenter = np.array([shape[0]/2., shape[1]/2.])

        # create aperture at center of image
        ap = utils2.EllipticalAperture( imgcenter, 1.5*self.rpet, 
                                        1.5*self.rpet/self.e, self.theta, clean_dat)
        # create galaxy and background masks from aperture
        apmask = ap.aper.astype('float')
        bkgmask = np.logical_not(apmask).astype('float')

        scale = utils2.scale(apmask, clean_dat.shape)
        delta = imgcenter - galcenter
 
        prior_points = np.zeros((9,2))
        asym = np.zeros((9,2))
        counter = 0
        while True:
            # generate shifts from initial pixel (delta) to 8 surrounding "pixels"
            deltas, points = utils2.generate_deltas(imgcenter, .2, delta)
            # find asymmetry in the new set of "pixels" -- need to optimize this!!!
            for idx, d in enumerate(deltas):
                # measure the asymmetry for 9 locations at and around the original delta
                asym[idx] = utils2.measure_asymmetry(clean_dat, apmask, bkgmask, d, scale)

            minloc = np.where(asym[:,0] == asym[:,0].min())[0]
            minasym = asym[minloc[0]]
            mindelta = deltas[minloc[0]]

            # if the asymmetry found at the original center is the minimum, we're done
            if asym[0,0] == asym[:,0].min():
                center = imgcenter - mindelta
                asymmetry = minasym
            
                self.asym_plot(clean_dat, mindelta, center)
                plt.savefig('asymfig_num'+str(counter)+'.png')

                exit()
                return asymmetry, center
            # if not, repeat the process until we find the asymmetry minimum
            else:
                delta = mindelta
                center = points[minloc[0]]
                self.asym_plot(clean_dat, mindelta, center)
                plt.savefig('asymfig_num'+str(counter)+'.png')
                counter += 1
                prior_points = points.copy()

        #return asymmetry, center
'''
