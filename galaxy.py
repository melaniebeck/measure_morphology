

import glob
import argparse
import os
import string
import pdb
import pyfits as fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.interpolation as sp_interp
from math import pi, ceil
from scipy.interpolate import interp1d
from astropy.table import Table
from collections import defaultdict
from random import gauss
from skimage import measure
from photutils import aperture_photometry, EllipticalAnnulus, \
                              EllipticalAperture, CircularAnnulus, \
                              CircularAperture
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
        self.med, self.rms = self.background(clean_dat, segmap)

        if not np.isnan(self.rpet):
            #aperture = EllipticalAperture(imgcenter, 1.5*self.rpet, 
            #                              1.5*self.rpet/self.e, self.theta)                            
            #self.A, self.ac = self.get_asymmetry(clean_dat, segmap)
            #self.C = self.get_concentration(clean_dat)
            #self.G, self.Gflag = self.get_gini(clean_dat)
            self.M20, self.mc = self.get_m20(clean_dat)
        else:
            #self.A, self.ac = np.nan, (self.x, self.y)
            #self.C = np.nan
            #self.G, self.Gflag = np.nan, 1
            self.M20, self.mc = np.nan, (self.x, self.y)
        self.elipt = cat['ELLIPTICITY']            

        #dir(self) in cmd line to see all the hidden shits
      
    def background(self, data, segmap):
        median = np.median(data[segmap==0])
        rms = np.sqrt(np.mean(np.square(data[segmap==0])))

        return median, rms
  
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
    
    def petro_plot(self, clean_dat, radii, sb, avgsb):

        # plot 1 petrosian radius on cleaned image
        an = EllipticalAnnulus((self.x, self.y), self.rpet, self.rpet/self.e, self.theta)
        plt.figure()
        an.plot(color='blue', lw=1.5, alpha=0.5)
        imgplot = plt.imshow(clean_dat, cmap='gray_r', origin='lower')
        imgplot.set_clim(-0.009, 0.022)
        plt.title(self.name)        
        plt.savefig('output/aperture/'+self.name+'_aper.png', bbox_inches='tight')
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
        
        
    def get_asymmetry(self, clean_dat, segmap):
        '''
        In this one, we're doing it Claudia's way
        1. make a smaller image of the galaxy -> 2*petrosian rad
        2. create a background image
        3. create an aperture 1.5*petrosian radius
        4. minimize asymmetry in the bakground img
        5. minimize asymmetry in the galaxy img

        #'''
        print "calculating Asymmetry..."

        galcenter = np.array([self.x, self.y])
        imgcenter = np.array([clean_dat.shape[0]/2., clean_dat.shape[1]/2.])
        delta = imgcenter - galcenter

        low = round(imgcenter[0] - 2*self.rpet)
        high = round(imgcenter[0] + 2*self.rpet)
        if not ((high-low) % 2 == 0.):
            high += 1
            delta = delta - 1
            #print delta
            
        smaller = clean_dat[low:high, low:high]
        imgcenter_sm = [smaller.shape[0]/2., smaller.shape[1]/2.]
        
        # create a new aperture at the new center of the new image!
        aperture = EllipticalAperture(imgcenter_sm, self.rpet, 
                                      self.rpet/self.e, self.theta)

        # create a SQUARE background/noise image approx same size as area of aperture
        # (we need to minimize calculations)
        size = ceil(np.sqrt(ceil(aperture.area()))) 
        bkg_img = np.zeros((size, size))
        mask = np.where(bkg_img == 0)
        # get the median and std of the background -- 
        # pixel values in clean_dat where segmap == 0
        med = np.median(clean_dat[segmap==0])
        rms = np.sqrt(np.mean(np.square(clean_dat[segmap==0])))
        #stddev = np.std(clean_dat[segmap==0])
        for pixel in zip(mask[0], mask[1]):
            bkg_img[pixel] = gauss(med, rms)
        
        # save the background image 
        bkg = fits.ImageHDU(data=bkg_img)
        bkg.writeto('output/bkgimgs/'+self.name+'.fits', clobber=True)
        #'''

        # minimize the background asymmetry --
        # doesn't matter WHERE the bkg min IS! Just that we find it somewhere!
        # also -- don't need to subpixel shift - it's just NOISE.
        ba = []
        count = 0
        for idx1 in range(bkg_img.shape[0]):
            for idx2 in range(bkg_img.shape[1]):
                #print idx-bkg_img.shape[0], idx
                shifted = bkg_img.take(range(idx1-bkg_img.shape[0], idx1), 
                                       mode='wrap', axis=0) \
                                 .take(range(idx2-bkg_img.shape[1], idx2), 
                                       mode='wrap', axis=1)
                rotated = np.rot90(shifted, 2) 
                #residual = np.abs(shifted - rotated)
                #numerator = aperture_photometry(residual, aperture)
                #ba.append(float(numerator['aperture_sum']))
                #bleh = np.sum(np.abs(shifted-rotated))
                ba.append(np.sum(np.abs(shifted-rotated)))

        # find the  minimum of all possible bkg asyms
        bkgasym = np.min(ba)*aperture.area()/(bkg_img.shape[0]**2)

        asyms = defaultdict(list)
        prior_points = []

        #  minimize the galaxy asymmetry
        while True:
            ga = []
            dd = []
            deltas, points = utils2.generate_deltas(imgcenter_sm, .3, delta)

            for d, p in zip(deltas, points): 
                # if the point already exists in the dictionary, 
                #don't run asym codes!
                if p not in asyms: 
                    newdata = sp_interp.shift(smaller, d)
                    rotdata = sp_interp.rotate(newdata, 180.)
                    residual = np.abs(newdata-rotdata)
                    numerator = aperture_photometry(residual, aperture)
                    denominator = aperture_photometry(np.abs(newdata), aperture)
                    num = float(numerator['aperture_sum']) 
                    den = float(denominator['aperture_sum'])
                    galasym = num/den

                    # create an array of asyms ... 
                    ga.append(galasym)
                    dd.append(den)
                    # ... and a dictionary that maps each asym to a 
                    #point on the image grid
                    asyms[p].append([galasym,den])

                # just take the value that's already in the dictionary 
                # for that point
                else:
                    ga.append(asyms[p][0][0])
                    dd.append(asyms[p][0][1])

            # if the asymmetry found at the original center 
            # (first delta in deltas) is the minimum, we're done!
            if ga[0] == np.min(ga):
                center = imgcenter - deltas[0]
                center_sm = imgcenter_sm - deltas[0]
                self.asym_plot(newdata, residual, aperture, 
                               imgcenter_sm, center_sm)
                return ga[0]-bkgasym/dd[0], center
            else:
                minloc = np.where(ga == np.min(ga))[0]
                delta = deltas[minloc[0]]
                prior_points = list(points)
                       

    def asym_plot(self, shifted, residual, aperture, imgcenter, galcenter):
        #low = imgcenter[0]-2*self.rpet
        #high = imgcenter[1]+2*self.rpet
        plt.figure()
        plt.imshow(shifted)
        plt.plot(galcenter[1], galcenter[0], 'k+', mew=2, ms=10)
        plt.plot(imgcenter[1], imgcenter[0], 'r+', mew=2, ms=10)
        aperture.plot()
        #plt.xlim(low, high)
        #plt.ylim(low, high)
        plt.title('Shifted Cleaned Image')
        plt.savefig('output/asyms/'+self.name+'_asym1_CS.png')
        plt.close()

        plt.figure()
        plt.imshow(residual)
        plt.plot(galcenter[1], galcenter[0], 'k+', mew=2, ms=10, 
                 label="Asymmetry Center")
        plt.plot(imgcenter[1], imgcenter[0], 'r+', mew=2, ms=10, 
                 label="Image Center")
        aperture.plot()
        #plt.xlim(low, high)
        #plt.ylim(low, high)
        plt.title('Asymmetry Residuals (I - I180)')
        plt.legend()
        plt.savefig('output/asyms/'+self.name+'_asym2_CS.png')
        plt.close()
        return 0.

    def get_concentration(self, clean_dat):#aper_sums, aper_radii, rpet
        '''
        To calculate the conctration we need to find the radius 
            -- which encloses 20% of the total light
            -- which encloses 80% of the total light
        So we need a running sum of the total pixel counts within various radii
        We also need to know the total flux -- 
           define as the total pixel counts within 
           an aperture of  one petrosian radius
        Divide the running sum by the fixed total flux value and 
           see where this ratio crosses .2 and .8
        #'''
        print "calculating Concentration..."

        # Calculate a running sum in CIRCULAR apertures instead of ELLIPTICAL
        position = [self.x, self.y]
        imgsize = clean_dat.shape[0]/2.
        radii = 10*np.logspace(-1.0, np.log10(imgsize/10.), num=20)

        flux = []
        for idx in range(1,len(radii)):          
            # for Avg SB, annuli from radius "before" to current radius in apix
            an = CircularAnnulus(position, radii[idx-1], radii[idx])
            #ans.append(an) 
            flux.append(aperture_photometry(clean_dat, an, method='exact'))

        phot_table = np.hstack(flux)
        counts = phot_table['aperture_sum']
        #annuli = np.hstack(ans)

        # calculate the average surface brightness
        cum_sum = [np.sum(counts[0:idx]) for idx in range(1,len(counts)+1)]
        cum_sum = np.array(cum_sum)      

        # Calculate the total flux in a Circular aperture of 1.5*rpet
        tot_ap = CircularAperture(position, self.rpet)
        tot_flux = aperture_photometry(clean_dat, tot_ap, method='exact')
        totflux = float(tot_flux['aperture_sum'])

        ratio = cum_sum/totflux
        x = radii[1:]

        #pdb.set_trace()

        plt.figure()
        plt.plot(x, ratio, 'ro')
        plt.axhline(y=0.2, linestyle='--', color='k')
        plt.axhline(y=0.8, linestyle='--', color='k')
        plt.ylim(0., 1.4)
        plt.title('Cumulative Flux / Total Flux vs. Radius')
        plt.xlabel('Radius [pixel]')
        plt.ylabel('Flux ratio [counts]')
        #plt.show()
        #pdb.set_trace()

        plt.close()

        # now we need to find the intersection of ratio with 0.2 and 0.8
        r = np.linspace(np.min(x), np.max(x), num=1000)
        f = interp1d(x, ratio, kind='cubic')        
        ratios = f(r)

        plt.plot(r, ratios)
        plt.savefig('output/conc/'+self.name+'_c_conc.png')

        r20 = np.nan
        r80 = np.nan

        if not np.any(np.isnan(ratios)):
            dirty20 = 0.2 - ratios
            # if any value in dirty is non negative, shit be rockin ... 
            if (dirty20 < 0).any():
                for idx,thing in enumerate(dirty20):
                    if thing < 0:
                        r20 = np.mean([r[idx-1],r[idx]])
                        break

            dirty80 = 0.8 - ratios
            # if any value in dirty is non negative, shit be rockin ... 
            if (dirty80 < 0).any():
                for idx,thing in enumerate(dirty80):
                    if thing < 0:
                        r80 = np.mean([r[idx-1],r[idx]])
                        break
        #pdb.set_trace()
        return  5*np.log10(r80/r20)
        
    def get_gini(self, clean_dat):
        '''
        Need all pixels associated with a galaxy -- use my aperture thing? 
        1. All pixels within 1 Rp
        2. All pixels above the mean SB at 1 Rp
        Gonna do 1. for  now but want to try 2. as well
        '''
        print "calculating Gini..."

        gflag = 0

        # create aperture at center of galaxy
        ap = utils2.EllipticalAperture( (self.x, self.y), self.rpet, \
                                            self.rpet/self.e, self.theta, \
                                            clean_dat)

        # create galaxy and background masks from aperture
        #apmask = ap.aper.astype('float')
        #bkgmask = np.logical_not(apmask).astype('float')
        
        pixels = ap.aper * clean_dat
        galpix = pixels[pixels > 0.]
        # sort galpixels
        galpix_sorted = sorted(galpix)
        # calculate the mean
        xbar = np.mean(galpix_sorted)
        n = len(galpix_sorted)

        #calculate G
        gsum = [2*i-n-1 for i, p in enumerate(galpix_sorted)]
        g = 1/(xbar*n*(n-1))*np.dot(gsum, galpix_sorted)

        #pdb.set_trace()

        plt.figure()
        imgplot = plt.imshow(pixels[200:300, 200:300], cmap='gray_r', origin='lower')
        imgplot.set_clim(-0.009, 0.022)
        plt.savefig('output/aperture/'+self.name+'_aperG.png')
        plt.close()

        #pdb.set_trace()

        if g > 1:
            gflag = 1

        return g, gflag
            
    def get_m20(self, clean_dat):

        print "Calculating M20..."
        
        shape = clean_dat.shape
        imgcenter = np.array([shape[0]/2., shape[1]/2.])
        galcenter = np.array([self.x, self.y])
        delta = imgcenter - galcenter

        # create an array that has the pixel elements as the array values
        x, y = np.ogrid[:clean_dat.shape[0], :clean_dat.shape[1]]
        dist = (x-imgcenter[0])*(x-imgcenter[0]) + (y-imgcenter[1])*(y-imgcenter[1])


        # create aperture at center of galaxy (mask)
        ap = utils2.EllipticalAperture(imgcenter, self.rpet, \
                                            self.rpet/self.e, self.theta, \
                                            clean_dat)

        mtots = defaultdict(list)
        prior_points = []
       
        #  minimize the galaxy asymmetry
        while True:
            mm = []
            deltas, points = utils2.generate_deltas(imgcenter, .3, delta)

            for d, p in zip(deltas, points): 
                # if the point already exists in the dictionary... 
                if p not in mtots: 
                    # shift the current galcenter to the center of the image
                    shifted = sp_interp.shift(clean_dat, d)

                    # apply the mask
                    galpix = ap.aper*shifted

                    # calculate mtot for this location
                    mtot = np.sum(galpix*dist)

                    # create an array of mtots ... 
                    mm.append(mtot)

                    # ... and a dictionary that maps each mtot to a 
                    # point on the image grid
                    mtots[p].append(mtot)

                # ...otherwise just take the value that's already in the dictionary 
                else:
                    mm.append(mtots[p][0])

            # if the mtot found at the original center 
            # (first delta in deltas) is the minimum, we're done!
            if mm[0] == np.min(mm):
                center = imgcenter - deltas[0]
                Mtot = mm[0] 
                break
            else:
                minloc = np.where(mm == np.min(mm))[0]
                delta = deltas[minloc[0]]
                prior_points = list(points)

        # Once we have a minimized mtot and the galcenter that minimizes it
        # we can calculate M20!

        # create aperture/dist map at center that minimizes Mtot
        m20_ap = utils2.EllipticalAperture(center, self.rpet, \
                                            self.rpet/self.e, self.theta, \
                                            clean_dat)
        m20dist = (x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1])

        # mask out the pixels belonging to the galaxy based on this aperture
        pixels = m20_ap.aper * clean_dat
        galpix = pixels[pixels > 0.]
        m20dist_pix = m20dist[pixels > 0.]

        # sort m20dist_pix based on the rank order of galpix
        m20dist_sorted = [x for y, x, in sorted(zip(galpix, m20dist_pix), reverse=True)]
        galpix_sorted = sorted(galpix, reverse=True)

        # compute 0.2 * ftot
        ftot20 = 0.2*np.sum(galpix_sorted)

        # compute the cumulative sum of galpix_sorted -- 
        fcumsum = np.cumsum(galpix_sorted)

        # find where fcumsum is less than ftot20 -- 
        # the elements in this array are those which will be used to compute M20
        m20_pixels = np.where(fcumsum < ftot20)
        
        m20_galpix = np.array(galpix_sorted[0:np.max(m20_pixels[0])])
        m20_distpix = np.array(m20dist_sorted[0:np.max(m20_pixels[0])])

        M20 = np.log10(np.sum(m20_galpix*m20_distpix)/Mtot)

        fig, ax = plt.subplots()
        
        imgplot = ax.imshow(clean_dat, cmap='gray_r')
        imgplot.set_clim(-0.01, 0.03)

        m20_ap2 = EllipticalAperture(center, self.rpet, \
                                            self.rpet/self.e, self.theta)
        contours = measure.find_contours(pixels,m20_galpix.min())
        for n, contour in enumerate(contours):
            ax.plot(contour[:,1], contour[:,0], linewidth=2)
        m20_ap2.plot()
        plt.plot(center[0], center[1], 'r+', mew=2, ms=10)
        ax.set_xlim(shape[0]/2.-3*self.rpet, shape[0]/2.+3*self.rpet)
        ax.set_ylim(shape[1]/2.-3*self.rpet, shape[1]/2.+3*self.rpet)
        plt.savefig('output/m20/'+self.name+'_m20.pdf')
        #plt.show()
        plt.close()
        
        return M20, center



####################### main ############################

def main():
    
    parser = argparse.ArgumentParser(description='Perform LLE/PCA/whatevs')
    parser.add_argument('directory', type=str, 
        help='Directory of fits images on which to run LLE.')
    parser.add_argument('output', type=str,
        help='Specify the desired name for output catalog.')
    args = parser.parse_args()
    
    fitsfiles = sorted(glob.glob(args.directory+'*.fits'))
    #fitsfiles = sorted(glob.glob(args.directory))

    outdir = 'output/datacube3/'

    galaxies = []
    t = Table(names=('name', 'Fidx', 'Fdist', 'Bdist', 
                     'F-B', 'Farea', 'Barea', 'Flag', 'aFlag'), 
              dtype=('S70', 'i', 'f4', 'f4', 'f4', 'f4', 'f4', 'i', 'i'))
    for f in fitsfiles: 
        basename = os.path.basename(f)
        filename = outdir+'f_'+basename
        if not os.path.isfile(filename):
            print "File not found! Running SExtractor before proceeding."
            print "Cleaning ", os.path.basename(f)
            fidx,fdist,bdist,fbdist,farea,barea,flag,aflag = utils2.clean_frame(f, outdir)
            t.add_row((basename, fidx, fdist, bdist, fbdist, farea, barea, flag, aflag))
        #else:
        # run everything else
        #print "Running", os.path.basename(f)
        #hdulist = fits.open(filename, memmap=True)
        #galaxies.append(Galaxy(hdulist,filename)) 
        #hdulist.close()


    t.write('data.txt', format='ascii.fixed_width', delimiter='')
    #info = Table(rows=[g.__dict__ for g in galaxies])
    #info.write(args.output, overwrite=True)

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
