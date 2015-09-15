
import re
import glob
import argparse
import os
import string
import math
import pdb
import pyfits as fits
import bisect
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.interpolation as sp_interp
from math import pi, ceil
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from collections import defaultdict
from random import gauss
from photutils import aperture_photometry, EllipticalAnnulus, \
                              EllipticalAperture, CircularAnnulus, \
                              CircularAperture
import utils
import clean
import galaxy_plot
import warnings

class Galaxy(object):

    def __init__(self, hdulist, filename, flags):
        # initialize fits image & catalog data
        try:
            image = hdulist['UCLN'].data
            catinfo = hdulist['UCLN'].header['SECATIDX']
        except:
            image = hdulist['CLN'].data
            catinfo = hdulist['CLN'].header['SECATIDX']

        cat = hdulist['CAT'].data[catinfo]
        segmap = hdulist['FSEG'].data
        
        # flags and naming attributes
        self.cat = flags[0]
        self.oflag, self.uflag, self.bflag = flags[1], flags[2], flags[3]
        self.name = os.path.basename(os.path.splitext(filename)[0])
        # The last index will change when I'm running on the regular SDSS sample
        self.objid = np.int64(os.path.splitext(\
                    os.path.basename(filename))[0].split('_')[3])

        # SExtractor attributes        
        self.e = cat['ELONGATION']
        self.x, self.y = cat['X_IMAGE'], cat['Y_IMAGE']
        self.kron = cat['KRON_RADIUS']
        self.a, self.b = cat['A_IMAGE'], cat['B_IMAGE']
        self.theta = cat['THETA_IMAGE']*pi/180. # in radians?
        self.ra, self.dec = cat['ALPHA_J2000'], cat['DELTA_J2000']
            
        # background attributes
        self.med, self.rms = self.background(image, segmap)
        
        # morphological attributes
        self.Rp, self.Rp_SB, self.Rpflag = self.get_petro(image)
        self.elipt = cat['ELLIPTICITY'] 
        
        if not np.isnan(self.Rp):
            self.A, self.Ax, self.Ay = self.get_asymmetry(image)
            self.r20, self.r80, self.C = self.get_concentration(image)
            self.G = self.get_gini1(image)
            self.G2 = self.get_gini2(image)
            self.M20, self.Mx, self.My = self.get_m20_1(image)
        else:
            self.A, self.Ax, self.Ay = np.nan, self.x, self.y
            self.r20, self.r80, self.C = np.nan, np.nan, np.nan
            self.G1, self.G2 = np.nan, np.nan 
            self.M20, self.Mx, self.My = np.nan, self.x, self.y
                
    def __enter__(self):
        return self

    def settonan(self):
         self.cat = flags[0]
         self.oflag, self.uflag, self.bflag = flags[1], flags[2], flags[3]
         self.name = os.path.basename(os.path.splitext(filename)[0])
         self.e = self.x= self.y = self.kron = self.a = self.b = np.nan
         self.theta = self.elipt = self.ra = self.dec= np.nan
         self.med =  self.rms = np.nan
         self.Rp = self.Rp_SB = self.Rpflag = np.nan
         self.A = self.Ax = self.Ay = np.nan
         self.r20 = self.r80 = self.C = np.nan
         self.M20 = self.Mx = self.My = self.G = np.nan 

    def background(self, data, segmap):
        mean, median, std = sigma_clipped_stats(data[segmap==0])
        #median = np.median(data[segmap==0])
        #rms = np.sqrt(np.mean(np.square(data[segmap==0])))
        return median, std

    def get_stn(self, mask):
        pix = np.where(mask.flatten() != 0.0)
        galpixels = mask.flatten()[pix]
        n = len(galpixels)
        return np.sum(galpixels/np.sqrt(self.rms**2+abs(galpixels)))
        
  
    def get_petro(self, image):
        r_flag = 0
       
        # condition of np.log10(imgsize/constant) ensures that the maximum
        # radius will never exceed the size of the image
        a = 10*np.logspace(-1.0, np.log10(image.shape[0]/2./10.), num=20)
        b = a/self.e
        position = [self.x, self.y]

        #annuli = np.hstack([EllipticalAnnulus(position, a[idx], a[idx+1], 
        #                                      b[idx], self.theta) \
        #                    for idx, radius in enumerate(a[:-1])])

        # test a small sample with CircularAnnulus to see if we get close
        # to SDSS Rp
        annuli = np.hstack([CircularAnnulus(position, a[idx], a[idx+1]) \
                            for idx, radius in enumerate(a[:-1])])
        
        counts = np.hstack([aperture_photometry(image, an, method='exact')\
                            for an in annuli])['aperture_sum']
        areas = [an.area() for an in annuli]

        sb_counts = np.array([c+counts[i+1] for i, c in enumerate(counts[:-1])])
        avgsb_counts = np.cumsum(counts)[:-1]

        sb_areas = np.array([ar+areas[i+1] for i, ar in enumerate(areas[:-1])])
        avgsb_areas = np.cumsum(areas)[:-1]

        # Local SBs averaged over an annulus at r (around r in log space)
        sb = self._sb = sb_counts/sb_areas
   
        # Mean SBs within r
        avgsb = self._avgsb = avgsb_counts/avgsb_areas

        # Petrosian Ratio -- find rp at which this ratio = 0.2
        self._ratio = sb/avgsb
        
        # need to test whether sb continues to decrease or if it's contaminated
        # by nearby light from other sources that wasn't fully cleaned
        # To do this: test for monotonicity of sb/<sb> beyond self.a 
        # (as given by SExtractor) 
        
        sb_avgsb = sb/avgsb
        tail = np.where(a[1:-1] >= .8*self.a)
        dx = np.diff(sb_avgsb[tail])
        loc = bisect.bisect(dx, 0.)
        
        if loc: 
            fitloc = tail[0][loc]
            fit = np.polyfit(a[fitloc+1:-1], sb[fitloc::], deg=0)

            if fit > 0.:
                newsb = np.concatenate((sb[0:fitloc], 
                                        sb[fitloc::]-fit[0]), axis=1)
                subtract = sb_counts-newsb*sb_areas
                newavgsb = (avgsb_counts - subtract)/avgsb_areas
                self._sb = newsb
                self._avgsb = newavgsb
                self._newratio = newsb/newavgsb

        '''
        # estimate error for SB and <SB>
        sberr = self._sberr = np.sqrt(areas1*self.rms**2)/areas1
        avgsberr = self._avgsberr = np.sqrt(asum*self.rms**2)

        # estimate error for the ratio of SB/<SB>
        ratio_err = self._ratio_err = np.sqrt((sberr/sb)**2 + 
                                              (avgsberr/avgsb)**2)
        #'''
        
        # now we need to find the intersection of sb/<sb> with 0.2:
        # define a finer spacing of radii to interpolate onto
        self._rads = a[1:-1]
        radii, ratios = utils.get_interp(self._rads, self._sb/self._avgsb)
        self._interprads, self._interpvals = radii, ratios

        if not np.any(np.isnan(ratios)):
            rp = utils.get_intersect(ratios, 0.2, radii, mono='dec')
            
            # if the ratio of rp/a is huge, something went wrong
            # this is a last ditch effort to catch strays
            if (rp > 0): 
                # Determine Surface Brightness at 1 Rp
                newsb = interp1d(self._rads, sb)
                rp_sb = newsb(rp)
                if rp_sb < 0:
                    rp_sb = np.nan
            else:
                rp_sb, r_flag = np.nan, 1
            return rp, rp_sb, r_flag
        else:
            print "Petrosian interpolation failed!"
            rp, rp_sb, r_flag = np.nan, np.nan, 2
            return rp, rp_sb, r_flag
 

    def bkg_asymmetry(self, aperture):

        # create a square background image approx same size 
        # as area of aperture (we need to minimize calculations)
        size = ceil(np.sqrt(ceil(aperture.area()))) 
        bkg_img = np.zeros((size, size))
        mask = np.where(bkg_img == 0)

        for pixel in zip(mask[0], mask[1]):
            bkg_img[pixel] = gauss(self.med, self.rms)
        
        # save the background image 
        bkg = fits.ImageHDU(data=bkg_img)
        utils.checkdir('output/asymings/')
        bkg.writeto('output/asymimgs/'+self.name+'.fits', clobber=True, 
                        output_verify='silentfix')          

        # minimize the background asymmetry
        ba = []
        for idx1 in range(bkg_img.shape[0]):
            for idx2 in range(bkg_img.shape[1]):
                shifted = bkg_img.take(range(idx1-bkg_img.shape[0], idx1), 
                                       mode='wrap', axis=0) \
                                 .take(range(idx2-bkg_img.shape[1], idx2), 
                                       mode='wrap', axis=1)
                rotated = np.rot90(shifted, 2) 
                ba.append(np.sum(np.abs(shifted-rotated)))

        # find the  minimum of all possible bkg asyms, normalized to the exact
        # area of the original aperture
        bkgasym = np.min(ba)*aperture.area()/(bkg_img.shape[0]*bkg_img.shape[1])
        return bkgasym
        
    def get_asymmetry(self, image):
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
        imgcenter = np.array([image.shape[0]/2., image.shape[1]/2.])
        delta = galcenter - imgcenter
        aper = EllipticalAperture(imgcenter, self.Rp, self.Rp/self.e,self.theta)
        bkg_asym = self.bkg_asymmetry(aper)

        asyms = defaultdict(list)
        prior_points = []

        #  minimize the galaxy asymmetry
        while True:
            ga = []
            dd = []
            deltas, points = utils.generate_deltas(imgcenter, .3, delta)
            #pdb.set_trace()
            for d, p in zip(deltas, points): 
                # if the point already exists in the dictionary, 
                #don't run asym codes!
                if p not in asyms: 
                    newdata = sp_interp.shift(image, d)
                    rotdata = sp_interp.rotate(newdata, 180.)
                    residual = newdata-rotdata
                    numerator = aperture_photometry(np.abs(residual), aper)
                    denominator = aperture_photometry(np.abs(newdata), aper)
                    num = float(numerator['aperture_sum']) 
                    den = float(denominator['aperture_sum'])
                    galasym = num/den

                    # create an array of asyms. 
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
                asym_center = imgcenter + deltas[0]
                
                # save the corresponding residual image
                new = sp_interp.shift(image, deltas[0])
                rot = sp_interp.rotate(new, 180.)
                resid = new - rot
                res = fits.ImageHDU(data=resid)
                utils.checkdir('output/asymings/')
                res.writeto('output/asymimgs/'+self.name+'_res.fits', 
                            clobber=True, output_verify='silentfix')

                # return A and A center (row, col)
                return ga[0]-bkg_asym/dd[0], asym_center[0], asym_center[1]
            else:
                minloc = np.where(ga == np.min(ga))[0]
                delta = deltas[minloc[0]]
                prior_points = list(points)

    def get_asymmetry2(self, image):
        
        center = [image.shape[0]/2., image.shape[1]/2.]
        galcenter = np.array([self.x, self.y])

        
        mxrange = [galcenter[0]-0.5*round(self.Rp), 
                   galcenter[0]+0.5*round(self.Rp)]
        myrange = [galcenter[1]-0.5*round(self.Rp), 
                   galcenter[1]+0.5*round(self.Rp)]

        aper = EllipticalAperture(center, self.Rp, self.Rp/self.e,self.theta)
        bkg_asym = self.bkg_asymmetry(aper)

        # Create a grid over which to calculate A at each point
        # This grid is finer than the original image
        xx, yy = np.ogrid[mxrange[0]:mxrange[1]:1, myrange[0]:myrange[1]:1]

        #pdb.set_trace()
        asyms = []
        for x in xx:
            print x
            for y in yy.transpose():
                dist = [center[0]-x[0], center[1]-y[0]]
                newdata = sp_interp.shift(image, dist)
                rotdata = sp_interp.rotate(newdata, 180.)
                residual = newdata-rotdata
                numerator = aperture_photometry(np.abs(residual), aper)
                denominator = aperture_photometry(np.abs(newdata), aper)
                num = float(numerator['aperture_sum']) 
                den = float(denominator['aperture_sum'])
                galasym = num/den
                asyms.append(galasym)
        
        A = np.min(asyms)-bkg_asym/den
        aa = np.reshape(asyms, (len(xx), len(yy.transpose())))
        #pdb.set_trace()
        loc = np.where(aa == np.min(aa))
        Ax, Ay = xx[loc[0][0]], yy.transpose()[loc[1][0]]
        print A, (Ax, Ay)
        return A

    def get_concentration(self, image):
        '''
        To calculate the conctration we need to find the radius 
            -- which encloses 20% of the total light
            -- which encloses 80% of the total light
        So we need a running sum of the total pixel counts in increasing radii
        We also need to know the total flux -- 
           define as the total pixel counts within 
           an aperture of  one petrosian radius
        Divide the running sum by the fixed total flux value and 
           see where this ratio crosses .2 and .8
        #'''
        print "calculating Concentration..."

        radii = 10*np.logspace(-1.0, np.log10(image.shape[0]/2./10.), num=20)

        # Build circular annuli centered on the ASYMMETRY CENTER of the galaxy
        annuli = [CircularAnnulus((self.Ax, self.Ay),radii[i-1],radii[i])\
                  for i in range(1,len(radii))]
        counts = np.hstack([aperture_photometry(image, an, method='center') \
                            for an in annuli])['aperture_sum']
        cum_sum = np.array([np.sum(counts[0:idx]) \
                            for idx in range(1,len(counts)+1)])

        # Calculate the total flux in a Circular aperture of 1*rpet
        tot_aper = CircularAperture((self.Ax, self.Ay), self.Rp)
        tot_flux = float(aperture_photometry(image, tot_aper, 
                                             method='center')['aperture_sum'])

        # ratio of the cumulative counts over the total counts in the galaxy
        ratio = cum_sum/tot_flux
        rads = radii[1:]

        # now we need to find the intersection of ratio with 0.2 and 0.8
        interp_radii, interp_ratio = utils.get_interp(rads, ratio)

        if not np.any(np.isnan(interp_ratio)):
            r20 = utils.get_intersect(interp_ratio, 0.2, interp_radii)
            r80 = utils.get_intersect(interp_ratio, 0.8, interp_radii)
        else:
            print "Concentration interpolation failed."
            r20 = r80 = np.nan
        
        return r20, r80, 5*np.log10(r80/r20)
        
    def get_gini1(self, image):
        print "calculating Gini..."

        # Mask 1: galaxy pixels defined exactly within 1 petrosian radius
        ap = utils.EllipticalAperture((self.x, self.y), self.Rp, self.Rp/self.e,
                                      self.theta, image)
        mask1 = ap.aper * image
        galpix = mask1[np.where(mask1 != 0.)].flatten()
        galpix_sorted = sorted(np.abs(galpix))
        xbar = np.mean(galpix_sorted)
        n = len(galpix_sorted)
        factor = 1/(xbar*n*(n-1))

        gsum = []
        for i in range(1,n+1):
            gsum.append(2*i-n-1)
        
        gini = factor*np.dot(gsum, galpix_sorted)

        return gini 

    def get_gini2(self, image):
        print "calculating Gini(2)..."

        # Mask 2: galaxy pixels defined as those with flux >= SB at 1 petro rad
        # This method is based on Lotz 2004
        mask2 = utils.get_SB_Mask(self.Rp, self.Rp_SB, image, self.name)*image
        if isinstance(mask2, int):
            return np.nan
        
        galpix = mask2[np.where(mask2 != 0.)].flatten()
        galpix_sorted = sorted(np.abs(galpix))
        xbar = np.mean(galpix_sorted)
        n = len(galpix_sorted)
        factor = 1/(xbar*n*(n-1))

        gsum = []
        for i in range(1,n+1):
            gsum.append(2*i-n-1)
        
        gini = factor*np.dot(gsum, galpix_sorted)

        return gini
        
    def get_m20_nomin(self, image):
        galcenter = np.array([self.Ax, self.Ay])

        # create aperture at center of galaxy (mask)
        gal_aper = utils.EllipticalAperture(galcenter, self.Rp, self.Rp/self.e,
                                            self.theta, image)
        galpix = gal_aper.aper*image


        # Create a "distance grid" - each element's value is it's distance
        # from the center of the image for the entire image
        x2, y2 = np.ogrid[:image.shape[0], :image.shape[1]]      
        dist_grid = (galcenter[0] - x2)**2 + (galcenter[1] - y2)**2

        mtot = np.sum(galpix*dist_grid)

        grid_sorted = np.array([i for j,i in sorted(zip(galpix.flatten(),\
                                                        dist_grid.flatten()),\
                                                    reverse=True)])
        galpix_sorted = np.array(sorted(galpix.flatten(), reverse=True))
        ftot20 = 0.2*np.sum(galpix_sorted)
        fcumsum = np.cumsum(galpix_sorted)
        m20_pix = np.where(fcumsum < ftot20)[0]
        if len(m20_pix) != 0:
            m20_galpix = galpix_sorted[m20_pix]
            m20_distpix = grid_sorted[m20_pix]
            M20 = np.log10(np.sum(m20_galpix*m20_distpix)/mtot) 
            self._Mlevel1 = np.min(m20_galpix)
            #pdb.set_trace()
            return M20
        else:
            self._Mlevel1 = np.nan
            return np.nan
        
    def get_m20_1(self, image):

        print "Calculating M20..."
        
        center = [image.shape[0]/2., image.shape[1]/2.]
        galcenter = np.array([self.x, self.y])

        # create .5*Rp 'box' centered on img center in which to calculate
        # Mtot at each pixel
        mxrange = [int(round(center[0]-0.5*self.Rp)), 
                   int(round(center[0]+0.5*self.Rp))]

        myrange = [int(round(center[1]-0.5*self.Rp)), 
                   int(round(center[1]+0.5*self.Rp))]

        # create grid which overlaps entire image
        x2, y2 = np.ogrid[:image.shape[0], :image.shape[1]]      

        # create 1Rp aperture at center of galaxy (mask)
        gal_aper = utils.EllipticalAperture(center, self.Rp, self.Rp/self.e,
                                            self.theta, image)
        mask1 = gal_aper.aper*image

        self.stn = self.get_stn(mask1)

        # create a 2d array for mtot values
        mtots = np.zeros_like(image, dtype='float32')

        for i in range(mxrange[0], mxrange[1]):
            for j in range(myrange[0], myrange[1]):

                # for each pixel in the box, create a distance mapping --
                # distance of each pixel from current "center"
                dist_grid = (i - x2)**2 + (j - y2)**2
    
                # calculate Mtot
                try:
                    mtots[i,j] = np.sum(mask1*dist_grid)
                except:
                    pdb.set_trace()

        # set all the zeros to nans so that we can find the true min
        mtots[np.where(mtots == 0)] = np.nan

        # find the minimum mtot value
        Mtot = np.nanmin(mtots)

        # find the coordinates of that minimum
        xc, yc = np.where(mtots == Mtot)

        # re-create the distance grid corresponding to those coordinates
        grid = (xc - x2)**2 + (yc - y2)**2
        
        # re-create a 1*Rp aperture centered on those coordinates
        m20_aper = utils.EllipticalAperture((xc, yc), self.Rp, self.Rp/self.e, 
                                            self.theta, image)
        
        # isolate the pixel flux within that aperture
        galpix = m20_aper.aper*image
        
        sortbyflux = sorted(zip(galpix.flatten(), grid.flatten()), reverse=True)
        sortedbyflux = [list(thing) for thing in zip(*sortbyflux)]
        sorted_galpix = np.array(sortedbyflux[0])
        sorted_gridpix = np.array(sortedbyflux[1])

        # 20% of the total galaxy flux ( .2*f_tot )
        ftot20 = 0.2*np.sum(sorted_galpix)

        # cumulative sum of the galaxy flux ( sum of f_i )
        fcumsum = np.cumsum(sorted_galpix)

        # Determine location where the cumulative sum is less than .2*f_tot
        m20_pix = np.where(fcumsum < ftot20)[0]

        if len(m20_pix) != 0:
            
            m20_galpix = sorted_galpix[m20_pix]
            m20_distpix = sorted_gridpix[m20_pix]

            M20 = np.log10(np.sum(m20_galpix*m20_distpix)/Mtot)

            self._Mlevel1 = np.min(m20_galpix)
            return M20, xc[0], yc[0]

        # if NO pixels satisfy the above condition, set M to NAN
        else:
            self._Mlevel1 = np.nan
            return np.nan, self.x, self.y

    def get_m20_2(self, image):

        print "Calculating M20(2)...(SB mask)"
        
        center = [image.shape[0]/2., image.shape[1]/2.]
        galcenter = np.array([self.x, self.y])

        mxrange = [center[0]-round(self.a), center[0]+round(self.a)]
        myrange = [center[1]-round(self.a), center[1]+round(self.a)]

        x, y = np.ogrid[mxrange[0]:mxrange[1], myrange[0]:myrange[1]]

        x2, y2 = np.ogrid[:image.shape[0], :image.shape[1]]      
        dist_grid = (center[0] - x2)**2 + (center[1] - y2)**2

        mask2 = utils.get_SB_Mask(self.Rp, self.Rp_SB, image, self.name)*image
        if isinstance(mask2, int):
            self._Mlevel2 = np.nan
            return np.nan, xc[0], yc[0]
        
        self.stn = self.get_stn(mask2)
        
        mtots = []
        for i in x:
            for j in y.transpose():
                indices0 = range(center[0]-i, center[0]-i + 
                                 dist_grid.shape[0])
                indices1 = range(center[1]-j, center[1]-j + 
                                 dist_grid.shape[1])
                shift_grid = dist_grid.take(indices0, axis=0, mode='wrap')\
                                      .take(indices1, axis=1, mode='wrap')
                mtots.append(np.sum(mask2*shift_grid))

        Mtot = np.min(mtots)
        mtot = np.array(mtots).reshape([len(x), len(y.transpose())])
        xc = x[np.where(mtot == Mtot)[0]][0]
        yc = y.transpose()[np.where(mtot == Mtot)[1]][0]
        
        indices0 = range(center[0]-xc,center[0]-xc+dist_grid.shape[0])
        indices1 = range(center[0]-yc,center[0]-yc+dist_grid.shape[0])
        grid = dist_grid.take(indices0, axis=0, mode='wrap')\
                        .take(indices1, axis=1, mode='wrap')
        galpix = mask2
        grid_sorted = np.array([i for j,i in sorted(zip(galpix.flatten(),\
                                                        grid.flatten()),\
                                                    reverse=True)])
        galpix_sorted = np.array(sorted(galpix.flatten(), reverse=True))
        ftot20 = 0.2*np.sum(galpix_sorted)
        fcumsum = np.cumsum(galpix_sorted)
        m20_pix = np.where(fcumsum < ftot20)
        if len(m20_pix) != 0:
            m20_galpix = galpix_sorted[m20_pix]
            m20_distpix = grid_sorted[m20_pix]
            M20 = np.log10(np.sum(m20_galpix*m20_distpix)/Mtot)
            self._Mlevel2 = np.min(m20_galpix)
            return M20, xc[0], yc[0]
        else:
            self._Mlevel2 = np.nan
            return np.nan, xc[0], yc[0]

    def table(self, init=False): 
          
        the_dict = self.__dict__
        r = re.compile(r"_.+")
        matching_keys = filter(r.match, the_dict.keys())
        
        for key in matching_keys:
            del the_dict[key]

        if init:
            names = ['name', 'objid', 'ra', 'dec',
                     'e', 'x', 'y', 'a', 'b', 'theta', 'elipt', 'kron', 
                     'Rp', 'Rpflag', 'Rp_SB', 'r20', 'r80',  'A', 'G','G2',
                     'C', 'M20', 'Ax', 'Ay', 'Mx', 'My', 'stn',
                     'med', 'rms', 'cat', 'oflag', 'uflag', 'bflag']
            dtypes = []
            for n in names:
                if n in ['name']:
                    dtypes.append('S80')
                elif n in ['cat', 'oflag', 'uflag', 'Rpflag','bflag']:
                    dtypes.append('i')
                elif n in ['objid']:
                    dtypes.append('int64')
                else: 
                    dtypes.append('f')
            t = Table(names=names, dtype=dtypes)
            return t, the_dict
        else:
            return the_dict

    def __exit__(self, type, value, traceback):
        self.stream.close()



####################### main ############################

def main():
    
    parser = argparse.ArgumentParser(description='Perform LLE/PCA/whatevs')
    parser.add_argument('directory', type=str, 
        help='Directory of fits images on which to run LLE.')
    parser.add_argument('catalog_name', type=str,
        help='Specify the desired name for output catalog.')
    parser.add_argument('--outdir', type=str, default='output/datacube/', 
        help='Specify the desired name for output directory.')
    args = parser.parse_args()

    fitsfiles = np.array(sorted(glob.glob(args.directory+'*.fits')))
    #fitsfiles = sorted(glob.glob(args.directory))

    warnings.filterwarnings('ignore', message='Overwriting existing file .*',
                            module='pyfits')

    utils.checkdir(args.outdir)

    try:
        t = Table.read(args.catalog_name)
        counter = len(t)
        fitsfiles = fitsfiles[counter::]
    except:
        counter = 0

    for idx, f in enumerate(fitsfiles): 
        basename = os.path.basename(f)
        filename = args.outdir+'f_'+basename

        #if not os.path.isfile(filename):
        #print "File not found! Running SExtractor before proceeding."
        print "Cleaning ", os.path.basename(f)
        flags = clean.clean_frame(f, args.outdir, sep=4., survey='SDSS')
        
        # check to see if SExtractor failed
        if np.any(np.array(flags)-9 < 0):

            print "Running", os.path.basename(f)
            hdulist = fits.open(filename, memmap=True)

            g = Galaxy(hdulist, filename, flags)
            if not np.isnan(g.Rp):
                galaxy_plot.plot(g, hdulist)
            if (idx == 0) and (counter == 0):
                t, gal_dict = g.table(init=True)
                t.add_row(gal_dict)
            else:
                t.add_row(g.table())
                    
            hdulist.close()
            t.write(args.catalog_name, overwrite=True)
            print counter," galaxies measured!"
            counter+=1
            del g
        else:
            print "SExtractor failed on "+basename

    #t.write('bigsample_mycat_flags.dat', format='ascii.fixed_width')
    print "Morphological parameter catalog complete.\n"
    #Saving catalog to file..."
    #t.write(args.output, format='ascii')
    exit()  


if __name__ == '__main__':
    main()
    
