
import re, glob, os, string, pdb
import argparse, warnings
import math, bisect
from math import pi, ceil
from collections import defaultdict
from random import gauss

import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt

import scipy.ndimage.interpolation as sp_interp
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry, EllipticalAnnulus, \
                              EllipticalAperture, CircularAnnulus, \
                              CircularAperture
import morph

class Morphology(object):

    def __init__(self, hdulist, filename, flags, outdir):

        # initialize fits image & catalog data
        try:
            image = hdulist['UCLN'].data
            catinfo = hdulist['UCLN'].header['SECATIDX']
        except:
            image = hdulist['CLN'].data
            catinfo = hdulist['CLN'].header['SECATIDX']

        cat = hdulist['CAT'].data[catinfo]
        segmap = hdulist['FSEG'].data

        self.xc, self.yc = image.shape[0]/2.,image.shape[1]/2.
        
        # FLAGS & NAMING ATTRIBUTES
        self.cat = flags[0]
        self.oflag, self.uflag, self.bflag = flags[1], flags[2], flags[3]
        self.name = os.path.basename(os.path.splitext(filename)[0])

        #pdb.set_trace()
        # The following line is only for SDSS cutouts specificy as I used 
        # the DR7 OBJID as the image filename and wanted to preserve the 
        # objid in the catalog
        #self.objid = np.int64(os.path.splitext(\
        #            os.path.basename(filename))[0].split('_')[1])
        self._outdir = outdir


        # SEXTRACTOR ATTRIBUTES        
        self.e = cat['ELONGATION']
        self.x, self.y = cat['X_IMAGE'], cat['Y_IMAGE']
        self.kron = cat['KRON_RADIUS']
        self.a, self.b = cat['A_IMAGE'], cat['B_IMAGE']
        self.theta = cat['THETA_IMAGE']*pi/180. # in radians?
        self.ra, self.dec = cat['ALPHA_J2000'], cat['DELTA_J2000']
        self.elipt = cat['ELLIPTICITY'] 

        # BACKGROUND VALUES
        self.med, self.rms = self.background(image, segmap)

        # PETROSIAN RADIUS & FRIENDS
        self.Rp, self.Rp_SB, self.Rpflag = self.get_petro_ell(image)
        self.Rp_c, self.Rp_SB_c, self.Rpflag_c = self.get_petro_circ(image)
        #self.Rp_c2, self.Rpflag_c2 = self.get_petro_circ2(image)
        #'''
        if not np.isnan(self.Rp) and not np.isnan(self.Rp_c):

            # CREATE SOME APERTURES IN WHICH TO MEASURE LIGHT DISTRIBUTIONS
            # get_asym requires apertures centered on image center
            ell_ap = EllipticalAperture((self.xc, self.yc), self.Rp, 
                                          self.Rp/self.e, self.theta)

            circ_ap = CircularAperture((self.xc, self.yc), self.Rp_c)

            # get_gini requires apertures centered on galaxy center
            gell_ap = morph.MyEllipticalAperture((self.x, self.y), self.Rp, 
                                            self.Rp/self.e, self.theta, image)

            gcirc_ap = morph.MyCircularAperture((self.x,self.y), self.Rp_c, image)

            self.stn = self.get_stn(gell_ap.aper*image)

            # MEASURE MORPHOLOGIES----------------------------------------
            self.A, self.Ax, self.Ay = self.get_asymmetry(image, ell_ap)
            self.A_c, self.Ax_c, self.Ay_c = self.get_asymmetry(image, circ_ap)

            self.r20, self.r50, self.r80, self.C = \
                self.get_concentration_ell(image)
            self.r20_c, self.r50_c, self.r80_c, self.C_c = \
                self.get_concentration_circ(image)

            [self.G, self.G_c] = self.get_gini1(image, [gell_ap, gcirc_ap])
            [self.G2, self.G2_c]= self.get_gini2(image)

            [self.M20, self.M20_c], self.Mx, self.My = \
                                    self.get_m20(image, gell_ap, gcirc_ap)

        else:
            print "Petrosian radius could not be calculated!!"
            pdb.set_trace()
            self.A, self.Ax, self.Ay = np.nan, self.x, self.y
            self.A_c, self.Ax_c, self.Ay_c = np.nan, self.x, self.y

            self.r20 = self.r50 = self.r80 = self.C = np.nan
            self.r20_c = self.r50_c = self.r80_c = self.C_c = np.nan
            self.G, self.G2 = np.nan, np.nan 
            self.G_c, self.G2_c = np.nan, np.nan 
            self.M20 = self.M20_c = np.nan 
            self.Mx, self.My = self.x, self.y
        #'''
                
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
        return median, std

    def get_stn(self, mask):
        pix = np.where(mask.flatten() != 0.0)
        galpixels = mask.flatten()[pix]
        n = len(galpixels)
        return np.sum(galpixels/np.sqrt(self.rms**2+abs(galpixels)))
        
    def get_petro_ell(self, image):
        r_flag = 0
    
        #'''
        # condition of np.log10(imgsize/constant) ensures that the maximum
        # radius will never exceed the size of the image
        a = 10*np.logspace(-1.0, np.log10(np.min([self.xc,self.yc])/10.),num=20)
        b = a/self.e
        #'''
		
        position = [self.x, self.y]
	
        annuli = np.hstack([EllipticalAnnulus(position, a[idx], a[idx+1], 
                                              b[idx], self.theta) \
                            for idx, radius in enumerate(a[:-1])])

        counts = np.hstack([aperture_photometry(image, an, method='exact')\
                            for an in annuli])['aperture_sum']
        
        areas = [an.area() for an in annuli]
        
        sb_counts = np.array([c+counts[i+1] for i, c in \
                              enumerate(counts[:-1])])
        avgsb_counts = np.cumsum(counts)[:-1]
        
        sb_areas = np.array([ar+areas[i+1] for i, ar in \
                             enumerate(areas[:-1])])
        
        avgsb_areas = np.cumsum(areas)[:-1]
        
        # Local SBs averaged over an annulus at r (around r in log space)
        sb = self._sb = sb_counts/sb_areas
        
        # Mean SBs within r
        avgsb = self._avgsb = avgsb_counts/avgsb_areas
        
        # Petrosian Ratio -- find rp at which this ratio = 0.2
        self._ratio = sb/avgsb
        
        # test for monotonicity of sb/<sb> beyond self.a to determine
        # if contaminated by nearly uncleaned source
        tail = np.where(a[1:-1] >= .8*self.a)
        dx = np.diff(self._ratio[tail])
        loc = bisect.bisect(dx, 0.)
        
        if loc: 
            fitloc = tail[0][loc]
            fit = np.polyfit(a[fitloc+1:-1], sb[fitloc::], deg=0)
            
            if fit > 0.:
                newsb = np.r_[sb[:fitloc], sb[fitloc::]-fit[0]]
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
        radii, ratios = morph.get_interp(self._rads, self._sb/self._avgsb)
        self._interprads, self._interpvals = radii, ratios
        
        if not np.any(np.isnan(ratios)):
            rp = morph.get_intersect(ratios, 0.2, radii, mono='dec')
            
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


    def get_petro_circ2(self, image):
        rflag = 0
        # Trying to match Rp from SDSS
        # minimize (ratio - 0.2)? Need an initial guess? 

        r0 = self.kron  #initial guess in pixels
        #r0 = 50
        eta = 0.2
        epsilon = 0.0001
        condition= True
        count = 1
        position = [self.x, self.y]

        while condition: 
            annulus = CircularAnnulus(position, 0.8*r0, 1.25*r0)
            
            an_counts = aperture_photometry(image, annulus, 
                                            method='exact')['aperture_sum']
            an_area = annulus.area()
            
            aperture = CircularAperture(position, r0)
            ap_counts = aperture_photometry(image, aperture, 
                                            method='exact')['aperture_sum']
            ap_area = aperture.area()
            
            ratio = (an_counts/an_area)/(ap_counts/ap_area)
            diff = (ratio-eta)[0]
            #print "ratio:", ratio[0]
            #print "diff:", diff
            
            #pdb.set_trace()
            if np.abs(ratio-eta) < epsilon:
                condition = False
                break

            print 'Before:', r0

            if diff < 0.:
                r0 = .95*r0
                print 'After:', r0
            else:
                r0 = 1.05*r0
                print 'After:', r0

            count+=1
            if count == 70:
                pdb.set_trace()

        if r0 > image.shape[0]/2.:
            rflag = 1
            #condition = False
        
        if r0 > 2*self.Rp_c1:
            pdb.set_trace()

        #print count
        return r0, rflag

    def get_petro_circ(self, image):
        r_flag = 0

        # condition of np.log10(imgsize/constant) ensures that the maximum
        # radius will never exceed the size of the image
        a = 10*np.logspace(-1.0, np.log10(np.min([self.xc,self.yc])/10.),num=20)
        b = a/self.e
        position = [self.x, self.y]
        
        annuli = np.hstack([CircularAnnulus(position, a[idx], a[idx+1]) \
                            for idx, radius in enumerate(a[:-1])])
        
        counts = np.hstack([aperture_photometry(image, an, method='exact')\
                            for an in annuli])['aperture_sum']
        
        areas = [an.area() for an in annuli]
        
        sb_counts = np.array([c+counts[i+1] for i, c in \
                              enumerate(counts[:-1])])
        avgsb_counts = np.cumsum(counts)[:-1]
        
        sb_areas = np.array([ar+areas[i+1] for i, ar in \
                             enumerate(areas[:-1])])
        
        avgsb_areas = np.cumsum(areas)[:-1]
        
        # Local SBs averaged over an annulus at r (around r in log space)
        sb = sb_counts/sb_areas
        
        # Mean SBs within r
        avgsb = avgsb_counts/avgsb_areas
        
        # Petrosian Ratio -- find rp at which this ratio = 0.2
        ratio = sb/avgsb
        
        # test for monotonicity of sb/<sb> beyond self.a to determine
        # if contaminated by nearly uncleaned source
        tail = np.where(a[1:-1] >= .8*self.a)
        dx = np.diff(ratio[tail])
        loc = bisect.bisect(dx, 0.)
        
        if loc: 
            fitloc = tail[0][loc]
            fit = np.polyfit(a[fitloc+1:-1], sb[fitloc::], deg=0)
            
            if fit > 0.:
                newsb = np.r_[sb[:fitloc],sb[fitloc::]-fit[0]]
                subtract = sb_counts-newsb*sb_areas
                newavgsb = (avgsb_counts - subtract)/avgsb_areas
                sb = newsb
                avgsb = newavgsb
                newratio = newsb/newavgsb

        # now we need to find the intersection of sb/<sb> with 0.2:
        # define a finer spacing of radii to interpolate onto
        rads = a[1:-1]
        radii, ratios = morph.get_interp(rads, sb/avgsb)
        interprads, interpvals = radii, ratios
        
        if not np.any(np.isnan(ratios)):
            rp = morph.get_intersect(ratios, 0.2, radii, mono='dec')
            
            if (rp > 0): 
                # Determine Surface Brightness at 1 Rp
                newsb = interp1d(rads, sb)
                rp_sb = newsb(rp)
                if rp_sb < 0:
                    rp_sb = np.nan
            else:
                rp_sb, r_flag = np.nan, 1
            return rp, rp_sb, r_flag

        else:
            rp, rp_sb, r_flag = np.nan, np.nan, 2
            return rp, rp_sb, r_flag


    def bkg_asymmetry(self, aperture, plot=False):

        # create a square background image approx same size 
        # as area of aperture (we need to minimize calculations)
        size = int(ceil(np.sqrt(ceil(aperture.area()))))
        bkg_img = np.zeros((size, size))
        mask = np.where(bkg_img == 0)

        for pixel in zip(mask[0], mask[1]):
            bkg_img[pixel] = gauss(self.med, self.rms)
        
        # save the background image 
        if plot:
            bkg = fits.ImageHDU(data=bkg_img)
            morph.checkdir(self._outdir+'asymimgs/')
            bkg.writeto(self._outdir+'asymimgs/'+self.name+'_bkg.fits', 
                        clobber=True, output_verify='silentfix')          


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
        
    def get_asymmetry(self, image, aper, save_residual=False):

        '''
        1. make a smaller image of the galaxy -> 2*petrosian rad
        2. create a background image
        3. create an aperture 1.5*petrosian radius
        4. minimize asymmetry in the background img
        5. minimize asymmetry in the galaxy img
        #'''

        print "calculating Asymmetry..."
       
        delta = np.array([self.x-self.xc, self.y-self.yc])
       
        bkg_asym = self.bkg_asymmetry(aper)

        asyms = defaultdict(list)
        prior_points = []

        while True:
            # These hold intermediary asym & denominator values
            ga, dd = [], []

            deltas, points = morph.generate_deltas([self.xc, self.yc], .3,delta)

            for d, p in zip(deltas, points):
 
                # if the point already exists in the dictionary, 
                #don't run asym codes!
                if p not in asyms: 
                    shifted = sp_interp.shift(image, d)
                    rotated = sp_interp.rotate(shifted, 180.)
                    residual = shifted - rotated

                    numerator = aperture_photometry(np.abs(residual), aper)
                    denominator = aperture_photometry(np.abs(shifted), aper)
                    num = float(numerator['aperture_sum']) 
                    den = float(denominator['aperture_sum'])
                    galasym = num/den
                    
                    ga.append(galasym)
                    dd.append(den)

                    # dict maps each asym to a point (p) on the image grid
                    asyms[p].append([galasym, den])

                # just take the value that's already in the dictionary 
                # for that point
                else:
                    ga.append(asyms[p][0][0])
                    dd.append(asyms[p][0][1])

            # If first value in asym is the minimum, we're done!
            if ga[0] == np.min(ga):
                
                asym_center = [self.xc, self.yc] + deltas[0]
                
                if save_residual:
                    # save the corresponding residual image
                    shift = sp_interp.shift(image, deltas[0])
                    rot = sp_interp.rotate(shift, 180.)
                    resid = shift - rot
                    res = fits.ImageHDU(data=resid)
                    morph.checkdir(self._outdir+'asymimgs/')
                    res.writeto(self._outdir+'asymimgs/'+self.name+'_res.fits', 
                                clobber=True, output_verify='silentfix')


                return ga[0]-bkg_asym/dd[0], asym_center[0], asym_center[1]

            else:
                minloc = np.where(ga == np.min(ga))[0]
                delta = deltas[minloc[0]]
                prior_points = list(points)


    def get_concentration_ell(self, image):

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

        a = 10*np.logspace(-1.0, np.log10(np.min([self.xc,self.yc])/10.),num=20)
        b = a/self.e
        position = [self.Ax, self.Ay]

        annuli = np.hstack([EllipticalAnnulus(position, a[idx], a[idx+1],
                                              b[idx], self.theta) \
                                for idx, radius in enumerate(a[:-1])])
        counts = np.hstack([aperture_photometry(image, an, method='exact') \
                            for an in annuli])['aperture_sum']
        cum_sum = np.cumsum(counts)[:-1]

        tot_aper = EllipticalAperture(position, 1.5*self.Rp, 
                                      1.5*self.Rp/self.e, self.theta)
        tot_flux = float(aperture_photometry(image, tot_aper, 
                                             method='center')['aperture_sum'])
        
        # ratio of the cumulative counts over the total counts in the galaxy
        ratio = cum_sum/tot_flux
        
        # now we need to find the intersection of ratio with 0.2 and 0.8
        interp_radii, interp_ratio = morph.get_interp(a[1:-1], ratio)
        
        if not np.any(np.isnan(interp_ratio)):
            r20 = morph.get_intersect(interp_ratio, 0.2, interp_radii)
            r50 = morph.get_intersect(interp_ratio, 0.5, interp_radii)
            r80 = morph.get_intersect(interp_ratio, 0.8, interp_radii)
        else:
            r20 = r50 = r80 = np.nan
            
        conc = 5*np.log10(np.divide(r80, r20))

        return r20, r50, r80, conc

    def get_concentration_circ(self, image):
        print "calculating Concentration..."

        radii = 10*np.logspace(-1.0, np.log10(np.min([self.xc, self.yc])/10.), 
                               num=20)

        # Build circular annuli centered on the ASYMMETRY CENTER of the galaxy
        annuli = np.hstack([CircularAnnulus((self.Ax_c, self.Ay_c), 
                                            radii[i-1], radii[i]) \
                            for i in range(1,len(radii))])
 
        counts = np.hstack([aperture_photometry(image, an, method='exact') \
                            for an in annuli])['aperture_sum']
        cum_sum = np.cumsum(counts)[:-1]

        # 2/24/16: Changed this to be 1.5*elliptical Rp (instead of circular)
        tot_aper = CircularAperture((self.Ax_c, self.Ay_c), 1.5*self.Rp)
        tot_flux = float(aperture_photometry(image, tot_aper, 
                                             method='center')['aperture_sum'])
        
        # ratio of the cumulative counts over the total counts in the galaxy
        ratio = cum_sum/tot_flux

        # now we need to find the intersection of ratio with 0.2 and 0.8
        interp_radii, interp_ratio = morph.get_interp(radii[1:-1], ratio)

        if not np.any(np.isnan(interp_ratio)):
            r20 = morph.get_intersect(interp_ratio, 0.2, interp_radii)
            r50 = morph.get_intersect(interp_ratio, 0.5, interp_radii)
            r80 = morph.get_intersect(interp_ratio, 0.8, interp_radii)
        else:
            r20 = r50 = r80 = np.nan
            
        conc = 5*np.log10(np.divide(r80, r20))

        return r20, r50, r80, conc
 
    def get_gini1(self, image, apertures):
        print "calculating Gini..."

        ginis = []
        for aper in apertures:
            mask = aper.aper * image

            galpix = mask[np.where(mask != 0.)].flatten()
            galpix_sorted = sorted(np.abs(galpix))
            xbar = np.mean(galpix_sorted)
            n = len(galpix_sorted)
            factor = 1/(xbar*n*(n-1))
            gsum = [2*i-n-1 for i in range(1,n+1)]
            
            ginis.append(factor*np.dot(gsum, galpix_sorted))
            
        return ginis

    def get_gini2(self, image):
        print "calculating Gini(2)..."
        ginis = []
        
        # Mask 2: galaxy pixels defined as those with flux >= SB at 1 petro rad
        # This method is based on Lotz 2004
        for rp, rp_sb in zip([self.Rp, self.Rp_c], [self.Rp_SB, self.Rp_SB_c]):
            
            morph.checkdir(self._outdir+'masks/')
            outname = self._outdir+'masks/'+self.name
            mask = morph.get_SB_Mask(rp, rp_sb, image, outname)*image

            if isinstance(mask, int):
                return np.nan
        
            galpix = mask[np.where(mask != 0.)].flatten()
            galpix_sorted = sorted(np.abs(galpix))
            xbar = np.mean(galpix_sorted)
            n = len(galpix_sorted)
            factor = 1/(xbar*n*(n-1))
            gsum = [2*i-n-1 for i in range(1,n+1)]
            ginis.append(factor*np.dot(gsum, galpix_sorted))

        return ginis
        
    def get_m20(self, image, ell_aper, circ_aper):

        print "Calculating M20..."

        # create .5*Rp 'box' centered on img center in which to calculate
        # Mtot at each pixel
        mxrange = [int(round(self.xc-0.5*self.Rp)),
                   int(round(self.xc+0.5*self.Rp))]

        myrange = [int(round(self.yc-0.5*self.Rp)), 
                   int(round(self.yc+0.5*self.Rp))]

        # create grid which overlaps entire image
        x2, y2 = np.ogrid[:2*self.xc, :2*self.yc]      

        # create aperture masks
        mask_ell = ell_aper.aper*image
        mask_circ = circ_aper.aper*image

        # create 2d array to store mtot values
        mtots_ell = np.zeros_like(image, dtype='float32')
        mtots_circ = np.zeros_like(image, dtype='float32')

        # calculate mtot at every pixel in our 'box'
        for i in range(mxrange[0], mxrange[1]):
            for j in range(myrange[0], myrange[1]):

                # for each pixel in the box, create a distance mapping --
                # distance of each pixel from current "center"
                dist_grid = (i - x2)**2 + (j - y2)**2

                #ell_ap = morph.EllipticalAperture((i,j), self.Rp,
                #                            self.Rp/self.e, self.theta, image)
                #mask_ell = ell_ap.aper*image
                
                #circ_ap = morph.CircularAperture((i, j), self.Rp_c, image)
                #mask_circ = circ_ap.aper*image

                # calculate Mtot
                try:
                    mtots_ell[i,j] = np.sum(mask_ell*dist_grid)
                    mtots_circ[i,j] = np.sum(mask_circ*dist_grid)
                except:
                    pdb.set_trace()

        M20s = []

        for idx, mtots in enumerate([mtots_ell, mtots_circ]):

            # set all the zeros to nans so that we can find the true min
            mtots[np.where(mtots == 0)] = np.nan
            
            # find the minimum mtot value
            Mtot = np.nanmin(mtots)
            
            # find the coordinates of that minimum
            xc, yc = np.where(mtots == Mtot)

            # re-create the distance grid corresponding to those coordinates
            grid = (xc - x2)**2 + (yc - y2)**2
            
            # re-create a 1*Rp aperture centered on those coordinates
            if idx == 0:
                m20_aper = morph.MyEllipticalAperture((xc, yc), self.Rp, 
                                            self.Rp/self.e, self.theta, image)
            else:
                m20_aper = morph.MyCircularAperture((xc, yc), self.Rp_c, image)

            # isolate the pixel flux within that aperture
            galpix = m20_aper.aper*image
            
            sortbyflux = sorted(zip(galpix.flatten(), grid.flatten()), 
                                reverse=True)
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
                
                M20s.append(M20)
                #return M20, xc[0], yc[0]

            # if NO pixels satisfy the above condition, set M to NAN
            else:
                self._Mlevel1 = np.nan
                M20s.append(np.nan)
                #return np.nan, self.x, self.y

        return M20s, xc[0], yc[0]


    def table(self, init=False): 
          
        the_dict = self.__dict__
        r = re.compile(r"_.+")
        matching_keys = filter(r.match, the_dict.keys())
        #pdb.set_trace()
        for key in matching_keys:
            del the_dict[key]

        if init:
            keys = [k for k in the_dict.keys()]
            dtypes = []
            for k in keys:
                if k in ['name']:
                    dtypes.append('S80')
                elif k in ['cat', 'oflag', 'uflag','Rpflag','Rpflag_c','bflag']:
                    dtypes.append('i')
                elif k in ['objid']:
                    dtypes.append('int64')
                else: 
                    dtypes.append('f')
            t = Table(names=keys, dtype=dtypes)
            return t, the_dict
        else:
            return the_dict

    def __exit__(self, type, value, traceback):
        self.stream.close()



####################### main ############################

def main():
    
    parser = argparse.ArgumentParser(description='Perform LLE/PCA/whatevs')
    parser.add_argument('-d', dest="directory", type=str, 
        help='Directory of fits images on which to run LLE.')
    parser.add_argument('-c', dest="catalog_name", type=str,
        help='Specify the desired name for output catalog.')
    parser.add_argument('--outdir', type=str, default='output/datacube/', 
        help='Specify the desired name for output directory.')
    args = parser.parse_args()

    # Select all FITS files in the given directory
    fitsfiles = np.array(sorted(glob.glob(args.directory+'*band_5*.fits')))
    #filename = "synthetic_image_0_band_5_camera_0_bg_0.fits"
    #fitsfiles = np.array([filename])
    #fitsfiles = sorted(glob.glob(args.directory))

    # There are a lot of useless warnings that pop up -- suppress them!
    warnings.filterwarnings('ignore', message='Overwriting existing file .*',
                            module='pyfits')

    #pdb.set_trace()
    # If the proper output directory doesn't exist, create it
    morph.checkdir(args.outdir+'datacube/')

    # The morphology catalog is built by going through each FITS image, 
    # cleaning it, processing it, measuring morphological parameters, and
    # making appropriate plots FOR EACH FITS IMAGE. The morphological 
    # parameters are appended to a huge table which is updated after each 
    # image. All the output (Cleaned FITS files, plots and figures) are 
    # stored in appropriately named directories which are created on the fly
    # and nestled within the specified subdirectory.
    #
    # Sometimes this shit crashes because I'm a terrible programmer.
    # In that case, I remove the offending FITS image, putting it in the 
    # "bad_cutouts" directory and then start up the code again.
    #
    # When the code is rerun, we don't want to start a new morph catalog - 
    # we want to keep appending to the old one so we test for this...


    # check to see if a current morphology catalog exists
    # if so -- find where cursor should go by getting the length
    try:
        t = Table.read(args.catalog_name)
        counter = len(t)
        # this "cursor" tells us which FITS file we should start at because, 
        # presumably, all previous images have already been processed.
        fitsfiles = fitsfiles[counter::]
    except:
        counter = 0

    tt = Table(names=('category', 'oFlag', 'uFlag', 'bFlag'))

    # Now that we have our list of FITS to process...
    for idx, f in enumerate(fitsfiles): 
        basename = os.path.basename(f)
        filename = args.outdir+'datacube/f_'+basename


        #if not os.path.isfile(filename):
        #print "File not found! Running SExtractor before proceeding."
        print "Cleaning ", os.path.basename(f)
        flags = morph.clean_frame(f, args.outdir+'datacube/', sep=4, 
                                  survey='SDSS')

        
        print flags
        tt.add_row((flags[0], flags[1], flags[2], flags[3]))
        """
        # check to see if SExtractor failed
        if np.any(np.array(flags)-9 < 0):
        
            print "Running", os.path.basename(f)
            hdulist = fits.open(filename, memmap=True)

            # Measure galaxy morphologies
            g = Morphology(hdulist, filename, flags, args.outdir)

            # Plot galaxy figures for quality control
            #if not np.isnan(g.Rp):
            #    galaxy_plot.plot(g, hdulist)


            # If this is the first galaxy; create a Table
            if (idx == 0) and (counter == 0):
                t, gal_dict = g.table(init=True)
                t.add_row(gal_dict)
            # Otherwise, append galaxy morphologies to existing Table
            else:
                t.add_row(g.table())
            
            # Close any remaining FITS files
            hdulist.close()
            
            # Write the file after each galaxy (costly but don't want to lose
            # data when my shit crashes, which it inevitably does
            t.write(args.catalog_name, overwrite=True)
            print counter+1," galaxies measured!"
            counter+=1
            del g
        else:
            print "SExtractor failed on "+basename
            #if (idx == 0) and (counter == 0):
            #t, gal_dict = g.table(init=True)
            #t.add_row(gal_dict)
     """

    tt.write('SEflags_subdir_000_band_5.fits')
    print "Morphological parameter catalog complete.\n"
    exit()  


if __name__ == '__main__':
    main()
    
