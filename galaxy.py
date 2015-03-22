
import re
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
        image = hdulist['CLN'].data
        segmap = hdulist['FSEG'].data
        cat = hdulist['CAT'].data[hdulist['CLN'].header['SECATIDX']]

        # flags and naming attributes
        self.cat = flags[0]
        self.oflag = flags[1]
        self.uflag = flags[2]
        self.name = os.path.basename(os.path.splitext(filename)[0])
  
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
        if self.Rp > 0.:
            self.A, self.Acx, self.Acy = self.get_asymmetry(image)
            self.r20, self.r80, self.C = self.get_concentration(image)
            self.G1, self.G2 = self.get_gini(image)
            self.M1, self.M2, self.Mcx1, \
            self.Mcy1, self.Mcx2, self.Mcy2 = self.get_m20(image)
        else:
            self.A, self.Acx, self.Acy = np.nan, self.x, self.y
            self.r20, self.r80, self.C = -1, -1, np.nan
            self.G1, self.G2 = np.nan, np.nan 
            self.M1, self.M2, self.Mcx1, self.Mcy1, self.Mcx2, \
                self.Mcy2 = np.nan,np.nan, self.x, self.y, self.x, self.y

        #dir(self) in cmd line to see all the hidden shits

    def __enter__(self):
        return self
      
    def background(self, data, segmap):
        mean, median, std = sigma_clipped_stats(data[segmap==0])
        #median = np.median(data[segmap==0])
        #rms = np.sqrt(np.mean(np.square(data[segmap==0])))
        return median, std
  
    def get_petro(self, image):
        r_flag = 0
       
        # condition of np.log10(imgsize/constant) ensures that the maximum
        # radius will never exceed the size of the image
        a = 10*np.logspace(-1.0, np.log10(image.shape[0]/2./10.), num=20)
        b = a/self.e
        position = [self.x, self.y]

        # determine flux (counts) in annuli at various radii
        annuli_atR = np.hstack([EllipticalAnnulus(position, a[idx-1], a[idx+1], 
                                              b[idx+1], self.theta) \
                                for idx in range(1,len(a[:-1]))])
        counts_atR = np.hstack([aperture_photometry(image, an,method='center')\
                                for an in annuli_atR])['aperture_sum']

        annuli_inR = np.hstack([EllipticalAnnulus(position, a[idx-1], a[idx], 
                                                  b[idx], self.theta) \
                                for idx in range(1,len(a[:-1]))])
        counts_inR = np.hstack([aperture_photometry(image, an,method='center')\
                                for an in annuli_inR])['aperture_sum']

        areas_atR = np.array([an.area() for an in annuli_atR])
        areas_inR = np.array([an.area() for an in annuli_inR])

        # surface brightness = total counts / area
        sb = self._sb = counts_atR/areas_atR   

        # average surface brightness is the cumulative sum of the counts within
        # the radius, R, divided by the cumulative area within in the same R 
        num = len(counts_inR)+1
        csum = np.array([np.sum(counts_inR[0:idx]) for idx in range(1, num)])
        asum = np.array([np.sum(areas_inR[0:idx]) for idx in range(1, num)])
        avgsb = self._avgsb = csum/asum  

        # need to test whether sb continues to decrease or if it's contaminated
        # by nearby light from other sources that wasn't fully cleaned
        # To do this: test for monotonicity of sb/avgsb beyond self.a (as given by 
        # SExtractor) 
        
        sb_avgsb = sb/avgsb
        tail = np.where(a[1:-1] >= self.a)
        test = sb_avgsb[tail]
        dx = np.diff(test)
        import bisect
        loc = bisect.bisect(dx, 0.)
        tofitfrom = tail[0][loc+1]

        fit = np.polyfit(a[tofitfrom+1:-1], sb[tofitfrom::], deg=0)
        p = np.poly1d(fit)
        
        newsb = np.concatenate((sb[0:tofitfrom], sb[tofitfrom::]-fit[0]), axis=1)

        plt.plot(a[1:-1], sb, 'k', a[1:-1], newsb, 'r')
        plt.axhline(0.)
        plt.show()
        #test_at = np.where(a[1:-1]>self.a)[0][test + 1]
        pdb.set_trace()

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
        radii, ratios = utils.get_interp(self._rads, sb/avgsb)
        self._interprads, self._interpvals = radii, ratios

        params=self.__dict__
        galaxy_plot.petro_SB2(params)
        
        pdb.set_trace()

        if not np.any(np.isnan(ratios)):
            rp = utils.get_intersect(ratios, 0.2, radii, mono='dec')
            if rp > -1:
                # Determine Surface Brightness at 1 Rp
                newsb = interp1d(self._rads, sb)
                rp_sb = newsb(rp)
            else:
                rp_sb, r_flag = -1, 1
            return rp, rp_sb, r_flag
        else:
            print "Petrosian interpolation failed!"
            rp, rp_sb, r_flag = -1, -1, 2
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
        delta = imgcenter - galcenter

        aper = EllipticalAperture(imgcenter, self.Rp, self.Rp/self.e,self.theta)
        bkg_asym = self.bkg_asymmetry(aper)

        asyms = defaultdict(list)
        prior_points = []

        #  minimize the galaxy asymmetry
        while True:
            ga = []
            dd = []
            deltas, points = utils.generate_deltas(imgcenter, .3, delta)

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

                    # create an array of asyms ... 
                    ga.append(galasym)
                    dd.append(den)
                    # ... and a dictionary that maps each asym to a 
                    #point on the image grid
                    asyms[p].append([galasym,den])

                    #galaxy_plot.asym_plot(newdata, residual, aperture, 
                    #                      galcenter, self.Rp)

                # just take the value that's already in the dictionary 
                # for that point
                else:
                    ga.append(asyms[p][0][0])
                    dd.append(asyms[p][0][1])

            # if the asymmetry found at the original center 
            # (first delta in deltas) is the minimum, we're done!
            if ga[0] == np.min(ga):
                center = imgcenter - deltas[0]
                # save the residual image 
                res = fits.ImageHDU(data=residual)
                res.writeto('output/asymimgs/'+self.name+'_res.fits', 
                            clobber=True, output_verify='silentfix')
                #galaxy_plot.asym_plot(newdata, residual, aperture, 
                #                      center, self.Rp, self.med, self.rms)
                return ga[0]-bkg_asym/dd[0], center[0], center[1]
            else:
                minloc = np.where(ga == np.min(ga))[0]
                delta = deltas[minloc[0]]
                prior_points = list(points)


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
        annuli = [CircularAnnulus((self.Acx, self.Acy),radii[i-1],radii[i])\
                  for i in range(1,len(radii))]
        counts = np.hstack([aperture_photometry(image, an, method='center') \
                            for an in annuli])['aperture_sum']
        cum_sum = np.array([np.sum(counts[0:idx]) \
                            for idx in range(1,len(counts)+1)])

        # Calculate the total flux in a Circular aperture of 1.5*rpet
        tot_aper = CircularAperture((self.Acx, self.Acy), self.Rp)
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
            r20 = r80 = -1
        
        return r20, r80, 5*np.log10(r80/r20)
        
    def get_gini(self, image):
        '''
        Need all pixels associated with a galaxy -- use my aperture thing? 
        1. All pixels within 1 Rp
        2. All pixels above the mean SB at 1 Rp
        Gonna do 1. for  now but want to try 2. as well
        '''
        print "calculating Gini..."

        # Create aperture at center of galaxy 
        # (Using my Aperture, not photutils, because I want access to the 
        # individual pixels/mask -- not just a sum of pixel values

        # Mask 1: galaxy pixels defined exactly within 1 petrosian radius
        ap = utils.EllipticalAperture((self.x, self.y), self.Rp, self.Rp/self.e,
                                      self.theta, image)
        mask1 = ap.aper * image

        # Mask 2: galaxy pixels defined as those with flux >= SB at 1 petro rad
        mask2 = utils.get_SB_Mask(self.Rp, self.Rp_SB, image, self.name)*image

        galpix = [mask1[mask1 != 0.], mask2[mask2 != 0.]]

        ginis = []
        for pix in galpix:
            galpix_sorted = sorted(np.abs(pix))
            xbar = np.mean(galpix_sorted)
            n = len(galpix_sorted)

            factor = 1/(xbar*n*(n-1))
            #print xbar, factor, n
            gsum = [2*i-n-1 for i, p in enumerate(galpix_sorted)]
            ginis.append(factor*np.dot(gsum, galpix_sorted))
        '''
        for pix in galpix:
            galpix_sorted= sorted(pix)
            xbar = np.mean(galpix_sorted)
            n = len(galpix_sorted)

            factor = 1/(xbar*n*(n-1))
            #print xbar, factor, n
            gsum = [2*i-n-1 for i, p in enumerate(galpix_sorted)]
            ginis.append(factor*np.dot(gsum, galpix_sorted))
        #'''
        #print "<S/N>:", stn
        print ginis

        return ginis 
            
    def get_m20(self, image):

        print "Calculating M20..."
        
        center = [image.shape[0]/2., image.shape[1]/2.]
        galcenter = np.array([self.x, self.y])

        mxrange = [center[0]-round(0.5*self.Rp), center[0]+round(0.5*self.Rp)]
        myrange = [center[1]-round(0.5*self.Rp), center[1]+round(0.5*self.Rp)]

        x, y = np.ogrid[mxrange[0]:mxrange[1], myrange[0]:myrange[1]]

        # Create a "distance grid" - each element's value is it's distance
        # from the center of the image for the entire image
        x2, y2 = np.ogrid[:image.shape[0], :image.shape[1]]      
        dist_grid = (center[0] - x2)**2 + (center[1] - y2)**2

        # create aperture at center of galaxy (mask)
        gal_aper = utils.EllipticalAperture(center, self.Rp, self.Rp/self.e,
                                            self.theta, image)
        mask1 = gal_aper.aper*image
        mask2 = utils.get_SB_Mask(self.Rp, self.Rp_SB, image, self.name)*image

        # We want the 0.0 element in dist_grid to correspond to each "test"
        # center that we calculate Mtot on so we have to shift the dist_grid
        # for each calculation of Mtot
        mtots1, mtots2 = [], []
        for i in x:
            for j in y.transpose():
                indices0 = range(center[0]-i, center[0]-i + 
                                     dist_grid.shape[0])
                indices1 = range(center[1]-j, center[1]-j + 
                                     dist_grid.shape[1])
                shift_grid = dist_grid.take(indices0, axis=0, mode='wrap')\
                                      .take(indices1, axis=1, mode='wrap')
                mtots1.append(np.sum(mask1*shift_grid))
                mtots2.append(np.sum(mask2*shift_grid))

        mtots = [mtots1, mtots2]
        Mtot = [np.min(mtots1), np.min(mtots2)]
        M20, xc, yc = [], [], []

        mlevel = []
        for idx, m in enumerate(mtots):

            mtot = np.array(m).reshape([len(x), len(y.transpose())])
            xc.append(x[np.where(mtot == Mtot[idx])[0]][0])
            yc.append(y.transpose()[np.where(mtot == Mtot[idx])[1]][0])

            indices0 = range(center[0]-xc[0],center[0]-xc[0]+dist_grid.shape[0])
            indices1 = range(center[0]-yc[0],center[0]-yc[0]+dist_grid.shape[0])
            grid = dist_grid.take(indices0, axis=0, mode='wrap')\
                            .take(indices1, axis=1, mode='wrap')
            
            if idx == 0: 
                # recreate the Aperture Mask now that we have the official 
                # center (the other mask doesn't depend on where the center 
                # is found since that mask is created based on SB cut
                m20_aper = utils.EllipticalAperture((xc[0], yc[0]), self.Rp, 
                                            self.Rp/self.e, self.theta, image)
                galpix = m20_aper.aper*image
            else:
                galpix = mask2

            # Sort the galaxy and grid pixels according to the descending order
            # of the galaxy pixels
            grid_sorted = np.array([i for j,i in sorted(zip(galpix.flatten(),\
                                                             grid.flatten()),\
                                                         reverse=True)])
            galpix_sorted = np.array(sorted(galpix.flatten(), reverse=True))

            # Calculate the 20% of the total flux of the galaxy
            ftot20 = 0.2*np.sum(galpix_sorted)
            
            # calculate the cumulative flux
            fcumsum = np.cumsum(galpix_sorted)
        
            # find where fcumsum is less than ftot20 -- 
            m20_pix_idx = np.where(fcumsum < ftot20)[0]
            
            # Calculate M20 from the brightest 20% of galaxy pixels
            m20_galpix = galpix_sorted[m20_pix_idx]
            m20_distpix = grid_sorted[m20_pix_idx]
            pdb.set_trace()
            M20.append(np.log10(np.sum(m20_galpix*m20_distpix)/Mtot[idx]))
            
            mlevel.append(np.min(m20_galpix))
            
        self._Mlevel1, self._Mlevel2 = mlevel[0], mlevel[1]
        return M20[0], M20[1], xc[0][0], yc[0][0], xc[1][0], yc[1][0]


        '''
        ##############################################################3
        # Old Method:
        x = y = np.arange(center[0]-round(0.5*self.Rp), 
                          center[0]+round(0.5*self.Rp))
        xx, yy = np.meshgrid(x, y)
        # distance grid: 
        x2, y2 = np.ogrid[:image.shape[0], :image.shape[1]]

        dist_grids = [(x2-xi)**2 + (y2-yi)**2 for xi, yi in \
                      zip(xx.flatten(), yy.flatten())]

        # create aperture at center of galaxy (mask)
        gal_aper = utils.EllipticalAperture(center, 
                            self.Rp, self.Rp/self.e, self.theta, image)

        mtots2 = [np.sum(gal_aper.aper*image*grid) for grid in dist_grids]
        Mtot2 = np.min(mtots2)
        xc2, yc2 = xx.flatten()[mtots2 == Mtot2], yy.flatten()[mtots2 == Mtot2]
        grid2 = dist_grids[np.where(mtots2 == Mtot2)[0]]

        print Mtot2, xc2, yc2
        pdb.set_trace()
        '''  
        
    def table(self, init=False): 
          
        the_dict = self.__dict__
        r = re.compile(r"_.+")
        matching_keys = filter(r.match, the_dict.keys())
        #matching_keys.append('name')
        
        for key in matching_keys:
            del the_dict[key]

        if init:
            names = ['name', 'cat', 'oflag', 'uflag', 'ra', 'dec', 'e', 'x', 
                     'y', 'a', 'b', 'theta', 'elipt', 'kron', 'Rp', 'Rpflag',
                     'Rp_SB', 'r20', 'r80', 'C', 'A', 'G1', 'G2', 'M1', 'M2', 
                     'Acx', 'Acy', 'Mcx1', 'Mcy1', 'Mcx2', 'Mcy2', 'med', 'rms']
            dtypes = []
            for n in names:
                if n in ['name']:
                    dtypes.append('S80')
                elif n in ['cat', 'oflag', 'uflag', 'Rpflag']:
                    dtypes.append('i')
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
    parser.add_argument('output', type=str,
        help='Specify the desired name for output catalog.')
    args = parser.parse_args()

    fitsfiles = np.array(sorted(glob.glob(args.directory+'*.fits')))
    fitsfiles = sorted(glob.glob(args.directory))

    #fitsfiles=np.full(10,fill_value=fitsfiles[0], dtype='|S80')
    #fitsfiles = fitsfiles[30:40]
    #pdb.set_trace()
    warnings.filterwarnings('ignore', message='Overwriting existing file .*',
                            module='pyfits')

    outdir = 'output/datacube/'

    for idx, f in enumerate(fitsfiles): 
        basename = os.path.basename(f)
        filename = outdir+'f_'+basename

        if not os.path.isfile(filename):
            print "File not found! Running SExtractor before proceeding."

        print "Cleaning ", os.path.basename(f)
        flags = clean.clean_frame(f, outdir)

        print "Running", os.path.basename(f)
        hdulist = fits.open(filename, memmap=True)
        g = Galaxy(hdulist, filename, flags)
        if g.Rp > 0.:
            galaxy_plot.plot(g, hdulist)
            if idx == 0:
                t, gal_dict = g.table(init=True)
                t.add_row(gal_dict)
            else:
                t.add_row(g.table())
        else:
            t.add_row(g.table())
        hdulist.close()
        t.write(args.output, format='ascii.fixed_width')
        del g

    print "Parameter catalog complete.\nSaving catalog to file..."
    t.write(args.output, format='ascii')
    exit()  


if __name__ == '__main__':
    main()
    
