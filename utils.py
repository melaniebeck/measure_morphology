import os
import string
import numpy as np
import pyfits as fits
from astropy.table import Table
from scipy.interpolate import interp1d
from collections import OrderedDict
from random import gauss
import pdb #"""for doing an IDL-like stop"""
from scipy.spatial import cKDTree
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import fast_ffts
import run_sextractor


def find_closest(point, listofpoints, k=1):
    
    # create a KDTree
    tree = cKDTree(listofpoints)
    # query the tree to find the set of coordinates closest to that point
    dist, index = tree.query(point, k)
    # return the value of the closest point
    return index, dist


def clean_pixels(data, mask, segmap):
    #mean = np.mean(data[segmap == 0])
    #std = np.std(data[segmap == 0])
    med = np.median(data[segmap == 0])
    rms = np.sqrt(np.mean(np.square(data[segmap==0])))
    for pixel in zip(mask[0], mask[1]):
        #pdb.set_trace()
        data[pixel] = gauss(med, rms)
    return data

def clean_image(image, SEseg, SEcat, idx, bkgseg):
    mask = np.where((SEseg != SEcat['NUMBER'][idx]) & (SEseg != 0)) 
    image = clean_pixels(image, mask, bkgseg)
    return image

def closest_above_thresh(SEcat, thing, center, coords, threshold=50., k=10):
    '''
    Find an object that is closest to the center of the image that also meets
    a minimum threshold value for some desired quantity (thing)
    Returns information about that object: 
       Distance to the center
       Coordinates of the object in the image
       Index of the object in the SE catalog
       Thing value that meets the set threshold

    For example: find galaxy closest to center that has a pixel area larger than 
                 fifty pixels**2

    Also returns a flag if the closest object does NOT meet the threshold
    '''
    Flag = 0

    indexes, distances = find_closest(center, coords, k)
    idxs, dists = indexes[distances != np.inf], distances[distances != np.inf]
    things = SEcat[thing][idxs]
    if np.any(np.array(things) > threshold):
        Dist = np.min(dists[things > threshold])
        Index = idxs[dists==Dist][0]
        Coord = coords[Index]
        Thing = SEcat[thing][Index]
        if Thing != things[0]:
            Flag = 1
    else:
        Dist = dists[0]
        Index = idxs[0]
        Coord = coords[Index]
        Thing = things[0]
        Flag = 1
    return Dist, Index, Coord, Thing, Flag

def savedata(mode, galnumber, data=[], names=[]):
    datacube=[]
    ihdr.set('SE_MODE', mode, 'Sextractor mode used to clean image')
    ihdr.set('SECATIDX', galnumber, 
             'Index in'+mode+'catalog denoting gal of interest')
    for datum, name in zip(data, names):
        datacube.append(fits.ImageHDU(data=datum, name=name))
    return datacube

def clean_frame(image, outdir, sep=17.):
    '''
    Gonna try a two-stage method of cleaning:
        1. run the cutout using the Faint SE parameters
        2. run again using the Bright SE parameters
        3. clean the image based on the Faint parameters
        4. check the seg maps -- if there is anything found in the Bright
            that wasn't found with the Faint, clean those too
        5. use (return) as the "official" parameters those that are found via
            the Bright run 

        sep is the minimum separation distance (in pixels) between: 
            1. the gal in BRIGHT and the center
            2. the gal in FAINT and the center
            3. the gal in FAINT and the gal in BRIGHT
    #'''

    import warnings
    warnings.filterwarnings('ignore')

    # initialize some shit for laterz
    Flag = 0
    Barea, B2dist = 0., 0.

    # commonly used SE keywords
    num = 'NUMBER'
    area = 'ISOAREA_IMAGE'
    x, y = 'X_IMAGE', 'Y_IMAGE'


    basename = os.path.basename(os.path.splitext(image)[0])
    outname = outdir+basename
    catnames = [outname+'_bright_cat.fits', outname+'_faint_cat.fits', 
                outname+'_smooth_cat.fits']
    segnames = [outname+'_bright_seg.fits', outname+'_faint_seg.fits',
                outname+'_smooth_seg.fits']

    # run SE in FAINT and BRIGHT modes
    run_sextractor.run_SE(image, 'BRIGHT', outdir)
    run_sextractor.run_SE(image, 'FAINT', outdir)

    # READ IN ORIG FITS-FILE AND both Bright/Faint SEG-MAPs
    img, ihdr = fits.getdata(image, header=True)
    cln = img.copy()

    bseg, fseg = fits.getdata(segnames[0]), fits.getdata(segnames[1]) 
    bcat, fcat = fits.getdata(catnames[0]), fits.getdata(catnames[1])

    center = [img.shape[0]/2., img.shape[1]/2.]
    Bdist = center[0]

    # test BRIGHT segmap for any objects whatsoever
    if np.any(np.nonzero(bseg)):
        # find object in BRIGHT closest to center
        BIndex,Bdist = find_closest(center,zip(bcat[x],bcat[y]))
        BCoord = zip(bcat[x], bcat[y])[BIndex]

        # if obj detected in BRIGHT -- find next closest
        idx, dst = find_closest(BCoord, zip(bcat[x], bcat[y]), k=2)

        if any(np.isinf(dst)):
            B2Index, B2dist = 0., 0.
        else:
            B2Index, B2dist = idx[1], dst[1]

            # find the combined area of these two obj (pixels**2)
            Barea = bcat[area][BIndex] + bcat[area][B2Index]

    coords = zip(fcat[x], fcat[y])
    Fdist, FIndex, FCoord, Farea, aFlag = closest_above_thresh(fcat, area, 
                                                               center, coords)

    DIST = abs(Fdist - Bdist)

    if (DIST < sep):
        # FLAG 1: MOST COMMON CATEGORY --> CLEAN IN FAINT MODE
        if (Bdist < sep) & (Fdist < sep):
            cln = clean_image(cln, fseg, fcat, FIndex, fseg)
            Flag, mode = 1, 'FAINT'

        # FLAG 2: CLEAN IN BRIGHT MODE & FLAG THESE!!
        if (Bdist < sep) & (Fdist > sep):
            ''' There is only one obj in here:
            Two super bright galaxies super close together
            Seen as distinct objs in BRIGHT (but with crazy square edges)
            Seen as one blended obj in FAINT
            JUST FLAG THIS BITCH
            '''
            cln = clean_image(cln, bseg, bcat, BIndex, fseg)
            Flag, mode = 2, 'BRIGHT'

        # FLAG 7: CLEAN IN FAINT MODE
        if (Bdist > sep) & (Fdist < sep):
            ''' There aren't many of these
            They're oddballs but most are well cleaned in FAINT
            '''
            cln = clean_image(cln, fseg, fcat, FIndex, fseg)
            Flag, mode = 7, 'FAINT'

        # FLAG 8: TWO STAGE CLEANING -- BRIGHT --> FAINT --> 
        # RUN SE AGAIN IN FAINT
        if (Bdist > sep) & (Fdist > sep):
            ''' If it's not detected in BRIGHT - should run SE in SMOOTH mode
            no we shouldn't. :(
            '''
            #run_sextractor.run_SE(image, 'SMOOTH', outdir)
            #sseg, scat = fits.getdata(segnames[2]), fits.getdata(catnames[2])
            #SIndex, Sdist = find_closest(center, zip(scat[x], scat[y]))
            #cln = clean_image(cln, sseg, scat, SIndex, sseg)
            cln = clean_image(cln, bseg, bcat, BIndex, fseg)

            cln_sv = cln.copy()
            cln_sv = fits.ImageHDU(data=cln_sv, name='MID_CLN')
            cln_sv.writeto(outname+'_mid_cln.fits', output_verify='silentfix', 
                           clobber=True)

            # run SE again in FAINT
            run_sextractor.run_SE(outname+'_mid_cln.fits', 'FAINT', 
                                  outdir, outstr2='run2')
            f2seg = fits.getdata(outname+'_mid_cln_faint_run2_seg.fits')
            f2cat = fits.getdata(outname+'_mid_cln_faint_run2_cat.fits')
            coords = zip(f2cat[x], f2cat[y])
            # find closest obj to center with area above threshold
            Fdist, FIndex, FCoord, Farea, aFlag = \
                        closest_above_thresh(f2cat, area, center, coords, k=5)
            cln = clean_image(cln, f2seg, f2cat, FIndex, f2seg)
            Flag, mode = 8, 'FAINT2'
            
    else:
        # FLAG 3: SHOULD NEVER HAPPEN
        if  (Bdist < sep) & (Fdist < sep):
            cln = clean_image(cln, bseg, bcat, BIndex, fseg)
            Flag, mode = 3, 'BRIGHT'

        # FLAG 4: TWO STAGE CLEANING - BRIGHT --> RUN SE AGAIN IN FAINT
        if  (Bdist < sep) & (Fdist > sep):
            '''
            If DIST large but Bdist small --> Fdist must be large
            This means that obj is likely detected in BRIGHT but is blended
            into a nearby bright object in FAINT.
            Try Two Stage Cleaning mode:
              1. Clean on BRIGHT then 
              2. Run SE again in FAINT mode and clean on that. 
            '''
            cln = clean_image(cln, bseg, bcat, BIndex, fseg)

            #save this image so that I can run SE on it
            cln_sv = cln.copy()
            cln_sv = fits.ImageHDU(data=cln_sv, name='MID_CLN')
            cln_sv.writeto(outname+'_mid_cln.fits', output_verify='silentfix',
                           clobber=True)
            
            # run SE again in FAINT
            run_sextractor.run_SE(outname+'_mid_cln.fits', 'FAINT', 
                                  outdir, outstr2='run2')
            f2seg = fits.getdata(outname+'_mid_cln_faint_run2_seg.fits')
            f2cat = fits.getdata(outname+'_mid_cln_faint_run2_cat.fits')
            coords = zip(f2cat[x], f2cat[y])
            # find closest obj to center with area above threshold
            Fdist, FIndex, FCoord, Farea, aFlag = \
                        closest_above_thresh(f2cat, area, center, coords, k=5)
            cln = clean_image(cln, f2seg, f2cat, FIndex, f2seg)
            Flag, mode = 4, 'FAINT2'
 
        # FLAG 5: CLEAN IN SMOOTH MODE
        if  (Bdist > sep) & (Fdist < sep):
            ''' These are mostly faint objects not detected in BRIGHT
            run SE in SMOOTH mode and then clean
            '''
            run_sextractor.run_SE(image, 'SMOOTH', outdir)
            sseg, scat = fits.getdata(segnames[2]), fits.getdata(catnames[2])
            SIndex, Sdist = find_closest(center, zip(scat[x], scat[y]))

            cln = clean_image(cln, sseg, scat, SIndex, sseg)
            Flag, mode = 5, 'SMOOTH'

        # FLAG 6: CLEAN IN FAINT MODE -- ALL GARBAGE ANYWAY
        if  (Bdist > sep) & (Fdist > sep):
            ''' this is mostly a garbage bin of crap 
            any object in here needs to be flagged and is likely not a true
            galaxy at all!
            '''
            #run_sextractor.run_SE(image, 'SMOOTH', outdir)
            #sseg, scat = fits.getdata(segnames[2]), fits.getdata(catnames[2])
            #SIndex, Sdist = find_closest(center, zip(scat[x], scat[y]))

            cln = clean_image(cln, fseg, fcat, FIndex, fseg)
            Flag, mode = 6, 'FAINT'

    # Now that we've done the cleaning -- Let's test it!    
    run_sextractor.run_SE(outdir+'f_'+basename+'.fits[1]', 'FAINT', 
                          outdir, outstr2='test')
    tseg = fits.getdata(outdir+'f_'+basename+'_faint_test_seg.fits')
    tcat = fits.getdata(outdir+'f_'+basename+'_faint_test_cat.fits')

    coords = zip(tcat[x], tcat[y])
    index, dist = find_closest(center, coords, k=10)
    tarea = tcat[area][index[np.where(dist != np.inf)]]

    # If we find obj near the center is too small then we overcleaned it
    ocln_flag = 0
    if (dist[0] < sep) & (tarea[0] < 50.): 
        print 'OVERCLEANED!!!'
        ocln_flag = 1

    # If we find large objs far from the center then we didn't clean enough
    if (np.any(dist[1::] > sep)) & (np.any(tarea[1::] > 100.)):
        print 'UNDER CLEANED!!'
        ocln_flag = 2


    # Save all major data products
    if mode == 'BRIGHT':
        datacube = savedata(mode, BIndex, data=[bseg, fseg, bcat], 
                            names=['BSEG', 'FSEG', 'BCAT'])

    elif mode == 'FAINT':
        datacube = savedata(mode, FIndex, data=[bseg, fseg, fcat], 
                            names=['BSEG', 'FSEG', 'FCAT'])

    elif mode == 'FAINT2':
        datacube = savedata('TWOSTAGE', FIndex,data=[bseg, fseg, f2seg, f2cat], 
                            names=['BSEG', 'FSEG', 'F2SEG', 'F2CAT'])
        datacube.insert(0, cln_sv)


    elif mode == 'SMOOTH':
        datacube = savedata(mode, SIndex, data=[bseg, fseg, sseg, scat], 
                            names=['BSEG', 'FSEG', 'SSEG', 'SCAT'])

        
    # SAVE ALL PRODUCTS TO DATA CUBE
    init0 = fits.ImageHDU(data=img, header=ihdr, name='ORG')
    init1 = fits.ImageHDU(data=cln, header=ihdr, name='CLN')
    datacube.insert(0, init1)
    datacube.insert(0, init0)

    newthing = fits.HDUList()
    for thing in datacube: 
        newthing.append(thing)
    newthing.update_extend()
    newthing.writeto(outdir+'f_'+basename+'.fits', output_verify='silentfix', 
                     clobber=True)


    # clean up directory
    #os.system("rm "+outdir+"*bright*.fits")
    #os.system("rm "+outdir+"*faint*.fits")
    #os.system("rm "+outdir+"*smooth*.fits")


    return FIndex, Fdist, Bdist, DIST, Farea, Barea, Flag, ocln_flag


def shift_image(data, deltax, deltay, phase=0, nthreads=1, use_numpy_fft=False,
        return_abs=True):
    """
    FFT-based sub-pixel image shift
    http://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation/content/html/efficient_subpixel_registration.html

    Will turn NaNs into zeros

    Adam Ginsberg code from agpy
    """

    fftn,ifftn = fast_ffts.get_ffts(nthreads=nthreads, use_numpy_fft=use_numpy_fft)

    if np.any(np.isnan(data)):
        data = np.nan_to_num(data)
    ny,nx = data.shape
    Nx = np.fft.ifftshift(np.linspace(-np.fix(nx/2),np.ceil(nx/2)-1,nx))
    Ny = np.fft.ifftshift(np.linspace(-np.fix(ny/2),np.ceil(ny/2)-1,ny))
    Nx,Ny = np.meshgrid(Nx,Ny)
    gg = ifftn( fftn(data)* np.exp(1j*2*np.pi*(-deltax*Nx/nx-deltay*Ny/ny)) * np.exp(-1j*phase) )
    if return_abs:
        return np.abs(gg)
    else:
        return gg

#----------------------------------------------------------------------------------#

class EllipticalAperture(object):

    def __init__(self, xycenter, a, b, theta, data):
        self.x, self.y = xycenter
        self.a, self.b = a, b
        self.theta = theta
        self.data = data
        #self.radius = radius
        #self.skyradius = skyradius

        self.aper = self.make_aperture()
        #self.skyaper = self.make_aperture(self.radius,self.skyradius)


    def make_aperture(self):
        u,v = [], [] 
        for x,y in np.ndindex(self.data.shape):
            point = (x-self.x, y-self.y)
            cosang = np.cos(self.theta)
            sinang = np.sin(self.theta)
            xtemp = point[0]*cosang + point[1]*sinang
            ytemp = -point[0]*sinang + point[1]*cosang
            u.append(xtemp**2)
            v.append(ytemp**2)
            
        ellipse1 = np.array([u0/self.a**2 + v0/self.b**2 for u0, v0 in zip(u, v)])

        mask = np.array([True if ellipse1[x] <= 1 else False \
                         for x in np.ndindex(ellipse1.shape)])
        mask = mask.reshape(self.data.shape)

        '''
        ellipse = Ellipse([self.xc, self.yc], self.a, self.b, self.theta)
        mask = np.array([True if ellipse.contains_point([x,y]) 
                         else False for x,y in np.ndindex(self.data.shape)])
        mask = mask.reshape(self.data.shape)
        #'''
        #pdb.set_trace()

        return mask.astype('float')
    
    @staticmethod
    def num_pix(weights):
        return np.sum(weights)    

    def run(self,longTable=False):
        '''Return:
        (skysub_aper_counts, aper_area, total_aper_counts,
         sky_area, total_sky_counts, med_sky_counts, std_sky_counts)
        '''
        print "here i am"
        aper = self.aper * self.data

        pdb.set_trace()
        #sky = self.skyaper * self.data
        #aper_area = self.radius*self.radius*np.pi
        #sky_area = self.skyradius**2*np.pi - aper_area

        return 0
        
def generate_deltas(center, shiftsize, shift):
    increments = (0.,-shiftsize,shiftsize)
    points = []
    deltas = []
    for s in increments:
        x = center[0]+shift[0]+s
        for h in increments:
            #pdb.set_trace()
            y = center[1]+shift[1]+h
            points.append((round(x,6),round(y,6)))
            deltas.append([round(s+shift[0],6),round(h+shift[1],6)])
    
    return np.array(deltas), points

def scale(galmask, imgsize):
    # determine the scale of the galaxy to the background
    galarea = np.sum(galmask)
    bkgarea = imgsize[0]*imgsize[1] - galarea
    return galarea/bkgarea    

def measure_asymmetry(data, galmask, bkgmask, shift, scale):
    
    # shift data such that the prescribed center corresponds to the center of the image
    newdata = shift_image(data, shift[1], shift[0])
    
    # mask out the galaxy and the background
    galpix = galmask * newdata
    galpixrot = galmask * np.rot90(newdata, 2)

    bkgpix = bkgmask * newdata
    bkgpixrot = bkgmask * np.rot90(newdata, 2)

    # calculate the asymmetry
    denom = np.sum(np.abs(galpix))
    galasym = np.sum(np.abs(galpix - galpixrot))/denom
    bkgasym = np.sum(np.abs(bkgpix - bkgpixrot))*scale/denom
    
    return galasym-bkgasym, bkgasym


def plot2d(data, cmap='gray_r', **kwargs):
    
    #if datax.ndim != datay.ndim: 
    #    print "data set sizes don't match!"
    #    return 0

    #num = np.arange(0, datax.ndim, 1)

    #f, axarr = plt.subplots(datax.ndim)
    #for n in num: 
    #    axarr[n].plot(datax[n], datay[n])
    #    axarr[n].set_title

    for name, value in kwargs.items(): 
        print "%s = %f" %(name, value)
    
    plt.figure()
    plt.imgshow(kwargs['img'], cmap=kwargs['cmap'] )
    

def hist2d(data, *args, **kwargs):
    if datax.ndim != datay.ndim: 
        print "data set sizes don't match!"
        return 0

    num = np.arange(0, datax.ndim, 1)

    f, axarr = plt.subplots(datax.ndim)
    for n in num: 
        axarr[n].plot(datax[n], datay[n])
        axarr[n].set_title    



def main():
    '''
    data = fits.getdata('output/f_0001_149.911728_2.70517_acs_I_mosaic_30mas_sci.fits', ext=1)
    catalog = Table.read('bigsample_mycatv4.fits')
    info = catalog[0]

    ap = EllipticalAperture( (info['x'], info['y']), info['rpet'], 
                             info['rpet']/info['e'], info['theta'], data)
    aper = ap.aper
    test = ap.run()

    
    pdb.set_trace()
    data[np.where(~aper)]=0.
    
    plt.figure()
    plt.imshow(data, cmap='gray_r', origin='lower')
    plt.show()

    exit()
    #'''

    center = [251,251]
    shiftsize = 1.
    deltas, points = generate_deltas(center,shiftsize,[0.,0.])
    print deltas

    pdb.set_trace()
    exit()

if __name__ == '__main__':
    main()
    








