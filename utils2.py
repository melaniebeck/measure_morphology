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
    mask = np.where((SEseg != SEcat[num][idx]) & (SEseg != 0)) 
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
    things = SEcat[thing][indexes]
    if any(thing > threshold):
        Dist = np.min(distances[np.where(things>threshold)])
        Index = indexes[distances==Dist][0]
        Coord = coords[Index]
        Thing = SEcat[thing][Index]
        if Thing != things[0]:
            Flag = 1
    else:
        Dist = distances[0]
        Index = indexes[0]
        Coord = coords[Index]
        Thing = thing[0]
        Flag = 1
    return Dist, Index, Coord, Thing, Flag

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

    # initialize some shit for laterz
    Flag = 0
    center = [img.shape[0]/2., img.shape[1]/2.]
    Bdist = center[0]
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

    # if distance between center of BRIGHT and FAINT objs is small
    # then it's most likely well detected in FAINT (not blended)
    # clean on FAINT only
    if (DIST < sep):
        if (Bdist < sep) & (Fdist < sep):
            cln = clean_image(cln, fseg, fcat, FIndex, fseg)
            Flag, mode = 1, 'FAINT'

        if (Bdist < sep) & (Fdist > sep):
            cln = clean_image(cln, bseg, bcat, BIndex, fseg)
            Flag, mode = 2, 'BRIGHT'

        if (Bdist > sep) & (Fdist < sep):
            cln = clean_image(cln, fseg, fcat, FIndex, fseg)
            Flag, mode = 7, 'FAINT'

        if (Bdist > sep) & (Fdist > sep):
            cln = clean_image(cln, bseg, bcat, BIndex, fseg)
            Flag, mode = 8, 'BRIGHT'

    else:
        # I don't think this one can ever happen....
        if  (Bdist < sep) & (Fdist < sep):
            cln = clean_image(cln, bseg, bcat, BIndex, fseg)
            Flag, mode = 3, 'BRIGHT'

        if  (Bdist < sep) & (Fdist > sep):
            '''
            If DIST large but Bdist small --> Fdist must be large
            This means that obj is likely detected in BRIGHT but is blended
            into a nearby bright object in FAINT.
            Try Two Stage Cleaning mode:
              1. Clean on BRIGHT then 
              2. Run SE again in FAINT mode and clean on that. 
            '''
            if (Barea > 0.6*Farea) & (Barea < Farea):
                cln = clean_image(cln, fseg, fcat, FIndex, fseg)

            else:
                cln = clean_image(cln, bseg, bcat, BIndex, fseg)
                
                # run SE again in FAINT
                run_sextractor.run_SE(image, 'FAINT', outdir)
                f2seg, f2cat = fits.getdata(segnames[2]),fits.getdata(catnames[2])
                F2Index, F2dist = find_closest(center, zip(scat[x], scat[y]))

                cln = clean_image(cln, f2seg, f2cat, F2Index, f2seg)

            Flag, mode = 4, 'FAINT'
 
        if  (Bdist > sep) & (Fdist < sep):
            # not detected in BRIGHT but it is in FAINT -- Smooth first?
            run_sextractor.run_SE(image, 'SMOOTH', outdir)
            sseg, scat = fits.getdata(segnames[2]), fits.getdata(catnames[2])
            SIndex, Sdist = find_closest(center, zip(scat[x], scat[y]))

            cln = clean_image(cln, sseg, scat, SIndex, sseg)
            Flag, mode = 5, 'SMOOTH'

        if  (Bdist > sep) & (Fdist > sep):
            # if central object not found in BRIGHT and not well 
            # detected in FAINT
                       
            run_sextractor.run_SE(image, 'SMOOTH', outdir)
            sseg, scat = fits.getdata(segnames[2]), fits.getdata(catnames[2])
            SIndex, Sdist = find_closest(center, zip(scat[x], scat[y]))

            cln = clean_image(cln, sseg, scat, SIndex, sseg)
            Flag, mode = 6, 'SMOOTH'

    # Now that we've done the cleaning -- Let's test it! 
    # Run SE again in FAINT. If nothing above the min area is found near the 
    # center of the image then we OVERCLEANED. FLAG THIS.
    run_sextractor.run_SE(cln, 'FAINT', outdir)
    tseg = fits.getdata(outname+'_faint_seg_faint.fits')
    tcat = fits.getdata(outname+'_faint_cat_faint.fits')

    coords = zip(tcat[x], tcat[y])
    Tdist, TIndex, TCoord, Tarea, tFlag = closest_above_thresh(tcat, area, 
                                                               center, coords)
    ocln_flag = 0
    if Tdist > sep: 
        print 'OVERCLEANED!!!'
        ocln_flag = 1


    if mode == 'BRIGHT':
        ihdr.set('SE_MODE', 'BRIGHT', 'Sextractor mode used to clean image')
        ihdr.set('SECATIDX', BIndex, 
                 'Index in BRIGHT catalog denoting gal of interest')
        init2 = fits.ImageHDU(data=bseg, name='BSEG')
        init3 = fits.ImageHDU(data=fseg, name='FSEG')
        init4 = fits.TableHDU(data=bcat, name='BCAT')

    elif mode == 'FAINT':
        ihdr.set('SE_MODE', 'FAINT', 'Sextractor mode used to clean image')
        ihdr.set('SECATIDX', FIndex, 
                  'Index in FAINT catalog denoting gal of interest')
        init2 = fits.ImageHDU(data=bseg, name='BSEG') 
        init3 = fits.ImageHDU(data=fseg, name='FSEG')
        init4 = fits.TableHDU(data=fcat, name='FCAT') 

    elif mode == 'SMOOTH':
        ihdr.set('SE_MODE', 'SMOOTH', 'Sextractor mode used to clean image')
        ihdr.set('SECATIDX', FIndex, 
                 'Index in FAINT catalog denoting gal of interest')
        init2 = fits.ImageHDU(data=fseg, name='FSEG')
        init3 = fits.ImageHDU(data=sseg, name='SSEG')
        init4 = fits.TableHDU(data=scat, name='SCAT') 

        
    # SAVE ALL PRODUCTS TO DATA CUBE
    init0 = fits.ImageHDU(data=img, header=ihdr, name='ORG')
    init1 = fits.ImageHDU(data=cln, header=ihdr, name='CLN')

    newthing = fits.HDUList()
    for thing in (init0, init1, init2, init3, init4):
        newthing.append(thing)
    newthing.update_extend()
    newthing.writeto(outdir+'f_'+basename+'.fits', output_verify='silentfix', 
                     clobber=True)
    
    # clean up directory
    #os.system("rm "+outdir+"*bright*.fits")
    os.system("rm "+outdir+"*faint*.fits")
    os.system("rm "+outdir+"*smooth*.fits")

    return FIndex, Fdist, Bdist, DIST, Farea, Barea, Flag, aFlag, ocln_flag


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
    








