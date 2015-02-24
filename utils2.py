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


def find_closest(point, listofpoints):
    
    # create a KDTree
    tree = cKDTree(listofpoints)
    # query the tree to find the set of coordinates closest to that point
    dist, index = tree.query(point)
    # return the value of the closest point
    return listofpoints[index], index, dist


def clean_pixels(data, mask, segmap):
    #mean = np.mean(data[segmap == 0])
    #std = np.std(data[segmap == 0])
    med = np.median(data[segmap == 0])
    rms = np.sqrt(np.mean(np.square(data[segmap==0])))
    for pixel in zip(mask[0], mask[1]):
        #pdb.set_trace()
        data[pixel] = gauss(mean, std)
    return data

def clean_frame(image, outdir):
    '''
    Gonna try a two-stage method of cleaning:
        1. run the cutout using the Faint SE parameters
        2. run again using the Bright SE parameters
        3. clean the image based on the Faint parameters
        4. check the seg maps -- if there is anything found in the Bright
            that wasn't found with the Faint, clean those too
        5. use (return) as the "official" parameters those that are found via
            the Bright run 
    #'''

    basename = os.path.basename(os.path.splitext(image)[0])
    outname = outdir+basename
    catnames = [outname+'_bright_cat.fits', outname+'_faint_cat.fits']
    segnames = [outname+'_bright_seg.fits', outname+'_faint_seg.fits']

    # run SE in FAINT and BRIGHT modes
    run_sextractor.run_SE(image, 'BRIGHT', outdir)
    run_sextractro.run_SE(image, 'FAINT', outdir)

    # READ IN ORIG FITS-FILE AND both Bright/Faint SEG-MAPs
    img, ihdr = fits.getdata(image, header=True)
    bseg, bshdr = fits.getdata(segnames[0], header=True)
    fseg, fshdr = fits.getdata(segnames[1], header=True)

    # READ IN THE SE CATALOGs (Bright & Faint)
    bcat, bchdr = fits.getdata(catnames[0], header=True)
    fcat, fchdr = fits.getdata(catnames[1], header=True)

    cln = img.copy()
    
    # FIND GALAXY OF INTEREST (it will be in the CENTER OF IMAGE)
    center = [img.shape[0]/2., img.shape[1]/2.]

    # test BRIGHT segmap for any objects whatsoever: if ANY nonzero then SOMETHING is detected
    if np.any(np.nonzero(bseg)):
        # find object in the Bright catalog closest to center
        myBCoord, myBIndex, Bdist = find_closest(center, zip(bcat['X_IMAGE'], bcat['Y_IMAGE']))

        # if dist is less than 10 pixels (0.3'') then this galaxy is our object
        if Bdist < 10.:
            myFCoord, myFIndex, Fdist = find_closest(center, zip(fcat['X_IMAGE'], fcat['Y_IMAGE']))

            # if this object is also found in Faint -- Only need to clean on FAINT
            if abs(Fdist - Bdist) < 7.:
                cleanmaskB = np.where((bseg != 0) & (fseg != fcat['NUMBER'][myFIndex]))
                cleanmaskF = np.where((fseg != 0) & (bseg == 0) & (fseg != fcat['NUMBER'][myFIndex]))
                #cln = clean_pixels(cln, cleanmask, fseg)
            
            # otherwise: clean on BRIGHT then on FAINT
            else:
                cleanmaskB = np.where((bseg != 0) & (bseg != bcat['NUMBER'][myBIndex]))
                cleanmaskF = np.where((fseg != 0) & (bseg == 0) & (bseg != bcat['NUMBER'][myBIndex]))

            cln = clean_pixels(cln, cleanmaskB, fseg)
            cln = clean_pixels(cln, cleanmaskF, fseg)

            # if obj is detected in BRIGHT, save BRIGHT catalog info
            bchdr.set('SECATIDX', myBIndex, 'Index in SE catalog denoting gal of interest')
            init4 = fits.TableHDU(data=bcat, header=bchdr, name='BCAT') 
    
            init3 = fits.ImageHDU(data=fseg, header=fshdr, name='FSEG')

        # if detection is > 10 pixels from center -- obj was NOT detected in BRIGHT
        # clean on FAINT only
        else:
            ''' SMOOTH ON THE ORIGINAL IMAGE
                RERUN SEXTRACTOR WITH FAINT PARAMETERS
                before cleaning '''
            run_sextractor.run_SE(img_name, 'output/cUF_'+fitsname, 'output/sUF_'+fitsname, 
                           'config_ultrafaint.sex')
            ufseg, ufshdr = fits.getdata('output/sUF_'+fitsname, header=True)
            ufcat, ufchdr = fits.getdata('output/cUF_'+fitsname, header=True)
            myUFCoord, myUFIndex, UFdist = find_closest(center, 
                                                    zip(ufcat['X_IMAGE'], ufcat['Y_IMAGE']))
            cleanmask = np.where((ufseg != 0) & (ufseg != ufcat['NUMBER'][myUFIndex]))
            cln = clean_pixels(cln, cleanmask, ufseg)

            # if only found in FAINT, save info from FAINT catalog
            ufchdr.set('SECATIDX', myUFIndex, 'Index in SE catalog denoting gal of interest')        
            init4 = fits.TableHDU(data=ufcat, header=ufchdr, name='UFCAT')
            
            init3 = fits.ImageHDU(data=ufseg, header=ufshdr, name='UFSEG')
            

    # if NOTHING is detected in BRIGHT clean with FAINT only
    else:
        ''' SMOOTH ON THE ORIGINAL IMAGE
        RERUN SEXTRACTOR WITH FAINT PARAMETERS
        before cleaning '''
        run_sextractor(img_name, 'output/cUF_'+fitsname, 'output/sUF_'+fitsname, 
                       'config_ultrafaint.sex')
        ufseg, ufshdr = fits.getdata('output/sUF_'+fitsname, header=True)
        ufcat, ufchdr = fits.getdata('output/cUF_'+fitsname, header=True)
        myUFCoord, myUFIndex, UFdist = find_closest(center, 
                                                    zip(ufcat['X_IMAGE'], ufcat['Y_IMAGE']))
        cleanmask = np.where((ufseg != 0) & (ufseg != ufcat['NUMBER'][myUFIndex]))
        cln = clean_pixels(cln, cleanmask, ufseg)
        
        # if only found in FAINT, save info from FAINT catalog
        ufchdr.set('SECATIDX', myUFIndex, 'Index in SE catalog denoting gal of interest')        
        init4 = fits.TableHDU(data=ufcat, header=ufchdr, name='UFCAT')
        
        init3 = fits.ImageHDU(data=ufseg, header=ufshdr, name='UFSEG')

    #pdb.set_trace()

    # SAVE ALL PRODUCTS TO DATA CUBE
    init0 = fits.ImageHDU(data=img, header=ihdr, name='ORG')
    init1 = fits.ImageHDU(data=cln, header=ihdr, name='CLN')
    init2 = fits.ImageHDU(data=bseg, header=bshdr, name='BSEG')
    #init3 = fits.ImageHDU(data=fseg, header=fshdr, name='FSEG')
    
    newthing = fits.HDUList()
    for thing in (init0, init1, init2, init3, init4):
        newthing.append(thing)
    newthing.update_extend()
    newthing.writeto(outdir+fitsname, output_verify='silentfix', clobber=True)
    
    # clean up directory
    os.system("rm output/datacube/[c,s,a]*.fits")

    return 0


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
    








