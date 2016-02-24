import os
import resource
import string
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from scipy.interpolate import interp1d
from collections import OrderedDict
from random import gauss
import pdb #"""for doing an IDL-like stop"""
from scipy.spatial import cKDTree
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
#import fast_ffts
import run_sextractor


def resource_getrusage():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    for name, desc in [
            ('ru_utime', 'User time'),
            ('ru_stime', 'System time'),
            ('ru_maxrss', 'Max. Resident Set Size'),
            ('ru_ixrss', 'Shared Memory Size'),
            ('ru_idrss', 'Unshared Memory Size'),
            ('ru_isrss', 'Stack Size'),
            ('ru_inblock', 'Block inputs'),
            ('ru_oublock', 'Block outputs')]:
        print '%-25s (%-10s) = %s' %(desc, name, getattr(usage, name))

def resource_getrlimits():
    
    for name, desc in [
            ('RLIMIT_CORE', 'core file size'),
            ('RLIMIT_CPU',  'CPU time'),
            ('RLIMIT_FSIZE', 'file size'),
            ('RLIMIT_DATA', 'heap size'),
            ('RLIMIT_STACK', 'stack size'),
            ('RLIMIT_RSS', 'resident set size'),
            ('RLIMIT_NPROC', 'number of processes'),
            ('RLIMIT_NOFILE', 'number of open files'),
            ('RLIMIT_MEMLOCK', 'lockable memory address')]:
        limit_num = getattr(resource, name)
        soft, hard = resource.getrlimit(limit_num)
        print 'Maximum %-25s (%-15s) : %20s %20s' % (desc, name, soft, hard)

def checkdir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)        

def get_interp(x, y, num=1000, kind='cubic'):
    newx = np.linspace(np.min(x), np.max(x), num=num)
    newy = interp1d(x, y, kind=kind)
    return newx, newy(newx)

def get_intersect(ydata, point, xvals, mono='inc'):
    '''
    I do this a lot: 
    data crosses a horizontal (y) line - at what x value does it intersect 
    that line?
    data: an array of y values
    point: y value defining the horizontal line that the data cross
    mono: data is assumed monotonic (it's really not) but is it 
          (inc)reasing or (dec)reasing monotonic?
    '''

    if mono=='inc':
        if np.any(point - ydata < 0):
            for idx, val in enumerate(point - ydata):
                if val < 0:
                    return np.mean([xvals[idx-1], xvals[idx]])
                    break
        else:
            print 'Data never cross the Horizontal Line.'
            print 'Cannot calculate the intersection.'
            return np.nan

    if mono=='dec':
        if np.any(ydata - point < 0):
            for idx, val in enumerate(ydata - point):
                if val < 0:
                    return np.mean([xvals[idx-1], xvals[idx]])
                    break
        else:
            print 'Data never cross the Horizontal Line.'
            print 'Cannot calculate the intersection.'
            return np.nan

def find_closest(point, listofpoints, k=1):
    
    # create a KDTree
    tree = cKDTree(listofpoints)
    # query the tree to find the set of coordinates closest to that point
    dist, index = tree.query(point, k)
    # return the value of the closest point
    return index, dist

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

def get_SB_Mask(Rp, Rp_SB, image, outname):
    '''
    Used to create a mask defining pixels belonging to a galaxy based on 
    the mean surface brightness at 1 Petrosian radius (Lotz 2004)

    Steps:
    1. convolve cleaned galaxy image with a Gaussian with sig=Rp/5
    2. measure the SB, mu, at Rp
    3. pixels in smoothed image with flux >= mu are assigned to the mask
    4. Return the mask
    '''
    im_center = [round(image.shape[0]/2), round(image.shape[1]/2)]

    conv = ndimage.gaussian_filter(image, sigma=Rp/5)
    if not np.any(conv > Rp_SB):
        return -1

    convimg = fits.ImageHDU(data=conv)
    convimg.writeto(outname+'_conv.fits', clobber=True)


    mask = np.array([True if conv[x] >= Rp_SB else False \
                     for x in np.ndindex(conv.shape)])
    mask = mask.reshape(conv.shape).astype('float')

    label_img, num_labels = ndimage.label(mask)
    
    # if there exists more than one object in the mask, we need to isolate the
    # correct one -- our object at the center
    if num_labels > 1:
        mask2 = np.array([True if label_img[x] == label_img[im_center[0], 
                                                            im_center[1]] \
                             else False for x in np.ndindex(conv.shape)])
        mask2 = mask2.reshape(conv.shape).astype('float')

        mm = fits.ImageHDU(data=mask2)
        mm.writeto(outname+'_mask.fits', clobber=True)
        return mask2
    else: 
        mm = fits.ImageHDU(data=mask)
        mm.writeto(outname+'_mask.fits', clobber=True)

        return mask


class MyEllipticalAperture(object):

    def __init__(self, xycenter, a, b, theta, data):
        self.x, self.y = xycenter
        self.a, self.b = a, b
        self.theta = theta
        self.data = data

        self.aper = self.make_aperture()
        self.phot = self.photometry()

    def make_aperture(self):
        cosang = np.cos(self.theta+np.pi/2.)
        sinang = np.sin(self.theta+np.pi/2.)
        xprime = np.array([(self.x-x)*cosang - (self.y-y)*sinang \
                           for x,y in np.ndindex(self.data.shape)])
        yprime = np.array([(self.x-x)*sinang + (self.y-y)*cosang \
                           for x,y in np.ndindex(self.data.shape)])

        ellipse = np.array([u**2/self.a**2 + v**2/self.b**2 \
                            for u,v in zip(xprime, yprime)])

        mask = np.array([True if ellipse[x] <= 1 else False \
                         for x in np.ndindex(ellipse.shape)])
        mask = mask.reshape(self.data.shape)

        return mask.astype('float')

    def plot(self, ax=None, fill=False,  **kwargs):
        
        kwargs['fill'] = fill

        if ax is None:
            ax = plt.gca()

        theta_deg = self.theta*180./np.pi
        patch = Ellipse((self.x, self.y), 2.*self.a, 2*self.b, 
                        theta_deg, **kwargs)
        ax.add_patch(patch)

    def photometry(self):
        #if self.asbval:
        return np.sum(np.abs(self.aper * self.data))
        #else:
        #    return np.sum(self.aper * self.data)

    def area(self):
        return len(self.aper[self.aper != 0.])

    @staticmethod
    def num_pix(weights):
        return np.sum(weights)    

    def run(self,longTable=False):
        '''Return:
        (skysub_aper_counts, aper_area, total_aper_counts,
         sky_area, total_sky_counts, med_sky_counts, std_sky_counts)
        '''
        aper = self.aper * self.data

        pdb.set_trace()
        #sky = self.skyaper * self.data
        #aper_area = self.radius*self.radius*np.pi
        #sky_area = self.skyradius**2*np.pi - aper_area

        return 0

class MyCircularAperture(object):

    def __init__(self, xycenter, r, data):
        self.x, self.y = xycenter
        self.r = r
        self.data = data

        self.aper = self.make_aperture()
        self.phot = self.photometry()

    def make_aperture(self):
        circle = np.array([(self.x-x)**2/self.r**2 + (self.y-y)**2/self.r**2 \
                           for x,y in np.ndindex(self.data.shape)])

        mask = np.array([True if circle[x] <= 1 else False \
                         for x in np.ndindex(circle.shape)])

        mask = mask.reshape(self.data.shape)

        return mask.astype('float')

    def plot(self, ax=None, fill=False,  **kwargs):
        
        kwargs['fill'] = fill

        if ax is None:
            ax = plt.gca()

        theta_deg = self.theta*180./np.pi
        patch = Circle((self.x, self.y), 2.*self.a, 2*self.b, 
                        theta_deg, **kwargs)
        ax.add_patch(patch)

    def photometry(self):
        return np.sum(np.abs(self.aper * self.data))

    def area(self):
        return len(self.aper[self.aper != 0.])



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
    








