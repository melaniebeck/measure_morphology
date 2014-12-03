#Asymmetry module
#measures the asymmetry of a galaxy on a "postage stamp" sized image

import os
import argparse
#import shutil
import string
import pyfits
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.spatial import cKDTree
from scipy.interpolate import interp1d
import sextutils        # eventually you need to write your own!
from collections import OrderedDict
from random import gauss
from photutils import aperture_photometry, EllipticalAnnulus
#from pylab import *
#import subprocess
#import matplotlib as mpl 
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from math import pi
import pdb #"""for doing an IDL-like stop"""


def write_config(cdict, outname):  
    ''' write config.sex out to file '''
    f = open(outname,'w')
    for k,v in cdict.iteritems():
        f.write('%s\t'%k)
        v = ','.join(v)
        f.write('%s\n'%v)
    f.close()

def find_closest(point, listofpoints):
    # create a KDTree
    tree = cKDTree(listofpoints)
    # query the tree to find the set of coordinates closest to that point
    dist, index = tree.query(point)
    # return the value of the closest point
    return listofpoints[index], index


def clean_frame(img_name, configs):
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

    fitsname = os.path.basename(img_name)

    bcatalogname = 'output/cB_'+fitsname     #bright cat name
    fcatalogname = 'output/cF_'+fitsname     #faint cat name
    bsegmapname = 'output/sB_'+fitsname
    fsegmapname = 'output/sF_'+fitsname    
#    apertname = 'output/a_'+fitsname

    catnames = [bcatalogname, fcatalogname] #
    segnames = [bsegmapname, fsegmapname] #
    
    outname = 'output/'+fitsname
    
    for config, cat, seg in zip(configs, catnames, segnames):
        # READ IN CONFIG.SEX FILE -- FAINT FIRST
        cfile = sextutils.parseconfig_se(config)
           
        # TAILOR CONFIG.SEX TO PRODUCE UNIQUE OUTPUT FOR EACH OBJECT
        cfile['CATALOG_NAME'] = [cat]
        cfile['CHECKIMAGE_NAME'] = [seg]
        
        # WRITE THE UPDATED CONFIG.SEX TO FILE
        write_config(cfile, 't_'+config)
             
        # RUN SEXTRACTOR WITH UPDATED CONFIG.SEX
        cmdline = "sextractor "+img_name+" -c t_"+config
        rcode = os.system(cmdline)
        if (rcode):
            "SExtractor command [%s] failed." % cmdline
    
    # READ IN ORIG FITS-FILE AND both Bright/Faint SEG-MAPs
    img, ihdr = pyfits.getdata(img_name, header=True)
    fseg, fshdr = pyfits.getdata(fsegmapname, header=True)
    bseg, bshdr = pyfits.getdata(bsegmapname, header=True)

    # READ IN THE SE CATALOGs (Bright & Faint)
    bcat, bchdr = pyfits.getdata(bcatalogname, header=True)
    fcat, fchdr = pyfits.getdata(fcatalogname, header=True)
      

    # FIND GALAXY OF INTEREST (it will be in the CENTER OF IMAGE)
    center = [img.shape[0]/2., img.shape[1]/2.]
    
    # find the coordinate closest to the center (Faint Catalog)
    myFCoord, myFIndex = find_closest(center, 
                            zip(fcat['X_IMAGE'],fcat['Y_IMAGE']))

    # find corresponding object in the Bright catalog
    myBCoord, myBIndex = find_closest(center, 
                            zip(bcat['X_IMAGE'], bcat['Y_IMAGE']))

    #compare objects from the Bright catalog to those in the Faint catalog:
    # for each obj in Bright, compare against everything in Faint
    # if, for the closest match, the separation is greater than 10 pixels (.3'')
    # flag this object (save its index in a list for cleaning later)
    
    toclean = []
    tree = cKDTree(zip(fcat['X_IMAGE'], fcat['Y_IMAGE']))
    for num in bcat['NUMBER']:
        dist,index = tree.query([bcat['X_IMAGE'][num-1], bcat['Y_IMAGE'][num-1]])
        if (dist > 10.) and (num != bcat['NUMBER'][myBIndex]): 
            toclean.append(num)
    

    # GENERATE VALUES FOR REPLACING OBJECTS WITH BKG NOISE
    bkg_mean = np.mean(img[fseg == 0])
    bkg_std = np.std(img[fseg == 0])
    
    cln = img.copy()
    # Clean using the FAINT catalog first!
    # CONSIDER EACH OBJECT IN THE Faint CATALOG
    for num in fcat['NUMBER']:
        # if the object is NOT our galaxy of interest then...
        if num != fcat['NUMBER'][myFIndex]:
            # select the pixels/array elements of object in the segmap
            # find corresponding pixels in image 
            # replace them with random bkg noise
            to_replace = np.where(fseg == num)
            for idx in zip(to_replace[0], to_replace[1]):
                cln[idx] = gauss(bkg_mean, bkg_std)
    
    for num in toclean:
        to_replace = np.where(bseg == num)
        for idx in zip(to_replace[0], to_replace[1]):
            cln[idx] = gauss(bkg_mean, bkg_std)
    

    # set header to reflect which index our galaxy will be found in
    bchdr.set('SECATIDX', myBIndex, 'Index in SE catalog denoting gal of interest')
    fchdr.set('SECATIDX', myFIndex, 'Index in SE catalog denoting gal of interest')
    
    # SAVE ALL PRODUCTS TO DATA CUBE (save info from BRIGHT catalog)
    init0 = pyfits.ImageHDU(data=img, header=ihdr, name='ORG')
    init1 = pyfits.ImageHDU(data=cln, header=ihdr, name='CLN')
    '''
    init2 = pyfits.ImageHDU(data=fseg, header=fshdr, name='SEG')
    init3 = pyfits.TableHDU(data=fcat, header=fchdr, name='CAT')
    '''
    init2 = pyfits.ImageHDU(data=bseg, header=bshdr, name='BSEG')
    init3 = pyfits.ImageHDU(data=fseg, header=fshdr, name='FSEG')
    init4 = pyfits.TableHDU(data=bcat, header=bchdr, name='BCAT')
    
    newthing = pyfits.HDUList()
    for thing in (init0, init1, init2, init3, init4):
        newthing.append(thing)
    newthing.update_extend()
    newthing.writeto(outname, output_verify='silentfix', clobber=True)
    
    # clean up directory
    os.system("rm output/[c,s,a]*.fits")
    
    #pdb.set_trace()

    # RETURN THE CLEANED IMAGE AND SE OUTPUT FOR OUR GAL
    return cln, fcat
    
def petro_radius(img, img_name, info):
    r_flag = False
    imgsize = img.shape[0]/2.
    # to measure the petrosian radius we need to find the surface brightness
    # profile
    position = [info['X_IMAGE'][0], info['Y_IMAGE'][0]]
    theta = info['THETA_IMAGE']*pi/180. # radians
    se_apsize = info['A_IMAGE']*info['KRON_RADIUS']
    
    a_out = se_apsize*np.logspace(-1.0, np.log10(imgsize/se_apsize), num=12)
    # check radii for appropriate values
    for idx, radius in enumerate(a_out):
        # radii should be no larger than 3*se_aperture
        #if radius > 3*se_apsize:
        #    a_out[idx] = 3*se_apsize
        # radii should not be larger than the image size -- flag if true
        if radius > imgsize:
            a_out[idx] = imgsize
            r_flag = True
            # print some sort of error

    a_in = np.zeros(len(a_out))
    a_in[1:len(a_out)] = a_out[0:len(a_out)-1]
    b_out = a_out/info['ELONGATION']

    #pdb.set_trace()    
    flux = []
    ans = []
    for a0_in, a0_out, b0_out in zip(a_in, a_out, b_out):
        an = EllipticalAnnulus(position,a0_in, a0_out, b0_out, theta)
        ans.append(an) 
        flux.append(aperture_photometry(img, an))
    
    phot_table = np.hstack(flux)
    counts = phot_table['aperture_sum']
    annuli = np.hstack(ans)
        
    # plot the annuli overlaid on cleaned image
    imgplot = plt.imshow(img, cmap='gray_r', origin='lower')
    imgplot.set_clim(-0.009, 0.022)
    areas = []
    for an in annuli:
        an.plot(color='blue', lw=1.5, alpha=0.5)
        areas.append(an.area())         
    #plt.show()
    fitsname = os.path.basename(img_name)
    name = string.split(fitsname)
    plt.savefig('output/sbaps/sbaps_'+name[0]+'.png', bbox_inches='tight')
    plt.clf()
    
    # surface brightness as a function of radius
    sb = counts/areas
    
    # calculate the average surface brightness
    c_sum = []
    a_sum = []
    for idx in range(1,len(counts)+1):
        c_sum.append(sum(counts[0:idx]))
    for idx in range(1,len(areas)+1):
        a_sum.append(sum(areas[0:idx]))
    c_sum = np.array(c_sum)
    a_sum = np.array(a_sum)

    avg_sb = c_sum/a_sum  
      
    # plot the surface brightness
    fig = plt.plot(a_out, counts/areas, 'ro') # counts/pixel**2
    # plot the average SB
    plt.plot(a_out, c_sum/a_sum, 'g^')
    # plot SB/avg_SB
    plt.plot(a_out,sb/avg_sb, 'b--')
    plt.xscale('log')
    plt.yscale('log')
    
    #plt.show()
    plt.clf()
    
    zeros = np.zeros(len(a_out))
    eta = []
    for x in range(0,len(a_out)):
        eta.append(0.2)   
   
    xlims = [a_out[0] - 0.1*a_out[0], a_out[len(a_out)-1]]
    ylims = [-.1*max(avg_sb), max(avg_sb) + 0.3*max(avg_sb)]

    plt.subplot(211)
    plt.title('Surface Brightness Profile')
    plt.plot(a_out, counts/areas, 'ro', label='Surface Brightness')
    plt.plot(a_out, c_sum/a_sum, 'g^', label='Average SB')
    plt.plot(a_out, zeros, 'k--')
    plt.xscale('log')
    plt.axis([xlims[0], xlims[1], ylims[0], ylims[1]])
    #plt.xlim(xlims)    
    legend = plt.legend(loc='upper right')
    legend.get_frame()
    
    plt.subplot(212)
    plt.title('u(R)/<u(R)>')
    plt.plot(a_out, sb/avg_sb, 'bo', label='SB/Avg_SB')
    plt.plot(a_out, eta, 'k--', label='eta=0.2')  
    plt.xscale('log')
    plt.axis([xlims[0], xlims[1], 0., 1.1])
    #plt.xlim(xlims)
    plt.legend(loc='upper right')

    plt.savefig('output/sbpro_'+name[0]+'.png', bbox_inches='tight')
    #plt.show()
    plt.clf()

    #pdb.set_trace()

    pet_ratio = sb/avg_sb
    # now we need to find the intersection of u(R)/<u(R)> with 0.2
    f1 = interp1d(pet_ratio,a_out, bounds_error=False, assume_sorted=False)
    f2 = interp1d(pet_ratio,a_out, bounds_error=False)
    #f3 = interp1d(pet_ratio,a_out)
    #pdb.set_trace()
    r1 = f1(0.2)
    r2 = f2(0.2)
    #r3 = f3(0.2)
    #pdb.set_trace()
    rad = r2.item()
     
    return rad
    
def find_center(im):


    return center


def asymmetry(img, se_output, petro_rad, angle):
    # need to define galaxy as those pixels within the petrosian radius
    # use the center of that to rotate ANGLE around
    
    

    return asym




def main():
    # this will provide basic documentation for this function as well as accept
    # a filename containing the list of galaxy fits files on which to compute 
    # the asymmetry
    parser = argparse.ArgumentParser(description='Measure asymmetry of galaxies')
    parser.add_argument('filename', type=str, 
        help='Input a text file containing list of fits files for computation of asymmetry')
    #parser.add_argument('SEconfig', type=str,
    #    help='Input the config.sex file for sextractor configuration parameters')
    args = parser.parse_args()
   
    #pdb.set_trace()
    
    # read in the list and store the names of the fits files into fnames
    fnames = np.genfromtxt(args.filename, dtype=(str))
    
    # define angle to rotate about (degrees)
    angle = 180

    prad = []
    for name in fnames:
        cln_img, se_output = clean_frame(name,['config_bright.sex','config_faint.sex']) #args.SEconfig
        #pdb.set_trace()
        
        #prad.append(petro_radius(cln_img, name, se_output))
        
    '''
    with open('prad.txt','w') as f:
        for thing in prad:
            f.write("%f\n" %thing) 
   
    test = Table.read('bright_gals.fits')
    rpet = test['rpet']
    line = np.arange(20, 300,dtype=float)
    
    
    fig = plt.plot(prad, rpet,'ro', line, line, 'k--', xlabel='asdfasd')
    plt.show()
    pdb.set_trace()
        
        
        #center = find_center(cln_img)
        
        #asym = asymmetry(cln_img, angle)
    #'''

if __name__ == '__main__':
    main()    














