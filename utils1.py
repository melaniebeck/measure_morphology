import os
import string
import numpy as np
import pyfits as fits
from astropy.table import Table
from scipy.interpolate import interp1d
import sextutils        # eventually you need to write your own!
from collections import OrderedDict
from random import gauss
from math import pi
import pdb #"""for doing an IDL-like stop"""
from scipy.spatial import cKDTree

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


def clean_frame(img_name): #, configs
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
    configs = ['config_bright.sex','config_faint.sex']
    
    outname = 'output/f_'+fitsname
    
    for config, cat, seg in zip(configs, catnames, segnames):
        # READ IN CONFIG.SEX FILE -- FAINT FIRST
        cfile = sextutils.parseconfig_se(config)
           
        # TAILOR CONFIG.SEX TO PRODUCE UNIQUE OUTPUT FOR EACH OBJECT
        cfile['CATALOG_NAME'], cfile['CHECKIMAGE_NAME'] = [cat], [seg]
         
        # WRITE THE UPDATED CONFIG.SEX TO FILE
        write_config(cfile, 't_'+config)
             
        # RUN SEXTRACTOR WITH UPDATED CONFIG.SEX
        cmdline = "sextractor "+img_name+" -c t_"+config
        rcode = os.system(cmdline)
        if (rcode):
            "SExtractor command [%s] failed." % cmdline

    # READ IN ORIG FITS-FILE AND both Bright/Faint SEG-MAPs
    img, ihdr = fits.getdata(img_name, header=True)
    fseg, fshdr = fits.getdata(fsegmapname, header=True)
    bseg, bshdr = fits.getdata(bsegmapname, header=True)

    # READ IN THE SE CATALOGs (Bright & Faint)
    bcat, bchdr = fits.getdata(bcatalogname, header=True)
    fcat, fchdr = fits.getdata(fcatalogname, header=True)
    
    '''
    Here is the problem: 
        nearly every galaxy is found in FAINT but not all are found in BRIGHT
        need to check to see if our galaxy is found in BRIGHT and in FAINT
        IF in both, proceed as below
        IF only in FAINT, then we need to smooth this galaxy before proceding
    #'''

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
    init0 = fits.ImageHDU(data=img, header=ihdr, name='ORG')
    init1 = fits.ImageHDU(data=cln, header=ihdr, name='CLN')
    '''
    init2 = fits.ImageHDU(data=fseg, header=fshdr, name='SEG')
    init3 = fits.TableHDU(data=fcat, header=fchdr, name='CAT')
    '''
    init2 = fits.ImageHDU(data=bseg, header=bshdr, name='BSEG')
    init3 = fits.ImageHDU(data=fseg, header=fshdr, name='FSEG')
    init4 = fits.TableHDU(data=bcat, header=bchdr, name='BCAT')
    
    newthing = fits.HDUList()
    for thing in (init0, init1, init2, init3, init4):
        newthing.append(thing)
    newthing.update_extend()
    newthing.writeto(outname, output_verify='silentfix', clobber=True)
    
    # clean up directory
    os.system("rm output/[c,s,a]*.fits")
    img.close()
    
    #pdb.set_trace()

    # RETURN THE CLEANED IMAGE AND SE OUTPUT FOR OUR GAL
    return cln, fcat










