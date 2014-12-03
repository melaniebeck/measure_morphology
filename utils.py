
import string
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
import sextutils        # eventually you need to write your own!
from collections import OrderedDict
from random import gauss
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
    from scipy.spatial import cKDTree
    
    # create a KDTree
    tree = cKDTree(listofpoints)
    # query the tree to find the set of coordinates closest to that point
    dist, index = tree.query(point)
    # return the value of the closest point
    return listofpoints[index], index


def clean_frame(img_name, config): #

    fitsname = os.path.basename(img_name)
    #basename = string.split(fitsname)
    catalogname = 'output/c_'+fitsname
    segmapname = 'output/s_'+fitsname
    apertname = 'output/a_'+fitsname
    
    outname = 'output/'+fitsname
            
    # READ IN CONFIG.SEX FILE
    cfile = sextutils.parseconfig_se(config)
    
    # TAILOR CONFIG.SEX TO PRODUCE UNIQUE OUTPUT FOR EACH OBJECT
    cfile['CATALOG_NAME'] = [catalogname]
    cfile['CHECKIMAGE_NAME'] = [segmapname, apertname]

    # WRITE THE UPDATED CONFIG.SEX TO FILE
    write_config(cfile, 't_'+config)
     
    # RUN SEXTRACTOR WITH UPDATED CONFIG.SEX
    cmdline = "sextractor "+img_name+" -c t_"+config
    rcode = os.system(cmdline)
    if (rcode):
        "SExtractor command [%s] failed." % cmdline
    
    # READ IN ORIG FITS-FILE AND SEG-MAP
    img, ihdr = pyfits.getdata(img_name, header=True)
    seg, shdr = pyfits.getdata(segmapname, header=True)

    # READ IN THE SE CATALOG
    cat, chdr = pyfits.getdata(catalogname, header=True)
    coords = zip(cat['X_IMAGE'],cat['Y_IMAGE'])
    
    # FIND GALAXY OF INTEREST (it will be in the CENTER OF IMAGE)
    center = [img.shape[0]/2., img.shape[1]/2.]
    
    # find the coordinate closest to the center
    myCoord = find_closest(center, coords)
    # find which index myCoord is at in the SE catalog
    myIndex = np.where(cat['X_IMAGE'] == myCoord[0])
    
    # GENERATE VALUES FOR REPLACING OBJECTS WITH BKG NOISE
    bkg_mean = np.mean(img[seg == 0])
    bkg_std = np.std(img[seg == 0])
    
    cln = img.copy()
    # CONSIDER EACH OBJECT IN THE CATALOG
    for num in cat['NUMBER']:
        # if the object is NOT our galaxy of interest then...
        if num != cat['NUMBER'][myIndex]:
            # select the pixels/array elements of object in the segmap
            # find corresponding pixels in image 
            # replace them with random bkg noise
            to_replace = np.where(seg == num)
            for idx in zip(to_replace[0], to_replace[1]):
                cln[idx] = gauss(bkg_mean, bkg_std)
        else:
            info = cat[myIndex]
    
    # SAVE ALL PRODUCTS TO DATA CUBE
    init0 = pyfits.ImageHDU(data=img, header=ihdr, name='ORG')
    init1 = pyfits.ImageHDU(data=cln, header=ihdr, name='CLN')
    init2 = pyfits.ImageHDU(data=seg, header=shdr, name='SEG')
    init3 = pyfits.TableHDU(data=cat, header=chdr, name='CAT')
    
    newthing = pyfits.HDUList()
    for thing in (init0, init1, init2, init3):
        newthing.append(thing)
    newthing.update_extend()
    newthing.writeto(outname, output_verify='silentfix', clobber=True)
    
    # clean up directory
    os.system("rm output/[c,s,a]_*.fits")
    
    #pdb.set_trace()

    # RETURN THE CLEANED IMAGE AND SE OUTPUT FOR OUR GAL
    return cln, info
    











