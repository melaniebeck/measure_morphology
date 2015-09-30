#! usr/bin/env python

'''
Various routines to create cutouts of the SDSS galaxies in the GZ2 sample. 
The filenames of the SDSS fields are different depending on whether the gal
is in a regular, shallow field or if it is in one of the deeper, stripe82
fields. 
'''


import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import bz2
import pdb
import warnings
import os

def makecutout(gal_catalog, hdulist, size):
    '''
    Function creates a stamp of a galaxy from a large field image.
    Stamp size is chosen as size*PETROR90_R radius from the SDSS catalog
    If this size is too large, default is set to 3*PETROR90_R
    
    Requires ancillary data including: 
        position of galaxy in pixel coordinates
        petrosian radius of the galaxy from SDSS
    Requires hdulist/FITS object of the entire SDSS field
        and it's header information
    '''
    img = hdulist[0].data
    hdr = hdulist[0].header

    rad = gal_catalog['PETROR90_R']/0.396 #[pixels]
    row, col = gal_catalog['ROWC_I'], gal_catalog['COLC_I']
    
    extents = [row-size*rad, row+size*rad, col-size*rad, col+size*rad]

    # test to see if extents are within img boundary
    if (extents[1] > hdr['NAXIS2']) or (extents[3] > hdr['NAXIS1']) or \
            np.any(extents < 0.):
        
        print "%i*Rp yields extents beyond image boundaries!"%size

        # if size*Rp doesn't work, default is 3*Rp 
        extents = [row-3*rad, row+3*rad, col-3*rad, col+3*rad]

        # which still might not work either ...  
        if (extents[1] > hdr['NAXIS2']) or (extents[3] > hdr['NAXIS1']) or \
            np.any(extents < 0.):
            
            # if that's the case, notify, flag, then exit the function
            print "3*Rp yields extents beyond image boundaries! \n"
            gal_catalog['bflag'] = 1

            return 0

    # otherwise, create the stamp!
    cutout = img[extents[0]:extents[1], extents[2]:extents[3]]

    # convert the array into a FITS image
    gal = fits.ImageHDU(data=cutout, header=hdr)

    return gal

def stripe82_only(gzdat, size, dirname):

    s82 = np.where((gzdat['RUN']==106) | (gzdat['RUN']==206))
    sub = gzdat[s82]
    print len(sub)

    outdir = 'SDSScutouts_4Rp_stripe82/'

    for obj in sub:
        rr, cc, ff = str(obj['RUN']), str(obj['CAMCOL']), str(obj['FIELD'])

        if not os.path.isfile(outdir+str(obj['OBJID'])+'.fits'):
            if rr =='106':
                rr = '100006'
            else:
                rr = '200006'
            
            name = dirname+'fpC-%s-i%s-%s.fit.gz'%(rr, cc, str(ff).zfill(4))
                
            try: 
                hdulist = fits.open(name)
                cutout = makecutout(obj, hdulist, size)
                if cutout:
                    outfile = outdir+str(obj['OBJID'])+'_4Rp.fits'
                    print "saving cutout "+outfile
                    cutout.writeto(outfile,output_verify='silentfix', 
                                   clobber=True)
            except:
                print name+" not found!"
                obj['imgflag_4Rp']=1

                sub.write('gz2sample_cutouts_stripe82.fits', overwrite=True)

def number_legit_stamps():
    # This little loop just goes through and looks to see which objects will
    # actually be found in the fields we have. 
    # It creates two columns in the catalog if an object's cutout is either not
    # created because its extents are outside the boundary of the field or if
    # the field isn't found (wasn't downloaded properly from SDSS)
    for obj in gzdat:
        rr, cc, ff = str(obj['RUN']), str(obj['CAMCOL']), str(obj['FIELD'])
        if rr not in ['106', '206']:
            name = dname+'frame-i-%s-%s-%s'%(str(rr).zfill(6),cc, 
                                             str(ff).zfill(4))
            if (os.path.isfile(name+'.fits.bz2') or \
                    os.path.isfile(name+'.fits')):
                rad = obj['PETROR90_R']/0.396 #[pixels]
                row, col = obj['ROWC_I'], obj['COLC_I']
                extents = [row-3*rad, row+3*rad, col-3*rad, col+3*rad]
                
                # test to see if extents are within img boundary
                if (extents[1] > 1498) or (extents[3] > 2048) or \
                        np.any(extents < 0.):
                    print "cutout extent beyond image boundaries!"
                    obj['bflag'] = 1
                else:
                    print name, "not found"
                    obj['imgflag']=1

    gzdat.write('gz2sample_cutouts.fits')


dname = 'SDSSimages/'
size = 4

gzdat = Table.read('gz2sample.fits')
gzdat['bflag_4Rp']=np.zeros_like(gzdat['OBJID'])
gzdat['imgflag_4Rp']=np.zeros_like(gzdat['OBJID'])

warnings.filterwarnings('ignore', message='Overwriting existing file .*',
                        module='fits')

stripe82_only(gzdat, size, 'SDSSimages_stripe82/')

# Instead of going through each object, you should go through each FRAME? 
# Ain't nobody got time to recode that though...

for obj in gzdat:
    # output directory for stamps
    outdir = 'SDSScutouts_4Rp/'

    # access RUN, CAMCOL, FIELD for this object
    rr, cc, ff = str(obj['RUN']), str(obj['CAMCOL']), str(obj['FIELD'])
    
    # why do I have this line? -- Now I remember: if script fails, we don't 
    # want to repeat stamps that have already been created
    if not os.path.isfile('%s%s_%iRp.fits'%(outdir, str(obj['OBJID']), size)):

        # if it's not in the deep stripe82...
        if rr not in ['106', '206']:

            # filenames for shallow fields follow this pattern...
            name = dname+'frame-i-%s-%s-%s'%(str(rr).zfill(6),cc, 
                                             str(ff).zfill(4))

            # but the way I downloaded them, they are two options for the suffix
            for suffix in '.fits.bz2', '.fits':
                try:
                    if suffix == '.fits.bz2':
                        imgfile = bz2.BZ2File(name+suffix)
                    else:
                        imgfile = name+suffix
                
                    # open the SDSS field in which our galaxy resides
                    hdulist = fits.open(imgfile)
            
                    # create/return the stamp
                    print "Creating stamp from", name+suffix
                    cutout = makecutout(obj, hdulist, size=size)
                
                    # if successful, cutout will be an array -- Save it
                    if cutout:
                        outfile = '%s%s_%iRp.fits'%(outdir, 
                                                    str(obj['OBJID']), size)
                        print "saving stamp "+outfile+'\n'
                        cutout.writeto(outfile, output_verify='silentfix', 
                                       clobber=True)
                    # if unsuccessful, cutout will be 0 -- Flag it
                    else:
                        obj['imgflag']=1    
                        print name+" not found!"

                # if one of the suffixes doesn't work -- don't sweat it.
                # The other one probably will. 
                except:
                    continue

        # if it is in the deep stripe82 field...
        else:
            # the RUN values are incorrect in the GZ2 catalog -- correct them
            if rr =='106':
                rr = '100006'
            else:
                rr = '200006'
            
            # filenames for deep stripe82 follow this pattern...
            name = dname+'fpC-%s-i%s-%s.fit.gz'%(rr, cc, str(ff).zfill(4))
    
            try: 
                # open the SDSS field in which the galaxy resides
                hdulist = fits.open(name)

                # create/return the cutout
                cutout = makecutout(obj, hdulist, size=size)

                if cutout:
                    outfile = '%s%s_%iRp.fits'%(outdir, str(obj['OBJID']), size)
                    print "saving cutout "+outfile
                    cutout.writeto(outfile,output_verify='silentfix', 
                                   clobber=True)
                else:
                    obj['imgflag_4Rp']=1    
                    print name+" not found!"                
            except:
                print name+" not found!"
                obj['imgflag_4Rp']=1

gzdat.write('gz2sample_cutouts.fits')
