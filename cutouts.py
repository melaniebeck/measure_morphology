#! usr/bin/env python

import numpy as np
import pyfits as fits
from astropy.table import Table
import matplotlib.pyplot as plt
import pdb
import warnings

dname = 'SDSS_images/'

gzdat = Table.read('gz2subsample.fits')
gzdat['bflag']=np.zeros_like(gzdat['OBJID'])
gzdat['imgflag']=np.zeros_like(gzdat['OBJID'])

warnings.filterwarnings('ignore', message='Overwriting existing file .*',
                        module='pyfits')

for obj in gzdat:
    # open the sdss image for this obj using RUN, CAMCOL, FIELD
    r, c, f = str(obj['RUN']), str(obj['CAMCOL']), str(obj['FIELD'])
    if r not in ['106', '206']:
        name = dname+'frame-i-%s-%s-%s.fits'%(str(r).zfill(6),c, 
                                              str(f).zfill(4))
        outdir = 'cutouts/'
    else:
        if r =='106':
            r = '100006'
        else:
            r = '200006'
        name = dname+'fpC-%s-i%s-%s.fit'%(r, c, str(f).zfill(4))
        outdir = 'cutouts/stripe82/'

    try:
        img, hdr = fits.getdata(name, header=True)

        rad = obj['PETROR90_R']/0.396 #[pixels]
        x, y = obj['ROWC_I'], obj['COLC_I']
        
        extents = [x-3*rad, x+3*rad, y-3*rad, y+3*rad]
        # test to see if extents are within img boundary
        if (extents[1] > hdr['NAXIS1']) or (extents[3] > hdr['NAXIS2']) or \
           np.any(extents < 0.):
            print "cutout extent beyond image boundaries!"
            obj['bflag'] = 1
        else:
            # cut out galaxy
            cutout = img[extents[0]:extents[1], extents[2]:extents[3]]
            
            gal = fits.ImageHDU(data=cutout)
            gal.writeto(outdir+str(obj['OBJID'])+'.fits', 
                        output_verify='silentfix', clobber=True)

    except:
        print name+" not found!"
        obj['imgflag']=1


        
