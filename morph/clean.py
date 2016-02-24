#! usr/bin/env python

import os, glob
import string
import numpy as  np
import astropy.io.fits as fits
from random import gauss
import run_sextractor
from utils import find_closest
import pdb
import matplotlib.pyplot as plt

def clean_pixels(data, mask, segmap):
    #mean = np.mean(data[segmap == 0])
    #std = np.std(data[segmap == 0])
    med = np.median(data[segmap == 0])
    rms = np.sqrt(np.mean(np.square(data[segmap==0])))
    for pixel in zip(mask[0], mask[1]):
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

    For example: find galaxy closest to center that has a pixel area 
                 larger than fifty pixels**2

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

def savedata(mode, hdr, galnumber, data=[], names=[]):
    datacube=[]
    hdr.set('SE_MODE', mode, 'Sextractor mode used to clean image')
    hdr.set('SECATIDX', galnumber, 
             'Index in '+mode+' catalog denoting gal of interest')
    for datum, name in zip(data, names):
        if name == 'CAT':
            datacube.append(fits.TableHDU(data=datum, name=name))
        else:
            datacube.append(fits.ImageHDU(data=datum, name=name))

    return datacube

def clean_directory(outdir):
    for desc in ['*bright*', '*faint*', '*smooth*']:
        files = glob.glob(outdir+desc+".fits")
        for f in files:
            os.remove(f)

def clean_frame(image, outdir, sep=17., survey='SDSS'):
    '''
    This is a multi-stage cleaning process for each galaxy cutout.

    ---------------------------------------------------------------
    INPUTS: 
    ---------------------------------------------------------------
    
    image -- Name of galaxy cutout on which to perform cleaning. 
          Galaxy must be in the center of the image! 
          Image does not have to be square.

    outdir -- Name of the directory to which SExtractor output should
          be directed. 

    sep -- Minimum separation distance (in pixels) between:
          1. obj in BRIGHT and center of image
          2. obj in FAINT and center of image
          3. obj in FAINT and obj in BRIGHT
    ---------------------------------------------------------------
    OUTPUTS:
    ---------------------------------------------------------------
    All kinds of SExtractor nonsense

    ---------------------------------------------------------------
    RETURNS:
    ---------------------------------------------------------------
    Multiple flags to assess the quality of the cleaning. After cleaning, 
    SExtractor is run again in SMOOTH mode and various flags are set 
    depending upon:
    
    oFlag -- Flag set if no obj is detected in the center of the image
    
    uFlag -- Flag set if large obj are detected far from the center

    bFlag -- Set if a bright source is detected 

    category -- This number designates which particular course image
          took through the cleaning process
    ---------------------------------------------------------------
    PROCESS:
    ---------------------------------------------------------------
    Initially, SExtractor runs with both BRIGHT and FAINT parameters.
    Based on the resulting segmentation maps, the distance of the object 
    closest to the center in each map is determined (Fdist and Bdist). 
    The distance between these two objects is also calculated (DIST). 
    Based on the possible combinations of these values, various cleanings
    are performed. 
    #'''

    # initialize some shit for laterz
    category = 0
    Barea, B2dist = 0., 0.

    # commonly used SE keywords
    num = 'NUMBER'
    area = 'ISOAREA_IMAGE'
    x, y = 'X_IMAGE', 'Y_IMAGE'
    flux = 'FLUX_AUTO'

    if survey == "SDSS":
        configfile = 'se_params_SDSS.cfg'
    elif survey == 'COSMOS':
        configfile = 'se_params_COSMOS.cfg'


    basename = os.path.basename(os.path.splitext(image)[0])
    outname = outdir+basename
    catnames = [outname+'_bright_cat.fits', outname+'_faint_cat.fits', 
                outname+'_smooth_cat.fits']
    segnames = [outname+'_bright_seg.fits', outname+'_faint_seg.fits',
                outname+'_smooth_seg.fits']
    pdb.set_trace()

    # run SE in FAINT and BRIGHT modes
    if not run_sextractor.run_SE(image, 'BRIGHT', cfg_filename=configfile, 
                                 outdir=outdir):
        return [9, 9, 9, 9]

    run_sextractor.run_SE(image, 'FAINT', cfg_filename=configfile, 
                          outdir=outdir)

    # READ IN ORIG FITS-FILE and both Bright/Faint SEGMAPs
    img, ihdr = fits.getdata(image, header=True)
    cln = img.copy()
    bseg, fseg = fits.getdata(segnames[0]), fits.getdata(segnames[1]) 
    bcat, fcat = fits.getdata(catnames[0]), fits.getdata(catnames[1])

    center = [img.shape[0]/2., img.shape[1]/2.]

    # check to see if ANYTHING is found ANYWHERE
    if len(bcat) == 0 and len(fcat) == 0:
        return [9,9,9,9]

    brightdist, bFlag = 0., 0
    # Test for at least ONE object in BRIGHT catalog
    if len(bcat) > 0:
        # find object in BRIGHT closest to center
        BIndex, Bdist = find_closest(center,zip(bcat[x],bcat[y]))
        BCoord = zip(bcat[x], bcat[y])[BIndex]

        # If more than one obj, determine if bright source nearby
        if (len(bcat) > 1) and (BIndex != 0):
            if np.any(bcat[flux] > 1000.):
                bx, by = bcat[x][0], bcat[y][0]
                brightdist = np.sqrt((bx-center[0])**2+(by-center[0])**2)
                if brightdist < 150.:
                    bFlag = 1        


    # If nothing found in BRIGHT, assign Bdist to edge of image (251.)
    else: 
        Bdist = center[0]

    # Find closest object in FAINT
    FIndex, Fdist = find_closest(center, zip(fcat[x], fcat[y]))

    DIST = np.abs(Fdist - Bdist)

    if (DIST < sep):
        # FLAG 1: MOST COMMON CATEGORY --> CLEAN IN FAINT MODE
        if (Bdist < sep) & (Fdist < sep):
            #cln = clean_image(cln, bseg, bcat, BIndex, fseg)
            cln = clean_image(cln, fseg, fcat, FIndex, fseg)
            category, mode = 1, 'FAINT'

        # FLAG 2: CLEAN IN BRIGHT MODE & FLAG THESE!!
        elif (Bdist < sep) & (Fdist > sep):
            ''' There is only one obj in here:
            Two super bright galaxies super close together
            Seen as distinct objs in BRIGHT (but with crazy square edges)
            Seen as one blended obj in FAINT
            JUST FLAG THIS BITCH
            9/6/15 -- These statements ^^^ appear only to be true for COSMOS 
            data. Going to clean in Faint mode instead for SDSS
            '''
            #cln = clean_image(cln, bseg, bcat, BIndex, fseg)
            cln = clean_image(cln, fseg, fcat, FIndex, fseg)
            category, mode = 2, 'FAINT'

        # FLAG 3: CLEAN IN FAINT MODE
        elif (Bdist > sep) & (Fdist < sep):
            ''' There aren't many of these
            They're oddballs but most are well cleaned in FAINT
            '''
            cln = clean_image(cln, fseg, fcat, FIndex, fseg)
            category, mode = 3, 'FAINT'

        # FLAG 4: TWO STAGE CLEANING -- BRIGHT --> RUN SE AGAIN IN FAINT
        elif (Bdist > sep) & (Fdist > sep): 

            cln = clean_image(cln, bseg, bcat, BIndex, fseg)

            cln_sv = cln.copy()
            cln_sv = fits.ImageHDU(data=cln_sv, name='MID_CLN')
            cln_sv.writeto(outname+'_mid_cln.fits', output_verify='silentfix', 
                           clobber=True)

            # run SE again in FAINT
            run_sextractor.run_SE(outname+'_mid_cln.fits', 'FAINT', 
                                  cfg_filename=configfile, 
                                  outdir=outdir, outstr2='run2')
            f2seg = fits.getdata(outname+'_mid_cln_faint_run2_seg.fits')
            f2cat = fits.getdata(outname+'_mid_cln_faint_run2_cat.fits')
            coords = zip(f2cat[x], f2cat[y])

            # find closest obj to center with area above threshold
            Fdist, FIndex, FCoord, Farea, aFlag = \
                        closest_above_thresh(f2cat, area, center, coords, k=5)
            cln = clean_image(cln, f2seg, f2cat, FIndex, f2seg)
            category, mode = 4, 'FAINT2'
            
    else:
        # FLAG 5: TWO STAGE CLEANING - BRIGHT --> RUN SE AGAIN IN FAINT
        if (Bdist < sep) & (Fdist > sep):

            cln = clean_image(cln, bseg, bcat, BIndex, fseg)

            #save this image so that I can run SE on it
            cln_sv = cln.copy()
            cln_sv = fits.ImageHDU(data=cln_sv, name='MID_CLN')
            cln_sv.writeto(outname+'_mid_cln.fits', output_verify='silentfix',
                           clobber=True)
            
            # run SE again in FAINT
            run_sextractor.run_SE(outname+'_mid_cln.fits', 'FAINT', 
                                  cfg_filename=configfile, 
                                  outdir=outdir, outstr2='run2')
            f2seg = fits.getdata(outname+'_mid_cln_faint_run2_seg.fits')
            f2cat = fits.getdata(outname+'_mid_cln_faint_run2_cat.fits')
            coords = zip(f2cat[x], f2cat[y])

            # find closest obj to center with area above threshold
            Fdist, FIndex, FCoord, Farea, aFlag = \
                        closest_above_thresh(f2cat, area, center, coords, k=5)
            cln = clean_image(cln, f2seg, f2cat, FIndex, f2seg)
            category, mode = 5, 'FAINT2'
 
        # FLAG 6: CLEAN IN SMOOTH MODE
        elif (Bdist > sep) & (Fdist < sep):
            ''' These are mostly faint objects not detected in BRIGHT
            run SE in SMOOTH mode and then clean
            '''
            run_sextractor.run_SE(image, 'SMOOTH', cfg_filename=configfile, 
                                  outdir=outdir)
            sseg, scat = fits.getdata(segnames[2]), fits.getdata(catnames[2])
            SIndex, Sdist = find_closest(center, zip(scat[x], scat[y]))

            cln = clean_image(cln, sseg, scat, SIndex, sseg)
            category, mode = 6, 'SMOOTH'

        # FLAG 7: CLEAN IN FAINT MODE -- ALL GARBAGE ANYWAY
        elif (Bdist > sep) & (Fdist > sep):
            ''' this is mostly a garbage bin of crap 
            any object in here needs to be flagged and is likely not a true
            galaxy at all!
            '''
            cln = clean_image(cln, fseg, fcat, FIndex, fseg)
            category, mode = 7, 'FAINT'

        # FLAG 8: 
        elif  (Bdist < sep) & (Fdist < sep):
            pdb.set_trace()
            cln = clean_image(cln, bseg, bcat, BIndex, fseg)
            category, mode = 8, 'BRIGHT'

    # Save all major data products
    if mode == 'BRIGHT':
        datacube = savedata(mode, ihdr, BIndex, data=[bseg, fseg, bcat], 
                            names=['BSEG', 'FSEG', 'CAT'])

    elif mode == 'FAINT':
        datacube = savedata(mode, ihdr, FIndex, data=[bseg, fseg, fcat], 
                            names=['BSEG', 'FSEG', 'CAT'])

    elif mode == 'FAINT2':
        datacube = savedata('TWOSTAGE', ihdr, FIndex, 
                            data=[bseg, fseg, f2seg, f2cat], 
                            names=['BSEG', 'FSEG', 'F2SEG', 'CAT'])
        datacube.insert(0, cln_sv)

    elif mode == 'SMOOTH':
        datacube = savedata(mode, ihdr, SIndex, data=[bseg, fseg, sseg, scat], 
                            names=['BSEG', 'FSEG', 'SSEG', 'CAT'])

        
    # SAVE ALL PRODUCTS TO DATA CUBE
    datacube.insert(0, fits.ImageHDU(data=cln, header=ihdr, name='CLN'))
    datacube.insert(1, fits.ImageHDU(data=img, header=ihdr, name='ORG'))

    newthing = fits.HDUList()
    for thing in datacube: 
        newthing.append(thing)
    newthing.update_extend()
    newthing.writeto(outdir+'f_'+basename+'.fits', output_verify='silentfix', 
                     clobber=True)


    # Now that we've done the cleaning -- Let's test it!    
    run_sextractor.run_SE(outdir+'f_'+basename+'.fits', 'SMOOTH', 
                          cfg_filename=configfile, outdir=outdir, 
                          outstr2='test')
    tseg = fits.getdata(outdir+'f_'+basename+'_smooth_test_seg.fits')
    tcat = fits.getdata(outdir+'f_'+basename+'_smooth_test_cat.fits')

    # If we find obj near the center is too small then we overcleaned it
    uFlag, oFlag = 0, 0
    if tcat.size:
        coords = zip(tcat[x], tcat[y])
        index, dist = find_closest(center, coords, k=10)
        objarea = tcat[area][index[0]]
        areas = tcat[area][np.where(tcat[area] != objarea)]

        if (dist[0] < sep): 
            if (objarea < 50.): 
                print 'OVERCLEANED!!!'
                oFlag = 1

            # If we find large objs far from the center- didn't clean enough
            if np.any(areas > 200.):       
                print 'UNDER CLEANED!!'
                uFlag = 1
                
                cln = clean_image(cln, tseg, tcat, index[0], tseg)
                data = fits.open(outdir+'f_'+basename+'.fits')
                data.insert(0,fits.ImageHDU(data=cln, name='UCLN'))
                data['UCLN'].header.set('SECATIDX', index[0], 'Index in SECAT')
                data['CAT'].data = tcat
                data['FSEG'].data = tseg
                data.writeto(outdir+'f_'+basename+'.fits', clobber=True,
                             output_verify='silentfix')
        else:
            oFlag = 1
    else:
        # if we don't find anything in the center anymore, something is wrong
        # either overcleaned or blended into nearby object
        oFlag = 1

    # I want to save the cleaned image separate -- 9/5/15 
    # get rid of this later, probably?
    cleanedup = fits.ImageHDU(data=cln, header=ihdr, name='CLN%i'%category)
    cleanedup.writeto(outdir+'f_'+basename+'_clnonly.fits', 
                      output_verify='silentfix', clobber=True)

    # clean up directory
    clean_directory(outdir) # doesn't work?! Harrumph.


    #FIndex, Fdist, Bdist, DIST, Farea, Barea,
    return  [category, oFlag, uFlag, bFlag]
