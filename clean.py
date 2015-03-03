#! usr/bin/env python

import os
import string
import numpy as  np
import pyfits as fits
from random import gauss
import run_sextractor
from utils import find_closest

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

