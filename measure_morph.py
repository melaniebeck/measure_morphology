import re, glob, os, string, pdb
import argparse, warnings

from astropy.table import Table
import astropy.io.fits as fits
import numpy as np

import morph


####################### main ############################

def main():
    
    parser = argparse.ArgumentParser(description='Perform LLE/PCA/whatevs')
    parser.add_argument('-d', dest="directory", type=str, 
        help='Directory of fits images on which to run LLE.')
    parser.add_argument('-c', dest="catalog_name", type=str,
        help='Specify the desired name for output catalog.')
    parser.add_argument('--outdir', type=str, default='output/datacube/', 
        help='Specify the desired name for output directory.')
    args = parser.parse_args()


    # Select all FITS files in the given directory
    fitsfiles = np.array(sorted(glob.glob(args.directory+"/*")))

    # There are a lot of useless warnings that pop up -- suppress them!
    warnings.filterwarnings('ignore', message='Overwriting existing file .*',
                            module='pyfits')

    # If the proper output directory doesn't exist, create it
    morph.checkdir(args.outdir+'datacube/')

    # The morphology catalog is built by going through each FITS image, 
    # cleaning it, processing it, measuring morphological parameters, and
    # making appropriate plots FOR EACH FITS IMAGE. The morphological 
    # parameters are appended to a huge table which is updated after each 
    # image. All the output (Cleaned FITS files, plots and figures) are 
    # stored in appropriately named directories which are created on the fly
    # and nestled within the specified subdirectory.
    #
    # Sometimes this shit crashes because I'm a terrible programmer.
    # In that case, I remove the offending FITS image, putting it in the 
    # "bad_cutouts" directory and then start up the code again.
    #
    # When the code is rerun, we don't want to start a new morph catalog - 
    # we want to keep appending to the old one so we test for this...


    # check to see if a current morphology catalog exists
    # if so -- find where cursor should go by getting the length
    try:
        t = Table.read(args.catalog_name)
        counter = len(t)
        # this "cursor" tells us which FITS file we should start at because, 
        # presumably, all previous images have already been processed.
        fitsfiles = fitsfiles[counter::]
    except:
        counter = 0

    #tt = Table(names=('category', 'oFlag', 'uFlag', 'bFlag'))

    # Now that we have our list of FITS to process...
    for idx, f in enumerate(fitsfiles): 
        basename = os.path.basename(f)

        # If we're cleaning the stamp, we tack an "f" in front of
        # the original filename
        filename = args.outdir+'f_'+basename

        # Otherwise, just pull up the regular filename
        filename = args.directory+basename

        #if not os.path.isfile(filename):
        #print "File not found! Running SExtractor before proceeding."
        #print "Cleaning ", os.path.basename(f)
        #flags = morph.clean_frame(f, args.outdir+'datacube/', sep=4, 
        #                          survey='SDSS')

        #tt.add_row((flags[0], flags[1], flags[2], flags[3]))

        #"""
        # check to see if SExtractor failed
        #if np.any(np.array(flags)-9 < 0):
        if True:

            print "Running", os.path.basename(f)
            hdulist = fits.open(filename, memmap=True)

            print hdulist
            # Measure galaxy morphologies
            g = morph.GalaxyMorphology(hdulist, filename, flags, args.outdir)

            # Plot galaxy figures for quality control
            if not np.isnan(g.Rp):
                morph.galaxyPlots(g, hdulist)


            # If this is the first galaxy; create a Table
            if (idx == 0) and (counter == 0):
                t, gal_dict = g.table(init=True)
                t.add_row(gal_dict)
            # Otherwise, append galaxy morphologies to existing Table
            else:
                t.add_row(g.table())
            
            # Close any remaining FITS files
            hdulist.close()
            
            # Write the file after each galaxy (costly but don't want to lose
            # data when my shit crashes, which it inevitably does
            t.write(args.catalog_name, overwrite=True)
            print counter+1," galaxies measured!"
            counter+=1
            del g
        else:
            print "SExtractor failed on "+basename
            #if (idx == 0) and (counter == 0):
            #t, gal_dict = g.table(init=True)
            #t.add_row(gal_dict)
        #"""

    t.write(args.directory+'/test.fits')
    print "Morphological parameter catalog complete.\n"
    exit()  


if __name__ == '__main__':
    main()
    
