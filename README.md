# How to use measure_morph.py

It measures galaxy morphology. 
It requires the following: `clean.py`, `galaxy_plot.py`, `run_sextractor.py`, and `utils.py`. All found in `morph/`. 
It has a BUNCH of dependences: `astropy`, SourceExtractor, other stuff, I'm sure. 
It is NOT parallelized so if you don't run this in chunks it will literally take FOREVER.

I should write a bash script which runs it in chunks FOR YOU... but I'm not going to do that just now.

To run it: `python measure_morph.py directory_to_stamps/ desired_catalog_name --outdir desired_name_for_output_directory/`

You can already run `python measure_morph.py -h` to see some options which probably don't make any sense anymore. 



Old shit: 

Various catalogs and what they contain:


bigsample_mycatv6.fits 
	contains: rpet, asymmetry (calculated my way with 1.5*rpet)

bigsample_mycatv7.fits
	contains: rpet, asymmetry (calculated CS's way with 1.5*rpet)

bigsample_mycat_1rpet_v8.fits
	contains: rpet, asymmetry (calculated CS's way with 1*rpet)

bigsample_mycat_v9.fits
	contains: rpet, asymmetry (as in 3), concentration (calculated with 
	elliptical annuli)

bigsample_mycat_v10.fits
	contains: rpet, asymmetry (as in 3), concentration (calculated with 
	circular annuli)

bigsample_mycat_v11_wG.fits
	contains: rpet and Gini only but Gini is all kinds of messed up

bigsample_mycat_v12.fits
	contains: rpet and Gini only and Gini is good!

bigsample_mycat_v13.fits
    contains: rpet and M20 only. M20 is not great -- lots of outliers; Mtot 
    centers need to be looked at again

i don't remember what 14-17 are...

bigsample_mycat_v18.dat
    contains all the updated parameters including two versions of Gini and M20
    (using both an elliptical aperture and the avg sb at 1 petrosian radius to 
    define the galaxy pixels)
    Except I think there was still a problem with Gini... one little bug

bigsample_mycat_v19.dat
    This is the latest and greatest. Everything in here is about as good as I 
    can get it! TWO versions of Gini and M20 and they both look great
    Asymmetry is still wonky but I don't know how to correct it and whether 
    that's where I should put my focus right now. 
    
Next Up:
1.  Test out the LLE using mycat_v19.dat. Need to see how it works and play 
    around with it a bit. Need to properly normalize all the parameters so that
    they are on equal footing.
2.  Make cutouts of the SDSS galaxies. These will be the training set for my codez
3.  Code up the MID statistics. 


