### Measuring galaxy morphological diagnostics 

Comes in two flavors: `MEASURE_MORPHOLOGY_parallel.py` and `measure_morph.py`. 
The former is parallized but does not call the image cleaning routine. The latter is not parallelized but does call the cleaning routine. 

Both script require the code found in the `morph/` directory, specifically `galaxyMorphology.py` which is an object that then calls a suite of morphological diagnostics to be measured on the galaxy image. It's not a very elegant design but it gets the job done. 

Other scripts of interest in the `morph/` directory include
* `clean.py`: takes a galaxy postage stamps, runs sextractor (`run_sextractory.py`) and cleans the image of light not belonging to the galaxy of interest, assumed to be at the center of the image. Returns a list of flags that are supposed to give an estimate of how well the cleaning was performed but they don't quite work as intended.
* `galaxyPlots.py`: if certain flags are set in `galaxyMorphology.py`, various diagnostic figures will be created for each galaxy that is processed. I usually turn this off and just call individual plotting functions after the fact. 
* `utils.py`: contains various functions needed for the other scripts to run. 



### Directory contents and structure

* `SDSSimages/` contains all the full SDSS fields
* `SDSScutouts_*/` contains cutouts for all individual galaxies using a box radius of 3Rp or 4Rp as designated. 
* `bad_cutouts/` contains postage stamps that failed various steps of postage stamp making, cleaning, or morphology measuring. I never went through them individually. 
* `output_*/` contains the output of `measure_morph.py` for either the 3Rp or 4Rp postage stamps. This output includes subdirectories:
	* `datacube/`: during cleaning I create a FITS cube containing the ORIG postage stamp, the BRIGHT and FAINT segmentation maps, and the resulting CLN postage stamp; also contains BRIGHT and FAINT catalogs for each postage stamp (sold separately). 
	* `asymimgs/`: while calculating asymmetry a difference image needs to be created and these are stored here
	* `masks/`: when I was experimenting with various ways to measure Gini I tried creating SB masks; these are stored here. 
	* `figures/`: if figures are created during morphology measurement, they are stored here
* `SDSSmorphology_catalogs/` contains all the resulting catalogs created during morphology measurement. I typically measured morphology on "chunks" of galaxies and combine the resulting catalogs into a master catalog in here. 
* `sexfiles/` contains various files necessary for SExtractor to run including filters which were occasionally used during the cleaning process. 
* `remeasure_ell_morph/` contains all the asymmetry difference images created when I re-did the morphology catalog to correct the elliptical apertures. 


To run `measure_morph.py`: `python measure_morph.py directory_to_stamps/ desired_catalog_name --outdir desired_name_for_output_directory/`

You can already run `python measure_morph.py -h` to see some options which probably don't make any sense anymore. 


### Morphology catalog descripion

'A': asymmetry
'Ax','Ay': asymmetry center
'C': concentration index
'G': gini coefficient
'G2': gini coefficient measured by using a SB mask (do not trust)
'G_c': gini coefficient measured using CIRCULAR aperture
'M20': M20
'Mlevel1': minimum pixel value defining brightest 20% of pixels
'Mx', 'My': M20 center
'Rp': Petrosian radius (elliptical aperture -- correct this time!)
'r20',
'r50', 
'r80', 
'Rp_SB': surface brightness as 1 Rp
'Rpflag': >0 indicates error in Rp measurement
'a': SExtractor semi major axis
'b': Sextractor semi minor axis
'theta': SExtractor position angle of ellipse
'e': SExtractor ellongation
'elipt':SExtractor ellipticity
'kron': SExtractor kron radius
'name': name of image datacube, of form "f_[SDSSid]_4Rp.fits"
'objid': SDSS object ID
'outdir': directory for all output files
'dec', 'ra': SExtractor coordinates 
'med': Median of background pixels as determined by SExtractor segmaps
'rms': RMS of background pixels as determined by SExtractor segmaps 
'stn': standard dev of background pixels as determined by SExtractor segmaps
'x','y': SDSS galaxy center
'xc', 'yc': rounded galaxy center
'Rp_corr': same as Rp (this is so that when matched with other catalogs it has a distinct column name for comparison with the uncorrected measurements)
'C_corr': same as C
'A_corr': same as A
'G_corr': same as G
'M20_corr': same as M20
'cat': "cleaning" category -- description in `clean.py`; not used for anything
'oflag': cleaning flag indicating image was over cleaned
'bflag': cleaning flag indicating... i forgot 
'uflag': cleaning flag indicating under cleaning -- None of the cleaning flags are very good
'SE_diff': difference between SDSS galaxy center and SExtractor's center
'A_diff': difference between SDSS galaxy center and asymmetry center
'M_diff': difference between SDSS galaxy center and M20 center



