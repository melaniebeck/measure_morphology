import re, glob, os, string, pdb, sys
import argparse, warnings

import astropy.visualization as vis #import ZScaleInterval
from astropy.table import Table
import astropy.io.fits as fits

import pandas as pd
import numpy as np

import morph

from joblib import Parallel, delayed

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
mpl.rcParams.update({
	'font.size': 24,
	'font.family': 'STIXGeneral',
	'mathtext.fontset':'stix',
	'xtick.bottom': False,
	'ytick.left': False})



def get_clean(hdulist):
	try:
		return hdulist['UCLN'].data
	except:
		return hdulist['CLN'].data


def do_the_things(args, filename, k):
	#Open their cleaned FITS iamges
	hdulist = fits.open(filename)

	flags = np.zeros(4)	
	g = morph.GalaxyMorphology(hdulist, filename, flags, args.outdir)

	dd = g.table(init=True)
	dfmini = pd.DataFrame(g.__dict__, index=[k])

	return dfmini


def main():  
	parser = argparse.ArgumentParser(description='Perform LLE/PCA/whatevs')
	parser.add_argument('-d', dest="directory", type=str, 
        help='Directory of fits images on which to run LLE.')
	parser.add_argument('-c', dest="catalog_name", type=str, default=None,
        help='Specify the desired name for output catalog.')
	parser.add_argument('--outdir', type=str, default='output/datacube/', 
        help='Specify the desired name for output directory.')
	parser.add_argument('--hard', dest='hardest', default=False, 
		help='run the "hardest" galaxies')
	args = parser.parse_args()


	# There are a lot of useless warnings that pop up -- suppress them!
	warnings.filterwarnings('ignore', message='Overwriting existing file .*',
                            module='pyfits')

	outdir = "/data/extragal/beck/gzcodez/SDSSmorphology_catalogs/110817"


	for chunk in range(26,30):
		#fitsfiles = "output_4Rp/chunk1/datacube/f_587725550679031929_4Rp.fits"
		#do_the_things(args, fitsfiles, 0)

		fitsfiles = glob.glob("output_4Rp/chunk{}/datacube/f*4Rp.fits".format(chunk))		

		result = Parallel(n_jobs=30, verbose=51)(delayed(do_the_things)(args,f, k) for k, f in enumerate(fitsfiles))

		#pdb.set_trace()

		df = pd.concat(result)
		df.to_csv("{}/SDSSmorphology_catalog_chunk{}.csv".format(outdir, chunk))

	pdb.set_trace()





if __name__ == "__main__": 
	main()



