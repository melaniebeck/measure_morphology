'''
Combines all morphology output catalogs into one comprehensive catalog
'''


from astropy.table import Table
import numpy as np
import pdb
import glob
import os


directory = 'SDSSmorphology_catalogs/'
rootname = 'SDSSmorphology_catalog_chunk'

files = np.array(sorted(glob.glob(directory+rootname+'*.fits')))	

cats = []
for f in files:
	temp = Table.read(f)
	outdirnum = int(os.path.splitext(os.path.basename(f))[0]\
				.split('_')[2].split('chunk')[1])
	temp['outdir']=np.full((len(temp),),fill_value=outdirnum,dtype='int32')
	cats.append(temp)

full_catalog = np.hstack([cat for cat in cats])
full_catalog = Table(full_catalog)
full_catalog.write(directory+'SDSSmorphology_full_catalog_92415.fits',
		   overwrite=True)

exit()
