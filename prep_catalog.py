
import os
import pdb
import argparse
import numpy as np
from astropy.table import Table, Column
import string


def colorcode_zest(cat, outname):
    '''
    This is an old version for ZEST data. I think.
    '''
    color = []
    for t,b in zip(cat['type'], cat['bulg']):
        if t == 1:
            color.append('red')
        elif t == 2:
            if b == 0.0:
                color.append('paleturquoise')
                #mine['color'][i]='paleturquoise'
            elif b == 1.0:
                color.append('dodgerblue')
                #mine['color'][i]='dodgerblue'
            elif b == 2.0:
                color.append('darkcyan')
                #mine['color'][i]='darkcyan'
            else:
                color.append('blue')
        elif t == 3:
            color.append('yellow')
        else:
            color.append('black')

    cat['color']=color
    cat.write(outname)


def color_data(data, classkey):
    '''
    For the Direct Training Samples
    ---------------------------------------------------------
    These already have the 'zclass' class numeric. 
    Based on this number, assign each galaxy a color
    '''
    colors = np.array(['yellow', 'darkred', 'red', 'lightsalmon', 'darkgreen', 
                       'lightgreen', 'lightseagreen', 'indigo', 
                       'darkviolet', 'plum'])
    classes = np.unique(data[classkey])
    data.add_column(Column(data=np.empty(len(data)), 
                           name='color', dtype='S15'))
    for clss, col in zip(classes, colors):
        colloc = np.where(data[classkey]==clss)
        data['color'][colloc[0]]=col

    #dat.write(args.catalog, overwrite=True)
    return data

def color_data2(data, classkey):
    '''
    For the Kaggle Training Samples:
    ------------------------------------------------------------
    For each galaxy, determines a class based on the 'gz2class' moniker
    Galaxy classes are described in my Google Doc.
    Creates two new columns: 
        'color' an actual color used in plotting
        'zclass' a numerical class descriptor used in some dimreduce algorithms
    '''
    gzclasses = np.unique(data[classkey])
    classes = {'M': ('yellow', 1.0), 'Er': ('darkred', 2.0),
               'Ei': ('red',2.1), 'Ec':('lightsalmon',2.2),  
               'Ser': ('darkgreen', 3.0), 'Seb':('lightgreen', 3.1), 
               'Sen': ('lightseagreen', 3.2), 'Sd': ('indigo', 4.2), 
               'Sc': ('darkviolet', 4.1), 'Sb':('plum', 4.0)}
    
    # remove stars. i mean seriously. wtf, guys.
    stars = np.where(data[classkey]==gzclasses[0])
    data.remove_rows(stars)
    data.add_column(Column(data=np.empty(len(data)), 
                           name='color', dtype='S15'))

    data.add_column(Column(data=np.zeros(len(data)), 
                           name='zclass', dtype=float))
    for obj in data:

        if obj[classkey][1] == 'B':
            key = obj[classkey][0:1]+obj[classkey][2:3]
        else:
            key = obj[classkey][0:2]    

        if key[1] == 'a':
            key = obj[classkey][0]+'b'
        if key[1] == 'e':
            key = obj[classkey][0:3]

        obj['color']=classes[key][0]
        obj['zclass']=classes[key][1]

    return data

def adjust_asym(data, key):
    '''
    Any value of Asymmetry that is less than zero is set equal to zero.
    '''
    bad = np.where(data[key] < 0.)
    data[key][bad]=0.
    return data

def whiten_data(data, keys):
    '''
    -- Stack all the morphological measurements in the proper fashion
    -- Remove entries (rows) which have non-numeric morphological values, 
       i.e. nan, inf
    -- Find the mean and standard deviation of each morphological feature
    -- "Whiten" each column by subtraction mean and dividing by std dev.
    '''
    bad, good = [], []
    subdat = data[keys]
    out = np.dstack(np.array([data[k] for k in subdat.columns]))[0]
    for i, params in enumerate(out):
        if np.any(np.isnan(params)) or np.any(np.isinf(params)):
            bad.append(i)
        else:
            good.append(i)

    means = [np.mean(subdat[good][k]) for k in keys] 
    stds = [np.std(subdat[good][k]) for k in keys]

    xx = np.array([(params-means)/stds for params in out])
    whitedata = xx[good]
    colors = np.array(data['color'][good])
    labels = np.array(data['zclass'][good])
    targets = np.dstack((labels, colors))[0]
    return whitedata, targets


def adjust_columnname(data, filename):
    '''
    My SDSS_morphology files have the dr7objid in string format with a 'f_'
    in front. Remove the 'f_' and cast the number as an int in order to match
    with other GZ catalogs
    '''
    objids = []
    i = 1
    for d in data:
        name = d['name'].strip()
        if name == 'f_mahshit':
            objids.append(000000000000000000+i)
            i+=1
        else:
            objids.append(int(string.split(name, '_')[1]))

    print i
    pdb.set_trace()
    #objids = np.array([int(string.split(n,'_')[1]) for n in data['name']])
    data['dr7objid'] = objids
    data.write(filename, overwrite=True)

def main():
    #'''
    #filename = 'testsample_direct_MAcut.fits'
    filename = 'SDSS_morphology_168K.fits'
    data = Table.read(filename)
    adjust_columnname(data, filename)

    #data = color_data(data, 'zclass')
    #data.write(filename, overwrite=True)
    #'''

    pdb.set_trace()

    filename = 'testsample_kaggle_MAcut.fits'
    data = Table.read(filename)
    #stars = np.where(data['gz2class']=='A      ')
    #pdb.set_trace()
    data = color_data2(data, 'gz2class')
    data.write(filename, overwrite=True)
    pdb.set_trace()


if __name__ == '__main__':
    main()
