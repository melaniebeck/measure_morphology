#! usr/bin/env python

import pdb
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
from sklearn.manifold import LocallyLinearEmbedding as LLE
from sklearn.manifold import SpectralEmbedding, Isomap, TSNE
from astropy.table import Table
from time import time

Axes3D

#dat = Table.read('bigsample_mycat_v19.dat', format='ascii.fixed_width')
dat = Table.read('catalogs/ZEST_catalog_colors.fits')
#dat = dat[0:10000]
#second = dat[10001:15000]
# We need to clean up the data -- can't have any galaxies with a nan or inf 
# in any of its parameters

#keys1 = ['elipt', 'C', 'A', 'G1', 'M1']
keys2 = ['elipt', 'C', 'A', 'G2', 'M2']

keys1 = ['cc', 'aa', 'm20', 'gg']

for keys in keys1, keys2:
    bad, good = [], []
    subdat = dat[keys]
    out = np.dstack(np.array([dat[k] for k in subdat.columns]))[0]
    for i, params in enumerate(out):
        if np.any(np.isnan(params)) or np.any(np.isinf(params)):
            bad.append(i)
        else:
            good.append(i)

    means = [np.mean(subdat[good][k]) for k in keys] 
    print "means done"
    stds = [np.std(subdat[good][k]) for k in keys]
    print "stds done"

    xx = np.array([(params-means)/stds for params in out])

    # this is the data that will go into the LLE! Yay!
    X = xx[good][0:10000]
    Xcolor = dat['color'][good][0:10000]
    X2 = xx[good][10001:15000]
    X2color = dat['color'][good][10001:15000]
    print "dataset ready"


    n_neighbors = [20]
    n_components = 3
    methods = ['hessian']#'hessian', 'ltsa', 'standard', 'modified'
    labels = ['LLE']#'Hessian LLE', 'LTSA', 'LLE', 'MLLE'
    method = 'modified'
    label = 'modified'
    
    fig = plt.figure(figsize=(20,12))
    plt.suptitle(label+" with %i points" % (len(X)), fontsize=14)
    
    #for i, method in enumerate(methods):
    for i, n_neigh in enumerate(n_neighbors):
        print "performing "+label+" with ",n_neigh
        t0 = time()
        A = LLE(n_neigh, n_components, eigen_solver='auto', method=method)
        error = A.fit(X).reconstruction_error_

        #A = Isomap(n_neigh, n_components, eigen_solver='auto')
        #error = A.fit(X).reconstruction_error()
        
        Y = A.fit_transform(X)
        Y2 = A.transform(X2)

        t1 = time()
        print "%s: %.2g sec" %(label, t1-t0)
        print "reconstruction error: ", error

        print "begin plotting"
        ax = fig.add_subplot(211, projection='3d')##, 
        ax.scatter(Y[:, 0], Y[:, 1],  Y[:, 2], 
                   color=Xcolor, 
                   marker='.', cmap=plt.cm.Spectral)#, 
        
        plt.title("%i neighbors (%.2g)" % (n_neigh, t1 - t0))
        ax.set_xlim(min(Y[:,0]), max(Y[:,0]))
        ax.set_ylim(-0.1,0.0 )
        ax.set_zlim(min(Y[:,2]), max(Y[:,2]))
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        #ax.xaxis.set_major_formatter(NullFormatter())
        #ax.yaxis.set_major_formatter(NullFormatter())
        #ax.zaxis.set_major_formatter(NullFormatter())

        ax = fig.add_subplot(212, projection='3d')##, 
        ax.scatter(Y2[:, 0], Y2[:, 1],  Y2[:, 2], 
                   color=X2color, 
                   marker='.', cmap=plt.cm.Spectral)#,
        plt.title('new points')
        ax.set_xlim(min(Y[:,0]), max(Y[:,0]))
        ax.set_ylim(-0.1,0.0 )
        ax.set_zlim(min(Y[:,2]), max(Y[:,2]))
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.axis('tight')
        plt.tight_layout()
    plt.savefig(label+'.pdf')
    plt.show()
    pdb.set_trace()
