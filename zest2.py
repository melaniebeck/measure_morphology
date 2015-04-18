#! usr/bin/env python

import pdb
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
from sklearn.manifold import LocallyLinearEmbedding as LLE
from astropy.table import Table
from time import time

Axes3D

dat = Table.read('bigsample_mycat_v19.dat', format='ascii.fixed_width')

# We need to clean up the data -- can't have any galaxies with a nan or inf 
# in any of its parameters
    
keys1 = ['elipt', 'C', 'A', 'G1', 'M1']
keys2 = ['elipt', 'C', 'A', 'G2', 'M2']

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
    stds = [np.std(subdat[good][k]) for k in keys]

    xx = np.array([(params-means)/stds for params in out])
    # this is the data that will go into the LLE! Yay!
    X = xx[good]
     
    n_neighbors = [5,7,10,12,15,20]
    n_components = 3
    methods = ['modified']#'hessian', 'ltsa', 'standard', 
    labels = ['Modified LLE']#'Hessian LLE', 'LTSA', 'LLE', 
    method = 'modified'
    label = 'MLLE'
    
    fig = plt.figure(figsize=(20,8))
    plt.suptitle("Modified LLE with %i points" % (len(X)), fontsize=14)
    
    #for i, method in enumerate(methods):
    for i, n_neigh in enumerate(n_neighbors):
        t0 = time()
        A = LLE(n_neigh, n_components, eigen_solver='auto', 
                 method=method)
        error = A.fit(X).reconstruction_error_
        Y = A.fit_transform(X)
        t1 = time()
        print "%s: %.2g sec" %(method, t1-t0)
        print "reconstruction error: ", error

        ax = fig.add_subplot(231 + i, projection='3d')
        ax.scatter(Y[:, 0], Y[:, 1], Y[:, 2], marker='^', cmap=plt.cm.Spectral)
        plt.title("%i neighbors (%.2g)" % (n_neigh, t1 - t0))
        ax.set_xlim(min(Y[:,0]), max(Y[:,0]))
        ax.set_ylim(min(Y[:,1]), max(Y[:,1]))
        ax.set_zlim(min(Y[:,2]), max(Y[:,2]))
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        #ax.xaxis.set_major_formatter(NullFormatter())
        #ax.yaxis.set_major_formatter(NullFormatter())
        #ax.zaxis.set_major_formatter(NullFormatter())
        plt.axis('tight')
        plt.tight_layout()
    plt.savefig(label+'.pdf')
    plt.show()
    pdb.set_trace()
