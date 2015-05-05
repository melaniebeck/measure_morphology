#! usr/bin/env python

import pdb
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
from sklearn.manifold import LocallyLinearEmbedding as LLE
from sklearn.manifold import SpectralEmbedding, Isomap, TSNE
from sklearn.lda import LDA
from astropy.table import Table
from time import time
import argparse

Axes3D


def whiten_data(data, keys):
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
    colors = data['color'][good]
    return whitedata, colors

def set_axes(axis):
    # determine along which axis to break into slices
    if axis == 0:
        xlab, ylab, lab = 'Y Label', 'Z Label', 'X'
        X, Y = 1, 2
    elif axis == 1:
        xlab, ylab, lab = 'X Label', 'Z Label', 'Y'
        X, Y = 0, 2
    elif axis == 2:
        xlab, ylab, lab = 'X Label', 'Y Label', 'Z'
        X, Y = 0, 1
    return xlab, ylab, lab, X, Y

def inner_ninety(data, targets):
    # Only plot central 90% of the data
    p95 = np.percentile(data, 95.)
    p5 = np.percentile(data, 5.)

    cut = np.where((data > p5) & (data<p95))
    #pdb.set_trace()

    dats = data[cut[0]]
    targets = targets[cut[0]]
    return dats, targets

def plot_dimreduce_3D(X, Xcol, y, ycol, method, n, t):

    fig = plt.figure(figsize=(20,12))
    plt.suptitle(method+" with training sample of %i points" 
                 %len(X), fontsize=14)
    
    ax = fig.add_subplot(121, projection='3d') 
    ax.scatter(X[:, 0], X[:, 1],  X[:, 2], color=Xcol, 
               marker='.', cmap=plt.cm.Spectral)
    
    if n:
        plt.title("%i neighbors (%.2g sec)" %(n, t))
    ax.set_xlim(min(X[:,0]), max(X[:,0]))
    ax.set_ylim(min(X[:,1]), max(X[:,1]))
    ax.set_zlim(min(X[:,2]), max(X[:,2]))
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    #ax.xaxis.set_major_formatter(NullFormatter())

    
    ax = fig.add_subplot(122, projection='3d') 
    ax.scatter(y[:, 0], y[:, 1], y[:, 2], color=ycol, 
               marker='.', cmap=plt.cm.Spectral)
    plt.title('Test sample of %i points' %len(y), fontsize=14)
    ax.set_xlim(min(y[:,0]), max(y[:,0]))
    ax.set_ylim(min(y[:,1]), max(y[:,1]))
    ax.set_zlim(min(y[:,2]), max(y[:,2]))
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.axis('tight')
    plt.tight_layout()
    plt.savefig(method+'.pdf')
    plt.show()

def plot_dimreduce(data, color, method, axis=0):
    xlab, ylab, lab, X, Y = set_axes(axis)

    fig = plt.figure(figsize=(20,12))
    steps = np.linspace(min(data[:,axis]), max(data[:,axis]), 9)

    for inc, step in enumerate(steps):
        if inc == 8:
            break
        else:
            slice = np.where((data[:,axis] > step) & 
                             (data[:,axis] < steps[inc+1]))
        #pdb.set_trace()
        dat = data[slice]
        col = color[slice]

        ax = fig.add_subplot(241+inc)
        ax.scatter(dat[:, X], dat[:, Y], c=col, marker='o')
        ax.set_title(str(step)+' < '+lab+' < '+str(steps[inc+1]))
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        plt.axis('tight')

    plt.tight_layout()
    plt.savefig(method+'_'+lab+'slices.pdf')    
    plt.show()

def plot_LDA(data, targets, classes, axis=0):
    
    xlab, ylab, lab, X, Y = set_axes(axis)

    fig = plt.figure(figsize=(20,12))
    #plt.suptitle('LDA method with 3 output dimensions')
    
    # Only plot central 90% of the data
    dat90, targ90 = inner_ninety(data, targets)

    # Plot planes of data through a individual slices through the third
    steps = np.linspace(min(dat90[:,axis]), max(dat90[:,axis]), 9)

    for inc, step in enumerate(steps):
        if inc == 8:
            break
        else:
            slice = np.where((dat90[:,axis] > step) & 
                             (dat90[:,axis] < steps[inc+1]))
        #pdb.set_trace()
        dat = dat90[slice]
        targs = targ90[slice]

        ax = fig.add_subplot(241+inc)
        for col, cls in zip("rbmgcyk", classes):
            ax.scatter(dat[targs==cls, X], dat[targs==cls, Y], 
                       c=col, marker='o')
        ax.set_title(str(step)+' < '+lab+' < '+str(steps[inc+1]))
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        plt.axis('tight')

    plt.tight_layout()
    plt.savefig('LDA_'+lab+'slices.pdf')    
    plt.show()

def plot_LDA_3D(data, targets, classes):

    fig = plt.figure(figsize=(20,12))
    #plt.suptitle('LDA method with 3 output dimensions')
    
    # Only plot central 90% of the data
    dat90, targ90 = inner_ninety(data, targets)

    pdb.set_trace()
    ax = fig.add_subplot(111, projection='3d')
    for col, cls in zip("rbmgcyk", classes):
        ax.scatter(dat90[targ90==cls, 0], dat90[targ90==cls, 1], 
                   dat90[targ90==cls, 2], c=col, marker='o')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.axis('tight')

    plt.tight_layout()
    plt.savefig('LDA.pdf')    
    plt.show()

def main():
    
    parser = argparse.ArgumentParser(description=
                                'Perform Dimensionality Reduction')
    parser.add_argument('--alg', type=str, default='MLLE',
        help='Algorithm to reduce dimensionality.')
    #parser.add_argument('catalog_name', type=str,
    #    help='Specify the desired name for output catalog.')
    #parser.add_argument('outdir_name', type=str,
    #    help='Specify the desired name for output directory.')
    args = parser.parse_args()

    dat = Table.read('catalogs/ZEST_catalog_colors.fits')

    mykeys = ['elipt', 'C', 'A', 'G', 'M20']
    zkeys = ['cc', 'aa', 'm20', 'gg']

    training_sample = dat[0:10000]
    testing_sample = dat[10001:20000]
    
    train, traincols = whiten_data(training_sample, zkeys)
    test, testcols = whiten_data(testing_sample, zkeys)

    n_neighbors = 20
    n_components = 3

    if args.alg in ['MLLE', 'LLE', 'LTSA', 'HLLE']:
        if args.alg == 'MLLE':
            method = 'modified'
        elif args.alg == 'LLE':
            method = 'standard'
        elif args.alg == 'LTSA':
            method = 'ltsa'
        elif args.alg == 'HLLE':
            method = 'hessian'
        #for i, n_neigh in enumerate(n_neighbors):

        print "performing "+method+" LLE with ",n_neighbors
        t0 = time()
        A = LLE(n_neighbors, n_components, eigen_solver='auto', method=method)
        error = A.fit(train).reconstruction_error_

        Y = A.fit_transform(train)
        Y2 = A.transform(test)

        t1 = time()
        print "%s: %.2g sec" %(args.alg, t1-t0)
        print "reconstruction error: ", error

        print "begin plotting"
        plot_dimreduce(Y2, testcols, method, axis=0)
        plot_dimreduce(Y2, testcols, method, axis=1)
        plot_dimreduce(Y2, testcols, method, axis=2)
        #plot_dimreduce_3D(Y, traincols, Y2, testcols, method, 
        #                  n_neighbors, (t1-t0))

    elif args.alg == 'ISO':
        method='IsoMap'
        print "performing IsoMap with ",n_neighbors

        t0 = time()
        A = Isomap(n_neighbors, n_components, eigen_solver='auto')
        error = A.fit(train).reconstruction_error()
      
        Y = A.fit_transform(train)
        Y2 = A.transform(test)

        t1 = time()
        print "%s: %.2g sec" %(args.alg, t1-t0)
        print "reconstruction error: ", error

        print "begin plotting"
        plot_dimreduce(Y2, testcols, method, axis=0)
        plot_dimreduce(Y2, testcols, method, axis=1)
        plot_dimreduce(Y2, testcols, method, axis=2)
        plot_dimreduce_3D(Y, traincols, Y2, testcols, method, 
                          n_neighbors, (t1-t0))

    elif args.alg == 'LDA':

        print "performing LDA"

        t0 = time()
        A = LDA(n_components, priors=None)

        y = []
        for c in traincols:
            if c == 'red':
                y.append(1)
            elif c == 'paleturquoise':
                y.append(2)
            elif c == 'dodgerblue':
                y.append(3)
            elif c == 'darkcyan':
                y.append(4)
            elif c == 'blue':
                y.append(5)
            elif c == 'yellow':
                y.append(6)
            elif c == 'black':
                y.append(7)
        y=np.array(y)

        Y = A.fit_transform(train, y)
        Y2 = A.fit(train, y).transform(train)

        t1 = time()
        print "%s: %.2g sec" %(args.alg, t1-t0)
        
        predict = A.predict(train)
        #print "Predicted classes:", predict
        #pdb.set_trace()
        
        #plot_LDA_3D(Y2, y,[1,2,3,4,5,6,7])
        plot_LDA(Y2, y, [1,2,3,4,5,6,7], axis=0)
        plot_LDA(Y2, y, [1,2,3,4,5,6,7], axis=1)
        plot_LDA(Y2, y, [1,2,3,4,5,6,7], axis=2)


if __name__ == '__main__':
    main()
