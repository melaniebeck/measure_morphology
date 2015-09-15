#! usr/bin/env python

import pdb
import string
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
from sklearn.manifold import LocallyLinearEmbedding as LLE
from sklearn.manifold import SpectralEmbedding, Isomap, TSNE
from sklearn.lda import LDA
from astropy.table import Table, Column
from time import time
import argparse
import matplotlib.gridspec as gridspec
from astroML.utils import split_samples
from astroML.utils import completeness_contamination

import prep_catalog

Axes3D


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

    
def reassign_colors(label, color):
    for i, l in enumerate(label):
        if l == 1:
            color[i] = 'yellow'
        elif l == 2:
            color[i] = 'red'
        elif l == 3:
            color[i] = 'green'
        elif l == 4:
            color[i] = 'purple' 
    return (label, color)

def simplify_classlabels(targets): 
    #pdb.set_trace()
    labels = [eval(i) for i in targets[:,0]]  
    labels = np.array(labels, dtype='int')
    simple_targets  = reassign_colors(labels, targets[:,1])
    return np.dstack(simple_targets)[0]

def classification_loss(predictions, truevals):
    predicted = np.atleast_2d(predictions)
    true = np.atleast_2d(truevals)

    loss_function = np.array((predicted == true))

def save_dimreduce(X, Y, y, meta, filename):

    X = Table(X)
    
    X['LLE_X'] = Y[:,0]
    X['LLE_Y'] = Y[:,1]
    X['LLE_Z'] = Y[:,2]
    X['target'] = y[:,0]
    X['t_color'] = y[:,1]

    X.meta = meta
    X.write('dimreduce/data/'+filename+'_%iN.fits'%meta['N'], overwrite=True)

def read_dimreduce(dataname, LLEcols, targcols):
    # read in a previous one that I saved. 
    data = Table.read(dataname)
    data = data.filled()
    LLE = np.dstack(np.array([data[k] for k in LLEcols]))[0]
    #targets = np.dstack((data[targcols[0]], data[targcols[1]]))[0]
    targets = data[targcols[0]]
    return data, LLE, targets

def open_previous_LLE(filename, directory, columns):
   
    X1, Y1 = read_dimreduce('%s%s_train.fits'%(directory,filename), columns)
    X2, Y2 = read_dimreduce('%s%s_test.fits'%(directory,filename), columns)

    return (X1, Y1), (X2, Y2)


def plot_dimreduce_3D(X, Xcol, y, ycol, method, n, t, err, sample, two=False):

    fig = plt.figure(figsize=(15,10))
    #plt.suptitle(method+" with training sample of %i points" 
    #             %len(X), fontsize=20)
    #plt.suptitle("LLE with training sample of %i points"%n, fontsize=24)

    #percentiles = np.percentile(X, (0.5, 95.))
    #cut = np.where((X > percentiles[0]) & (X < percentiles[1]))
    #X90 = X[cut]
    #Xcol90 = Xcol[cut]
    #Xcol = Xcol.tolist()
    #Xcol = [col.strip() for col in Xcol]

    if not two:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X[:, 0], X[:, 1],  X[:, 2], color=Xcol, 
                   marker='^', cmap=plt.cm.Spectral, alpha=0.65)
        #ax.set_title('Subsample of ~6100 galaxies in LLE space', fontsize=20)
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        #ax.set_zticklabels([])
        plt.text(0, 0, 'Ellipticals', fontsize=16, color='red')
        ax.grid(True)

        plt.axis('tight')
        plt.tight_layout()
        plt.savefig('%s.png'%sample)

    else:
        ax = fig.add_subplot(121, projection='3d') 
        ax.scatter(X[:, 0], X[:, 1],  X[:, 2], color=Xcol, 
                   marker='.', cmap=plt.cm.Spectral, alpha=0.5)
        #ax.set_title('Subsample of ~6100 galaxies in LLE space', fontsize=20)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.grid(True)

        ax = fig.add_subplot(122, projection='3d') 
        ax.scatter(y[:, 0], y[:, 1], y[:, 2], color=ycol, 
                   marker='^', cmap=plt.cm.Spectral, s=20)
        #ax.set_title('1330 galaxies re-classified via PANOPTES', fontsize=20)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.grid(True)

        plt.axis('tight')
        plt.tight_layout()
        plt.savefig('%s_2figs.png'%sample)

    plt.show()

    #textstr = '$\mu=%.2f$\n$\sigma=%.4f$'%(stats[0], stats[1])
    #ax.text(0.5, 0.5,'matplotlib', horizontalalignment='center',
    #        verticalalignment='center',
    #        transform=ax.transAxes)
    #ax.text(0.05, 0.95, 'ellipticals', color='red', 
    #        horizontalalignment='center', verticalalignment='center',
    #        transform=ax.transAxes, fontsize=14)
    #ax.text(0.04, 0.95, 'edge on', color='green',
    #        horizontalalignment='center', verticalalignment='center',
    #        transform=ax.transAxes, fontsize=14)
    #ax.text(0.03, 0.95, 'disks', color='purple',
    #        horizontalalignment='center', verticalalignment='center',
    #        transform=ax.transAxes, fontsize=14)
    #props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    #ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    #        verticalalignment='top', bbox=props)
    #ax.xaxis.set_major_formatter(NullFormatter())


def plot_dimreduce(data, color, method, n, sample, axis=0):
    xlab, ylab, lab, X, Y = set_axes(axis)

    percentiles = np.percentile(data, (0.5, 95.))
    cut = np.where((data > percentiles[0]) & (data < percentiles[1]))
    dat90 = data[cut]
    col90 = color[cut]

    #pdb.set_trace()
    fig = plt.figure(figsize=(20,12))
    steps = np.linspace(min(dat90[:,axis]), max(dat90[:,axis]), 9)

    for inc, step in enumerate(steps):
        if inc == 8:
            break
        else:
            slice = np.where((dat90[:,axis] > step) & 
                             (dat90[:,axis] < steps[inc+1]))
        #pdb.set_trace()
        dat = dat90[slice]
        col = col90[slice]

        ax = fig.add_subplot(241+inc)
        #pdb.set_trace()
        ax.scatter(np.array(dat[:, X]), np.array(dat[:, Y]), c=col, marker='o')
        ax.set_title(str(step)+' < '+lab+' < '+str(steps[inc+1]))
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        plt.axis('tight')

    plt.tight_layout()
    plt.savefig(method+'_%iN_%sslices_%s.png'%(n,lab,sample))    
    #plt.show()


def plot_LDA(data, targets, classes, colors, sample, axis=0):
    
    xlab, ylab, lab, X, Y = set_axes(axis)

    fig = plt.figure(figsize=(20,12))
    
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
        for col, cls in zip(colors, classes):
            ax.scatter(dat[targs==cls, X], dat[targs==cls, Y], 
                       c=col, marker='o')
        ax.set_title(str(step)+' < '+lab+' < '+str(steps[inc+1]))
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        plt.axis('tight')

    plt.tight_layout()
    plt.savefig('LDA_%sslices_%s.png'%(lab,sample))    
    plt.show()
    plt.close()

def plot_complete_contam(Nparams, completeness, contamination, method, sample):
    # plot completeness vs Nparams
    ax = fig.add_subplot(222)
    ax.plot(Nparams, completeness, 'o-k', ms=6, label='unweighted')
    
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    
    ax.set_ylabel('completeness')
    ax.set_xlim(0.5, 4.5)
    ax.set_ylim(-0.1, 1.1)
    ax.grid(True)
    
    # plot contamination vs Nparams
    ax = fig.add_subplot(224)
    ax.plot(Nparams, contamination, 'o-k', ms=6, label='unweighted')
    
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%i'))
    
    ax.set_xlabel('N morphological parameters')
    ax.set_ylabel('contamination')
    
    ax.set_xlim(0.5, 4.5)
    ax.set_ylim(-0.1, 1.1)
    ax.grid(True)
    plt.savefig('complete_contam_%s_%s.fits'%(method, sample))
    plt.show()

def plot_LDA_3D(data, targets, classes, colors, sample):

    fig = plt.figure(figsize=(20,12))
    #plt.suptitle('LDA method with 3 output dimensions')
    
    # Only plot central 90% of the data
    dat90, targ90 = inner_ninety(data, targets)

    #pdb.set_trace()
    ax = fig.add_subplot(111, projection='3d')#
    for col, cls in zip(colors, classes):
        ax.scatter(dat90[targ90==cls, 0], dat90[targ90==cls, 1], 
                   dat90[targ90==cls, 2], c=col, marker='o')# 
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.axis('tight')

    plt.tight_layout()
    plt.savefig('LDA_%s.png'%sample)    
    plt.show()
    plt.close()

def cut_out_mixedup_region(dat, Y_train):

    cut = np.where((Y_train[:,0] > -0.012 ) & (Y_train[:,0] < -0.005) &
                   (Y_train[:,1] < 0.009 ) & (Y_train[:,1] > -0.009) &
                   (Y_train[:,2] < 0.01) & (Y_train[:,2] > -0.03) )
    sub = dat[cut]

    newdat = dat
    for idx in cut[0]:
        newdat.remove_row(idx)
        print idx, len(newdat)
    pdb.set_trace()
    for idx, thing in enumerate(dat):
        pdb.set_trace()
        if idx in cut[0]:
            print idx, cut[0][idx]

    pdb.set_trace()
    unique, idx, counts =np.unique(sub['dr7objid'], 
                                   return_index=True, 
                                   return_counts=True)
    subY = Y_train[cut]
    suby = y_train[cut]
    subxc = xc_train[cut]
    subcol = np.array(color_train)[cut]

    #subsample = subsample[idx]
            
    subsample = sub[idx]
    subsample = Table(subsample)
    subsample.write('panoptes_to_test.fits', overwrite=True)
    
    f = open('expert_list_2.txt', 'w+')
    f.write('# dr7objid, imgname \n')
    for thing in subsample:
        f.write('%i,%s\n'%(thing['dr7objid'], thing['image_url']))
    f.close()

            
    # match with those that actually got re-classified
    reclass = Table.read('panoptes_run2_classifications_mergers_final.fits')
    matches = []
    for fuck in reclass:
        matches.append(np.where(fuck['dr7objid']==sub['dr7objid'])[0])
        matched = [m[0] for m in matches]
        
        sub_subY = subY[matched]
        
        panclass = reclass['panclass']
        sub_subcol = []
        for p in panclass:
            if p == 1.0:
                sub_subcol.append('yellow')
            elif p == 2.0:
                sub_subcol.append('red')
            elif p == 3.0:
                sub_subcol.append('purple')
            elif p == 4.0:
                sub_subcol.append('purple')
                
    # plot subsection that I selected (~6100 gals)
    # and those that got reclassified
    plot_dimreduce_3D(subY, subcol, 
                      sub_subY, sub_subcol, 
                      method, len(Y_train), 0.0, 0.0, sample)

def QC_plots(dat):

    fig = plt.figure()
    gs1 = gridspec.GridSpec(1,2)
    
    #fig.suptitle('"Direct" sample used in LLE (%i objs)'%len(dat), 
    #             fontsize=20)
    
    ax1 = fig.add_subplot(gs1[0])
    ax1.plot(dat['REDSHIFT'], dat['PETROMAG_R'], 'k.')
    ax1.set_ylabel('Apparent Magnitude [r-band]', fontsize=16)
    ax1.set_xlabel('z', fontsize=16)
    ax1.set_ylim(max(dat['PETROMAG_R']), min(dat['PETROMAG_R']))
    
    ax2 = fig.add_subplot(gs1[1], sharey=ax1)
    ax2.plot(dat['Rp']*0.396, dat['PETROMAG_R'], 'k.')
    ax2.set_xlabel('Petrosian Radius [arc sec]', fontsize=16)
    #ax2.set_ylabel('PETROMAG_R')
    
    gs1.tight_layout(fig, rect=[0, 0, 1, 0.95])
    plt.savefig('QC_directbig_arc.png')
    #plt.show()
    #'''

def replace_panoptes(data):
    # replace original GZ classifications with those from PANOPTES
    count=0
    panoptes = Table.read('panoptes_run2_classifications_mergers_final.fits')
    for p in panoptes:
        obj = np.where(data['dr7objid']==p['dr7objid'])
        if obj[0]:
            count+=1
        data['zclass'][obj] = p['panclass']
        
    pdb.set_trace()
    data.write('zoo2_directbig_MAcut_panoptesclass.fits', overwrite=True)
    #pdb.set_trace()
        
def main():
    
    parser = argparse.ArgumentParser(description=
                                'Perform Dimensionality Reduction')
    parser.add_argument('--alg', type=str, default='MLLE',
        help='Algorithm to reduce dimensionality.')
    parser.add_argument('catalog', type=str,
        help='Specify the catalog on which to perform DimReduce.')
    args = parser.parse_args()

    #dat = Table.read('catalogs/ZEST_catalog_colors.fits')
    #training_sample = dat[0:10000]
    #testing_sample = dat[10001:20000]
    #zkeys = ['cc', 'aa', 'm20', 'gg']

    base = os.path.basename(args.catalog)
    filename = os.path.splitext(base)[0]

    dat = Table.read(args.catalog)
    mkeys = ['elipt', 'C', 'A_1a', 'G', 'M20']#

    #dat.remove_column('color')
    if 'color' not in dat.colnames:
        if 'kaggle' in sample:
            dat = prep_catalog.color_data2(dat, 'gz2class')
        if 'direct' in sample:
            dat = prep_catalog.color_data(dat, 'zclass')
        dat.write(args.catalog, overwrite=True)

    #dat = prep_catalog.adjust_asym(dat, mkeys[2])
    #train, traincols, targets = prep_catalog.whiten_data(dat, mkeys)

    n_neighbors = [10,12,15,20]
    #n_neighbors = [7]
    n_components = 3

    for i, n_neigh in enumerate(n_neighbors):
        
        if args.alg in ['MLLE', 'LLE', 'LTSA', 'HLLE']:
            if args.alg == 'MLLE':
                method = 'modified'
            elif args.alg == 'LLE':
                method = 'standard'
            elif args.alg == 'LTSA':
                method = 'ltsa'
            elif args.alg == 'HLLE':
                method = 'hessian'
                           
            #replace_panoptes(dat)
            #pdb.set_trace()
            #sample = 'directbig_panoptes'

            X, y = prep_catalog.whiten_data(dat, mkeys)

            (dat1, dat2),(thing1,thing2) = split_samples(dat, dat,[0.75, 0.35], 
                                                       random_state=0)
            
            (X_train, X_test), (y_train, y_test) = split_samples(X, y, 
                                                [0.75, 0.35], random_state=0)

            y_train = simplify_classlabels(y_train)
            y_test = simplify_classlabels(y_test)

            #filename = 'modified_7_directbig_new'

            X_train = X
            y_train = simplify_classlabels(y)

            #'''
            #sample ='direct_zcut'

            #Y_train, Y_test = open_previous_LLE(filename)

            #cut = np.where(X1['REDSHIFT'] <= 0.05)
            #X1_cut = X1[cut]
            #QC_plots(X1_cut)
            #Y_train = np.array(Y_train)[cut]
            #col_train = np.array(col_train)[cut]
            #X = Table(X)
            #cut_out_mixedup_region(X, np.array(Y_train))

            #'''
            print "performing "+method+" LLE with",n_neigh,\
                "nearest neighbors"
            print "on training sample of",len(X_train),"objects"

            t0 = time()
            A = LLE(n_neigh, n_components, eigen_solver='auto', method=method)
            error = A.fit(X_train).reconstruction_error_
            
            Y_train = A.fit_transform(X_train)
            Y_test = A.transform(X_train)
            t1 = time()
            #'''        

            metadata = {'method':method, 'N':n_neigh, 'd':n_components, 
                        'error':error, 'time':t1-t0, 'sample':filename+'_total'}
            save_dimreduce(dat, Y_train, y_train, metadata, filename+'_total')

            #metadata = {'method':method, 'N':n_neigh, 'd':n_components, 
            #            'error':error, 'time':t1-t0, 'sample':filename+'_test'}
            #save_dimreduce(X2, Y_test, y_test, metadata, filename+'_test')

            # plot in 3D
            plot_dimreduce_3D(Y_train, y_train[:,1], Y_test, y_test[:,1], 
                              method, n_neigh, error, t1-t0, filename, two=False)

        #====================================================================#

        elif args.alg == 'ISO':
            method='IsoMap'
                
            print "performing IsoMap with",n_neigh,"nearest neighbors"
            print "on training sample of",len(dat),"objects"
            
            t0 = time()
            A = Isomap(n_neigh, n_components, eigen_solver='dense')
            error = A.fit(train).reconstruction_error()
            
            Y = A.fit_transform(train)
            #Y2 = A.transform(test)
            
            t1 = time()
            print "%s: %.2g sec" %(args.alg, t1-t0)
            print "reconstruction error: ", error
            
            print "begin plotting"
            plot_dimreduce(Y, traincols, method, n_neigh, sample, axis=0)
            plot_dimreduce(Y, traincols, method, n_neigh, sample, axis=1)
            plot_dimreduce(Y, traincols, method, n_neigh, sample, axis=2)
            plot_dimreduce_3D(Y, traincols, Y, traincols, method, 
                              n_neigh, (t1-t0), error, sample)
            
        elif args.alg == 'LDA':
            
            print "performing LDA"
            
            X, Xc, y = prep_catalog.whiten_data(dat, mkeys)

            (X_train, X_test), (y_train, y_test) = split_samples(X, y, 
                                                [0.75, 0.25], random_state=0)

            DRclf = LDA(3, priors=None)
            #DRclf.fit(X_train, y_train)
            DRtrain = DRclf.fit(X_train, y_train).transform(X_train)
            DRtest = DRclf.fit(X_train, y_train).transform(X_test)

            classes = np.unique(y_train)
            colors = np.array(['darkred', 'red', 'lightsalmon', 
                               'darkgreen', 'lightgreen', 'lightseagreen', 
                               'indigo', 'darkviolet', 'plum'])
            plot_LDA_3D(DRtrain, y_train, classes, colors, sample)

            pdb.set_trace()

            #classifiers = []
            #predictions = []
            #Nparams = np.arange(1, X.shape[1]+1)
            #for nc in Nparams:
            clf = LDA()
            clf.fit(DRtrain, y_train)
            y_pred = clf.predict(DRtest)
            
            matchesLDA = (y_pred == y_test)
            print np.sum(matchesLDA)

            pdb.set_trace()

            #------------------------------------------

            from sklearn.neighbors import KNeighborsClassifier
            knc = KNeighborsClassifier(5)
            knc.fit(DRtrain, y_train)
            y_pred = knc.predict(DRtest)

            matchesKNN = (y_pred == y_test)
            print np.sum(matchesKNN)

            pdb.set_trace()
            #------------------------------------------

            from astroML.classification import GMMBayes
            gmmb = GMMBayes(9)
            gmmb.fit(DRtrain, y_train)
            y_pred = gmmb.predict(DRtest)

            matchesGMMB = (y_pred == y_test)
            print np.sum(matchesGMMB)

            pdb.set_trace()
            #------------------------------------------

            # plot the results
            fig = plt.figure(figsize=(5, 2.5))
            fig.subplots_adjust(bottom=0.15, top=0.95, hspace=0.0,
                                left=0.1, right=0.95, wspace=0.2)

            # left plot: data and decision boundary
            ax = fig.add_subplot(121)
            pdb.set_trace()
            im = ax.scatter(X[:, 3], X[:, 4], color=Xc, cmap=plt.cm.Spectral, 
                            s=4, lw=0) #cmap=plt.cm.binary,, zorder=2
            im.set_clim(-0.5, 1)
            
            #im = ax.imshow(Z, origin='lower', aspect='auto',
            #               cmap=plt.cm.binary, zorder=1,
            #               extent=xlim + ylim)
            #im.set_clim(0, 1.5)
            
            #ax.contour(xx, yy, Z, [0.5], colors='k')
            
            #ax.set_xlim(xlim)
            #ax.set_ylim(ylim)
            
            ax.set_xlabel('$G$')
            ax.set_ylabel('$M20$')

            #pred, true = classification_loss(predictions, y_test)
            #completeness, contamination = completeness_contamination(pred, true)

            pdb.set_trace()


            #'''
            #t0 = time()
            #A = LDA(n_components, priors=None)
            #Y = A.fit_transform(train, targets)
            #Y2 = A.fit(train, targets).transform(train)
                
            #t1 = time()
            #print "%s: %.2g sec" %(args.alg, t1-t0)
            
            predict = A.predict(train)
            #print "Predicted classes:", predict
            #pdb.set_trace()
            

            #pdb.set_trace()
            #'''
            
            plot_LDA_3D(Y2, targets, classes, colors, sample)
            plot_LDA(Y2, targets, classes, colors, sample, axis=0)
            plot_LDA(Y2, targets, classes, colors, sample, axis=1)
            plot_LDA(Y2, targets, classes, colors, sample, axis=2)
            
            pdb.set_trace()

if __name__ == '__main__':
    main()


'''        
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
'''

'''
# Classify
#------------------------------------------

clf = LDA()
clf.fit(Y_train, y_train)
y_pred = clf.predict(Y_test)

matchesLDA = (y_pred == y_test)
print np.sum(matchesLDA)

pdb.set_trace()

#------------------------------------------

from sklearn.neighbors import KNeighborsClassifier
knc = KNeighborsClassifier(5)
knc.fit(Y_train, y_train)
y_pred = knc.predict(Y_test)

matchesKNN = (y_pred == y_test)
print np.sum(matchesKNN)

pdb.set_trace()
#------------------------------------------

from astroML.classification import GMMBayes
gmmb = GMMBayes(9)
gmmb.fit(Y_train, y_train)
y_pred = gmmb.predict(Y_test)

matchesGMMB = (y_pred == y_test)
print np.sum(matchesGMMB)

pdb.set_trace()
#------------------------------------------
#'''

'''
print "begin plotting"
plot_dimreduce_3D(Y, traincols, Y, traincols, method, 
n_neigh, (t1-t0), error, sample)
plot_dimreduce(Y, traincols, method, n_neigh, sample, axis=0)
plot_dimreduce(Y, traincols, method, n_neigh, sample, axis=1)
plot_dimreduce(Y, traincols, method, n_neigh, sample, axis=2)
'''
        
'''
for color in xc_train:
if color in ['red', 'lightsalmon', 'darkred']:
color_train.append('red')
elif color in ['darkgreen', 'lightgreen', 'lightseagreen']:
color_train.append('green')
elif color in ['indigo', 'darkviolet', 'plum']:
color_train.append('purple')
elif color == 'yellow':
color_train.append('yellow')

for color in xc_test:
if color in ['red', 'lightsalmon', 'darkred']:
color_test.append('red')
elif color in ['darkgreen', 'lightgreen', 'lightseagrean']:
color_test.append('green')
elif color in ['indigo', 'darkviolet', 'plum']:
color_test.append('purple')
elif color == 'yellow':
color_test.append('yellow') 
'''
