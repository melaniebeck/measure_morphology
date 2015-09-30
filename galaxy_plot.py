'''
Collection of plotting functions used by the galaxy.py class
Has a couple methods that can be accessed outside of the the class
'''

#! usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import numpy as np
import pyfits as fits
from photutils import EllipticalAperture, EllipticalAnnulus
from skimage import measure

import galaxy
import pdb
import utils

    
'''
things we might want to plot:

1. asymmetry residual image with asym center marked
2. m20 images with brightest 20% contours and mc center marked 
(also include asym center?)
3. concentration circles of r80 and r20
4. surface brightnes profiles used to calcualte petrosian radius
5. 1 petrosian radius overplotted onto galaxy image


#'''


def petro_radius(gal, image):
    '''
    Show the galaxy with a 1 Petrosian Radius elliptical aperture
    '''
    shape = [image.shape[0]/2., image.shape[1]/2.]
    size = 2*gal.Rp
    aper = EllipticalAperture((gal.x, gal.y), gal.Rp, gal.Rp/gal.e, gal.theta)

    fig = plt.figure(figsize=(8,8))
    gs = plt.GridSpec(3,3)
    ax = fig.add_subplot(gs[:,:])
    plt.setp(ax, xlim=(shape[0]-size, shape[0]+size),
             ylim=(shape[1]-size, shape[1]+size))

    imgplot = plt.imshow(image, cmap='gray_r', origin='lower')
    imgplot.set_clim(gal.med - 2*gal.rms, np.max(image.flatten())/2.)
    aper.plot(color='blue', lw=1.5, alpha=0.5)

    plt.title('1 Petrosian Radius') 
    gs.tight_layout(fig)
    plt.savefig(gal._outdir+'figures/'+gal.name+'_RpAper.png')
    plt.close()

def petro_radius2(params, image):
    '''
    Show the galaxy with a 1 Petrosian Radius elliptical aperture
    '''
    shape = [image.shape[0]/2., image.shape[1]/2.]
    size = 2*params['Rp']
    aper = EllipticalAperture((params['x'], params['y']), params['Rp'], 
                              params['Rp']/params['e'], params['theta'])
    
    #iqr = np.percentile(imgflat, 75) - np.percentile(imgflat, 25)
    #binsize = 2*iqr/len(imgflat)**(1./3.)
    #print iqr, binsize

    fig = plt.figure(figsize=(8,8))
    gs = plt.GridSpec(3,3)
    ax = fig.add_subplot(gs[:,:])
    plt.setp(ax, xlim=(shape[0], shape[0]+size),
             ylim=(shape[1]-size, shape[1]+size))

    imgplot = plt.imshow(image, cmap='gray_r', origin='lower')
    imgplot.set_clim(np.min(image), np.max(image)/2.)
    aper.plot(color='blue', lw=1.5, alpha=0.5)

    plt.title('1 Petrosian Radius') 
    gs.tight_layout(fig)
    plt.savefig('output/figures/'+params['name']+'_RpAper.png')
    plt.show()
    plt.close()
    
def petro_SB(gal):
    '''
    Plot the SB as a function of radius
    Plot the ratio SB/<SB> as a function of radius
    Include the interpolation
    '''

    sb = gal._sb
    avgsb = gal._avgsb
    radii = gal._rads
    interp_v = gal._interpvals
    interp_r = gal._interprads
    
    xlims = [np.min(radii), np.max(radii)]

    fig = plt.figure(figsize=(9,7))
    gs = plt.GridSpec(2,1)
    
    ax2 = fig.add_subplot(gs[1])
    ax1 = fig.add_subplot(gs[0], sharex=ax2)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2.set_xlabel('Radius (pixels)', fontsize=16)
    ax2.set_xlim(xlims[0], xlims[1])
    ax1.set_xlim(xlims[0], xlims[1])

    ax1.semilogx(radii, sb, 'ro', label='SB')
    #ax1.errorbar(radii, sb, yerr=gal._sberr, fmt=None)
    ax1.semilogx(radii, avgsb, 'go', label='<SB>')
    ax1.hlines(0., 1., np.max(radii), linestyle='--')
    ax1.set_ylim(-0.005, np.max(sb)+0.2*np.max(sb))
    ax1.set_ylabel(r'$\mu$(R)', fontsize=16)

    try:
        ax2.semilogx(radii, gal._newratio, 'bo')
    except:
        pass
        
    ax2.semilogx(radii, gal._ratio, 'ro', label='SB/<SB>')
    #ax2.errorbar(radii, sb/avgsb, yerr=gal._ratio_err, fmt=None)
    #ax2.semilogx(interp_r, interp_v, 'k', label='Interpolation')
    ax2.hlines(0.2, 1., np.max(radii), linestyle='--')
    ax2.vlines(gal.Rp,-0.5, 0.2, linestyle='-.')
    ax2.text(.03, 0.05, "Rp = %3.2f"%gal.Rp, fontsize=16, color='k', 
             transform=ax2.transAxes)
    ax2.set_ylim(-0.5, 1.05)
    ax2.set_ylabel(r'$\mu$(R)/<$\mu$(<R)>', fontsize=16)
    #ax2.legend()

    gs.tight_layout(fig)
    plt.savefig(gal._outdir+'figures/'+gal.name+'_SBprofile.png')
    #plt.show()
    plt.close()
    
    #plt.show()

def petro_SB2(gal):
    '''
    Plot the SB as a function of radius
    Plot the ratio SB/<SB> as a function of radius
    Include the interpolation
    '''

    sb = gal['_sb']
    avgsb = gal['_avgsb']
    radii = gal['_rads']
    
    xlims = [np.min(radii), np.max(radii)]

    fig = plt.figure(figsize=(9,7))
    gs = plt.GridSpec(2,1)
    
    ax2 = fig.add_subplot(gs[1])
    ax1 = fig.add_subplot(gs[0], sharex=ax2)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2.set_xlabel('Radius (pixels)', fontsize=16)
    ax2.set_xlim(xlims[0], xlims[1])
    ax1.set_xlim(xlims[0], xlims[1])

    ax1.semilogx(radii, sb, 'ro', label='SB')
    #ax1.errorbar(radii, sb, yerr=gal._sberr, fmt=None)
    ax1.semilogx(radii, avgsb, 'go', label='<SB>')
    ax1.hlines(0., 1., np.max(radii), linestyle='--')
    ax1.set_ylim(-0.005, np.max(sb)+0.2*np.max(sb))
    ax1.set_ylabel(r'$\mu$(R)', fontsize=16)

    try:
        ax2.semilogx(radii, gal['_newratio'], 'r^')
    except:
        pass
    ax2.semilogx(radii, gal['_ratio'], 'ro', label='SB/<SB>')
    #ax2.errorbar(radii, sb/avgsb, yerr=gal._ratio_err, fmt=None)
    ax2.semilogx(gal['_interprads'], gal['_interpvals'], 'k')
    ax2.hlines(0.2, 1., np.max(radii), linestyle='--')
    ax2.set_ylim(-0.5, 1.05)
    ax2.set_ylabel(r'$\mu$(R)/<$\mu$(<R)>', fontsize=16)
    #ax2.legend()

    gs.tight_layout(fig)
    plt.savefig('output/figures/'+gal['name']+'_SBprofile.png')
    plt.show()
    plt.close()
    

def asym_plot(gal, image):
    shape = [image.shape[0]/2., image.shape[1]/2.]
    size = 2*gal.Rp

    aper = EllipticalAperture((gal.Ax, gal.Ay), gal.Rp, gal.Rp/gal.e, 
                              gal.theta)

    # read in residual figure created during asymmetry calculation
    residual = fits.getdata(gal._outdir+'asymimgs/'+gal.name+'_res.fits')
   
    hist, bins = np.histogram(residual[shape[0]-size:shape[0]+size, 
                                       shape[1]-size:shape[1]+size])
    '''
    pdb.set_trace()
    plt.bar(bins, hist, width=np.mean(np.diff(bins)))
    plt.show()
    pdb.set_trace()
    '''

    font = {'fontname':'Helvetica'}
    fig = plt.figure(figsize=(8,8))
    gs = plt.GridSpec(3,3)
    ax = fig.add_subplot(gs[:,:])
    plt.setp(ax, xlim=(shape[0]-size, shape[0]+size),
             ylim=(shape[1]-size, shape[1]+size))
    plt.rc('font', family='serif')

    imgplot = ax.imshow(residual, cmap='gray_r', origin='center')
    imgplot.set_clim(bins[0], bins[-1])
    ax.plot(gal.Ax, gal.Ay, 'b+', mew=2, ms=20)
    aper.plot(lw=2)
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.text(.05, .05, 'A = %1.3f'%gal.A, fontsize=30, color='w', 
            transform=ax.transAxes)
    gs.tight_layout(fig)
    plt.savefig(gal._outdir+'figures/'+gal.name+'_asym.png')
    #plt.show()   
    plt.close()

def conc_plot(gal, image):
    shape = [image.shape[0]/2., image.shape[1]/2.]
    size = 1.5*gal.Rp

    fig = plt.figure(figsize=(8,8))
    gs = plt.GridSpec(3,3)
    ax = fig.add_subplot(gs[:,:])
    plt.setp(ax, xlim=(shape[0]-size, shape[0]+size),
             ylim=(shape[1]-size, shape[1]+size))

    imgplot = ax.imshow(image, cmap='gray_r')
    imgplot.set_clim(gal.med-.5*gal.rms, gal.med+10*gal.rms)
    patches = [mpatches.Circle((gal.Ax, gal.Ay), gal.r20, fill=None,
                               color='black', lw=2), 
               mpatches.Circle((gal.Ax, gal.Ay), gal.r80, fill=None,
                               color='black', lw=2)]
    for patch in patches:
        ax.add_patch(patch)
    ax.text(.05, .05, 'C = %1.3f'%gal.C, fontsize=20, color='k', 
            transform=ax.transAxes)

    gs.tight_layout(fig)
    plt.savefig(gal._outdir+'figures/'+gal.name+'_conc.png')
    #plt.show()
    plt.close()

def m20_plot(gal, image):
    '''
    Plot the 2 different methods for calculating M20:
    1. Galaxy pixels defined by 1 petrosian radius mask
    2. Galaxy pixels defined by SB cut
    '''
    aperture = utils.EllipticalAperture((gal.Mx, gal.My), 
                                         gal.Rp, gal.Rp/gal.e, 
                                         gal.theta, image)
    mask1 = aperture.aper*image

    contours1 = measure.find_contours(mask1, gal._Mlevel1)

    #seg = fits.getdata('output/gini/'+gal.name+'_mask.fits')
    #contours3 = measure.find_contours(seg, gal.Rp_SB)
    
    shape = [image.shape[0]/2., image.shape[1]/2.]
    size = 2*gal.Rp
    
    fig = plt.figure(figsize=(10,6))
    gs = plt.GridSpec(3,3)
    
    ax1 = fig.add_subplot(gs[:,:])
    plt.setp(ax1, xlim=(shape[0]-size, shape[0]+size),
             ylim=(shape[1]-size, shape[1]+size))

    imgplot = ax1.imshow(image, cmap='gray_r', origin='center')
    imgplot.set_clim(gal.Rp_SB, gal.Rp_SB+10*gal.Rp_SB)
    for n, contour in enumerate(contours1):
        ax1.plot(contour[:,1], contour[:,0], color='blue', linewidth=2)
    aperture.plot(linewidth=2, color='black')
    ax1.plot(gal.Mx, gal.My, 'r+', mew=2, ms=20)
    ax1.text(.05, .05, 'M = %s'%"{0:.3f}".format(gal.M20), fontsize=20, 
             color='k', transform=ax1.transAxes)
    ax1.text(.05, .12, 'G = %s'%"{0:.3f}".format(gal.G), fontsize=20, 
             color='k', transform=ax1.transAxes)
    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax1.yaxis.set_major_formatter(plt.NullFormatter())
    ax1.set_title('M20: Elliptical Aperture')

    '''
    if not isinstance(mask2, int):
        ax2 = fig.add_subplot(gs[:,2:4])
        plt.setp(ax2, xlim=(shape[0]-size, shape[0]+size),
                 ylim=(shape[1]-size, shape[1]+size))
        imgplot = ax2.imshow(image, cmap='gray_r', origin='center')
        imgplot.set_clim(gal.Rp_SB, gal.Rp_SB+10*gal.Rp_SB)
        for n, contour in enumerate(contours2):
            ax2.plot(contour[:,1], contour[:,0], color='blue', linewidth=2)
        for n, contour in enumerate(contours3):
            ax2.plot(contour[:,1], contour[:,0], color='black', linewidth=2)
        ax2.plot(gal.Mcx2, gal.Mcy2, 'r+', mew=2, ms=20)
        ax2.text(.05, .05, 'M = %s'%"{0:.3f}".format(gal.M2), fontsize=20, 
                 color='k', transform=ax2.transAxes)
        ax2.text(.05, .12, 'G = %s'%"{0:.3f}".format(gal.G2), fontsize=20, 
                 color='k', transform=ax2.transAxes)
        ax2.xaxis.set_major_formatter(plt.NullFormatter())
        ax2.yaxis.set_major_formatter(plt.NullFormatter())
        ax2.set_title('M20: SB @ 1Rp on Gaussian-smoothed Image')
    #'''
    gs.tight_layout(fig)

    #pdb.set_trace()
    plt.savefig(gal._outdir+'figures/'+gal.name+'_M20.png')
    #plt.show()
    plt.close()

def plot(galaxy, hdulist):
    try:
        clean_img = hdulist['UCLN'].data
    except:
        clean_img = hdulist['CLN'].data

    utils.checkdir(galaxy._outdir+'figures/')
    petro_radius(galaxy, clean_img)
    petro_SB(galaxy)
    asym_plot(galaxy, clean_img)
    conc_plot(galaxy, clean_img)
    m20_plot(galaxy, clean_img)
