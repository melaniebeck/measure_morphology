#! usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pyfits as fits
from photutils import EllipticalAperture, EllipticalAnnulus
import galaxy
import pdb

    
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
    an = EllipticalAperture((gal.x, gal.y), gal.rp, gal.rp/gal.e, gal.theta)
    plt.figure()
    an.plot(color='blue', lw=1.5, alpha=0.5)
    imgplot = plt.imshow(image, cmap='gray_r', origin='lower')
    imgplot.set_clim(-0.009, 0.022)
    plt.title('1 Petrosian Radius')        
    plt.savefig('output/figures/'+gal.name+'_RpAper.png')
    plt.close()
    
    
#def petro_SB(sb, avgsb, radii, interp_vals, interp_radii):
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

    ax2.semilogx(radii, sb/avgsb, 'ro', label='SB/<SB>')
    #ax2.errorbar(radii, sb/avgsb, yerr=gal._ratio_err, fmt=None)
    ax2.semilogx(interp_r, interp_v, 'k', label='Interpolation')
    ax2.hlines(0.2, 1., np.max(radii), linestyle='--')
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_ylabel(r'$\mu$(R)/<$\mu$(<R)>', fontsize=16)
    #ax2.legend()

    gs.tight_layout(fig)

    plt.savefig('output/figures/'+gal.name+'_SBprofile.png')
    plt.close()
    
    plt.show()



def asym_plot(gal, image):

    #an = EllipticalAperture((gal.x, gal.y), gal.rp, gal.rp/gal.e, gal.theta)

    # read in residual figure created during asymmetry calculation
    residual = fits.getdata('output/asymimgs/'+gal.name+'_res.fits')

    shape = [image.shape[0]/2., image.shape[1]/2.]
    size = 2*gal.rp

    fig = plt.figure(figsize=(8,8))
    gs = plt.GridSpec(3,3)
    ax = fig.add_subplot(gs[:,:])
    plt.setp(ax, xlim=(shape[0]-size, shape[0]+size),
             ylim=(shape[1]-size, shape[1]+size))

    imgplot = ax.imshow(residual, cmap='gray_r', origin='center')
    imgplot.set_clim(gal.med-gal.rms, gal.med+10*gal.rms)
    ax.plot(gal.Ac[0], gal.Ac[1], 'b+', mew=2, ms=20)
    ax.plot(gal.Ac[1], gal.Ac[0], 'r+', mew=2, ms=10)
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.text(left, 0.1*bottom, 'A = %f3'%gal.A, fontsize=16, color='k', 
            transform=ax.transAxes)
    gs.tight_layout(fig)
    plt.savefig('output/figures/'+gal.name+'_asym.png')

    plt.show()

    pdb.set_trace()


    #plt.plot(galcenter[1], galcenter[0], 'k+', mew=2, ms=10)
    # plt.plot(imgcenter[1], imgcenter[0], 'r+', mew=2, ms=10)
    #plt.xlim(low, high)
    #plt.ylim(low, high)
    plt.title('Shifted Cleaned Image')
    plt.savefig('output/asyms/'+self.name+'_asym1_CS.png')
    plt.close()
    
    plt.figure()
    plt.imshow(residual)
    plt.plot(galcenter[1], galcenter[0], 'k+', mew=2, ms=10, 
             label="Asymmetry Center")
    plt.plot(imgcenter[1], imgcenter[0], 'r+', mew=2, ms=10, 
             label="Image Center")
    aperture.plot()
    #plt.xlim(low, high)
    #plt.ylim(low, high)
    plt.title('Asymmetry Residuals (I - I180)')
    plt.legend()
    plt.close()
    return 0.
    

def conc():
    stuff

    #pdb.set_trace()
    
    plt.figure()
    plt.plot(x, ratio, 'ro')
    plt.axhline(y=0.2, linestyle='--', color='k')
    plt.axhline(y=0.8, linestyle='--', color='k')
    plt.ylim(0., 1.4)
    plt.title('Cumulative Flux / Total Flux vs. Radius')
    plt.xlabel('Radius [pixel]')
    plt.ylabel('Flux ratio [counts]')
    #plt.show()
    #pdb.set_trace()
    plt.plot(r, ratios)
    #plt.savefig('output/conc/'+self.name+'_c_conc.png')
    plt.close()

def m20_plot():

    m20_ap2 = EllipticalAperture((xc, yc), self.rp, self.rp/self.e, 
                                 self.theta)
    contours = measure.find_contours(clean_dat, m20_galpix.min())
    for n, contour in enumerate(contours):
        ax.plot(contour[:,1], contour[:,0], linewidth=2)
    pdb.set_trace()
    m20_ap2.plot()
    plt.plot(center[0], center[1], 'r+', mew=2, ms=10)
    ax.set_xlim(shape[0]/2.-3*self.rp, shape[0]/2.+3*self.rp)
    ax.set_ylim(shape[1]/2.-3*self.rp, shape[1]/2.+3*self.rp)
    plt.savefig('output/m20/'+self.name+'_m20.pdf')
    plt.show()

    fig, ax = plt.subplots()
    
    imgplot = ax.imshow(image, cmap='gray_r')
    imgplot.set_clim(-0.01, 0.03)
    
    plt.close()
    

def plot(galaxy, hdulist):
    clean_img = hdulist['CLN'].data

    petro_radius(galaxy, clean_img)
    petro_SB(galaxy)
    asym_plot(galaxy, clean_img)
