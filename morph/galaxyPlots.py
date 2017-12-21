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
from astropy.visualization import ZScaleInterval
from skimage import measure

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
    plt.savefig(gal.outdir+'figures/'+gal.name+'_RpAper.png')

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

    try:
        sb = gal._sb
        avgsb = gal._avgsb
        radii = gal._rads
        interp_v = gal._interpvals
        interp_r = gal._interprads
        ratio = gal._ratio

        try:
            newratio = gal._newratio
        except:
            pass
    except:
        sb = gal['sb']
        avgsb = gal['avgsb']
        radii = gal['rads']
        interp_v = gal['interpVals']
        interp_r = gal['interpRads'] 
        ratio = gal['ratio'] 

        try:
            newratio = gal['newratio']
        except:
            pass

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
        ax2.semilogx(radii, newratio, 'bo')
    except:
        pass
        
    ax2.semilogx(radii, ratio, 'ro', label='SB/<SB>')
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
    plt.savefig(gal.outdir+'figures/'+gal.name+'_SBprofile.png')

    #plt.show()
    plt.close()
    


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
    

def asym_plot(gal, image, ax=None):
	shape = [image.shape[0]/2., image.shape[1]/2.]
	size = 2*gal.Rp
	
	aper = EllipticalAperture((gal.Ax, gal.Ay), gal.Rp, gal.Rp/gal.e, 
                              gal.theta)

	# read in residual figure created during asymmetry calculation
	residual = fits.getdata('../'+gal.outdir+'asymimgs/'+gal.name+'_res.fits')
   
	#hist, bins = np.histogram(residual[shape[0]-size:shape[0]+size, 
    #                                   shape[1]-size:shape[1]+size])
    
	if not ax:
		fig = plt.figure(figsize=(8,8))
		gs = plt.GridSpec(3,3)
		ax = fig.add_subplot(gs[:,:])
		plt.setp(ax, xlim=(shape[0]-size, shape[0]+size),
                 ylim=(shape[1]-size, shape[1]+size))
		plt.rc('font', family='serif')

	else:
		interval = ZScaleInterval()
		imgplot = ax.imshow(interval(residual), cmap='gray_r', origin='center')

		ax.plot(gal.Ax, gal.Ay, 'r+', mew=.5, ms=10)
		aper.plot(lw=1)
		ax.xaxis.set_major_formatter(plt.NullFormatter())
		ax.yaxis.set_major_formatter(plt.NullFormatter())
		ax.text(.05, .05, 'A = %1.3f'%gal.A, fontsize=18, color='yellow', 
            	transform=ax.transAxes)

	if not ax:	
		gs.tight_layout(fig)
		plt.savefig(gal.outdir+'figures/'+gal.name+'_asym.png')
		#plt.show()   
		plt.close()

def conc_plot(gal, image, ax=None):

	if not ax:
		fig = plt.figure(figsize=(8,8))
		gs = plt.GridSpec(3,3)
		ax = fig.add_subplot(gs[:,:])

	else:
		interval = ZScaleInterval()
		ax.imshow(interval(image), cmap='Greys_r', interpolation='nearest', origin='lower')


		apr20 = utils.MyEllipticalAperture((gal.Ax, gal.Ay), 
                                   		   	gal.r20, gal.r20/gal.e, 
                                   		   	gal.theta, image)

		apr80 = utils.MyEllipticalAperture((gal.Ax, gal.Ay), 
                                   		   	gal.r80, gal.r80/gal.e, 
                                   		   	gal.theta, image)

		aptot = utils.MyEllipticalAperture((gal.Ax, gal.Ay), 
                                   		   	1.5*gal.Rp, 1.5*gal.Rp/gal.e, 
                                   		   	gal.theta, image)

		apr20.plot(linewidth=1, color='black')
		apr80.plot(linewidth=1, color='k')
		aptot.plot(linestyle='--', color='k')

		ax.text(.05, .05, 'C = %1.3f'%gal.C, fontsize=18, 
				color='yellow', 
            	transform=ax.transAxes)
		ax.xaxis.set_major_formatter(plt.NullFormatter())
		ax.yaxis.set_major_formatter(plt.NullFormatter())

	if not ax:
		gs.tight_layout(fig)
		plt.savefig(gal.outdir+'figures/'+gal.name+'_conc.png')
		#plt.show()
		plt.close()

def m20_plot(gal, image, ax=None):
	aperture = utils.MyEllipticalAperture((gal.xc, gal.yc), 
                                   		   gal.Rp, gal.Rp/gal.e, 
                                   		   gal.theta, image)
	mask1 = aperture.aper*image
	contours1 = measure.find_contours(mask1, gal.Mlevel1)

    #seg = fits.getdata('output/gini/'+gal.name+'_mask.fits')
    #contours3 = measure.find_contours(seg, gal.Rp_SB)
    
	if not ax:
		fig = plt.figure(figsize=(10,6))
		gs = plt.GridSpec(3,3)
		ax = fig.add_subplot(gs[:,:])
	else:
		interval = ZScaleInterval()
		ax.imshow(interval(image), cmap='Greys_r', origin='lower', interpolation='nearest')

		for n, contour in enumerate(contours1):
			ax.plot(contour[:,1], contour[:,0], color='blue', linewidth=2)
    	
		aperture.plot(linewidth=1, color='black')
		ax.plot(gal.Mx, gal.My, 'r+', mew=.5, ms=10)
		ax.text(.05, .05, r"M$_{20} = $"+"{0:.2f}".format(gal.M20), fontsize=18, 
             	color='yellow', transform=ax.transAxes)
		ax.text(.05, .15, r"$G = ${0:.2f}".format(gal.G), fontsize=18, 
             	color='yellow', transform=ax.transAxes)
		ax.xaxis.set_major_formatter(plt.NullFormatter())
		ax.yaxis.set_major_formatter(plt.NullFormatter())

	if not ax:
		plt.savefig(gal.outdir+'figures/'+gal.name+'_M20.png')
		#plt.show()
		plt.close()
	else:
		return ax


#def main(galMorph, hdulist):

def plot(galMorph, hdulist):
    try:
        clean_img = hdulist['UCLN'].data
    except:
        clean_img = hdulist['CLN'].data

    utils.checkdir(galMorph.outdir+'figures/')

    petro_radius(galMorph, clean_img)
    petro_SB(galMorph)
    asym_plot(galMorph, clean_img)
    conc_plot(galMorph, clean_img)
    m20_plot(galMorph, clean_img)


def main():

	df = pd.read_csv("gz2sample.csv")
	


#if __name__ == "__main__":
#    pyplot()

    #gs.tight_layout(fig)
