#! usr/bin/env python

import matplotlib.pyplot as plt
from photutils import EllipticalAperture, EllipticalAnnulus
import galaxy


class GalPlots(galaxy.Galaxy):
    
    '''
    things we might want to plot:
    
    1. asymmetry residual image with asym center marked
    2. m20 images with brightest 20% contours and mc center marked 
       (also include asym center?)
    3. concentration circles of r80 and r20
    4. surface brightnes profiles used to calcualte petrosian radius
    5. 1 petrosian radius overplotted onto galaxy image

    
    #'''

    def petro_plotradius(self, clean_dat):
        '''
        Show the galaxy with a 1 Petrosian Radius elliptical aperture
        '''
        an = EllipticalAnnulus((self.x, self.y), self.rpet, \
                               self.rpet/self.e, self.theta)
        plt.figure()
        an.plot(color='blue', lw=1.5, alpha=0.5)
        imgplot = plt.imshow(clean_dat, cmap='gray_r', origin='lower')
        imgplot.set_clim(-0.009, 0.022)
        plt.title('1 Petrosian Radius')        
        plt.savefig('output/figures/'+self.name+'_RpAper.png')
        plt.close()


    def petro_plotSB(self, sb, avgsb, radii, interp_vals=[], interp_radii=[]):
        '''
        Plot the SB as a function of radius
        Plot the ratio SB/<SB> as a function of radius
        Include the interpolation
        '''

        xlims = [radii[0], radii[:-1]]
        ylims = [-.2*max(avgsb), max(avgsb) + 0.2*max(avgsb)] 
        
        fig, axes = plt.subplots(nrows=2, sharex=True)
        plt.setp(axes[0], title='Surface Brightness')
        plt.setp(axes[1], title='u(R)/<u(R)>')

        fig.tight_layout()
        fig.subplots_adjust(top=0.9)

        axes[0].plot(radii, sb, 'ro', label='SB')
        axes[0].plot(radii, avgsb, 'go', label='<SB>')
        axes[0].hlines(0., 0, np.max(radii))
        axes[0].set_xlabel('Radius (pixels)')
        axes[0].set_xscale('log')
        plt.legend()

        axes[1].plot(radii, sb/avgsb, 'ro', label='SB/<SB>')
        axes[1].plot(interp_radii, interp_vals, label='Interpolation')
        axes[1].hlines(0.2, 0, np.max(radii))
        axes[1].set_xscale('log')
        plt.legend()

        plt.savefig('output/figures/'+self.name+'_SBprofile.png')
        plt.close()

        '''
        plt.figure()
        plt.subplot(211)
        plt.title('Surface Brightness Profile')
        plt.plot(rpix, sb, 'ro', label='Surface Brightness')
        plt.plot(rpix, avgsb, 'g^', label='Average SB')
        plt.plot(rpix, zeros, 'k--')
        plt.xscale('log')
        plt.axis([xlims[0], xlims[1], ylims[0], ylims[1]])
        plt.legend(loc='upper right')
        
        plt.subplot(212)
        plt.title('u(R)/<u(R)>')
        plt.plot(rpix, sb/avgsb, 'bo', label='SB/Avg_SB')
        plt.plot(rpix, eta, 'k--', label='eta=0.2')  
        plt.xscale('log')
        plt.axis([xlims[0], xlims[1], 0., 1.1])
        plt.legend(loc='upper right')    
        plt.savefig('output/profiles/'+self.name+'_prof.png', \
                    bbox_inches='tight')
        plt.close() 
        '''

    def asym_plot():
        stuff

    def conc_plot():
        stuff
    

