import glob
import argparse
import os
import string
import pdb
import pyfits as fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from math import pi
from scipy.interpolate import interp1d
from astropy.table import Table

import utils

#f = open('prad2.txt','r')
#f.readlines()
'''
with open('prad2.txt','r') as f:
    rpet = map(float, f)
rpet = np.asarray(rpet, dtype=float)
'''

cdat = Table.read('bright_gals.fits')
mnew = Table.read('bright_gals_newpetro.fits')
mdat = Table.read('bright_gals_config1.fits')
mdatb = Table.read('bright_gals_config_bright.fits')
mdatf = Table.read('bright_gals_config_faint.fits')
mdatc = Table.read('bright_gals_combined1.fits')
mdatct = Table.read('bright_gals_cstheta_test.fits')
mdatcp = Table.read('bright_gals_csparams_test.fits')


datasets = [mdat, mdatb, mdatf, mdatc, mdatct, mdatcp]
filenames = ['mdat','mdatb','mdatf','mdatc','mdatct','mdatcp']
formats = ['ro', 'bo', 'yo', 'go', 'co', 'mo']

testfile = np.genfromtxt(r'testingshits/test_mel_prof_IDL_gal96.dat', dtype=float)
testfile

#pdb.set_trace()


cut = np.where(cdat['rpet'] > 1000.)
cdat['rpet'][cut] = 0.

#pdb.set_trace()

img = fits.getdata('acs_mosaic_2.0/0001_149.935944_2.666823_acs_I_mosaic_30mas_sci.fits')

#mydat = fits.getdata('output/0001_149.935944_2.666823_acs_I_mosaic_30mas_sci.fits',3)

plt.figure()
plt.subplot(2,2,1)  # 4 plots in a grid
plt.plot(#cdat['acs_elongation'], mdat['ELONGATION'],'ro',
         #cdat['acs_elongation'], mdatb['ELONGATION'],'bo', 
         #cdat['acs_elongation'], mdatct['ELONGATION'],'go',
         cdat['acs_elongation'], mnew['ELONGATION'],'go',
         #cdat['acs_elongation'], mdatf['ELONGATION'],'yo',
         [0,7], [0,7], 'k--')
plt.title('elongation comparison')

plt.subplot(2,2,2)
plt.plot(#ccdat['acs_a_image'], mdat['A_IMAGE'], 'ro',
         #cdat['acs_a_image'], mdatb['A_IMAGE'], 'bo',   
         #cdat['acs_a_image'], mdatct['A_IMAGE'], 'go', 
         cdat['acs_a_image'], mnew['A_IMAGE'], 'go', 
         #cdat['acs_a_image'], mdatf['A_IMAGE'], 'yo',
         [0,60], [0,60], 'k--')
plt.title('semi-major axis comparison')

plt.subplot(2,2,3)
plt.plot(#cdat['acs_theta_image'], mdat['THETA']*180./pi+90, 'ro',
         #cdat['acs_theta_image'], mdatb['THETA']*180./pi+90, 'bo',
         #cdat['acs_theta_image'], mdatct['THETA']*180./pi+90, 'go',
         cdat['acs_theta_image'], mnew['THETA']*180./pi+90, 'go',
         #cdat['acs_theta_image'], mdatf['THETA']*180./pi+90, 'yo', 
         [-10, 80], [0,90], 'g--',
         [0,90], [0,90], 'k--',
         [-100, -10], [90,180], 'g--',
         [-90,0], [90,180], 'k--')
plt.title('theta comparison (degrees)')

plt.subplot(2,2,4)
plt.plot(cdat['rpet'], mnew['PETRO_RAD'], 'ro',
         #cdat['rpet'], mdat['PETRO_RAD'], 'ro',
         #cdat['rpet'], mdatb['PETRO_RAD'], 'bo', 
         #cdat['rpet'], mdatct['PETRO_RAD'], 'go',
         #cdat['rpet'], mdatcp['PETRO_RAD'], 'mo',
         #cdat['rpet'], mdatf['PETRO_RAD'], 'yo',
         [0,300], [0,300], 'k--')
plt.title('petrosian radius comparison')

plt.savefig('newpetro_compare.png')
#plt.show()

plt.clf()


plt.plot(cdat['rpet'], mdat['PETRO_RAD'], 'ro', label='Original')
#plt.plot(cdat['rpet'], mdatb['PETRO_RAD'], 'bo', label='bright SE')
#plt.plot(cdat['rpet'], mdatf['PETRO_RAD'], 'yo', label='faint SE')
plt.plot(cdat['rpet'], mdatc['PETRO_RAD'], 'co', label='Combo SE')
#plt.plot(cdat['rpet'], mdatct['PETRO_RAD'], 'go', label='CS theta')
#plt.plot(cdat['rpet'], mdatcp['PETRO_RAD'], 'mo', label='CS params')
plt.plot(cdat['rpet'], mnew['PETRO_RAD'], 'go', label='New Method')
plt.plot([0,250], [0,250], 'k--')
#'''
plt.title('Petrosian Radius Comparison')
plt.xlabel('Claudias')
plt.ylabel('Mine')
plt.legend(loc='lower right')
plt.savefig('NewPetroComparison.png')
plt.show()

pdb.set_trace()

plt.subplot(2,1,1)
plt.plot(cdat['rpet'], mnew['PETRO_RAD'], 'ro')
plt.plot([0,250], [0,250], 'k--')

plt.title('petrosian radius comparison')
plt.xlabel('Claudias')
plt.ylabel('Mine')

plt.subplot(2,1,2)
plt.title('MyPetro/CPetro vs CPetro')
ratio = mnew['PETRO_RAD']/cdat['rpet']

'''
subset1 = np.where(ratio > 1.1)
subset2 = np.where(ratio < .9)
idx = 0
for xy in zip(cdat['rpet'][subset1],ratio[subset1]):
    plt.annotate(str(subset1[0][idx]+1), xy=xy)
    idx += 1
idx = 0
for xy in zip(cdat['rpet'][subset2],ratio[subset2]):
    plt.annotate(str(subset2[0][idx]+1), xy=xy)
    idx += 1
#'''

# 34% OF GALAXIES ARE OFF BY MORE THAN 10% COMPARED TO CLAUDIA'S ORIGINAL CATALOG!
        
plt.plot(cdat['rpet'], mnew['PETRO_RAD']/cdat['rpet'], 'ro', [0,250], [1.0, 1.0], 'k--')

plt.savefig('petro_compare_wratio.png')
plt.show()
plt.clf()

pdb.set_trace()

for idx, dataset in enumerate(datasets):
    plt.subplot(2,2,1)  # 4 plots in a grid
    plt.plot(cdat['acs_elongation'], dataset['ELONGATION'], formats[idx], 
             [0,7], [0,7], 'k--')
    plt.title('elongation')
    plt.ylabel('Mine')
    plt.subplot(2,2,2)
    plt.plot(cdat['acs_a_image'], dataset['A_IMAGE'], formats[idx],
             [0,60], [0,60], 'k--')
    plt.title('semi-major axis')
    plt.subplot(2,2,3)
    plt.plot(cdat['acs_theta_image'], dataset['THETA']*180./pi+90., formats[idx],
             [0,90], [0,90], 'k--')
    plt.ylabel('Mine')
    plt.xlabel('Claudias')
    plt.title('theta (degrees)')
    plt.subplot(2,2,4)
    plt.plot(cdat['rpet'], dataset['PETRO_RAD'], formats[idx],
             [0,300], [0,300], 'k--')
    plt.xlabel('Claudias')
    plt.title('petrosian radius')
    plt.savefig(filenames[idx]+'.png')
    plt.show()

pdb.set_trace()

fig = plt.figure()
ax = fig.add_subplot(111)
imgplot = plt.imshow(img,cmap='gray_r', origin='lower')
imgplot.set_clim(-0.009, 0.022)

a = dat['acs_a_image'][0]
b = dat['acs_b_image'][0]
t = dat['acs_theta_image'][0]
print t

el = Ellipse(xy=(250.,250.), width=2*a, height=2*b, angle=0.0)

ax.add_artist(el)
el.set_alpha(0.5)

ax.autoscale(False)

plt.show()

pdb.set_trace()



subset = np.where(rpet > 175.)
idx = 0
for xy in zip(dat['rpet'][subset],rpet[subset]):
    ax.annotate(str(subset[0][idx]+1), xy=xy)
    idx += 1

plt.plot(dat['rpet'], rpet,'ro',[0,300],[25,325],'k--')
plt.xlabel('claudias')
plt.ylabel('mine')   

plt.show()
pdb.set_trace()

#gridspec - multiplot
