#! /usr/bin/env python

import os
import subprocess
import pyfits
import ConfigParser
import pdb
import argparse

def single_SE(image, outstr, outdir='', params={}, outstr2=0):
    basename = os.path.basename(os.path.splitext(image)[0])
    if isinstance(outstr2, int):
        cat = '%s%s_%s_cat.fits' %(outdir, basename, outstr)
        seg = '%s%s_%s_seg.fits' %(outdir, basename, outstr)
    else:
        cat = '%s%s_%s_%s_cat.fits' %(outdir, basename, outstr, outstr2)
        seg = '%s%s_%s_%s_seg.fits' %(outdir, basename, outstr, outstr2)
    
    params['-catalog_name'] = cat
    params['-checkimage_type'] = 'segmentation'
    params['-checkimage_name'] = seg

    args = ['/usr/bin/sextractor', image]
    for key, value in params.iteritems():
        args.append(key)
        args.append(value)

    subprocess.check_call(args)

def run_SE(image, section, outdir='', outstr2=0):
    ''' Run SExtractor on COSMOS/ZEST cutouts using the parameters
        in se_param.cfg
         
        If section = 'BRIGHT', the parameters are geared toward finding
        the brightest objects -- MINCONT = 0.04, THRESH = 2.2

        If section = 'FAINT', parametrs are set to find the faintest 
        objects -- MINCONT = 0.065, THRESH = 1.0

        If section = 'SMOOTH', the parameters are identical to FAINT
        except with the addition of gaussian smoothing
    '''

    config = ConfigParser.ConfigParser()
    config.read('se_params.cfg')
    options = config.options(section)
    
    params = {}
    for option in options:
        params[option] = config.get(section, option)

    if section == 'BRIGHT':
        outstr = 'bright'
    if section == 'FAINT':
        outstr = 'faint'
    if section == 'SMOOTH':
        outstr = 'smooth'
    
    if isinstance(outstr2, int):
        single_SE(image, outstr, outdir, params)
    else:
        single_SE(image, outstr, outdir, params, outstr2)
        
def main():
    
    parser = argparse.ArgumentParser(description='Run SExtractor')
    parser.add_argument('Imagename', type=str, 
        help='Name of fits image to run SExtractor on')
    parser.add_argument('--mode', default='BRIGHT',
        help='Specify which config file to use: Bright, Smooth, or Faint')
    args = parser.parse_args()

    run_SE(args.Imagename, args.mode)


if __name__ == '__main__':
    main()
