
import pdb
import string
import numpy as np
from astropy.table import Table, vstack, hstack, Column
import matplotlib.pyplot as plt




def select_sample():
    '''
    Put together a "pure" sample of GZ2 galaxies by using the thresholds stipulated
    by Willett et al 2013 for the Ellipticals (T07), Edge On disks (T09) and 
    Regular Disks (T05)

    Pull these particular galaxies from the main sample of galaxies contained in
    zoo2MainSpecz.fits 

    Pull merger classification from GZ1 cross-matched with GZ2.

    Define a new column called "zest_class" that assigns each galaxy a letter based
    on it's p-fraction for it's given GZ2 category. 

    M -- Mergers
    E[1-3] -- Ellipticals of varying shapes
    D[1-3] -- Edge on disks with various bulge shapes
    B[1-3] -- Regular disks with varying degrees of buldge size

    Save this subsample of galaxies with a new catalog name (zoo2_puresample.fits)
    
    '''
    gzclass = Table.read('zoo2MainSpecz.fits')
    gzmerg = Table.read('zoo2Mergers.fits')
    
    gzclass.add_column(Column(data=np.empty((len(gzclass),1)), 
                              name='zclass', dtype='S3'))

    # Mark mergers from gzmerg in gzclass with 'M'
    for thing in gzclass:
        if thing['dr7objid'] in gzmerg['dr7objid']:
            thing['zclass']='M'

    mergers = np.where(gzclass['zclass']=='M')
    mer = gzclass[mergers[0]]
    mer.write('mergers.fits', overwrite=True)
    
    gzclass.remove_rows(mergers[0])

   
    ellipts = np.where(
        (gzclass['t01_smooth_or_features_a01_smooth_count']>=20) &
        (gzclass['t01_smooth_or_features_a01_smooth_debiased']>=0.469) )
    
    ell = gzclass[ellipts]
    
    for thing in ell:
        megamax = np.max([thing['t07_rounded_a17_in_between_debiased'],
                          thing['t07_rounded_a18_cigar_shaped_debiased'], 
                          thing['t07_rounded_a16_completely_round_debiased']])
        
        if thing['t07_rounded_a16_completely_round_debiased'] == megamax:
            thing['zclass']='E1'
            
        elif thing['t07_rounded_a17_in_between_debiased'] == megamax:
            thing['zclass']='E2'
            
        elif thing['t07_rounded_a18_cigar_shaped_debiased'] == megamax:
            thing['zclass']='E3'

    ell.write('ellipticals.fits', overwrite=True)

            ###------------------------------------------------------------###

    edgeons = np.where(
        (gzclass['t01_smooth_or_features_a02_features_or_disk_debiased']>=0.430) &
        (gzclass['t02_edgeon_a04_yes_debiased']>=0.715) &
        (gzclass['t02_edgeon_a04_yes_count']>=20) )
    
    edg = gzclass[edgeons]

    for thing in edg:
        megamax = np.max([thing['t09_bulge_shape_a25_rounded_debiased'],
                          thing['t09_bulge_shape_a26_boxy_debiased'],
                          thing['t09_bulge_shape_a27_no_bulge_debiased']])
        
        if thing['t09_bulge_shape_a25_rounded_debiased'] == megamax:
            thing['zclass']='D1'
            
        elif thing['t09_bulge_shape_a26_boxy_debiased'] == megamax:
            thing['zclass']='D2'
            
        elif thing['t09_bulge_shape_a27_no_bulge_debiased'] == megamax:
            thing['zclass']='D3'
    
    edg.write('edgeons.fits', overwrite=True)

            ###-------------------------------------------------------------###

    bulgies = np.where(
        (gzclass['t01_smooth_or_features_a02_features_or_disk_debiased']>=0.430) &
        (gzclass['t02_edgeon_a05_no_debiased']>=0.715) &
        (gzclass['t02_edgeon_a05_no_count']>=20) )
    
    bul = gzclass[bulgies]
    
    for thing in bul:
        megamax = np.max([thing['t05_bulge_prominence_a10_no_bulge_debiased'],
                          thing['t05_bulge_prominence_a11_just_noticeable_debiased'],
                          thing['t05_bulge_prominence_a12_obvious_debiased']+
                          thing['t05_bulge_prominence_a13_dominant_debiased']])

        if thing['t05_bulge_prominence_a10_no_bulge_debiased'] == megamax:
            thing['zclass']='B1'
            
        elif thing['t05_bulge_prominence_a11_just_noticeable_debiased'] == megamax:
            thing['zclass']='B2'
            
        elif (thing['t05_bulge_prominence_a12_obvious_debiased'] + 
              thing['t05_bulge_prominence_a13_dominant_debiased']) == megamax:
            thing['zclass']='B3'

    bul.write('bulgedisks.fits', overwrite=True)

    pure_sample = vstack([mer,ell,edg,bul])
    pure_sample.write('zoo2_pure_sample.fits', overwrite=True)

def match_morphology():

    gzmorph = Table.read('SDSS_morphology_v1.fits')
    gzclass = Table.read('zoo2_puresample.fits')

    idx = []
    for gal in gzmorph:
        name = int(string.split(gal['name'],'_')[1])
        loc = np.where(name == gzclass['dr7objid'])
        if loc[0]:
            idx.append(loc[0][0])


    #training_sample = hstack([])
    pdb.set_trace()


def adjust_columnname():
    gzmorph = Table.read('SDSS_morphology_v2.fits')
    objids = np.array([int(string.split(n,'_')[1]) for n in gzmorph['name']])
    gzmorph['dr7objid'] = objids
    gzmorph.write('SDSS_morphology_v1.fits', overwrite=True)
    pdb.set_trace()
    
def main():

    # Step 1: determine GZ2 class categories/select sample
    select_sample()

    # Step 2: From what morphologies we've already measured, see
    # which galaxies overlap with the selected sample
    #adjust_columnname()
    
    #match_morphology()

    # Step 3: Use this sample as the training set for the LLE. 
    # RUN. DAT. BITCH.


if __name__ == '__main__':
    main()
