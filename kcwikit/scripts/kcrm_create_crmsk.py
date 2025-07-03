#!/usr/bin/env python

"""
Create crmsk files from grouped json files.
"""

import argparse
import json
import numpy as np
import os
from astropy.table import Table
from astropy.io import fits
import re
import warnings
from astropy.stats import sigma_clip,sigma_clipped_stats

def parser_init():
    """Create command-line argument parser for this script."""
    parser = argparse.ArgumentParser(description="Create crmsk files from grouped json files.")
    parser.add_argument(
        'filename',
        type=str,
        help='Input json file that contains the file grouping.',
        nargs=1
        )

    parser.add_argument(
        '-l', '--label',
        help='List of grouping labels to operate.',
        type=str,
        nargs='*',
        default = []
        )
    
    parser.add_argument(
        '-s', '--sigma',
        help='Sigma value for sigma clipping.',
        type=float,
        default = 3
    )
    
    parser.add_argument(
        '--base_dir',
        help='Location of the raw directory',
        nargs=1,
        default='./'
        )
    parser.add_argument(
        '--redux_dir',
        help='Location of the redux directory',
        nargs=1,
        default='redux'
        )
    return parser


#### Need to move to utils
def read_proctab(tfil='kcwi.proc'):
    # from official DRP
    proctab = Table.read(tfil, format='ascii.fixed_width')
    proctab.dtypes = ('int32', 'S24', 'int64', 'S9', 'S12',
                            'float64', 'S4', 'S6', 'S5', 'float64',
                            'float64', 'S4', 'S5', 'float64', 'int32',
                            'S5', 'S25', 'S25', 'S25')

    # format columns
    proctab['GANG'].format = '7.2f'
    proctab['CWAVE'].format = '8.2f'
    proctab['MJD'].format = '15.6f'
    # prevent string column truncation
    for col in proctab.itercols():
        if col.dtype.kind in 'SU':
            proctab.replace_column(col.name, col.astype('object'))

    return proctab

#### Need to move to utils
def search_proctab(proctab, frame, target_type=None, target_group=None,
        nearest=False):
    # copied from official DRP

    if target_type is not None and proctab is not None:
        tab = proctab[(proctab['CAM'] ==
                            frame.header['CAMERA'].upper().strip())]
        # get target type images
        tab = tab[(proctab['TYPE'] == target_type)]
        # BIASES must have the same CCDCFG
        if 'BIAS' in target_type:
            tab = tab[(tab['DID'] == int(frame.header['CCDCFG']))]
            if target_group is not None:
                tab = tab[(tab['GRPID'] == target_group)]
        # raw DARKS must have the same CCDCFG and TTIME
        elif target_type == 'DARK':
            tab = tab[tab['DID'] == int(frame.header['CCDCFG'])]
            tab = tab[tab['TTIME'] == float(frame.header['TTIME'])]
            if target_group is not None:
                tab = tab[tab['GRPID'] == target_group]
        # MDARKS must have the same CCDCFG, will be scaled to match TTIME
        elif target_type == 'MDARK':
            tab = tab[(tab['DID'] == int(frame.header['CCDCFG']))]
        elif target_type == 'OBJECT':
            tab = tab[tab['GRPID'] == target_group]
        else:
            tab = tab[(tab['CID'] == frame.header['STATEID'])]
        # Check if nearest entry is requested
        if nearest and len(tab) > 1:
            tfno = frame.header['MJD']
            minoff = 99999
            trow = None
            for row in tab:
                off = abs(row['MJD'] - tfno)
                if off < minoff:
                    minoff = off
                    trow = row
            if trow is not None:
                tab = tab[(tab['MJD'] == trow['MJD'])]
    else:
        tab = None
    return tab


def create_crmsk(filename, label=[], sigma=3, base_dir='./', redux_dir='redux'):
    
    if isinstance(filename, list):
        filename = filename[0]

    if isinstance(redux_dir, list):
        redux_dir = redux_dir[0]

    if isinstance(base_dir, list):
        base_dir = base_dir[0]
    
    redux = os.path.join(base_dir, redux_dir)

    obj_dict = json.load(open(filename, 'r'))

    if not os.path.isdir(os.path.join(redux, 'int_median')):
        os.mkdir(os.path.join(redux, 'int_median'))

    for obj in list(obj_dict.keys()):

        if len(label) > 0:
            if obj not in label:
                continue

        if len(obj_dict[obj]) >= 3:
            msci = np.median([fits.getdata(os.path.join(redux, re.sub('.fits', '_int.fits', f))) for f in obj_dict[obj]], axis = 0)
        else:
            msci = np.min([fits.getdata(os.path.join(redux, re.sub('.fits', '_int.fits', f))) for f in obj_dict[obj]], axis = 0)
        mhdu = fits.open(os.path.join(redux, re.sub('.fits', '_int.fits', obj_dict[obj][0])))
        mhdu[0].data = msci
        mhdu.writeto(os.path.join(redux, 'int_median', '%s_int.fits'%obj), overwrite = True)
        print('######################%s##############'%obj)
        for f in obj_dict[obj]:
            print(f)
            img_hdu = fits.open(os.path.join(redux, re.sub('.fits', '_int.fits', f)))
            img = img_hdu[0].data

            # get wavemap
            proctab = read_proctab(os.path.join(base_dir, 'kcwir.proc'))
            tab = search_proctab(proctab,
                    frame=img_hdu[0], target_type='MARC',
                    nearest=True)
            
            ofn = tab['filename'][0].replace('.fits', "_wavemap.fits")
            wavemap_file = os.path.join(redux, ofn)
            wavemap = fits.open(wavemap_file)[0].data

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                ratio = np.nanmedian(img[wavemap > 0] / msci[wavemap > 0])
            #    _,ratio,__=sigma_clipped_stats(img[wavemap > 0] / msci[wavemap > 0], sigma = 3, maxiters=1,grow=2)

            crs = img - msci#*ratio
            #print("rescaling msci by a factor of %f"%ratio)
            crs_ratio=np.abs((img - msci)/msci)
            crmsk = (sigma_clip(crs, sigma = 3, maxiters=1, grow=2).mask & (crs_ratio>0.6)).astype(float)
            crmsk_surr = sigma_clip(crs, sigma = 3, maxiters=1, grow=10).mask.astype(float)


            idx = np.where((crmsk < 1e-4) & (crmsk_surr > 1e-4) & (wavemap > 0) )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                ratio2 = np.nanmedian(img[idx]/msci[idx])
            
            print(ratio, ratio2)

            crmsk_med = crmsk * msci * ratio2
            crmsk_med[crmsk_med < 1e-4] = 0

            primary = fits.PrimaryHDU(crmsk, header=img_hdu[0].header)
            med = fits.ImageHDU(crmsk_med, name='MEDSCI')
            sur = fits.ImageHDU(crmsk_surr, name = 'CRMSK_SUR')
            hdul = fits.HDUList([primary, med, sur])
            hdul.writeto(os.path.join(redux, re.sub('.fits', '_crmsk.fits', f)), overwrite=True)



def main():
    arg_parser = parser_init()
    args = arg_parser.parse_args()
    create_crmsk(**vars(args))

#Call if run from command-line
if __name__ == "__main__":
    main()
