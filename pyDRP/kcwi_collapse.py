import argparse
import os

#Third-party Imports
from astropy.io import fits
import numpy as np

import pdb

def parser_init():
    description = 'Create 2D pseudo-white-light image.'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('file', type=str, nargs='+', help='File name')
    parser.add_argument('--wavebin', dest='wavebin', type=float,
            nargs=2, help='Lower and upper boundaries of the wavelength bin')
    parser.add_argument('--usepix', dest='use_pix',
            action='store_true', help='Use pixel index for the wavelength direction')

    args = parser.parse_args()
    return args


def main(file, wavebin=None, use_pix=False):

    pre = 'kcwi_collapse.py'

    # use_pix cannot be True if wavebin is not set
    if wavebin is None:
        use_pix = False

    # safety check
    if isinstance(file, str):
        file = [file]

    for fn in file:
        # check file exist
        if os.path.isfile(fn):
            print(pre + ': Operating - ' + fn)

            # Read cube
            hdl = fits.open(fn)
            ohdl = hdl.copy()
            # number of extensions
            n_ext = len(hdl)

            # Loop through all extensions
            for i_ext, hdu in enumerate(hdl):
                # make a hard copy
                ohdl[i_ext] = hdu.copy()
                ohdl[i_ext].header = hdu.header

                # read cube
                cube = hdu.data
                hdr = hdu.header
                mask = hdl['MASK'].data
                flag = hdl['FLAGS'].data
                mask = mask | (flag != 0)

                # check if header contains WCS
                if ('CRVAL3' not in hdr) and (i_ext==0):
                    print(pre + ': File has no WCS - ' + fn)
                    return

                if 'CRVAL3' in hdr:
                    # get wavelength
                    zz = np.arange(hdr['NAXIS3'])
                    wave = (np.arange(hdr['NAXIS3']) - hdr['CRPIX3'] + 1) * hdr['CD3_3'] + hdr['CRVAL3']
                    
                    # Find wavelength range
                    if wavebin is None:
                        wrange = [hdr['WAVGOOD0'], np.min((hdr['WAVGOOD1'], 5500))]
                    else:
                        wrange = wavebin

                    # wavelength index
                    if use_pix:
                        qwave = (zz > wrange[0]) & (zz < wrange[1])
                    else:
                        qwave = (wave > wrange[0]) & (wave < wrange[1])

                    
                # Collapse
                tmp_cube = np.array(cube, dtype=float)
                tmp_cube[mask != 0] = np.nan
                if 'EXTNAME' in hdr:
                    if hdr['EXTNAME']=='MASK':
                        # binary mask
                        img = np.nansum(tmp_cube[qwave, :, :], axis=0)
                    elif hdr['EXTNAME']=='UNCERT':
                        # error
                        img = np.sqrt(np.nansum(tmp_cube[qwave, :, :]**2, axis=0)) /\
                                np.sum((mask[qwave, :, :] == 0), axis=0)
                    elif hdr['EXTNAME']=='FLAGS':
                        # Flag
                        img = np.nanmax(tmp_cube[qwave, :, :], axis=0)
                else:
                    # Intensity
                    img = np.nanmean(tmp_cube[qwave, :, :], axis=0)


                # header
                ohdr = ohdl[i_ext].header
                if 'CRVAL3' in ohdr:
                    del hdr['NAXIS3']
                    del hdr['CD3_3']
                    del hdr['CTYPE3']
                    del hdr['CUNIT3']
                    del hdr['CNAME3']
                    del hdr['CRVAL3']
                    del hdr['CRPIX3']

                ohdl[i_ext].data = img

            ofn = fn.replace('.fits', '.thum.fits')
            ohdl.writeto(ofn, overwrite=True)


        else:
            print(pre + 'File not found - ' + fn)


if __name__ == '__main__':
    args = parser_init()
    main(**vars(args))





