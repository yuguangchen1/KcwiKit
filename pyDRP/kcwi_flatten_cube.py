import argparse
import os
from datetime import datetime

import pdb

#Third-party Imports
from astropy.io import fits
import numpy as np

def parser_init():
    description = "Flatten or de-flatten cube. This is a replica of kcwi_flatten_cube.pro in the IDL DRP."

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        'file',
        type=str,
        nargs='+',
        help='File name')
    parser.add_argument('-r', '--reverse', dest='reverse',
            action='store_true',
            help='Reverse the flattening')
    parser.add_argument('-t', '--trim', dest='trim',
            action='store_const', const=0, default=0,
            help='Pixels to be trimmed')

    args = parser.parse_args()
    return args




def main(file, reverse=False, trim=0):

    pre = 'kcwi_flatten_cube.py'

    # safty to convert string to list
    if isinstance(file, str):
        file = [file]

    for cfile in file:
        # Check file exist
        if os.path.isfile(cfile):

            print(pre + 'Operating - ' + cfile)
            # Read cube
            hdl_cub = fits.open(cfile)
            # Number of extensions
            n_ext = len(hdl_cub)
            
            # Loop through all extensions
            for i_ext, hdu_cub in enumerate(hdl_cub):

                # 3D?
                if len(hdu_cub.shape)==3:

                    # Constructing a 3d from 2d
                    if reverse:
                        #
                        # Read in 2d file
                        tfil = cfile.replace('.fits', '_2d.fits')
                        # check 2d
                        if os.path.exists(tfil):
                            # read in 2d file
                            hdu_img = fits.open(tfil)[i_ext]
                            img = hdu_img.data
                            thdr = hdu_img.header

                            # unpack slices
                            for i in range(hdu_cub.shape[2]):
                                ix0 = int(trim)
                                ix1 = int(hdu_cub.shape[1] - trim)
                                ox0 = int(i * hdu_cub.shape[1] + trim)
                                ox1 = int((i + 1) * hdu_cub.shape[1] - trim)
                                hdu_cub.data[:, ix0:ix1, i] = img[:, ox0:ox1]
                        else:
                            print(pre + ': 2D version not found - ' + tfil + ' Ext [{0:d}]'.format(i_ext))
                            return
                    # flattening a 3d to 2d
                    else:
                        # output image
                        oim = np.zeros((hdu_cub.shape[0], hdu_cub.shape[2] * hdu_cub.shape[1]))
                        
                        # pack slices
                        for i in range(hdu_cub.shape[2]):
                            ix0 = trim
                            ix1 = hdu_cub.shape[1] - trim
                            ox0 = i * hdu_cub.shape[1] + trim
                            ox1 = (i + 1) * hdu_cub.shape[1] - trim
                            oim[:, ox0:ox1] = hdu_cub.data[:, ix0:ix1, i]

                        # get wavelength values
                        hdr = hdu_cub.header
                        if 'CRVAL3' in hdr:
                            w0 = hdr['CRVAL3']
                            dw = hdr['CD3_3']
                            crpixw = hdr['CRPIX3']
                            ctypew = hdr['CTYPE3']
                            # set spatial scale
                            s0 = 0.
                            ds = 24. / (hdu_cub.shape[2] * hdu_cub.shape[1])
                            crpixs = 1.
                            # update header
                            hdr['HISTORY'] = '  ' + pre + datetime.now().ctime()
                            hdr['INCUBEF'] = (cfile, 'Input cube filename')
                            # remove old WCS
                            del hdr['NAXIS3']
                            del hdr['CTYPE1']
                            del hdr['CTYPE2']
                            del hdr['CTYPE3']
                            del hdr['CUNIT1']
                            del hdr['CUNIT2']
                            del hdr['CUNIT3']
                            del hdr['CNAME1']
                            del hdr['CNAME2']
                            del hdr['CNAME3']
                            del hdr['CRVAL1']
                            del hdr['CRVAL2']
                            del hdr['CRVAL3']
                            del hdr['CRPIX1']
                            del hdr['CRPIX2']
                            del hdr['CRPIX3']
                            del hdr['CD1_1']
                            del hdr['CD1_2']
                            del hdr['CD2_1']
                            del hdr['CD2_2']
                            del hdr['CD3_3']

                            # set wavelength axis WCS values
                            hdr['WCSDIM'] = 2
                            hdr['CTYPE1'] = ('SPATIAL', 'SLICE')
                            hdr['CUNIT1'] = ('slu', 'SLICE units')
                            hdr['CNAME1'] = ('KCWI SLICE', 'SLICE name')
                            hdr['CRVAL1'] = (s0, 'SLICE zeropoint')
                            hdr['CRPIX1'] = (crpixs, 'SLICE reference pixel')
                            hdr['CDELT1'] = (ds, 'SLICE per pixel')
                            hdr['CTYPE2'] = (ctypew, 'Wavelengh')
                            hdr['CUNIT2'] = ('Angstrom', 'Wavelength units')
                            hdr['CNAME2'] = ('KCWI 2D Wavelength', 'Wavelength name')
                            hdr['CRVAL2'] = (w0, 'Wavelength zeropoint')
                            hdr['CRPIX2'] = (crpixw, 'Wavelength reference pixel')
                            hdr['CDELT2'] = (dw, 'Wavelength Anstroms per pixel')

                        else:
                            del hdr['NAXIS3']

                        hdu_cub.data = oim

                    #hdl_cub[i_ext] = hdu_cub    

                else:
                    print(pre + ': Data not in 3D - ' + cfile + ' Ext [{0:d}]'.format(i_ext))



            # Write
            if reverse:
                # get 3D file name
                ofil = tfil.replace('_2d.fits', '_3d.fits')
            else:
                # get 2d file name
                ofil = cfile.replace('.fits', '_2d.fits')

            
            hdl_cub.writeto(ofil, overwrite=True)
            
    return
            





if __name__ == '__main__':
    args = parser_init()
    main(**(vars(args)))
    





