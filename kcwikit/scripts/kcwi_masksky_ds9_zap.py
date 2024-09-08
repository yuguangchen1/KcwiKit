#!/usr/bin/env python

"""Creates mask image from ds9 region file.

Args:
    imagename (string): The name of a *_intf.fits image
    regionname (string): The name of a ds9 region file

Returns:
    None

Note:
    To use this routine, process your data with kcwi_stage4flat.pro.
    Then display the target *_intf.fits file in ds9.
    Use region shapes to indicate non-sky pixels in image (box, circle, etc.).
    Write out ds9 region file (*.reg).
    Run this routine:

    python ~/kderp/devel/kcwi_masksky_ds9_zap.py kb180101_00111_wlimg.fits kb180101_00111.reg

    (replace paths/filenames with your local paths/filenames)

    This should create kb180101_00111_smsk.fits, which will be used when you
    run kcwi_stage5sky.pro.
"""
try:
    import astropy
except ImportError:
    print("Please install astropy: required for image I/O")
    quit()
try:
    import pyregion
except ImportError:
    print("Please install pyregion: required for DS9 region I/O")
    quit()
import numpy as np
import sys
import os
import pdb

def main():

    # check args
    narg=len(sys.argv)

    # should be three (including routine name)
    if narg != 3 and narg != 2:
        print("Usage: python kcwi_masksky_ds9_zap.py <imagename> <regionname>")
        print("imagename : used for array dimensions and filename purposes, ")
        print("            must be an _wlimg image.")
        print("regionname: name of region file containing ds9 mask regions")
        print("            (typically a .reg)")
        exit()

    # read arg values
    imfname=sys.argv[1]
    if narg == 3:
        regfname=sys.argv[2]
    else:
        regfname = imfname.replace('_wlimg.fits', '.reg')


    # make sure it's an _intf image
    #if '_intf.fits' not in imfname:
    #    print("imagename must be _intf.fits image")
    #    exit()

    # do inputs exist?
    if os.path.exists(imfname) == False:
        print("File does not exist: "+imfname)
        exit()

    if os.path.exists(regfname) == False:
        print("File does not exist: "+regfname)
        exit()

    # create output mask image name
    outfile=imfname.replace("_wlimg.fits", "_zap_mask.fits")
    print("Creating: "+outfile)

    # load the header from the pointed-to image.
    hdu_list = astropy.io.fits.open(imfname)
    mhdu = hdu_list[1]
    edgemask = mhdu.data > 1
    header= hdu_list[0].header

    # determine the size of the array
    shape = (header["NAXIS1"], header["NAXIS2"])

    try:
        # load in the region file
        r = pyregion.open(regfname).as_imagecoord(header)
    except ValueError:
        print('Region File is Empty (contains no objects)')
        sys.exit(1)

    region_mask  = r.get_mask(hdu=mhdu)
    allmask = np.zeros_like(mhdu.data)
    allmask[region_mask] = 2
    allmask[edgemask] = 1

    # write out the mask
    hdu = astropy.io.fits.PrimaryHDU(allmask, header = mhdu.header)
    hdu.writeto(outfile, overwrite=True)

    print("Done.")

#Call if run from command-line
if __name__ == "__main__":
    main()
