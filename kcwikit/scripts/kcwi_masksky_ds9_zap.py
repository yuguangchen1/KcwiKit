#!/usr/bin/env python

from astropy.io import fits as pf
from kcwikit.kcwi.kcwi_regions import regions_to_mask

import numpy as np
import sys
import os

"""Creates mask image from ds9 region file.

Args:
    imagename (string): The name of a *_wlimg.fits image
    regionname (string): The name of a ds9 region file

Returns:
    None

Note:
    To use this routine, process your data with KSkyWizard.
    Then display the target *_wlimg.fits file in ds9.
    Use region shapes to indicate non-sky pixels in image (box, circle, etc.).
    Write out ds9 region file (*.reg).
    Run this routine:

    kcwi_masksky_ds9_zap kb180101_00111_wlimg.fits kb180101_00111.reg

    (replace paths/filenames with your local paths/filenames)

    This should create kb180101_00111_zap_mask.fits, which will be used when you
    run KSkyWizard.
"""

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
    hdu_list = pf.open(imfname)
    mhdu = hdu_list[1]
    edgemask = mhdu.data > 1
    header= hdu_list[0].header
    shape = hdu_list[0].shape

    # rasterise the region file (handles image/physical/sky coordinates,
    # '-' exclusions and shapes outside the image)
    mask = regions_to_mask(regfname, header=header, shape=shape)
    if not mask.any():
        print("[Warning] no region covers the image; writing an empty mask")

    allmask = np.zeros_like(mhdu.data)
    allmask[mask] = 2
    allmask[edgemask] = 1

    # write out the mask
    hdu = pf.PrimaryHDU(allmask, header = mhdu.header)
    hdu.writeto(outfile, overwrite=True)

    print("Done.")

#Call if run from command-line
if __name__ == "__main__":
    main()
