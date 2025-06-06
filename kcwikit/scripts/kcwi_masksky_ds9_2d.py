#!/usr/bin/env python

from astropy.io import fits as pf
from regions import Regions

import numpy as np
import sys
import os

def main():

    """Creates mask image from ds9 region file.

    Args:
        imagename (string): The name of a *_icube_2d.fits image
        regionname (string): The name of a ds9 region file

    Returns:
        None

    Note:
        To use this routine, create flattened data cubes with `kcwi_flatten_cube`.
        Then display the target *_icube_2d.fits file in ds9.
        Use region shapes to indicate non-sky pixels in image (box, circle, etc.).
        Write out ds9 region file (*.reg).
        Run this routine:

        kcwi_masksky_ds9_2d.py kb180101_00111_icube_2d.fits ds9.reg

        (replace paths/filenames with your local paths/filenames)

        This should create kb180101_00111_icube_2d.mask.fits, which will be used 
        for median fitlering
    """

    # check args
    narg=len(sys.argv)

    # should be three (including routine name)
    if narg != 3:
        print("Usage: python kcwi_masksky_ds9_2d.py <imagename> <regionname>")
        print("imagename : used for array dimensions and filename purposes, ")
        print("            must be an _icube_2d image.")
        print("regionname: name of region file containing ds9 mask regions")
        print("            (typically a .reg)")
        exit()

    # read arg values
    imfname=sys.argv[1]
    regfname=sys.argv[2]

    # do inputs exist?
    if os.path.exists(imfname) == False:
        print("File does not exist: "+imfname)
        exit()

    if os.path.exists(regfname) == False:
        print("File does not exist: "+regfname)
        exit()

    # create output mask image name
    outfile=imfname.replace(".fits", ".mask.fits")
    print("Creating: "+outfile)

    # load the header from the pointed-to image.
    hdu_list = pf.open(imfname)
    header= hdu_list[0].header
    shape = hdu_list[0].shape

    # load the region file
    with open(regfname, 'r') as f:
        # Read it out as a string
        regstr = f.read()
        
        # Check if the region file is in physical coordinates
        if 'physical' in regstr:
            print("[Warning] 'physical' coordinates not supported by regions. Replacing with 'image'")
            regstr = regstr.replace('physical', 'image')

        r = Regions.parse(regstr, format='ds9')
        mask = None
        for region in r.regions:
            if mask is None:
                mask = region.to_mask().to_image(shape).astype(bool)
            else:
                mask = mask | region.to_mask().to_image(shape).astype(bool)

    # write out the mask
    hdu = pf.PrimaryHDU(np.uint8(mask))
    hdu.writeto(outfile, overwrite=True)

    print("Done.")

if __name__ == "__main__":
    main()
