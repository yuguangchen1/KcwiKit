#!/usr/bin/env python

from astropy.io import fits as pf
from kcwikit.kcwi.kcwi_regions import regions_to_mask

import numpy as np
import sys
import os

def main():
    """Creates an astroscrappy *protection* mask from a ds9 region file.

    Display the *_intf.fits frame in ds9 and draw regions over the bright,
    variable sources (QSO core, strong emission lines) that astroscrappy
    mistakes for cosmic rays. Write out the region file, then run:

        * ``kcwi_masksky_ds9_bcrr kb180101_00111_intf.fits kb180101_00111_bcrmsk.reg``

    This creates kb180101_00111_intf_bcrmsk.fits (1 = protect), which the
    KCWI_DRP RemoveCosmicRays primitive passes to astroscrappy as ``inmask``
    when ``bcrmsk = True`` is set in kcwi.cfg. ``kcwi_automask --bcrmsk`` makes
    the same two files automatically.

    Args:
        imagename (string): The name of a *_intf.fits image
        regionname (string): The name of a ds9 region file

    Returns:
        None

    """
    # check args
    narg = len(sys.argv)

    # should be three (including routine name)
    if narg != 3:
        print("Usage: kcwi_masksky_ds9 <imagename> <regionname>")
        print("imagename : used for array dimensions and filename purposes, ")
        print("            must be an kb*.fit image.")
        print("regionname: name of region file containing ds9 mask regions")
        print("            (typically a .reg)")
        exit()

    # read arg values
    imfname = sys.argv[1]
    regfname = sys.argv[2]

    # do inputs exist?
    if not os.path.exists(imfname):
        print("Image file does not exist: "+imfname)
        exit()

    if not os.path.exists(regfname):
        print("Region file does not exist: "+regfname)
        exit()

    # create output mask image name
    outfile = imfname.replace(".fits", "_bcrmsk.fits")
    print("Creating: "+outfile)

    # load the header from the pointed-to image.
    hdu_list = pf.open(imfname)
    header = hdu_list[0].header
    shape = hdu_list[0].shape

    # rasterise the region file (handles image/physical/sky coordinates,
    # '-' exclusions and shapes outside the image)
    mask = regions_to_mask(regfname, header=header, shape=shape)
    if not mask.any():
        print("[Warning] no region covers the image; writing an empty mask")

    # write out the mask
    hdu = pf.PrimaryHDU(np.uint8(mask))
    hdu.writeto(outfile, overwrite=True)

    print("Done.")


if __name__ == "__main__":
    main()
