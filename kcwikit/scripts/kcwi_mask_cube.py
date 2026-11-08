#!/usr/bin/env python
"""
Flag the pixels inside hand-drawn DS9 regions (<frame>_rm.reg) in the FLAGS
extension of a data cube, so that they are excluded from the final stack.

The region file is drawn on a white-light image of the cube; the same
2-D region is applied to every wavelength channel by setting bit 32 of
FLAGS. The cube itself is left unchanged and a copy is written:

    <frame>_zap_icubes.fits -> <frame>_zap_new_icubes.fits   (--red)
    <frame>_icubes.fits     -> <frame>_new_icubes.fits       (--blue)
"""

import os
import argparse
import numpy as np
from astropy.io import fits

from kcwikit.kcwi.kcwi_regions import regions_to_mask, header_to_2d

FLAG_REGION = 32   # FLAGS bit set for pixels inside the region file


def parser_init():
    parser = argparse.ArgumentParser(
        description='Set FLAGS bit %d inside the regions of <frame>_rm.reg for '
                    'every wavelength channel of the matching cube.' % FLAG_REGION)
    parser.add_argument('files', type=str, nargs='+',
                        help='Region file name(s); must end with _rm.reg')
    chan = parser.add_mutually_exclusive_group(required=True)
    chan.add_argument('-r', '--red', dest='red', action='store_true',
                      help='Red channel: apply to <frame>_zap_icubes.fits')
    chan.add_argument('-b', '--blue', dest='blue', action='store_true',
                      help='Blue channel: apply to <frame>_icubes.fits')
    return parser


def makemask(files, blue=False, red=False):
    """Apply each <frame>_rm.reg to the FLAGS extension of its cube.

    Args:
        files (list of str): region file names ending with _rm.reg.
        blue, red (bool): which channel's cube naming convention to use.
    """
    if red:
        in_suffix, out_suffix = '_zap_icubes.fits', '_zap_new_icubes.fits'
    else:
        in_suffix, out_suffix = '_icubes.fits', '_new_icubes.fits'

    for regfn in files:
        if not regfn.endswith('_rm.reg'):
            print('Skipping %s: region file name must end with _rm.reg' % regfn)
            continue
        if not os.path.isfile(regfn):
            print('Skipping %s: file not found' % regfn)
            continue
        cubefn = regfn.replace('_rm.reg', in_suffix)
        if not os.path.isfile(cubefn):
            print('Skipping %s: cube %s not found' % (regfn, cubefn))
            continue

        with fits.open(cubefn) as hdl:
            scihdr = hdl[0].header
            # only the spatial WCS matters for the region; no white-light image
            # (and no grating-specific wavelength window) is needed
            hdr2d = header_to_2d(scihdr)
            shape2d = hdl[0].data.shape[1:]
            mask2d = regions_to_mask(regfn, header=hdr2d, shape=shape2d)

            flags = hdl['FLAGS'].data
            # OR the bit in, so running twice (or on a _new cube) is idempotent
            flags[:, mask2d] |= np.array(FLAG_REGION, dtype=flags.dtype)
            hdl['FLAGS'].header['HISTORY'] = 'kcwi_mask_cube: bit %d set from %s' % (
                FLAG_REGION, os.path.basename(regfn))

            fdest = cubefn.replace(in_suffix, out_suffix)
            hdl.writeto(fdest, overwrite=True)
        print('%s: %d spaxel(s) flagged from %s -> %s' % (
            os.path.basename(cubefn), int(mask2d.sum()), os.path.basename(regfn),
            os.path.basename(fdest)))


def main():
    args = parser_init().parse_args()
    makemask(**vars(args))


if __name__ == '__main__':
    main()
