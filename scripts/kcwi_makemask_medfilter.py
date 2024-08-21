#!/usr/bin/env python

"""
Converting region files to sky masks specifically for median filtering. Similar to kcwi_makemask.pro.
"""

import glob, os, sys
import argparse

import pathlib
fileloc = str(pathlib.Path(__file__).parent.resolve())

def parser_init():
    """Create command-line argument parser for this script."""
    parser = argparse.ArgumentParser(description="Convert region masks to binary masks.")
    parser.add_argument(
        'dir',
        type=str,
        help='Directory name in which region files are located.',
        default='./',
        nargs='?'
        )
    parser.add_argument(
        '-d', '--cubed',
        help='Act on DAR corrected cubed files.',
        action='store_true'
        )
    return parser

def makemask(dir, cubed=False):
    """Converting region files to sky masks specifically for median filtering.

    Args:
        dir (str): directory of which the region files are stored.

    Returns:
        None. Region files are converted to binary FITS files. 

    """

    if not cubed:
        search_for = dir + "/*_icube.thum.reg"
    else:
        search_for = dir + '/*_icubed.thum.reg'

    for regfn in glob.glob(search_for):
        imgfn=regfn.replace('.reg','.fits')
        os.system(fileloc+"/kcwi_masksky_ds9_thum.py "+imgfn+" "+regfn)

    if not cubed:
        search_for = dir + "/*_icube_2d.reg"
    else:
        search_for = dir + "/*_icubed_2d.reg"

    for regfn in glob.glob(search_for):
        imgfn=regfn.replace('.reg','.fits')
        os.system(fileloc+"/kcwi_masksky_ds9_2d.py "+imgfn+" "+regfn)

    print('Sky Masks Complete')

def main():
    arg_parser = parser_init()
    args = arg_parser.parse_args()
    makemask(**vars(args))

#Call if run from command-line
if __name__ == "__main__":
    main()
