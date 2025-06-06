import glob, os, sys
import argparse

import pathlib
fileloc = str(pathlib.Path(__file__).parent.resolve())

def parser_init():
    """Create command-line argument parser for this script."""
    parser = argparse.ArgumentParser(description="Convert region masks to binary masks for sky subtraction.")
    parser.add_argument(
        'dir',
        type=str,
        help='Directory name in which region files are located.',
        default='./',
        nargs='?'
        )
    return parser

def makemask(dir):
    """Converting region files to sky masks.

    Args:
        dir (str): directory of which the region files are stored.

    Returns:
        None. Region files are converted to binary FITS files. 

    """

    for regfn in glob.glob(dir+"/*.reg"):
        imgfn=regfn.replace('.reg','_intf.fits')
        os.system(fileloc+"/kcwi_masksky_ds9.py "+imgfn+" "+regfn)

    print('Sky Masks Complete')

def main():
    arg_parser = parser_init()
    args = arg_parser.parse_args()
    makemask(**vars(args))

if __name__ == "__main__":
    main()