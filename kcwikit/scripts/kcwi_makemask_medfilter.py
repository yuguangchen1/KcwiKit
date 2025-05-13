#!/usr/bin/env python

"""
Converting region files to sky masks specifically for median filtering. Similar to kcwi_makemask.pro.
"""

import glob, os, sys, re
import argparse
import json

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
    parser.add_argument(
        '-r', '--red',
        help='Act on red channel exposures.',
        action='store_true'
        )
    parser.add_argument(
        '-f','--filename',
        type=str,
        help='Input json file that contains the file grouping',
        default=None,
        nargs='?'
        )
    return parser
def makemask(dir,filename=None, cubed=False, red=False):
    """Converting region files to sky masks specifically for median filtering.

    Args:
        dir (str): directory of which the region files are stored.

    Returns:
        None. Region files are converted to binary FITS files. 

    """
    if red:
        obj_dict = json.load(open(filename, 'r'))
        for obj in list(obj_dict.keys()):
            for f in obj_dict[obj]:
                print("producing mask for",f,"with",obj_dict[obj][0])
                regfn1=os.path.join(dir, re.sub('.fits', '_icube.thum.reg', obj_dict[obj][0]))
                regfn2=os.path.join(dir, re.sub('.fits', '_icube_2d.reg', obj_dict[obj][0]))
                
                if os.path.isfile(regfn1):
                    imgfn=os.path.join(dir, re.sub(".fits","_icube.thum.fits",f))
                    os.system(fileloc+"/kcwi_masksky_ds9_thum.py "+imgfn+" "+regfn1)
                else:
                    print(regfn1+" not found")
                
                if os.path.isfile(regfn2):
                    imgfn=os.path.join(dir, re.sub(".fits","_icube_2d.fits",f))
                    os.system(fileloc+"/kcwi_masksky_ds9_2d.py "+imgfn+" "+regfn2)
                else:
                    print(regfn2+" not found")
    else:
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
