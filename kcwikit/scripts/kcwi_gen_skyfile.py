#!/usr/bin/env python

"""
Generate log files based on FITS headers. 
"""

import argparse
import glob, os, sys, re
from astropy.io import fits
from astropy.table import Table
import json
def parser_init():
    """Create command-line argument parser for this script."""
    parser = argparse.ArgumentParser(description="Create log files based on FITS headers")
    parser.add_argument(
        'outfile',
        type=str,
        help='Filename of the output CSV table',
        nargs=1
        )
    parser.add_argument(
        '-f', '--filename',
        type=str,
        help='log file',
        nargs='*',
        required=False
        )
    parser.add_argument(
        '-l', '--referencelog',
        type=str,
        help='json file documenting the frame grouping',
        nargs='*',
        required=False
        )
    parser.add_argument(
        '-r', '--RED',
        help='Red channel only.',
        action='store_true'
        )
    parser.add_argument(
        '-b', '--BLUE',
        help='Blue channel only.',
        action='store_true'
        )
    parser.add_argument(
        '-d', '--display',
        help='Print table',
        action='store_true'
        )
    return parser

def create_log(outfile, filename=[], RED=False, BLUE=False, display=False,referencelog=[]):
    
    if not isinstance(filename, list):
        filename = [filename]

    if isinstance(outfile, list):
        outfile = outfile[0]

    if isinstance(referencelog, list):
        referencelog = referencelog[0]
    
    nos = []
    skyframe = []
    skymask = []
    #if RED:
    obj_dict = json.load(open(referencelog, 'r'))
    selected_keys = [k for k in obj_dict.keys() if "sky" in k]
    for obj in list(obj_dict.keys()):
        if obj+"_sky" in selected_keys:
            sky_group=obj_dict[obj+"_sky"]
        else:
            sky_group=None
        for fn in obj_dict[obj]:
            nos.append(fn)
            if sky_group is not None:
                skyframe.append(sky_group[0])
                skymask.append('')
            else:
                skyframe.append(fn)
                skymask.append(re.sub(".fits","_smsk.fits",obj_dict[obj][0]))
    
    t = Table( (nos, skyframe,skymask), \
        names=('frame', 'skyframe', 'skymask'))
    import pandas as pd
    if display:
        print(t)
    import pickle; pickle.dump(t, open('tmp.pickle', 'wb'))

    df = t.to_pandas()
    df.replace('', pd.NA, inplace=True)
    df.to_csv(outfile, index=False, na_rep='',sep=' ')
    #t.write(outfile, overwrite=True, format='ascii.basic',fill_values=[('', '')])


def main():
    arg_parser = parser_init()
    args = arg_parser.parse_args()
    create_log(**vars(args))

#Call if run from command-line
if __name__ == "__main__":
    main()
