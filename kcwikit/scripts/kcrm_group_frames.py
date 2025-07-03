#!/usr/bin/env python

"""
Group red channel frames for median cosmic ray rejection. 
"""

import argparse
import json
from astropy.io import ascii
import numpy as np
import os

def parser_init():
    """Create command-line argument parser for this script."""
    parser = argparse.ArgumentParser(description="Group red channel frames for median cosmic ray rejection. ")
    parser.add_argument(
        'logfile',
        type=str,
        help='Filename of the log CSV table',
        nargs=1
        )
    parser.add_argument(
        'outfile',
        type=str,
        help='Output JSON dictionary file',
        nargs=1
        )
    parser.add_argument(
        '-c', '--continuous',
        help='Split non continuous file numbers',
        action='store_true'
        )
    parser.add_argument(
        '-d', '--display',
        help='Print json',
        action='store_true'
        )
    parser.add_argument(
        '-i', '--check_int',
        help='Check if the *.int file exists',
        action='store_true'
        )
    parser.add_argument(
        '-b', '--blue',
        help='process blue channel',
        action='store_true'
        )
    parser.add_argument(
        '--redux_dir',
        help='Location of the redux directory',
        nargs=1,
        default='./redux'
        )
    return parser


def group_frames(logfile, outfile, continuous=False, display=False, check_int=False,blue=False, redux_dir = './redux'):
    # group frames based on PA

    if isinstance(logfile, list):
        logfile = logfile[0]

    if isinstance(outfile, list):
        outfile = outfile[0]

    log = ascii.read(logfile, format = 'csv')
    if blue:
        kcrm = log[log['Camera'] == 'BLUE']
    else:
        kcrm = log[log['Camera'] == 'RED']

    kcrm.show_in_notebook()

    obj_dict = {}

    for i, row in enumerate(kcrm):

        #if row['Camera'] != 'RED':
        #    continue

        if row['ImType'] != 'OBJECT':
            continue

        if check_int:
            #the redux_dir is a list if --redux_dir is specified
            if isinstance(redux_dir, list):
                redux_dir = redux_dir[0]
            intfn = os.path.join(redux_dir, row['No.'].replace('.fits', '_int.fits'))
            if not os.path.isfile(intfn):
                continue

        key = '{0}_{1:d}_{2}_{3:d}'.format(row['Object'], int(np.round(row['Sky PA'])), row['State Name'], int(np.round(row['Exptime'])))

        if key not in obj_dict.keys():
            obj_dict[key] = []
        obj_dict[key].append(row['No.'])


    if continuous:  # requires file numbers to be continuous
        obj_dict_cont = {}
        for key in obj_dict.keys():
            fns = obj_dict[key]
            fnos = [int(fn.split('_')[1][:5]) for fn in obj_dict[key]]

            if len(fnos) <= 1:
                obj_dict_cont[key] = obj_dict[key]
                continue

            results_no = []
            tmp_no = [fnos[0]]
            results_fn = []
            tmp_fn = [fns[0]]
            for i in range(1, len(fnos)):
                if fnos[i] == fnos[i - 1] + 1:
                    tmp_no.append(fnos[i])
                    tmp_fn.append(fns[i])
                else:
                    results_no.append(tmp_no)
                    results_fn.append(tmp_fn)
                    tmp_no = [fnos[i]]
                    tmp_fn = [fns[i]]
            results_no.append(tmp_no)
            results_fn.append(tmp_fn)

            if len(results_no) == 1:
                # do nothing
                obj_dict_cont[key] = obj_dict[key]
            else:
                for i, result_fn in enumerate(results_fn):
                    obj_dict_cont[key+'_{0:d}'.format(i)] = result_fn

        obj_dict = obj_dict_cont

    if display:
        print(obj_dict)

    json.dump(obj_dict, open(outfile, 'w'), indent=4)

    return


def main():
    arg_parser = parser_init()
    args = arg_parser.parse_args()
    group_frames(**vars(args))

#Call if run from command-line
if __name__ == "__main__":
    main()


