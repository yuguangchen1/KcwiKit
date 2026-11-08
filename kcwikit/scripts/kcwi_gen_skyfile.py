#!/usr/bin/env python

"""
Generate the KCWI DRP sky-frame table (kcwi.sky) from the frame-grouping JSON
written by kcrm_group_frames.

Each science frame gets a row "frame skyframe skymask":
  * a group <obj> with a companion group <obj>_sky uses the first frame of
    <obj>_sky as an external sky frame and no sky mask (per the KCWI DRP
    instructions);
  * otherwise the frame is its own sky frame and the group's first-frame
    _smsk.fits mask is used.
"""

import argparse
import re
import json
from astropy.table import Table


def parser_init():
    """Create command-line argument parser for this script."""
    parser = argparse.ArgumentParser(
        description="Create the DRP sky-frame table (kcwi.sky) from the frame-grouping JSON")
    parser.add_argument(
        'outfile',
        type=str,
        help='Name of the output table (usually kcwi.sky)'
        )
    parser.add_argument(
        '-l', '--referencelog',
        type=str,
        help='JSON file documenting the frame grouping (from kcrm_group_frames)',
        required=True
        )
    parser.add_argument(
        '-d', '--display',
        help='Print the table',
        action='store_true'
        )
    return parser


def create_skyfile(outfile, referencelog, display=False):
    """Build the sky table and write it as a space-separated file."""
    with open(referencelog, 'r') as f:
        obj_dict = json.load(f)

    rows = []
    for obj, frames in obj_dict.items():
        if obj.endswith('_sky'):
            continue                       # sky groups only serve their object group
        sky_group = obj_dict.get(obj + '_sky')
        for fn in frames:
            if sky_group:
                # external sky frame: no sky mask, per KCWI pipeline instructions
                rows.append((fn, sky_group[0], ''))
            else:
                rows.append((fn, fn, re.sub(r'\.fits$', '_smsk.fits', frames[0])))

    t = Table(rows=rows, names=('frame', 'skyframe', 'skymask'))
    if display:
        print(t)

    # MakeMasterSky reads this line by line with whitespace splitting and skips
    # '#' lines; a frame with an external sky frame simply has no third column
    with open(outfile, 'w') as f:
        f.write('# frame skyframe skymask\n')
        for frame, skyframe, skymask in rows:
            f.write('%s %s %s\n' % (frame, skyframe, skymask) if skymask else '%s %s\n' % (frame, skyframe))
    print('%s: %d frame(s) written' % (outfile, len(rows)))


def main():
    args = parser_init().parse_args()
    create_skyfile(**vars(args))


#Call if run from command-line
if __name__ == "__main__":
    main()
