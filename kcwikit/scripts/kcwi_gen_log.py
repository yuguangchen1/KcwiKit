#!/usr/bin/env python

"""
Generate log files based on FITS headers. 
"""

import argparse
from astropy.io import fits
from astropy.table import Table

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
        help='FITS files',
        nargs='*',
        required=True
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

def create_log(outfile, filename=[], RED=False, BLUE=False, display=False):
    
    if not isinstance(filename, list):
        filename = [filename]

    if isinstance(outfile, list):
        outfile = outfile[0]

    
    nos = []
    cameras = []
    objects = []
    imtypes = []
    groupids = []
    statenames = []
    targets = []
    lsts = []
    uts = []
    ras = []
    decs = []
    els = []
    airmasses = []
    exptimes = []
    binnings = []
    ampmodes = []
    ccdspeeds = []
    ccdmodes = []
    gratings = []
    filters = []
    cwaves = []
    pwaves = []
    nsmasks = []
    ifus = []
    rotmodes = []
    skypas = []

    
    for i, fn in enumerate(filename):
        hdr = fits.getheader(fn)

        if RED:
            if hdr['CAMERA'] != 'RED':
                continue
        if BLUE:
            if hdr['CAMERA'] != 'BLUE':
                continue

        nos.append(fn)
        cameras.append(hdr['CAMERA'])
        objects.append(hdr['OBJECT'])
        imtypes.append(hdr['IMTYPE'])
        groupids.append(hdr['GROUPID'])
        statenames.append(hdr['STATENAM'])
        targets.append(hdr['TARGNAME'])
        lsts.append(hdr['LST'])
        uts.append(hdr['UT'])
        ras.append(hdr['RA'])
        decs.append(hdr['DEC'])
        els.append(hdr['EL'])
        airmasses.append(hdr['AIRMASS'])
        exptimes.append(hdr['TTIME'])
        binnings.append(hdr['BINNING'])
        ampmodes.append(hdr['AMPMODE'])
        ccdspeeds.append(hdr['CCDSPEED'])
        ccdmodes.append(hdr['CCDMODE'])

        if cameras[-1] == 'BLUE':
            gratings.append(hdr['BGRATNAM'])
            filters.append(hdr['BFILTNAM'])
            cwaves.append(hdr['BCWAVE'])
            pwaves.append(hdr['BPWAVE'])
            nsmasks.append(hdr['BNASNAM'])
        else:
            gratings.append(hdr['RGRATNAM'])
            filters.append('')
            cwaves.append(hdr['RCWAVE'])
            pwaves.append(hdr['RPWAVE'])
            nsmasks.append(hdr['RNASNAM'])

        ifus.append(hdr['IFUNAM'])
        rotmodes.append(hdr['ROTMODE'])
        skypas.append(hdr['ROTPOSN'] + hdr['ROTREFAN'])

    t = Table( (nos, cameras, objects, imtypes, groupids, statenames, targets, lsts, uts, ras, decs, els, airmasses, exptimes, binnings, ampmodes, ccdspeeds, ccdmodes, gratings, filters, cwaves, nsmasks, ifus, rotmodes, skypas), \
        names=('No.', 'Camera', 'Object', 'ImType', 'GroupID', 'State Name', 'Target', 'LST', 'UT', 'RA', 'DEC', 'ELS', 'AIRMASS', 'Exptime', 'Binning', 'Ampmode', 'CCD Speed', 'CCD Mode', 'Grating', 'Filter', 'Central wave', 'NS mask', 'IFU', 'Rot Mode', 'Sky PA'))
    
    if display:
        print(t)

    t.write(outfile, overwrite=True, format='csv')


def main():
    arg_parser = parser_init()
    args = arg_parser.parse_args()
    create_log(**vars(args))

#Call if run from command-line
if __name__ == "__main__":
    main()
