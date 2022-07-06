#!/usr/bin/env python

"""
Converting region files to sky masks. Same as kcwi_makemask.pro.
"""

import glob, os, sys
import pdb

narg=len(sys.argv)
if narg==1:
    dir='./'
if narg==2:
    dir=sys.argv[1]
if narg>2:
    print('Incorrect number of arguments.')
    exit()

for regfn in glob.glob(dir+"/*_icube.thum.reg"):
    #pdb.set_trace()

    imgfn=regfn.replace('.reg','.fits')
    os.system("kcwi_masksky_ds9_thum.py "+imgfn+" "+regfn)

for regfn in glob.glob(dir+"/*_icube_2d.reg"):
    #pdb.set_trace()

    imgfn=regfn.replace('.reg','.fits')
    os.system("kcwi_masksky_ds9_2d.py "+imgfn+" "+regfn)

print('Sky Masks Complete')
