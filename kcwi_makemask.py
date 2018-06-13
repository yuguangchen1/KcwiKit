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

for regfn in glob.glob(dir+"/*.reg"):
    #pdb.set_trace()
    
    imgfn=regfn.replace('.reg','_intf.fits')
    os.system("/scr/yuguangchen/Soft/idl/kcwi/kcwi_masksky_ds9.py "+imgfn+" "+regfn) 



