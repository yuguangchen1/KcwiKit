# Cosmic Ray Mask Construction

## Instructions for cosmic ray mask (`*_crmsk.fits`) construction for CorrectDefects stage of KCWI\_DRP

1. Find CR in QFitsView or by some other means and note approximate wavelength

2. Find CR in `*_intk.fits` file (assuming we're in the redux directory)

`ds9 kb210811_00031_intk.fits kb210811_00019_wavemap.fits -lock frame image -zscale -lock scale -zoom to fit -lock crosshair image -mode crosshair &`

3. Find CR in `*crr.fits` file

`ds9 kb210811_00031_crr.fits -zscale kb210811_00031_crr.fits\[FLAGS\] -minmax -lock frame image -zoom to fit -mode region -regions shape box &`

4. Draw region around CR

5. Save region mask with `*_crr.reg` extension (e.g. `kb210811_00031_crr.reg`)

6. Convert regions file to fits file

`python ~/Software/kcwi/pyDRP/kcwi_reg2fits.py kb210811_00031_crr.fits kb210811_00031_crr.reg kb210811_00031_crmsk.fits`

Reminder: `python kcwi_reg2fits.py <imagename> <regionname> <outfname>`

7. Re-reduce image

`reduce_kcwi -f kb210811_00031.fits`
