# Complete Reference of the Parameter File

Here is the complete reference on the parameter file. However, not all of the parameters are commonly used. An example for common usage is provided [here](../examples/q2343-BX610.par). 

Some parameters can be overridden by manually setting them when calling the functions. For example, `kcwi.kcwi_align(..., search_size=15)` will manually increase the alignment box to 15 arcsec. 

```
dimension 100 100       # set up the number of pixels in the X and Y directions
orientation 0           # PA of the frame, starting north and increasing counter-clockwise
xpix 0.3                # pixel scale (arcsec) in the X direction
ypix 0.3                # pixel scale (arcsec) in the Y direction
wavebin                 # wavelength range in which the white-light images are created, most useful for the red channel
drizzle                 # drizzle factor, default 0.7
align_box 45 37 65 65   # x0, y0, x1, y1 of the alignment box
align_dimension         # number of pixels specifically for the alignment frame
align_xpix              # pixel scale (arcsec) in the X direction, specifically for the alignment frame
align_ypix              # pixel scale (arcsec) in the Y direction, specifically for the alignment frame
align_orientation       # PA of the frame, specifically for the alignment frame
align_search_size       # half width of the search size, default 10 pixels
align_conv_filter       # size of the convolution filter when searching for the local maximum in cross correlation, default 2 pixels
align_upfactor          # upscaling factor for finer local maximum searching, default 10 times.
align_ad                # RA, Dec of the center of the alignment frame, default center of the first frame
background_subtraction  # conducting background subtraction for alignment? Default false
background_level        # background level for background subtraction, default median of the frame
stack_dimension         # number of pixels, specifically for the stacking frame
stack_xpix              # X pixel scale, specifically for the stacking frame
stack_ypix              # Y pixel scale, specifically for the stacking frame
stack_orientation       # PA of the frmae, specifically for the stacking frame
stack_ad                # RA, Dec of the center of the stacking frame, default center of the first frame
wave_ref                # CRPIX3 and CRVAL3 of the wavelength grid of the final cube. Default equal to the first frame
nwave                   # number of wavelength pixels of the final cube. Default equal to the first frame
dwave                   # pixel scale (Angstrom) in the wavelength direction of the final cube. Default equal to the first frame
ref_xy 51 52            # the X, Y pixel location of the reference object for astrometry
ref_ad 356.53929 12.822 # the RA, Dec (in degrees) of the reference object for astrometry
ref_fn /path/to/image   # an external reference FITS image for astrometry
ref_search_size         # search size for astrometry, see align_search_size
ref_conv_filter         # convolution filter size for astrometry, see align_conv_filter
ref_upfactor            # upscaling factor for astrometry, see align_upfactor
ref_nocrl               # if set to nonzero, astrometry is skipped
```