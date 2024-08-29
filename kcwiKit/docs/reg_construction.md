# Region File Creation
## A guide for making regions files for the KCWI DRP.
### Last updated: 9/19/2022

For some quick examples see `mizar:/scr/kbss-kcwi/tier1drp/kcwi_jul21/2021jul05/redux`

## Sky Subtraction Regions (`.reg`)

These regions files are constructed using the `*intk.fits` images for a second round of sky subtraction with continuum objects masked.

Use `_intk.fits` files because this is where sky subtraction failed and mark regions to *exclude* i.e. lines/continuum.
Save in **Physical** coordinates (that's all you have at this point!).

**Sample command:**
  - `for i in {71..92}; do ds9 kb220528_000${i}_intk.fits -height 1000 -width 1100 -cmap cool -zscale -regions save kb220528_000${i}.reg -mode region -zoom to fit -regions shape box &; done`
    - Be careful with the `-save` command - this overwrites existing regions files!
    - Some other colormaps that might work: sls, green

**To go individually:**
  - `touch kb210705_00030.reg`
  - `ds9 kb210705_00030_intk.fits -height 1000 -width 1100 -cmap sls -zscale -regions load kb210705_00030.reg -mode region -zoom to fit -regions shape box &`
  - Save region file
  - `for i in {31..53}; do cp kb210705_00030.reg kb210705_000$i.reg; done`
  - `ds9 kb210705_00031_intk.fits -height 1000 -width 1100 -cmap sls -zscale -regions load kb210705_00031.reg -mode region -zoom to fit -regions shape box &`
  - `^31^32^:G` (for same command, global replace)
  - ...


## Median Filtering Regions

### Continuum Objects (`_icube.thum.reg`)

The recently constructed datacubes are flattened to form a whitelight image and continuum objects are masked for median filtering.

For standard stars and science targets, draw regions (ellipses) around continuum objects and save in **Physical** coordinates.

**Sample command:**
  - `for i in {50..55}; do ds9 kb210705_000${i}_icube.thum.fits -regions save kb210705_000${i}_icube.thum.reg -height 1200 -width 800 -zscale -mode region -regions shape ellipse -zoom to fit &; done`

**To go individually:**
  - `touch kb210811_00026_icube.thum.reg`
  - `ds9 kb210811_00026_icube.thum.fits -regions load kb210811_00026_icube.thum.reg -height 1200 -width 800 -zscale -mode region -regions shape ellipse -zoom to fit &`
  - Save region file
  - `for i in {35..40}; do cp kb210811_00026_icube.thum.reg kb210811_000${i}_icube.thum.reg; done`
  - `ds9 kb210811_00035_icube.thum.fits -regions load kb210811_00035_icube.thum.reg -height 1200 -width 800 -zscale -mode region -regions shape ellipse -zoom to fit &`
  - `^35^36^:G` (for same command, global replace)
  - ...


### Line Emission (`_icube_2d.reg`)

The science datacubes are deprojected back into a 2D image and line emission (typically Ly&alpha;) is marked.

For standard stars and science targets, find emission lines and draw boxes around emission. Don't worry about the continuum emission since the `_icube.thum.reg` files will take care of that. Save regions in **Physical** coordinates as well.

**Sample command:**
  - `for i in {50..55}; do ds9 kb210705_000${i}_icube_2d.fits -regions save kb210705_000${i}_icube_2d.reg -height 1200 -width 2000 -zscale -mode region -smooth function boxcar -smooth radius 3 -smooth yes -regions shape box -zoom to fit &; done`

**To go individually:**
  - `touch kb210811_00026_icube_2d.reg`
  - `ds9 kb210811_00026_icube_2d.fits -regions load kb210811_00026_icube_2d.reg -height 1200 -width 2000 -zscale -mode region -smooth function boxcar -smooth radius 3 -smooth yes -regions shape box -zoom to fit &`
  - Save region file
  - `for i in {35..40}; do cp kb210811_00026_icube_2d.reg kb210811_000${i}_icube_2d.reg; done`
  - `ds9 kb210811_00035_icube_2d.fits -regions load kb210811_00035_icube_2d.reg -height 1200 -width 2000 -zscale -mode region -smooth function boxcar -smooth radius 3 -smooth yes -regions shape box -zoom to fit &`
  - `^35^36^:G` (for same command, global replace)
  - ...

## Miscellaneous

If images are not sequential, use nested brace expansion: `for i in {{099..102},{131,132,138}}; ...`
