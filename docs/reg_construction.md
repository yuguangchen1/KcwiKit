# Region File Creation
## A guide for making regions files for the KCWI DRP.
### Last updated: 9/19/2022

For some quick examples see `mizar:/scr/kbss-kcwi/tier1drp/kcwi_jul21/2021jul05/redux`

## Sky Subtraction Regions (`.reg`)

These regions files are constructed using the `*intk.fits` images for a second round of sky subtraction with continuum objects masked.

## Median Filtering Regions

### Continuum Objects (`_icube.thum.reg`)

The recently constructed datacubes are flattened to form a whitelight image and continuum objects are masked for median filtering.


### Line Emission (`_icube_2d.reg`)

The science datacubes are deprojected back into a 2D image and line emission (typically Ly$\alpha$) is marked.
