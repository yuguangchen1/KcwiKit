# KCWI Post-DRP Stacking Instructions
## Initial Version: 12/29/2023

# Prerequisites

Incomplete list: `Montage` (or `MontagePy`), `kcwi`, `astropy`, `fpdf`, `tqdm`
# Inputs

1. `KCWI_DRP` science outputs all centered on a particular target
2. An `<object_name>.list` file specifying the input cubes to be used to make the stack
3. An `<object_name>.par` file with a list of arguments and parameters for the particular stack
4. A `<object_name>.ipynb` `jupyter` notebook in which to run the stacking routine. This could be substituted for a `python` script in one so desires.

The input files nominally come in two flavors:
- `kb*_icubes.fits` and `kr*_icubes.fits` if flux calibrated **OR**
- `kb*_icubed.fits` and `kr*_icubed.fits` if making a quick stack while observing
but in principle, one could use cubes with any naming convention (e.g. KOA) as long as the `*_icube?.fits` tail is retained. Note that the `cubed = True` argument will be required in the stacking function calls if the cubes are not flux calibrated.

## `.list`  File

A sample `.list` file is provided at `kcwi/py/examples/q2343-BX610.list` and the first few lines are displayed here for convenience.

```
/scr/yuguangchen/obs/kcwi/kcwi_nov19/2019nov29/redux/kb191129_00041 3 3
/scr/yuguangchen/obs/kcwi/kcwi_nov19/2019nov29/redux/kb191129_00040 3 3
/scr/yuguangchen/obs/kcwi/kcwi_nov19/2019nov29/redux/kb191129_00042 3 3
/scr/yuguangchen/obs/kcwi/kcwi_nov19/2019nov29/redux/kb191129_00043 3 3
/scr/yuguangchen/obs/kcwi/kcwi_nov19/2019nov29/redux/kb191129_00044 3 3
```

Note that the full path of each file is required and no cube type tail (e.g. `icubes.fits`) is provided here. The second and third columns list the number of pixels to be trimmed from the bottom and top of each input cube, thus removing some bad edge pixels. 

In other words, the line `kb230105_00056 2 5` would remove two pixels from the bottom of the cube and five pixels from the top (for all wavelengths). These pixels correspond to slice edge rolloff in the 2D (detector) format. This clipping only applies to the "data-region" of the cube and does not trim the DAR padding.

We have found three pixel top and bottom clipping to work well for the Medium slicer, low resolution grating, 2x2 binning setup but your mileage may vary. You can check if the trimming needs to be adjust by opening the `kcwi_align/*thum0.fits` file and inspecting each file/extension.

Here's a simple command that can be used to produce sequential lists of cubes in the correct format: `printf "%s\n" /path/to/dir/2021aug13/redux/kb210811_000{26..34}\ 3\ 3`

```
/path/to/dir/2021aug13/redux/kb210811_00026 3 3
/path/to/dir/2021aug13/redux/kb210811_00027 3 3
.
.
.
```

## `.par`  File

The parameter file takes in a list of options used while stacking and to define the output stack (e.g. final stack dimensions, pixel scales, reference astrometry images, etc.).

A sample `.par` file is again provided at `kcwi/py/examples/q2343-BX610.par` and displayed here.

```
dimension 100 100
align_box 45 37 65 65 
orientation 0
xpix 0.3
ypix 0.3
ref_xy 51 52
ref_ad 356.53929 12.822
ref_fn /scr/yuguangchen/obs/kcwi/plan/q2343/q2343Rs.fits
```

For a basic stack without a final WCS adjustment (skipping `kcwi_astrometry`), one can omit the `ref_*` lines.

Each of these parameters correspond to arguments in the `kcwi_align`, `kcwi_stack`, or `kcwi_astrometry` functions and additional parameters can be added on new lines. We encourage the reader to check out each of the relevant docstrings (LINK!) for more details and implementation strategies. We outline the basic parameters below:

**General Parameters:**
`dimension x y`: final cube x and y dimensions
`orientation deg`: final cube position angle (degrees east of north; 0 corresponds to north up, east left)
`xpix` and `ypix` are the pixel scales in units of arcsec/pix

**`kcwi_align` Parameters**:
`align_box x_lower_left y_lower_left x_upper_right y_upper_right`: sets up an alignment box in the `kcwi_align` stage used to cross-correlate each frame and align on a point source common to all frames. Details on constructing this box are discussed in the `kcwi_align` section.

**`kcwi_astrometry` Parameters**:
`ref_xy x y`: the coordinates of the reference object in the stacked cube
`ref_ad RA Dec`: the physical coordinates (in degrees) of the reference object in the reference image (`ref_fn`). 
`ref_fn path`: the full path to the reference (finder) image used to correct the cube's WCS.

## The `jupyter` notebook

The `jupyter` notebook is used to display the stacking (text) outputs and several diagnostic plots for each stage. A sample notebook is provided at `kcwi/py/examples/q2343-BX610.ipynb`. Note that the `list`, `par`, and notebook should all have the same file name (abbreviated as `fn`). 

The notebook is quite simple and only has a few cells. In general, all one needs is the following (split up into individual cells) to produce a WCS corrected mosaic:

```
import kcwi
fn = '<object_name>.list'
kcwi.kcwi_align(fn, noalign=True)
kcwi.kcwi_align(fn)
kcwi.kcwi_stack(fn)
kcwi.kcwi_astrometry(fn)
```

We discuss each section individually below.

# `kcwi_align`

The `kcwi_align` stage is typically run twice: the first time with `noalign = True`, and the second with `noalign = False` (the default). The `noalign` argument skips the cross-correlation of frames. 

To align each of the constituent cubes, we start by running `kcwi.kcwi_align(fn, noalign=True)`. This makes a trivial `<object_name>.shift.list` file, a `kcwi_align` subdirectory in the working directory, and writes a `object_name.thum0.fits` diagnostic cube into that subdirectory. 

The `thum0` file is a multi-extension cube with a whitelight (wavelength averaged) image of each constituent cube on each extension. The wavelength range used to make each whitelight image can be adjusted using the `wavebin` argument. At this stage, one can see if the trimming has done its job or if further edge pixel mitigation might be required, and one can (hopefully!) see their target in each one of the images. If the trimming needs to be adjusted, simply adjust the values in the `.list` file and rerun `kcwi_align`.

If there is a nearby point source that can be picked up in all of the cubes, that is the target around which to draw the alignment box; otherwise, one could use the cube target (assumed here to be an extended source) as a backup at the cost of some alignment accuracy.

The `thum0` file places each frame on the same final grid and the goal is to find a box large enough that the alignment source (or at least the brightest part of it) is within the box in all images. Once such a box is constructed, the coordinates should be fed into the `align_box` argument of the `.par` file. Then one can run `kcwi.kcwi_align(fn, noalign=False)`, or more simply `kcwi.kcwi_align(fn)`.

For each cube, a diagnostic plot should be displayed roughly centered on the alignment source, plus a smaller subgrid aligmnet box and cross, all centered on the brightest part of the target. If this is not the case, the alignment has failed for that object. Some causes and possible remedies are listed below:
- The alignment box is misaligned or not well positioned $\rightarrow$ reposition the box, expanding/contracting as required
- If the object is at the edge of the search box for a few of the frames, one can adjust the `search_size` parameter from 10 (default) to 20 pixels (i.e. `kcwi_align(fn, search_size=20)`)

## `.pre.list` file

If the WCS information for a few frames is significantly inaccurate or the frames are considerably offset from the group, one can also construct a `.pre.list` file. An example is shown below:

```
kb210415_00060 -1.8 0
kb210415_00061 -1.8 0
kb210705_00030 -1.5 -1.5
kb210705_00031 -4.0 0.0
kb210705_00032 -4.0 -3.0
kb210706_00049 1 -1
```

Similar in format to the `.list` file, in this case only the affected frames are listed and the latter two columns show the x and y shift _in arcseconds_ respectively. For example, if frame 60 needs to be moved 10 pixels up with respect to the first (fiducial) frame, and the pixel scale (`ypix`) was 0.3"/pixel, the appropriate line item would be `kb210415_00060 0 3`. A similar calculation would be required for any horizontal translation.

Finally, one can look at the output `kcwi_align/<object_name>.thum.fits` multi-extension cube (note `thum` not `thum0`) to see that all cubes have been aligned successfully.

# Optional: `kcwi_check_flux` and `kcwi_norm_flux`

If a few of your exposures were through clouds or underwent other throughput variations, you can normalize the overall flux to the first (fiducial) frame. 

Running `kcwi.kcwi_check_flux(fn)` will output a plot of each frame's flux relative to the first. The user can take note of which frames appear significantly different from one another and input that list into `kcwi_norm_flux`.

For instance, if frames 2, 5, and 7 had a significant loss in flux due to cloud cover, the user can run `kcwi.kcwi_norm_flux(fn, frame=[2,5,7])` and this function will generate a `<object_name>.flx.list` file that will show the necessary multiplicative flux correction required for each frame. The `kcwi_stack` function will read this file in by default and make the appropriate corrections during stacking.

# `kcwi_stack`

The `kcwi_stack` function converts the aligned cubes into a format readable by `Montage` and stacks them. If alignment was successful, this stage should generally run on rails and a `tqdm` progress bar provides some eye candy. This step generally takes several minutes on a decent system, obviously dependent on the number of exposures.

Some pertinent flags for this stage are enumerated below.

`cubed`: set to `True` to use `*icubed.fits` files instead of `*icubes.fits` (default)
`overwrite`: overwrite any previous files
`montagepy`: set to `True` if your system has `MontagePy` installed instead of `Montage`
`multiple_grangles`: set to `True` if you are combining exposures from the same grating with different central wavelengths. Currently under development - talk to Nik!
`npix_trim`: the number of pixels to trim from the side of a cube. Default is 3, may require more or less depending on slicer/binning.

The **output** of this stage should be four cubes:
- `<object_name>_icubes.fits` - stacked cube
- `<object_name>_vcubes.fits` - stacked variance cube
- `<object_name>_mcubes.fits` - stacked `MASK` and `FLAGS` cube (concatenated)
- `<object_name>_ecubes.fits` - stacked exposure time cube; masked pixels not included in the total exposure time

# Optional: `kcwi_astrometry`

The `kcwi_astrometry` function updates the final stacked cube's WCS information based on an external reference image. For extragalactic targets, NED optical band finder charts tend to work well for this sort of thing. 

To find the center of your object in the cube, one can use `QFitsView` and the `g` key for a rough estimate. Those coordinates get fed into `ref_xy`, and the physical coordinates (in degrees) of the same target in the reference image (`ref_fn`) become `ref_ad` in the parameter file. 

With `ref_xy`, `ref_ad`, and `ref_fn` set we can proceed with running `kcwi_astrometry(fn)`. Two rounds of cross-correlations are run and the output diagnostic plot should appear similar to those of `kcwi_align`. 

Note: if you'd prefer to just change the WCS information without doing a full cross-correlation (this is sometimes useful for QSOs), you can add `ref_nocrl 1` to the parameter file. 

The **output** WCS-corrected cube will appear in the working directory as `<object_name>_icubes_wcs.fits`.