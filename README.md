# KCWI Post-Processing and Improvements

`KCWIKit` extends the official KCWI DRP with a variety of stacking tools and DRP improvements. The software offers masking and median filtering scripts to be used while running the KCWI DRP, and a step-by-step KCWI_DRP implementation for finer control over the reduction process. Once the DRP has finished, `KCWIKit` can be used to stack the output cubes via the `Montage` package. Various functions cross-correlate and mosaic the constituent cubes and the final stacked cubes are WCS corrected. Helper functions can then be used to deproject the stacked cube into lower-dimensional representations should the user desire.

This repo is organized as follows:
- docs/ Documentation and Instructions
- pro/ Improvinng the [IDL pipeline](https://github.com/Keck-DataReductionPipelines/KcwiDRP).
- py/ Post-DRP alignment, stacking.
- pyDRP/ Improving the [Python pipeline](https://kcwi-drp.readthedocs.io/en/latest/).
- gui/ An interactive GUI for the post-DRP reduction, including performing/refining sky subtraction, flux calibration, and telluric correction. 

Check the subdirectories for additional instructions and prerequisites. 

## Citing KCWIKit

If you use `KCWIKit`, please cite [Chen et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508...19C).

Developed and maintained by [Nikolaus Prusinski](mailto:nik@astro.caltech.edu), [Yuguang Chen](mailto:yugchen@ucdavis.edu) and [Zhuyun Zhuang](mailto:zzhuang@astro.caltech.edu). If you encounter bugs, have questions on how to run the code, or have feature requests, drop us a line! 


