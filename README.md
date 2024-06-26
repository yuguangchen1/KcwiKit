# KCWI Post-Processing and Improvements

<a href="https://ascl.net/2404.003"><img src="https://img.shields.io/badge/ascl-2404.003-blue.svg?colorB=262255" alt="ascl:2404.003" /></a>

`KCWIKit` extends the official KCWI DRP with a variety of stacking tools and DRP improvements. The software offers masking and median filtering scripts to be used while running the KCWI DRP, and a step-by-step KCWI_DRP implementation for finer control over the reduction process. Once the DRP has finished, `KCWIKit` can be used to stack the output cubes via the `Montage` package. Various functions cross-correlate and mosaic the constituent cubes and the final stacked cubes are WCS corrected. Helper functions can then be used to deproject the stacked cube into lower-dimensional representations should the user desire.

This repo is organized as follows:
- docs/ Documentation and Instructions
- pro/ Improvinng the [IDL pipeline](https://github.com/Keck-DataReductionPipelines/KcwiDRP).
- py/ Post-DRP alignment, stacking.
- pyDRP/ Improving the [Python pipeline](https://kcwi-drp.readthedocs.io/en/latest/). 

Check the subdirectories for additional instructions and prerequisites. 

## Citing KCWIKit

If you use `KCWIKit`, please cite the [ASCL entry](https://ascl.net/2404.003) (BibTeX format on [ADS](https://ui.adsabs.harvard.edu/abs/2024ascl.soft04003P)) and [Chen et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508...19C) where the code is implemented.

Developed and maintained by [Nikolaus Prusinski](mailto:nik@astro.caltech.edu) and [Yuguang Chen](mailto:yugchen@ucdavis.edu). If you encounter bugs, have questions on how to run the code, or have feature requests, drop us a line! 


