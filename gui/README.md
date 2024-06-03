# KCRM Post-DRP GUI instruction
## Version: 06/03/2024

## About The GUI

This GUI is designed to refine the sky subtraction, flux calibration, and telluric correction of the KCRM (KCWI-red) data. The users would still have to reduce the data with the [`KCWI DRP`](https://kcwi-drp.readthedocs.io/en/latest/). 
The sky subtraction is performed using the [`ZAP`](https://zap.readthedocs.io/en/latest/) package with a PCA approach. 
The telluric correction makes use of the tellfit function of the [`Pypeit`](https://pypeit.readthedocs.io/en/release/telluric.html) package to derive the best-fit telluric model of the standard star.

## Prerequisites

1. Install and run all the way through the [`KCWI DRP`](https://kcwi-drp.readthedocs.io/en/latest/). You can skip the sky subtraction as the GUI will handle it later.
   
   `reduce -r -f kr*fits -k`
   
   Please note that the released branch of the KCWI DRP only contains the standard star spectra up to 9200A. If your configuration goes beyond that, please replace the files in kcwidrp/data/stds/ with the ones in their develop branch.

2. Install the [`Pypeit`](https://pypeit.readthedocs.io/en/release/telluric.html) package and the modified version of the [`ZAP`](https://github.com/jasonpeng17/zap_for_kcwi) package.
   The major difference between the official ZAP and the modified one is described in the Continuum Filter Widths and Wavelength Segments Section [here](https://github.com/jasonpeng17/zap_for_kcwi/blob/master/doc/index.rst).

3. Download the telluric model grid via `pypeit_install_telluric TelFit_MaunaKea_3100_26100_R20000.fits`. Update the `telgridfile` in the kcrm_viewer.py as
   `telgridfile = your_path_of_telfit_file`

4. Set the `initial_dir` in the kcrm_viewer.py as your favorite data directory, or set it as None.

## Usage
   
 
   
