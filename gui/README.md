# KCWI Post-DRP GUI instruction
## Version: 06/03/2024

## About The GUI

This GUI is designed to refine the sky subtraction, flux calibration, and telluric correction of the KCWI/KCRM data. The users would still have to reduce the data with the [`KCWI DRP`](https://kcwi-drp.readthedocs.io/en/latest/). 
The sky subtraction is performed using the [`ZAP`](https://zap.readthedocs.io/en/latest/) package with a PCA approach. 
The telluric correction makes use of the tellfit function of the [`Pypeit`](https://pypeit.readthedocs.io/en/release/telluric.html) package to derive the best-fit telluric model of the standard star.

## Prerequisites

1. Install and run all the way through the [`KCWI DRP`](https://kcwi-drp.readthedocs.io/en/latest/). You can skip the sky subtraction as the GUI will handle it later.
   
   `reduce -r -f kr*fits -k`
   
   Please note that the released branch of the KCWI DRP only contains the standard star spectra up to 9200A. If your configuration goes beyond that, please replace the files in kcwidrp/data/stds/ with the ones in their develop branch.

2. Install the [`Pypeit`](https://pypeit.readthedocs.io/en/release/telluric.html) package and the modified version of the [`ZAP`](https://github.com/jasonpeng17/zap_for_kcwi) package.
   The major difference between the official ZAP and the modified one is described in the Continuum Filter Widths and Wavelength Segments Section [here](https://github.com/jasonpeng17/zap_for_kcwi/blob/master/doc/index.rst).

3. Download the telluric model grid via `pypeit_install_telluric TelFit_MaunaKea_3100_26100_R20000.fits`. Update the `telgridfile` in the kcwi_viewer.py as
   `telgridfile = your_path_of_telfit_file`. If you can't find the downloaded file, it is usually stored at ~/.pypeit/cache/download/url/*. The file `contents` is TelFit_MaunaKea_3100_26100_R20000.fits. Alternatively, you can directly download the file using the link in `url` in the same directory.  

4. Set the `initial_dir` in the kcwi_viewer.py as your favorite data directory, or set it as None.

5. Install the entire KcwiKit package. To make it available in your Python file directly, please add its path to the PYTHONPATH environmental variable:

   `export PYTHONPATH=${PYTHONPATH}:YOUR_PATH_OF_KCWIKIT/py`

## Usage

### 1. Run the GUI
   
   `python kcwi_viewer.py`

   A GUI window will pop up.
   <img width="1208" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/1802ca72-4767-47ef-a876-2a96ec4aa216">

### 2. Set the input and output directory
Browse and select the input directory to be the `redux` directory where all the KCWI DRP outputs are stored. Also browse and select the output directory (please create a separate working directory to avoid overwriting the DRP output).

### 3. Flux calibration and telluric correction.
  
   Browse and select the *_invsens.fits of the DRP with the `Browse DRP invsens` button.
   <img width="1206" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/177802ee-f669-4470-b844-f1fda3531ff7">
   
   If you are not happy with the default region, you can select the regions interactively. To active this function, please first click the plotting canvas first.
   - Put the mouse and press 'e' on each side of a region to exclude it from the fit
   - Put the mouse and press 'i' on each side of a region to include it from the fit
   - Put the mouse and press 'a' on a single data point to add it into the fit (highly recommend for regions beyond 9000A where most regions are contaminated by telluric)
   - Put the mouse and press 'd' to delete the single data point mistakenly added earlier.
   - Note: Please do not remove the start and end point of the spectrum to avoid running into problems.
   - Press 'f' to re-fit the sensitivity function. You can iteratively choose the regions and re-fit the sensitivity function until you are happy with the model (the cyan line)
     <img width="1193" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/a0125d01-9fe0-4bc0-807f-d857b89b0463">
   - Press 't' to derive the telluric model. This step will take a while and the GUI would become frozen. Be patient until you see the output.
     <img width="1180" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/14e5ab89-2ac4-42b4-9c01-01be296a35b4">
   - Press the `Save updated invsens` button if you are happy with the results.

### 4. Sky subtraction

   a. Enter the science image number. For instance, enter 85 and press the `return ` button if the frame is kr230923_00085. If you want to use the off-field sky frame, please check the `Use Off-field Sky Frame No.`, also enter the sky image frame number and press return. You can also update the science image number by pressing `Previous` and `Next` button.
   
   b. If this is the first time that you reduce this frame, please press `Load Raw Cube`, `Save Cropped Cube` and `Load Cropped Cube`. Otherwise, simply press `Load Cropped Cube` to load the data. This step is simply to crop the datacube outside of the good wavelength region to avoid running into edge problem with ZAP.

   c. If you plan to use the in-field sky, you need to mask out the science target. You can check the `_wlimg.fits` in the output directory,  mask the source using DS9, and save the DS9 region as kr230923_00085.reg for frame kr230923_00085. After that, enter the mask frame number, press enter, and press `Update ZAP mask` button. You can update ZAP mask multiple times whenever you press the button

   d. Now it's time to run ZAP. 
   - If you want to use the default ZAP setup, please uncheck `Use multiple skyseg in ZAP` and `Additional Sky Seg near Halpha`, and press `Run ZAP`. In this case, it adopts one sky segment and fit the overall variation across the entire wavelength. By default, it adopts `cfwidth=300` to remove the science signal. For more information related to the `cfwidth`, please check the github [page](https://github.com/jasonpeng17/zap_for_kcwi/blob/master/doc/index.rst) of the modified version and the original ZAP paper by [Soto et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.3210S/abstract). To proceed, please type 'q' in the text widget like the example below and press enter to start ZAP.
     <img width="1160" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/e320095f-625b-48dc-8bc8-301f4f439044">
     
     Alternatively, you can update the cfwidth by changing the text directly (the last value in the tuple), and type 'u' to update. The screenshot below shows an example to update the cfwidth to 100.
     <img width="1160" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/48e7def7-e149-4c6f-8e71-4ae3c8ea8db4">
     
     The text widget will then print out the updated cfwidth. If you are happy with it, type 'q' in the text wdiget to proceed. If you want to reset the default setting, you can type 'r' instead.
   - We find that using multiple sky segment sometimes can reduce the sky residual. This approach was adopted in the original ZAP by [Soto et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.3210S/abstract) when it first came out, but the updated version switched to one sky segment. For how the sky segment is chosen, please refer to [Soto et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.3210S/abstract) for more details. If you want to use this option, please check `Use multiple skyseg in ZAP`. Besides, we find that if a very strong emission line happens to fall in a sky segment with relatively sparse sky lines (e.g., Halpha emission), a large `cfwidth` would aritificially introduce more noise in the sky-subtracted spectrum. By checking `Additionally Sky Seg near Halpha` and setting the redshift, the code would generate an additional segment in the region of Halpha with `cfwidth=5`.
     <img width="1149" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/fd3cb55d-6935-47e2-88e0-267ee55428ed">

     If you are happy with this segment, type 'q' in the text widget to proceed. Otherwise, you can update the segment in the text widget directly to whatever regions you want to use, and type 'u' to update the change. Every tuple indicates (the start of the segment, the end of the segment, cfwdith). You can either delete a tuple or add a new one, or simply changing some of the values in it. The screenshot below shows an example when you want to combine segment (7200, 7700, 300) and (7700, 8265, 300) into a new segment with cfwidth = 50.
     <img width="1157" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/51f3fc57-47d7-425c-8dd0-768468c5718b">

     Now the text widget would print out the updated sky segment. Type 'q' if you wanna proceed, or 'r' if you want to reset the sky segment.
     <img width="1149" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/116ee095-48f3-46a1-8c3b-237b24494e8e">

   - If you uncheck or check the box mistakenly, you can simply update the check box and press `Run ZAP` button again, which will reactivate the text widget.
   - It will take a while for ZAP to run, please be patient. The output will be in the terminal window.

### 5. Examine the output
   Now you can examine the sky-subtracted spectrum to see if that looks good to you. By default, it will extract the spectrum from a 10x10 box centered on the center of the FoV (i.e., a box with lower left at (x1, y1) = (13,44) and upper right at (x2, y2) = (23, 54). You can update the region used to extract the spectrum simply by updating the DS9 pixel indices in the box below the plotting canvas as `x1, y1, x2, y2` (lower left + upper right). 
   <img width="1178" alt="image" src="https://github.com/zhuyunz/KcwiKit/assets/33030986/527c7005-6f49-4b5e-be97-4bfa37c6728c">

   The pre-ZAP spec is scaled to have the same median as the sky-subtracted spectrum for better display purpose.

   If you set the redshift, it will also indicate where the emission lines are expected to be. You can use the navigation bar to examine the plot, and go back to the full plot with right click. Besides, you can also examine the plot with following options (click the canvas first to activate the keyboard interaction):
   
   - Press 'n' to turn off the pre-ZAP spec.
   - Press 's' to turn on the pre-ZAP spec so that you would know where the sky lines are.
   - Press 'o' to move the spectrum to the observed frame. It is useful when you want to update the sky segment (which requires observed frame) and re-run ZAP.
   - Press 'r' to move the spectrum to the rest frame.

   Go to the next frame if you are happy with the spectrum; otherwise, press `Run ZAP` button and run the sky subtraction again with the updated sky segment. It's always a good idea to check the output `_zap_icubes.fits`! For the frame using the in-field sky, you can also check the white-lighted image of the clean datacube `_zapclean_wlimg.fits` to make sure the source is properly masked. You can update the DS9 region mask and re-run the ZAP again!


## Contact

If you encounter any problems or find bugs, please email [Zhuyun Zhuang](mailto:zzhuang@astro.caltech.edu).
