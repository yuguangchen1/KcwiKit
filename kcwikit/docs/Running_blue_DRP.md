# Running the DRP for the blue channel

Here summarizes the typical routine with optional improvements to reduce the blue-channel KCWI data using the [updated DRP](../docs/install_DRP.md).

**Before you start**: Running the default DRP does not require to follow this instruction. Simply run ```reduce_kcwi -b -f kb*.fits -g```. 

1. Copy the [configuration file `kcwi.cfg`](../pyDRP/configs/kcwi.cfg) to the working directory. Running the DRP until after the sky subtraction.

    ```bash
    reduce_kcwi -b -f kb*.fits -g -c kcwi.cfg -st 1
    ```

2. (optional) Avoid oversubtraction in sky model. 

    Create sky masks to remove astronomical sources from the sky model. See [here](../docs/reg_construction.md) for details to create `region` masks and covert them to binary FITS files. 

    Remove all the entries (`OBJECT` and `SKY`) in the `kcwib.proc` table for the science observations. For example, for ```kb240211_00089``` frame, remove the following lines. 

    ```
    |      89 | 65c6d022d14508b9bd6870ea | 2201009 |   OBJECT |  2024-02-11-89 |   15.0 | BLUE | Medium |   BL |   28.52 |  4499.90 | 2,2 | KBlue |    60351.382006 |     5 |  intk | kb240211_00089.fits |      sdss08579 | kb240211_00089.fits |
    |      89 | 65c6d022d14508b9bd6870ea | 2201009 |      SKY |  2024-02-11-98 |   35.0 | BLUE | Medium |   BL |   28.52 |  4499.90 | 2,2 | KBlue |    60351.382006 |     5 |  sky | kb240211_00089.fits |      sdss08579 | kb240211_00089.fits |
    ```

    Remove the `*_sky.fits` in the `redux/` directory for the frames of which the sky model needs update.

    Rerun the DRP until after the cube creation
    ```bash
    reduce_kcwi -b -f kb*.fits -g -c kcwi.cfg -st 2
    ```
3. (optional) Check for possible CR rejection false positives
   It is very common for the current DRP to mistakenly consider bright continuum sources or strong emission lines as cosmic rays, if the count rate surpasses $\simeq$2000. Fortunately, it is easy to identify the fake positives from icube.fits.

   Open the icube file and inspect the FLAGS extension. If there are values greater than 3 clustering around the QSO, then it is most likely CR false positives.

   To fix this problem, open the corresponding intf file, create region masks covering the continuum/strong lines, and use kcwi_masksky_ds9_bcrr to convert them into FITS masks:
   ```bash
   kcwi_masksky_ds9_bcrr kb250421_00030_intf.fits kb250421_00030.reg
   ```
   You also have to set bcrmsk=True in the .cfg file
4. (optional) Running the median filtering.

    Create white-light and 2D-flattened images for `icube.fits` files. 

    ```bash
    kcwi_collapse kb*_icube.fits
    kcwi_flatten_cube kb*_icube.fits
    ```

    Create region masks on the 2D-flattened images and white-light images. See [here](../docs/reg_construction.md) for instructions. 

    Run median filtering. 
    ```bash
    kcwi_medfilter kb*_icube.fits
    ```
    See [kcwi_medfilter](../docs/scripts_instruction.md) instructions for details to run and configure the parameters. 

5. Once again remove all the science entries in `kcwib.proc`. Run the DRP to the end.

    ```bash
    reduce_kcwi -b -f kb*.fits -g -c kcwi.cfg -st 3
    ```

6. (optional) Improving the flux-calibration accuracy by combining multiple standard star frames. 

    Create a `*.list` file of the full path of the `*_invsens.fits` files to be combined. Run the `kcwi_combinestd` command. 

    ```bash
    kcwi_combinestd <name>.list
    ``` 
    This will create a new `<name>_invsens.fits` that stores the combined inverse sensitivity curves. It will also create a `bokeh` figure comparing the individual curves. If some frames are outliers, most likely because of clouds, remove or comment out (using `#`) that frame in the `*.list` file, and rerun the command again. 

    Replace the flux calibration entries in the `kcwib.proc` file. For example, replace this line, 
    ```
    |      80 | 5b906d4dce9e2434e0316846 | 2201009 |  INVSENS |  NONE |   10.0 | BLUE | Medium |   BL |   27.75 |  4499.87 | 2,2 | KBlue |    58433.645553 |     8 | invsens |  kb181111_00080.fits |               feige34 |   kb181111_00080.fits |
    ```
    with this line
    ```
    |      80 | 5b906d4dce9e2434e0316846 | 2201009 |  INVSENS |  NONE |   10.0 | BLUE | Medium |   BL |   27.75 |  4499.87 | 2,2 | KBlue |    58433.645553 |     8 | invsens |         med_bl.fits |               feige34 |         med_bl.fits |
    ```
    Note that one either need to remove or replace all the other `invsens` entries, because the DRP automatically use the nearest `invsens` in UT to create flux-calibrated science frames. 

    Once again, remove entries of the science frames in the `kcwib.proc` table to ensure the DRP can update the `icubes.fits` files. 

    Rerun the last stage of the DRP. 
    ```bash
    reduce_kcwi -b -f kb*.fits -g -c kcwi.cfg -st 3
    ```

There is an example bash code to run all the steps above in batch. See [here](./master_reduce_instructions.md) for details. 

Now, you can continue to the [instruction on performing post-DRP stacking](./KCWI_post-DRP_stacking.md) to combine your frames. 

