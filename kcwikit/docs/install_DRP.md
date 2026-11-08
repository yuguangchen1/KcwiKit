# Installing the modified DRP

**Note**: A new DRP version that supports `Python >= 3.11` is currently being prepared by the observatory. Some of the improvements on this page will be included in the new DRP. 

**Disclaimer**: Some features introduced in this modified version has not been sufficiently tested (which is why it is not yet merged into the official version). We are happy to receive bug reports, but we cannot guarantee that this version would work for your observations. So use with caution!

There are two methods to install the modified DRP. 

1. If you have installed the official DRP with `Python ~= 3.8`, you can download the modified files from the [pyDRP](../pyDRP/) directory in this repository and replace the older ones with them. For example, copy all files in the [`pyDRP/primitives`](../pyDRP/primitives/) to the `KCWI_DRP/kcwidrp/primitives` directory where your local version of the DRP is installed.

    Then, rerun the installation script (`python setup.py install`). 

2. Alternatively, we keep a separate version of DRP at `https://github.com/yuguangchen1/KCWI_DRP/tree/KCWIKit`. Download this updated version.

    ```bash
    git clone https://github.com/yuguangchen1/KCWI_DRP.git
    git checkout KCWIKit
    ```

    And install this updated version. 

    ```bash
    python setup.py install
    ```


## Summary of the improvements

The fork adds, relative to the official 1.1.1 release (details and dates in the fork's `changelog.txt`):

- Sine clock-pattern subtraction (`subtract_sine`), a flat-field priority list (`flat_order`), configurable interpolation orders (`warp_order`, `DAR_shift_order`, `wave_interp_order`) and DAR reference wavelength (`DAR_wref`), all set in `kcwi.cfg`.
- `DAR_apply = False` to skip the DAR shift of science frames (apply it when drizzling in KCWIKit instead); standard stars are always corrected.
- Cosmic rays: red-channel masks from a comparison of frames (`kcrm_create_crmsk`, `crmsk = True`), an optional extra astroscrappy pass (`extra_crr`), and a protection mask (`bcrmsk = True`, `<frame>_intf_bcrmsk.fits` from `kcwi_automask --bcrmsk` or `kcwi_masksky_ds9_bcrr`) that keeps astroscrappy from cleaning bright, variable sources such as a QSO core.
- Camera-specific saturation levels; MASK/FLAGS warped with constant fill outside the field.
- Partial runs: `reduce_kcwi -st ff` (stages 1-4, until `_intf`) and `-st cr` (stage 1 until cosmic-ray rejection).

The `pyDRP/` directory in this repository is an older snapshot of these primitives and has drifted from the fork; prefer method 2 above.
