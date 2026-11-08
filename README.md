# Keck Cosmic Web Imager: Post-Processing and Improvements

<a href="https://ascl.net/2404.003"><img src="https://img.shields.io/badge/ascl-2404.003-blue.svg?colorB=262255" alt="ascl:2404.003" /></a>

Welcome to `KCWIKit`!

`KCWIKit` is a community developed package to extend the official KCWI DRP with a variety of stacking tools and DRP improvements. The software offers masking and median filtering scripts to be used while running the KCWI DRP, and a step-by-step KCWI_DRP implementation for finer control over the reduction process. Once the DRP has finished, `KCWIKit` can be used to stack the output cubes via the [`Montage`](http://montage.ipac.caltech.edu/) package. Various functions cross-correlate and mosaic the constituent cubes and the final stacked cubes are WCS corrected. Helper functions can then be used to deproject the stacked cube into lower-dimensional representations should the user desire.

This repo is organized as follows (all under `kcwikit/`):
- `docs/` Documentation and instructions.
- `kcwi/` Post-DRP stacking (`kcwi.py`), combination of stacks (`kcwi_combine.py`), region-file handling (`kcwi_regions.py`) and helpers (`kcwi_tools.py`).
- `scripts/` Command line tools for quick operations on KCWI data (masking, cosmic rays, alignment, sky-level matching, ...). See [Command Line Tools](./kcwikit/docs/scripts_instruction.md).
- `pyDRP/` Modified primitives for the [Python pipeline](https://kcwi-drp.readthedocs.io/en/latest/); the maintained fork lives in [KCWI_DRP](https://github.com/yuanzeding/KCWI_DRP).
- `arch/` Archival programs for the [IDL pipeline](https://github.com/Keck-DataReductionPipelines/KcwiDRP) (no longer supported).

What changed recently is listed in the [CHANGELOG](./CHANGELOG.md).

We maintain a separate tool [`KCWI Sky Wizard`](https://github.com/zhuyunz/KSkyWizard) for the telluric correction and sky subtraction. 

## Installation

KCWIKit is recommended to be installed with `python=~3.12`. 

1. Activate the KCWI DRP environment.
    ```bash
    conda activate kcwidrp
    ```
    
    Or create a new Python environment. 

    ```bash
    conda create -n kcwikit python=3.12
    ```

2. Install `MontagePy` or the `Montage` command line tools. The two forms of the `Montage` package perform the same functions. You can choose either to install.

    `MontagePy` is generally easier to install if you have the correct `Python` environment:
    ```bash
    pip install MontagePy
    ```
    Note: the `MontagePy` 2.2.0 wheel on Apple silicon reads values with |value| < 1 as 0 ([Montage issue #80](https://github.com/Caltech-IPAC/Montage/issues/80)). `kcwi_stack` works around this automatically when it uses `MontagePy` (`ISS6FIX` header keyword in the trimmed cubes); the command-line `mProjectCube` is not affected.

    <details>
    <summary> Alternative: Install the command line version </summary>

    If the above installation fails, alternatively, you can choose to install the command line tools by compiling from the [source code](http://montage.ipac.caltech.edu/docs/download2.html).
    </details>

3. Git clone this repository and run the setup tools. 

    ```bash
    git clone https://github.com/yuguangchen1/KcwiKit.git 
    cd KCWIKit
    pip install .
    ```
    DS9 region files are read with the astropy [`regions`](https://astropy-regions.readthedocs.io/) package (installed by `pip`); `pyregion` is used as a fallback if it is installed instead.

4. <details>
    <summary>(optional) Install KCWI Sky Wizard. </summary>

    `KSkyWizard` is a standalone tool to perform telluric correction and advanced sky subtraction based on PCA models for the red channel (KCRM). It is recommended to be run on `Python >= 3.11` environments. We refer to the [KSkyWizard documentation](https://github.com/zhuyunz/KSkyWizard) for the installation and instruction. 
    </details>

## Instructions

1. Installing and running the improved DRP. 

    1.1. [Installing the modified DRP](./kcwikit/docs/install_DRP.md)

    1.2. [Running the DRP for the blue channel](./kcwikit/docs/Running_blue_DRP.md)

    1.3. [Running the DRP for the red channel](./kcwikit/docs/Running_red_DRP.md)

2. [Running the post-DRP stacking code.](./kcwikit/docs/KCWI_post-DRP_stacking.md)

3. [Command Line Tools.](./kcwikit/docs/scripts_instruction.md)

4. Newer tools (2025-2026), all documented in the [Command Line Tools](./kcwikit/docs/scripts_instruction.md) page and in each script's `--help`:

    - `kcwi_automask` - automatic sky, median-filter, ZAP and cosmic-ray protection masks instead of drawing regions in DS9.
    - `kcwi_xalign` / `kcwi_shiftcheck` - blue-red frame registration without a point source, and a check of the resulting shift lists.
    - `kcwi_skylevel`, `kcwi_dichroic_step`, `kcwi_cube_windows` - measure and correct residual sky levels across the dichroic, anchored to HST / Legacy Survey / Pan-STARRS images or another IFU cube.
    - `kcwi_stack(crr_split=True)` and `kcwi_combine_stacks` (in `kcwikit.kcwi`) - cosmic-ray rejection on the stacking grid and co-addition of several stacks.
    - `kcwi_masksky_ds9_bcrr`, `kcwi_mask_cube`, `kcwi_gen_skyfile`, `kcrm_group_frames --blue` - red/blue channel helpers for the modified DRP.

## Examples
[Here](./kcwikit/examples/) are some examples for the stacking code which are referenced in the instructions. 

Some KCWI-red commissioning data are shared [here](https://drive.google.com/drive/folders/1j94cx0IZTZLTpLdm6QRKcpDQm-RvMWOJ) if you would like to try running the code. 

## KCWIKit at KSM24

KCWIKit is hosting a breakout session at the [2024 Keck Science Meeting](https://kecksciencemeeting.org/breakout-sessions/)! 


## Citing KCWIKit

If you use `KCWIKit`, please cite the [ASCL entry](https://ascl.net/2404.003) (BibTeX format on [ADS](https://ui.adsabs.harvard.edu/abs/2024ascl.soft04003P)) and [Chen et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508...19C) where the code is implemented.

Developed and maintained by [Nikolaus Prusinski](mailto:nik@astro.caltech.edu), [Yuguang Chen](mailto:yugchen@ucdavis.edu), [Zhuyun Zhuang](mailto:zzhuang@astro.caltech.edu) and [Yuanze Ding](mailto:yding@caltech.edu) (2025-2026 additions: automatic masking, frame alignment, sky-level matching, split-frame cosmic-ray rejection). If you encounter bugs, have questions on how to run the code, or have feature requests, drop us a line! 


