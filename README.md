# Keck Cosmic Web Imager: Post-Processing and Improvements

<a href="https://ascl.net/2404.003"><img src="https://img.shields.io/badge/ascl-2404.003-blue.svg?colorB=262255" alt="ascl:2404.003" /></a>

Welcome to `KCWIKit`!

`KCWIKit` is a community developed package to extend the official KCWI DRP with a variety of stacking tools and DRP improvements. The software offers masking and median filtering scripts to be used while running the KCWI DRP, and a step-by-step KCWI_DRP implementation for finer control over the reduction process. Once the DRP has finished, `KCWIKit` can be used to stack the output cubes via the [`Montage`](http://montage.ipac.caltech.edu/) package. Various functions cross-correlate and mosaic the constituent cubes and the final stacked cubes are WCS corrected. Helper functions can then be used to deproject the stacked cube into lower-dimensional representations should the user desire.

This repo is organized as follows:
- docs/ Documentation and Instructions
- kcwi/ Post-DRP stacking.
- scripts/ Covenient command line tools to conduct quick operations for the KCWI data. 
- DRP/ Improving the [Python pipeline](https://kcwi-drp.readthedocs.io/en/latest/).
- archs/ archival programs for the [IDL pipeline](https://github.com/Keck-DataReductionPipelines/KcwiDRP).(No longer supported.)

Check the subdirectories for additional instructions and prerequisites. 

## Installation

KCWIKit is recommended to be installed with `python <= 3.9`. 

1. Activate the KCWI DRP environment.
    
    Or create a new Python environment. 

    ```bash
    conda create -n kskywizard python=3.8
    ```

2. Install `MontagePy` or the `Montage` command line tools. The two forms of the `Montage` package performs the same functions. You can choose either to install.

    `MontagePy` is generally easier to install if you have the correct `Python` environment:
    ```bash
    pip install MontagePy
    ```

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

4. <details>
    <summary>(optional) Install KCWI Sky Wizard. </summary>

    `KSkyWizard` is a standalone tool to perform telluric correction and advanced sky subtraction based on PCA models for the red channel (KCRM). It is recommended to be run on `Python >= 3.11` environments. We refer to the [KSkyWizard documentation](https://github.com/zhuyunz/KSkyWizard) for the installation and instruction. 
    </details>

## Instuctions

1. Installing and running the improved DRP. 

    1.1. [Installing the modified DRP](./kcwikit/docs/install_DRP.md)

    1.2. Running the DRP for the blue channel

    1.3. Running the DRP for the red channel (under construction)

2. [Running the post-DRP stacking code.](./kcwikit/docs/KCWI_post-DRP_stacking.md)

3. [Command Line Tools.](./kcwikit/docs/scripts_instruction.md)

## Examples

## KCWIKit at KSM24

KCWIKit is hosting a breakout session at the [2024 Keck Science Meeting](https://kecksciencemeeting.org/breakout-sessions/)! 


## Citing KCWIKit

If you use `KCWIKit`, please cite the [ASCL entry](https://ascl.net/2404.003) (BibTeX format on [ADS](https://ui.adsabs.harvard.edu/abs/2024ascl.soft04003P)) and [Chen et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508...19C) where the code is implemented.

Developed and maintained by [Nikolaus Prusinski](mailto:nik@astro.caltech.edu), [Yuguang Chen](mailto:yugchen@ucdavis.edu) and [Zhuyun Zhuang](mailto:zzhuang@astro.caltech.edu). If you encounter bugs, have questions on how to run the code, or have feature requests, drop us a line! 


