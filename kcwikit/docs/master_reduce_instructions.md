# KCWI DRP Instructions
## A guide to [`master_reduce`](../examples/master_reduce) aka the complete and thorough KCWI reduction
### Last updated: 9/19/2022

Quick reference: `master_reduce -h` or `head ~/Software/kcwi/pyDRP/master_reduce`

Requirements: Forked [KCWI\_DRP](https://github.com/prusinski/KCWI_DRP) and [kcwi](https://github.com/yuguangchen1/kcwi) repositories
[KCWI\_DRP installation instructions](https://kcwi-drp.readthedocs.io/en/latest/installing.html)
Create the `kcwidrp` environment using the `environment.yml` file but run `python setup.py install` with the forked KCWI\_DRP

Dependencies: `zsh`, `astropy`, `gsed` (Mac), `python 3`, emoji terminal support :sunglasses:
Additional Dependencies (mostly for post-DRP): `tqdm`, `colorcet`, `bokeh 2.0.0`
Location: `/path/to/dir/kcwi/pyDRP/mater_reduce` i.e. kcwi directory structure is intact!

## Prerequisites

1. `cd` to the raw data directory (you can specify this with the `-d` flag otherwise the code will assume you’re in the raw directory)

2. Construct a list of science target base names - this helps the code differentiate standard stars and twilights from science exposures
  - See a sample text file at [`sample_sci_objs.txt`](../examples/sample_sci_objs.txt) for reference
  - Example: `echo kb210115_000{53..71} | xargs -n 1 > sci.txt` if the science images were numbered 53 through 71.
  - File name is arbitrary (`sci.txt` is default)

## Running Straight Through

### Regions files (`.reg`, `_icube_2d.reg`, and `_icube.thum.reg`) are needed for this

Simply run `master_reduce` specifying the file with science object names.
If `master_reduce` is located in `~/Software/kcwi/pyDRP`, then run `~/Software/kcwi/pyDRP/master_reduce -o sci.txt`

Or, if you want to modify the default configuration and supply your own config file (`kcwi_qso_saturated.cfg`): `~/Software/kcwi/pyDRP/master_reduce -o sci.txt -c kcwi_qso_saturated.cfg`

See `kcwi/pyDRP` and `KCWI_DRP/kcwidrp/configs` for pre-made config files based on individual use cases.


## Making regions files as you go along

If the regions files do not exist (or only some do), use the stage flag (`-s`) and provide a comma separated list of stages that need to be run.
For instance we’d normally do something like `stage1drp`, make the `.reg` files, `stage1plusToStage2pre`, make the second set of regions files, and finally `stage2postregToEnd`.

The latter two stages are functions which compose several shorter functions - these can be mixed and matched as required and the list can be as long as needed (e.g. `~/Software/kcwi/pyDRP/master_reduce -o sci.txt -s stage3drp,stage3plus,stage4drp`).

Essentially we run a stage of the DRP and then need to cleanup the proc table and files before proceeding to the next stage. By and large, the `*drp` functions run the DRP and the `*plus` functions run the cleanup and prep for the next stage.

## Miscellaneous
See `reg_construction.md` for information and instructions on how to construct `.reg` files.

Time saving option for your `.zshrc` file: `alias mr="~/Software/kcwi/pyDRP/master_reduce"`
