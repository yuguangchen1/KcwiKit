# KCWI Post-DRP Stacking Installation Instructions
## Initial Version: 12/29/2023

# Montage

This KCWI stacking code is based on the [`Montage`](http://montage.ipac.caltech.edu/docs/download2.html) package. As such, this is the first major package that needs to get installed.

The current GitHub repository is located [here](https://github.com/Caltech-IPAC/Montage) and one can install either the command line version or `MontagePy` with `pip install MontagePy`. Be advised that it may be necessary to use the `montagepy = True` flag in `kcwi_stack` if you choose to go this route.

It may be necessary to add a line like `export PATH=$PATH:/Users/nik/Software/Montage/bin` to your `.zshrc` (or `.bashrc`) file in order to run `Montage` directly from the command line.

To test your installation you can try running `mProjectCube` in a terminal. Doing so should provide the following output: 

```
$ mProjectCube

[struct stat="ERROR", msg="Usage: mProjectCube [-z factor][-d level][-s statusfile][-h hdu][-x scale][-w weightfile][-t threshold][-X(expand)][-e(nergy-mode)][-f(ull-region)] in.fits out.fits hdr.template"]
```

# `kcwi` stacking code

The KCWI stacking code can be downloaded from the main GitHub [page](https://github.com/yuguangchen1/kcwi/tree/master). One can run `git clone https://github.com/yuguangchen1/kcwi.git` in whatever directory they'd like to house the code.

The next step is to add that directory to your `PYTHONPATH` with a line like `export PYTHONPATH=/Users/nik/Software/kcwi/py` in your `.zshrc` file. 

If this was successful, you should be able to run `import kcwi` in an `ipython` shell without error.

# Next steps

The stacking code requires `Montage` and the `kcwi` packages to be installed and functioning. Some additional necessary packages include `astropy`, `fpdf`, and `tqdm`.

You should now be able to stack cubes - see the instructions [here](https://github.com/yuguangchen1/kcwi/blob/master/docs/KCWI_post-DRP_stacking.md)!
