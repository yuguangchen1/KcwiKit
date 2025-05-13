import glob, os, sys, re
import argparse
import json
import numpy as np
import pyregion
import argparse
from astropy.io import fits
import pathlib
fileloc = str(pathlib.Path(__file__).parent.resolve())


def parser_init():
    description = 'Conduct median filtering and subtract the median-filtered background.'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('files', type=str, nargs='+',
            help='Region file name; must end with _rm.reg')
    parser.add_argument('-r','--red', dest='red',
            action='store_true',
            help='Red channel ZAPPED cube.')
    parser.add_argument('-b','--blue', dest='blue',
            action='store_true',
            help='Blue channel flux calibrated cube.')
    args = parser.parse_args()
    return args

def collapse_header(hdr):
    """
    Quick wrapper to collapse a 3-D header into a 2-D one.
    Copied from KCWIKit

    Parameters
    ----------
    hdr: header

    Returns
    -------
    hdr_img: collapsed header

    """

    hdr_img=hdr.copy()
    hdr_img['NAXIS']=2
    del hdr_img['NAXIS3']
    del hdr_img['CD3_3']
    del hdr_img['CTYPE3']
    del hdr_img['CUNIT3']
    del hdr_img['CNAME3']
    del hdr_img['CRVAL3']
    del hdr_img['CRPIX3']

    return hdr_img

def makemask(files,blue=False, red=False):
    """Converting region files to binary masks specifically for removal in the final stack.
    Args:
        dir (str): directory of which the region files are stored.
    Returns:
        None. Region files are converted to binary FITS files. 
    """

    #search_for = dir + "/*_rm.reg"

    #for regfn in glob.glob(search_for):
    #        imgfn=regfn.replace('_rm.reg','.fits')
    #        os.system(fileloc+"/kcwi_masksky_ds9_thum.py "+imgfn+" "+regfn)
    for file in files:
        if os.path.isfile(file):
            regfn=file
            if red:
                cubefn=regfn.replace('_rm.reg','_zap_icubes.fits')
                #imgfn=regfn.replace('_rm.reg','_wlimg.fits')
                hdl=fits.open(cubefn)
                scihdr=hdl[0].header
                obswave= (np.arange(scihdr['NAXIS3']) + 1 - scihdr['CRPIX3']) * scihdr['CD3_3'] + scihdr['CRVAL3']
                wlimg_index = np.where((obswave >= 5650.) & (obswave <= 8930.))[0]
                wlimg = np.sum(hdl[0].data[wlimg_index], axis = 0)
                mhdu=fits.PrimaryHDU(wlimg, header = collapse_header(scihdr))
                #mhdu=fits.open(imgfn)[0]
            elif blue:
                cubefn=regfn.replace('_rm.reg','_icubes.fits')
                hdl=fits.open(cubefn)
                scihdr=hdl[0].header
                obswave= (np.arange(scihdr['NAXIS3']) + 1 - scihdr['CRPIX3']) * scihdr['CD3_3'] + scihdr['CRVAL3']
                wlimg_index = np.where((obswave >= 3550.) & (obswave <= 5500))[0]
                wlimg = np.sum(hdl[0].data[wlimg_index], axis = 0)
                mhdu=fits.PrimaryHDU(wlimg, header = collapse_header(scihdr))
                #imgfn=regfn.replace('_rm.reg','_thum.fits')
            
            
            
            
            
            r = pyregion.open(regfn)
            region_mask  = r.get_mask(hdu=mhdu)
            mask2d=np.zeros_like(mhdu.data)
            mask2d[region_mask] = 1
            
            expanded_mask2d = np.repeat(mask2d[np.newaxis, :, :], hdl[0].data.shape[0], axis=0)
            #maskhdu = fits.PrimaryHDU(expanded_mask2d, header = hdl[0].header)
            flagcube=hdl['FLAGS'].data.copy()
            flagcube[expanded_mask2d > 0] += 32
            hdl['FLAGS'].data = flagcube
            print("Updating cubes with mask from", regfn)
            if red:
                fdest=cubefn.replace('_zap_icubes.fits','_zap_new_icubes.fits')
            elif blue:
                fdest=cubefn.replace('_icubes.fits','_new_icubes.fits')
            hdl.writeto(fdest, overwrite = True)
            print(fdest, "created")
            # Create the mask FITS file 
            #filename = regfn.replace('_rm.reg','_rm_mask.fits')
            #maskhdu.writeto(filename, overwrite = True)

def main():
    args = parser_init()
    #args = arg_parser.parse_args()

    # Generate a master parameter table
    makemask(**vars(args))

if __name__ == '__main__':

    main()
