import numpy as np
import os.path as path
import os
import shutil
import glob
from astropy.io import ascii
from astropy.io import fits
from astropy import wcs
from astropy import coordinates
from astropy import units as u
from astropy import time
from astropy import stats
from astropy import table
from .kcwi_regions import regions_to_mask
from reproject import reproject_interp
from reproject import reproject_exact
#from MontagePy.main import mProjectCube
# from PyAstronomy import pyasl
from scipy import interpolate
from scipy import signal
from scipy import ndimage
#from image_registration import chi2_shift
#from image_registration.fft_tools import shift
import matplotlib
import matplotlib.pyplot as plt
from fpdf import FPDF #fpdf2
from tqdm import tqdm
import pdb
import time as ostime
import warnings
# from specutils.utils import wcs_utils
import pathlib

# MontagePy
try:
    from MontagePy.main import *
    flag_montagepy = True
except:
    print('MontagePy not installed, using command-line version.')
    flag_montagepy = False

def get_wav_axis(header):
    """Returns a NumPy array representing the wavelength axis of a cube.

    Args:
        header (astropy.io.fits.Header): header that contains wavelength
            or velocity axis that is specified in 'CTYPE' keywords in any
            dimension.

    Returns:
        numpy.ndarray: Wavelength axis for this data.

    Examples:

        If you wanted to plot your spectrum vs. wavelength in matplotlib:

        >>> import matplotlib.pyplot as plt
        >>> from cwitools import cubes
        >>> from astropy.io import fits
        >>> spec,header = fits.getdata("myspectrum.fits",header=True)
        >>> wav_axis = cubes.get_wav_axis(header)
        >>> fig,ax = plt.subplots(1,1)
        >>> ax.plot(wav_axis,spec)
        >>> fig.show()

    Note: copied from CWITools
    """

    #Select the appropriate axis.
    naxis = header['NAXIS']
    flag = False
    for i in range(naxis):
        #Keyword entry
        card = "CTYPE{0}".format(i+1)
        if not card in header:
            raise ValueError("Header must contain 'CTYPE' keywords.")

        #Possible wave types.
        if header[card] in ['AWAV', 'WAVE', 'VELO']:
            axis = i+1
            flag = True
            break

    #No wavelength axis
    if flag == False:
        raise ValueError("Header must contain a wavelength/velocity axis.")

    #Get keywords defining wavelength axis
    nwav = header["NAXIS{0}".format(axis)]
    wav0 = header["CRVAL{0}".format(axis)]
    dwav = header["CD{0}_{0}".format(axis)]
    pix0 = header["CRPIX{0}".format(axis)]

    #Calculate and return
    return np.array([wav0 + (i - pix0 + 1) * dwav for i in range(nwav)])


# Read parameter files for alignment, stacking, and astrometry
def kcwi_stack_readpar(parname='q0100-bx172.par'):
    """
    Parsing parameter files.

    Args:
        parname (string): file name of the parameter file.

    Returns:
        dict: loaded parameters.

    """

    #parname='q0100-bx172.par'
    par={"align_box":np.array([-1,-1,-1,-1]),
        "align_dimension":np.array([-1,-1]),
        "align_xpix":-1.,
        "align_ypix":-1.,
        "align_orientation":-1000.,
        "align_search_size":-1000,
        "align_conv_filter":-1000,
        "align_upfactor":-1000.,
        "stack_dimension":np.array([-1,-1]),
        "stack_xpix":-1.,
        "stack_ypix":-1.,
        "stack_orientation":-1000.,
        "wavebin":np.array([-1.,-1.]),
        "ref_fn":'',
        "ref_xy":np.array([-1.,-1.]),
        "ref_ad":np.array([-1.,-1.]),
        "ref_search_size":-1000,
        "ref_conv_filter":-1000,
        "ref_upfactor":-1000,
        "ref_box":np.array([-1,-1,-1,-1]),
        "ref_nocrl":0,
        "align_ad":np.array([-1.,-1.]),
        "stack_ad":np.array([-1.,-1.]),
        "stepsig":0,
        "drizzle":0.,
        "med_x":0.,
        "med_y":0.,
        "med_z":0.,
        "background_subtraction":False,
        "background_level":-1000,
        "wave_ref":np.array([0., 0.]),
        "nwave":0,
        "dwave":0.
    }

    with open(parname,'r') as file:
        lins=file.readlines()
        keys=[]
        for i in range(len(lins)):
            if not lins[i].isspace():
                keys.append(lins[i].split()[0].casefold())


    # Alignment keywords
    # one line per box: several align_box lines give several alignment boxes,
    # which kcwi_align cross-correlates together (see _align_crls)
    q=np.where(np.array(keys)=="align_box")[0]
    if len(q) > 0:
        bxs=[]
        for qq in q:
            ele=lins[qq].split()
            tmp=np.array(ele[1:5]).astype(float)
            bxs.append([tmp[0],tmp[2],tmp[1],tmp[3]])
        bxs=np.array(bxs).astype(int)
        par['align_box']=bxs[0] if len(bxs)==1 else bxs

    q=np.where(np.array(keys)=="align_dimension")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["align_dimension"]=np.array(ele[1:3]).astype(int)

    q=np.where(np.array(keys)=="align_xpix")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["align_xpix"]=float(ele[1])

    q=np.where(np.array(keys)=="align_ypix")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["align_ypix"]=float(ele[1])

    q=np.where(np.array(keys)=="align_orientation")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["align_orientation"]=float(ele[1])

    q=np.where(np.array(keys)=="align_ad")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["align_ad"]=np.array(ele[1:3]).astype(float)

    q=np.where(np.array(keys)=='align_search_size')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['align_search_size']=float(ele[1])

    q=np.where(np.array(keys)=='align_conv_filter')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['align_conv_filter']=float(ele[1])

    q=np.where(np.array(keys)=='align_upfactor')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['align_upfactor']=float(ele[1])


    # Stacking keywords
    q=np.where(np.array(keys)=="stack_dimension")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_dimension"]=np.array(ele[1:3]).astype(int)

    q=np.where(np.array(keys)=="stack_xpix")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_xpix"]=float(ele[1])

    q=np.where(np.array(keys)=="stack_ypix")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_ypix"]=float(ele[1])

    q=np.where(np.array(keys)=="stack_orientation")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_orientation"]=float(ele[1])

    q=np.where(np.array(keys)=="stack_ad")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_ad"]=np.array(ele[1:3]).astype(float)


    # Global keywords
    q=np.where(np.array(keys)=="dimension")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_dimension"]=np.array(ele[1:3]).astype(int)
        par["align_dimension"]=np.array(ele[1:3]).astype(int)

    q=np.where(np.array(keys)=="xpix")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_xpix"]=float(ele[1])
        par["align_xpix"]=float(ele[1])

    q=np.where(np.array(keys)=="ypix")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_ypix"]=float(ele[1])
        par["align_ypix"]=float(ele[1])

    q=np.where(np.array(keys)=="orientation")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_orientation"]=float(ele[1])
        par["align_orientation"]=float(ele[1])

    q=np.where(np.array(keys)=="ad")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_ad"]=np.array(ele[1:3]).astype(float)
        par["align_ad"]=np.array(ele[1:3]).astype(float)

    q = np.where(np.array(keys)=='wave_ref')[0]
    if len(q) > 0:
        q = q[0]
        ele = lins[q].split()
        par['wave_ref'] = np.array(ele[1:3]).astype(float)

    q = np.where(np.array(keys)=='nwave')[0]
    if len(q) > 0:
        q = q[0]
        ele = lins[q].split()
        par['nwave'] = int(ele[1])

    q = np.where(np.array(keys)=='dwave')[0]
    if len(q) > 0:
        q = q[0]
        ele = lins[q].split()
        par['dwave'] = float(ele[1])

    q=np.where(np.array(keys)=='search_size')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['align_search_size']=float(ele[1])
        par['ref_search_size']=float(ele[1])

    q=np.where(np.array(keys)=='conv_filter')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['align_conv_filter']=float(ele[1])
        par['ref_conv_filter']=float(ele[1])

    q=np.where(np.array(keys)=='upfactor')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['align_upfactor']=float(ele[1])
        par['ref_upfactor']=float(ele[1])


    # wavebin
    q=np.where(np.array(keys)=="wavebin")[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par["wavebin"]=np.array(ele[1:3]).astype(float)

    # astrometry
    q=np.where(np.array(keys)=="ref_xy")[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par["ref_xy"]=np.array(ele[1:3]).astype(float)

    q=np.where(np.array(keys)=="ref_ad")[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par["ref_ad"]=np.array(ele[1:3]).astype(float)

    q=np.where(np.array(keys)=='ref_fn')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['ref_fn']=ele[1]

    q=np.where(np.array(keys)=='ref_search_size')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['ref_search_size']=float(ele[1])

    q=np.where(np.array(keys)=='ref_conv_filter')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['ref_conv_filter']=float(ele[1])

    q=np.where(np.array(keys)=='ref_upfactor')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['ref_upfactor']=float(ele[1])

    q=np.where(np.array(keys)=='ref_box')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        tmp=np.array(ele[1:5]).astype(float)
        par['ref_box']=np.array([tmp[0],tmp[2],tmp[1],tmp[3]]).astype(int)

    q=np.where(np.array(keys)=='ref_nocrl')[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['ref_nocrl']=float(ele[1])


    q=np.where(np.array(keys)=="stepsig")[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par["stepsig"]=float(ele[1])


    # drizzle
    q=np.where(np.array(keys)=="drizzle")[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par["drizzle"]=float(ele[1])


    # median filtering
    q=np.where(np.array(keys)=="med_x")[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par['med_x']=float(ele[1])

    q=np.where(np.array(keys)=='med_y')[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par['med_y']=float(ele[1])

    q=np.where(np.array(keys)=='med_x')[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par['med_z']=float(ele[1])

    # background subtraction
    q=np.where(np.array(keys)=="background_subtraction")[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['background_subtraction']=bool(float(ele[1]))

    q=np.where(np.array(keys)=="background_level")[0]
    if len(q)>0:
        q=q[0]
        ele=lins[q].split()
        par['background_level']=float(ele[1])


    return par


def a2v_conversion(wave):
    """ Convert air-based wavelengths to vacuum

    Adapted from wave.py in: https://github.com/pypeit/PypeIt/
    Formula from https://ui.adsabs.harvard.edu/abs/1996ApOpt..35.1566C/
    This is the Goddard IDL Air to Vacuum Conversion. Stolen from KCWI_DRP.

    Parameters
    ----------
    wave: Quantity array
        Wavelengths
    Returns
    -------
    wave: Quantity array
        Wavelength array corrected to vacuum wavelengths
    """
    # Convert to AA
    wave = wave.to(u.AA)
    wavelength = wave.value

    # Standard conversion format
    sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
    factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
    factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

    # Convert
    wavelength = wavelength*factor
    # Units
    new_wave = wavelength*u.AA
    new_wave.to(wave.unit)

    return new_wave

# Perform air-to-vac, and barycentric correction in wavelength

def empty_hdul_like(hdul):
    """
    Create an empty HDUList with the same structure as an existing HDUList.
    Data arrays are replaced with zeros of the same shape/dtype.
    This is to solve the problem that hdl.copy() is a shallow copy and still links the data arrays.
    """
    new_hdus = []
    for i, hdu in enumerate(hdul):
        if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)) and \
                not isinstance(hdu, fits.CompImageHDU):
            # Tables are copied as they are (deep copy of the record array).
            new_hdu = hdu.__class__(data=hdu.data.copy() if hdu.data is not None else None,
                                    header=hdu.header.copy(), name=hdu.name)
        elif hdu.data is None:
            # Header-only HDU: keep the header, no data.
            cls = fits.PrimaryHDU if i == 0 else fits.ImageHDU
            new_hdu = cls(data=None, header=hdu.header.copy())
        else:
            empty_data = np.zeros_like(hdu.data)
            cls = fits.PrimaryHDU if i == 0 else fits.ImageHDU
            new_hdu = cls(data=empty_data, header=hdu.header.copy())
        new_hdus.append(new_hdu)
    return fits.HDUList(new_hdus)

def kcwi_vachelio(hdl, hdr_ref=None, mask=False, method='barycentric'):
    """
    Convert air wavelength axis to vacuum wavelength axis. Correct heliocentric
        velocity.

    Args:
        hdl (astropy.io.fits.HDUList or astropy.io.fits.PrimaryHDU): HDUList/HDU
            of the input data cube.
        hdr_ref (astropy.io.fits.Header): reference header if different to
            the original header.
        mask (bool): set if 'hdu' is a mask cube. Only act when input is hdu.
        method (bool): 'heliocentric' or 'barycentric'

    Returns:
        astropy.io.fits.Primary: new HDU with vacuum, and heliocentric-corrected
            wavelength axis.
        float: heliocentric velocity.

    """

    if isinstance(hdl, fits.HDUList):
        inputflag = 'hdl'
        hdl_in = hdl#.copy() this is useless since it's a shallow copy
    else:
        inputflag = 'hdu'
        hdl_in = fits.HDUList(hdl)

    # Check if wavelength has already been corrected
    flag_vac = False
    flag_vcorr = False

    hdr_in = hdl_in[0].header
    if 'CTYPE3' not in hdr_in:
        raise ValueError('Wavelength info must be in header.')
        return

    if hdr_in['CTYPE3']=='WAVE':
        flag_vac = True
        print("Skipping air-to-vac conversion, already done.")
    if ('VCORR' in hdr_in) and (np.abs(hdr_in['VCORR'])>=1e-5):
        flag_vcorr = True
        print("Skipping heliocentric/barycentric correction, already done.")

    if flag_vac and flag_vcorr:
        # nothing to be done
        print('kcwi_vachelio - Already processed. Returning.')
        if inputflag=='hdl':
            return hdl_in, hdr_in['VCORR']
        elif inputflag=='hdu':
            return hdl_in[0], hdr_in['VCORR']

    if hdr_ref is None:
        hdr_new=hdl_in[0].header.copy()
    else:
        hdr_new=hdl_in[0].header.copy()
        hdr_new['NAXIS3']=hdr_ref['NAXIS3']
        hdr_new['CTYPE3']=hdr_ref['CTYPE3']
        hdr_new['CUNIT3']=hdr_ref['CUNIT3']
        hdr_new['CNAME3']=hdr_ref['CNAME3']
        hdr_new['CRVAL3']=hdr_ref['CRVAL3']
        hdr_new['CRPIX3']=hdr_ref['CRPIX3']
        hdr_new['CD3_3']=hdr_ref['CD3_3']

    hdr_new['CTYPE3']='WAVE'

    # setup output
    hdl_out = empty_hdul_like(hdl_in)
    hdl_out[0].header = hdr_new


    wave_old = get_wav_axis(hdr_in)
    wave_new = get_wav_axis(hdr_new)

    # air -> vac
    if flag_vac==False:
        # wave_vac=pyasl.airtovac2(wave_old) # PyAstronomy doesn't work on Mac M1 8/14/22?
        # wave_vac=wcs_utils.air_to_vac(wave_old*u.Angstrom).value # spec_utils, not using Goddard conversion
        wave_vac=a2v_conversion(wave_old*u.AA).value # Goddard Conversion, matches KCWI_DRP
        print("Converting air to vacuum wavelengths.")
    else:
        wave_vac = wave_old.copy()

    if flag_vcorr==False:
        targ=coordinates.SkyCoord(hdr_in['TARGRA'],hdr_in['TARGDEC'],unit='deg',obstime=hdr_in['DATE-BEG'])
        keck=coordinates.EarthLocation.of_site('Keck Observatory')
        vcorr=targ.radial_velocity_correction(kind=method,location=keck)
        vcorr = vcorr.to('km/s').value

        wave_hel=wave_vac*(1+vcorr/2.99792458e5)
        print("Applying %s correction of %.2f km/s." %(method, vcorr))
    else:
        wave_hel = wave_vac.copy()
        vcorr = 0.

    #resample
    for i_ext, hdu_in in enumerate(hdl_in):
        skipflag=False
        #cube_old = np.nan_to_num(hdu_in.data)
        cube_old = np.nan_to_num(hdu_in.data,posinf=0.,neginf=0.)

        shape_new = (hdr_new['NAXIS3'], hdr_new['NAXIS2'], hdr_new['NAXIS1'])
        # How to resample this extension:
        #   'science'  - cubic interpolation of the values
        #   'uncert'   - linear interpolation of the variance (sigma**2), then sqrt;
        #                cubic can overshoot and nearest-neighbour misaligns sigma
        #                and flux by up to half a channel
        #   'mask'     - nearest-neighbour ('nearest-up') so integer flags stay integers
        if inputflag=='hdu':
            kind = 'mask' if mask else 'science'
            fill_value = 128 if mask else np.nan
        else:
            if 'EXTNAME' not in hdu_in.header:
                kind = 'science'
                fill_value = np.nan
                print("Stacking primary HDU.")
            else:
                if hdu_in.header['EXTNAME']=='MASK':
                    print(hdu_in.header['EXTNAME'])
                    kind = 'mask'
                    fill_value = 1
                elif hdu_in.header['EXTNAME']=='FLAGS':
                    print(hdu_in.header['EXTNAME'])
                    kind = 'mask'
                    fill_value = 128
                elif hdu_in.header['EXTNAME']=='UNCERT':
                    print(hdu_in.header['EXTNAME'])
                    kind = 'uncert'
                    fill_value = np.nan
                else:
                    print("Warning: EXTNAME not recognized, skipping HDU.")
                    skipflag=True
                    #print("Warning: EXTNAME not recognized, assuming science.")
                    #mask_flag = False
                    #fill_value = np.nan

        # interpolation
        if skipflag:
            hdl_out[i_ext].data = cube_old
            continue
        cube_new=np.zeros(shape_new)
        for i in range(shape_new[2]):
            for j in range(shape_new[1]):
                spec = cube_old[:, j, i]
                if kind == 'science':
                    f_cubic=interpolate.interp1d(wave_hel,spec,kind='cubic',fill_value=fill_value,bounds_error=False)
                    spec_new=f_cubic(wave_new)
                elif kind == 'uncert':
                    f_var=interpolate.interp1d(wave_hel,spec**2,kind='linear',fill_value=fill_value,bounds_error=False)
                    spec_new=np.sqrt(np.clip(f_var(wave_new), 0, None))
                else:
                    f_nearest=interpolate.interp1d(wave_hel,spec,kind='nearest-up',bounds_error=False,fill_value=fill_value)
                    spec_new=f_nearest(wave_new)
                cube_new[:, j, i] = spec_new

        hdl_out[i_ext].data = cube_new

#   plt.clf()
#   plt.plot(wave_old,cube_old[:,45,15],drawstyle='steps-mid')
#   plt.plot(wave_new,cube_new.T[:,45,15],drawstyle='steps-mid')
#   plt.xlim(3300,3500)
#   plt.ylim(-0.01,0.01)

    if inputflag=='hdl':
        return hdl_out, vcorr
    elif inputflag=='hdu':
        return hdl_out[0], vcorr


def kcwi_checkexptime(dir='./',redux=False):

    if type(dir) is fits.hdu.image.PrimaryHDU:
        hdu=dir
        hdu_all=[hdu]
    else:
        if redux==True:
            dir='./redux/'

        print('Checking EXPTIME...')

        fn_all=glob.glob(dir+'/kb*.fits')

        hdu_all=[(fits.open(fn))[0] for fn in fn_all]

    copyflag=np.zeros(len(hdu_all),dtype=bool)
    for i in range(len(hdu_all)):
        exptime=time.TimeDelta(hdu_all[i].header['XPOSURE'],format='sec')
        expend=time.Time(hdu_all[i].header['DATE-END'],format='isot')
        rdend=time.Time(hdu_all[i].header['DATEREND'],format='isot')
        if rdend <= expend:
            copyflag[i]=True
            if not type(dir) is fits.hdu.image.PrimaryHDU:
                print(os.path.basename(fn_all[i]))
            expbeg=time.Time(hdu_all[i].header['DATE-BEG'],format='isot')
            rdtime=time.TimeDelta(53.64,format='sec')
            expend=rdend-rdtime
            exptime=expend-expbeg

            hdu_all[i].header['XPOSURE']=exptime.sec
            hdu_all[i].header['TELAPSE']=exptime.sec+0.005

            print('    Setting EXPTIME = '+str(exptime.sec))

    if type(dir) is fits.hdu.image.PrimaryHDU:
        return hdu_all[0]
    else:
        if not os.path.exists(dir+'/old'):
            os.makedirs(dir+'/old')
        for i in range(len(fn_all)):
            if copyflag[i]==True:
                shutil.copyfile(fn_all[i],dir+'/old/')
                hdu_all[i].writeto(fn_all[i])

        return True




def kcwi_check_flux(fnlist, thumfn=None, nsig=1.5, cubed=False):
    """
    Check the relative flux difference of overlapping sources in a list of FITS files
        for better flux calibration.

    Args:
        fnlist (str): file name of a list file in which the individual FITS file names
            are stored in 'col1'.
        thumfn (array of str): set if thum files are stored in different directories.
        nsig (float): number of sigmas in sigma-clipping.
        cubed (bool): cubed files?

    Returns:
        None

    """

    tab = ascii.read(fnlist)
    fn = tab['col1']

    if thumfn is None:
        thumfn = 'kcwi_align/'+fnlist.replace('.list','.thum.fits')
    hdu_thum = fits.open(thumfn)[0]

    # convert to SB units
    if cubed==False:
        for i in range(hdu_thum.shape[0]):
            hdu_i = fits.open(fn[i]+'_icubes.fits')[0]
            dx=np.sqrt(hdu_i.header['CD1_1']**2+hdu_i.header['CD2_1']**2)*3600.
            dy=np.sqrt(hdu_i.header['CD1_2']**2+hdu_i.header['CD2_2']**2)*3600.
            area=dx*dy
            hdu_thum.data[i,:,:] = hdu_thum.data[i,:,:] / area

    # calculate relative flux
    sig = np.nanstd(hdu_thum.data, axis=(1,2))
    med = np.nanmedian(hdu_thum.data, axis=(1,2))

    frame_all = []
    flux_rel_all = []
    flux_rel = np.zeros(hdu_thum.shape[0])
    for i in range(hdu_thum.shape[0]):
        index = (hdu_thum.data[i, :, :] >= nsig * sig[i] + med[i])
        if i==0:
            index0 = index.copy()

        tmp0 = hdu_thum.data[0, (index & index0)].flatten()
        tmp = hdu_thum.data[i, (index & index0)].flatten()

        frame_all = np.append(frame_all, np.repeat(i, len(tmp)))
        flux_rel_all = np.append(flux_rel_all, tmp/tmp0)

        if np.sum(index & index0) <= 5:
            continue
        flux_rel[i] = np.median(tmp/tmp0)


    fig, ax = plt.subplots(figsize=(10,6))
    ax.scatter(frame_all, flux_rel_all, s=10, color='C0', alpha=0.5)
    ax.plot(np.arange(len(flux_rel)), flux_rel, 'o-', color='C1')
    ax.set_xlabel('Frame #')
    ax.set_ylabel('Relative Flux')
    xlim = ax.get_xlim()
    ax.plot(xlim, [1,1], '--', color='black')
    ax.set_yscale('log')
    ax.set_xlim(xlim)
    plt.show()

    return


def kcwi_norm_flux(fnlist, frame=[], thumfn=None, nsig=1.5, cubed=False):
    """
    Generate a table of flux correction factor.

    Args:
        fnlist (str): file name of the list that contains individual FITS files.
        frame (list-like): list of frame numbers that need to be corrected.
        thumfn (list): list of thum files if stored in non-default locations.
        nsig (float): number of sigmas for sigma-clip.
        cubed (bool): using cubed files?

    Returns:
        astropy.ascii.table.Table: table containing the normalizing factor to be
            applied. This table is also written as '*.flx.list'.


    """

    if cubed:
        suffix="cubed"
    else:
        suffix="cubes"

    tab = ascii.read(fnlist)
    fn = tab['col1']

    if thumfn is None:
        thumfn = 'kcwi_align/'+fnlist.replace('.list','.thum.fits')
    hdu_thum = fits.open(thumfn)[0]

    # convert to SB units
    for i in range(hdu_thum.shape[0]):
        hdu_i = fits.open(fn[i]+'_i'+suffix+'.fits')[0]
        dx=np.sqrt(hdu_i.header['CD1_1']**2+hdu_i.header['CD2_1']**2)*3600.
        dy=np.sqrt(hdu_i.header['CD1_2']**2+hdu_i.header['CD2_2']**2)*3600.
        area=dx*dy
        hdu_thum.data[i,:,:] = hdu_thum.data[i,:,:] / area


    sig = np.nanstd(hdu_thum.data, axis=(1,2))
    med = np.nanmedian(hdu_thum.data, axis=(1,2))

    # relative flux
    flux_rel = np.zeros(hdu_thum.shape[0]) + np.nan
    flag = np.zeros(hdu_thum.shape[0])
    for i in range(hdu_thum.shape[0]):

        if i in frame:
            flag[i] = 1

        index = (hdu_thum.data[i, :, :] >= nsig * sig[i] + med[i])
        if i==0:
            index0 = index.copy()

        tmp0 = hdu_thum.data[0, (index & index0)].flatten()
        tmp = hdu_thum.data[i, (index & index0)].flatten()

        if np.sum(index & index0) <= 5:
            continue
        flux_rel[i] = np.median(tmp/tmp0)

    flux_mean = np.nanmean(flux_rel[flag == 0])

    # Correction factor
    flux_corr = np.zeros(hdu_thum.shape[0])+1
    for i in frame:
        flux_corr[i] = flux_mean / flux_rel[i]


    fluxtable=table.Table([np.array(fn), flux_corr])
    ascii.write(fluxtable, fnlist.replace('.list','.flx.list'),overwrite=True,format='no_header')

    return fluxtable


def kcwi_check_samewave(hdr0, hdr1):
    """
    Check if the wavelength axes are the same in two headers.

    Args:
        hdr0 (astropy.io.fits.header) - input header #0
        hdr1 (astropy.io.fits.header) - input header #1

    Returns:
        boolean: whether wave axes are the same
    """

    if hdr0['NAXIS3'] != hdr1['NAXIS3']:
        print("different NAXIS3: %d vs %d" % (hdr0['NAXIS3'], hdr1['NAXIS3']))
        # Not the same amount of pixels
        return False

    wave0 = (np.arange(hdr0['NAXIS3']) - hdr0['CRPIX3'] + 1) * hdr0['CD3_3'] + hdr0['CRVAL3']
    wave1 = (np.arange(hdr1['NAXIS3']) - hdr1['CRPIX3'] + 1) * hdr1['CD3_3'] + hdr1['CRVAL3']

    if not np.isclose(wave0[0], wave1[0]):
        print("different starting point: %.3f vs %.3f" % (wave0[0], wave1[0]))
        # Starting point different
        return False
    if not np.isclose(wave0[1] - wave0[0], wave1[1] - wave1[0]):
        print("different delta w: %.3f vs %.3f" % (wave0[1] - wave0[0], wave1[1] - wave1[0]))
        # delta w different
        return False

    return True

def kcwi_resample_wave(hdu, newhdr, method='cubic',plot=False):
    """
    Resample a cube to match the wavelength direction to a different header.

    Args:
        hdu (astropy HDU): input hdu
        newhdr (astropy header): header containing the new wavelength grid
        order (str): interpolation method or 'mask'

    Returns:
        astropy.io.fits.PrimaryHDU: resampled data cube in the form of HDU
    """

    hdr = hdu.header
    wave = (np.arange(hdr['NAXIS3']) - hdr['CRPIX3'] + 1) * hdr['CD3_3'] + hdr['CRVAL3']
    data = hdu.data.copy()
    #print("oldwave with length %d" % hdr['NAXIS3'],"starting at %.2f" % hdr['CRVAL3'], " delta %.4f" % hdr['CD3_3'])
    #print("newwave with length %d" % newhdr['NAXIS3'],"starting at %.2f" % newhdr['CRVAL3'], " delta %.4f" % newhdr['CD3_3'])
    newwave = (np.arange(newhdr['NAXIS3']) - newhdr['CRPIX3'] + 1) * newhdr['CD3_3'] + newhdr['CRVAL3']
    newdata = np.zeros((len(newwave), hdu.shape[1], hdu.shape[2])) + np.nan
    #print("newdata with shape", newdata.shape)
    #print("olddata with shape", data.shape)

    data = data.reshape(len(wave), -1)
    newdata = newdata.reshape(len(newwave), -1)
    #print(newdata.shape)
    for i in range(newdata.shape[1]):
        spec = data[:, i]

        if method != 'mask':
            mask = ~np.isfinite(spec)
            spec = np.nan_to_num(spec)

            # No good data, skip
            if np.sum(spec)==0:
                continue

            ci = interpolate.interp1d(wave, spec, kind=method, bounds_error=False, fill_value=np.nan)
            #print(wave)
            #print(newwave)
            
            newspec = ci(newwave)

            mi = interpolate.interp1d(wave, mask, kind='linear', bounds_error=False, fill_value=1)
            newmask = mi(newwave)

            newspec[newmask != 0] = np.nan

        else:
            # mask cube
            mask = spec

            mi = interpolate.interp1d(wave, mask, kind='linear', bounds_error=False, fill_value=128)
            newmask = mi(newwave)

            newspec = newmask

        newdata[:, i] = newspec

    newdata = newdata.reshape((len(newwave), hdu.shape[1], hdu.shape[2]))

    if plot:
        plt.plot(wave,hdu.data[:,68,15],"g")
        plt.plot(newwave,newdata[:,68,15],"r")
    newhdu = hdu.copy()
    newhdu.header['NAXIS3'] = newhdr['NAXIS3']
    newhdu.header['CRPIX3'] = newhdr['CRPIX3']
    newhdu.header['CRVAL3'] = newhdr['CRVAL3']
    newhdu.header['CD3_3'] = newhdr['CD3_3']
    newhdu.data = newdata

    return newhdu



def _crr_zsmooth(a, size, step=1):
    """Running median along wavelength, evaluated every `step` channels.

    Used for the per-spaxel scale between two cubes, which is smooth along
    wavelength by construction, so there is no point evaluating it at every
    channel: scipy's median_filter costs O(size) per output voxel and is ~80%
    of the cosmic-ray rejection's run time on a long stack. Evaluating it on a
    decimated grid and interpolating is about eight times faster and agrees to
    ~1%, which then largely cancels in the ratio of two such medians.

    Args:
        a (ndarray): [nz, ny, nx].
        size (int): width of the running median [channels].
        step (int): evaluate every this many channels; 1 reproduces
            ndimage.median_filter(a, (size, 1, 1), mode='nearest') exactly.

    Returns:
        ndarray: same shape as `a`.
    """
    nz = a.shape[0]
    step = max(int(step), 1)
    if step == 1 or nz < 2 * step:
        return ndimage.median_filter(a, size=(int(size), 1, 1), mode='nearest')
    h = int(size) // 2
    zc = np.arange(0, nz, step)
    if zc[-1] != nz - 1:
        zc = np.append(zc, nz - 1)
    sm = np.empty((len(zc),) + a.shape[1:], dtype=a.dtype)
    for i, z in enumerate(zc):
        sm[i] = np.median(a[max(z - h, 0):min(z + h + 1, nz)], axis=0)
    w = np.interp(np.arange(nz), zc, np.arange(len(zc)))
    lo = np.floor(w).astype(int)
    hi = np.minimum(lo + 1, len(zc) - 1)
    fr = (w - lo).reshape((-1,) + (1,) * (a.ndim - 1))
    return sm[lo] * (1 - fr) + sm[hi] * fr


def _crr_flag_outliers(I, V, good, nsig=5., rel=0.5, window=5, zwindow=3, zmed=15,
                       grow_xy=2, grow_z=2, grow_nsig=2., pad_xy=0, loo=None,
                       zratio=51, zstep=8, fsys=0.05, fgrad=0.5, grow_rel=0.1,
                       smax=5.):
    """
    Flag positive outliers (cosmic rays) in each of N cubes on a common grid by
    comparing each cube with the other ones (used by kcwi_stack with crr_split).

    For cube k, the comparison estimate `est` is the other cube (N=2) or the
    median of the others (N>=3). The excess D = I_k - est has a running median
    along wavelength (zmed pixels) subtracted, so smooth differences from
    seeing, alignment or flux calibration are removed; cosmic rays are narrower
    than the instrumental resolution. A voxel is flagged only if BOTH
        absolute:  D > nsig * sqrt(V_k + var_est)
        relative:  D > rel * L,  L = local max of |est| in a (zwindow, window, window) box.
    Flags are then grown, first spectrally (up to grow_z pixels) and then
    spatially (up to grow_xy pixels, 8-connected), only into connected
    neighbors with D > grow_nsig * sigma. pad_xy adds an unconditional
    spatial margin.

    With N >= 5 cubes the median of ALL of them is used as the comparison
    estimate (one outlier does not move it), which is computed once instead of
    once per cube; with fewer, the median of the others (leave-one-out) is used.
    `loo` overrides that choice.


    Three guards keep a bright quasar from being clipped, where frame-to-frame
    seeing and transparency differences are large in absolute terms:
      - `zratio`: each cube is scaled to the comparison estimate per spaxel,
        smoothly along wavelength, before the difference is taken;
      - `fgrad`: a term proportional to the spatial gradient is added to the
        noise, so that a sub-pixel registration error is not significant on
        the steep flanks of a point source;
      - `fsys`: a fractional term is added as well, but see below -- at the
        default thresholds it never binds;
      - `grow_rel`: the relative test is applied while growing too.
    Set zratio=0, fsys=0, fgrad=0, grow_rel=0 for the previous behaviour.

    Args:
        I, V (ndarray [N, nz, ny, nx]): intensity and variance.
        good (bool ndarray [N, nz, ny, nx]): valid voxels.
        zratio (int): channels over which the per-spaxel scale between a cube
            and the comparison estimate is measured; 0 disables the matching.
        zstep (int): channels between evaluations of that scale, which is
            smooth along wavelength. 1 evaluates it everywhere and is the
            slowest part of this function; 8 is ~8x faster, agrees to ~1% and
            largely cancels in the ratio.
        fsys (float): fractional systematic added in quadrature to the noise.
            Inert at the default thresholds, and deliberately kept anyway.
            Because L is the local MAXIMUM of |scale*est|, L >= |scale*est|
            always, so the part of the absolute test this term controls,
                D > nsig * fsys * |scale*est|,
            is weaker than the relative test D > rel * L whenever
                nsig * fsys <= rel,
            which at the defaults is 5 * 0.05 = 0.25 against 0.5. Measured on
            a synthetic bright quasar it exceeded the relative threshold in
            0.000% of the bright voxels (median ratio 0.199, maximum 0.500,
            the bound). The relative test is itself a fractional floor, ten
            times stronger, so this one does nothing while rel is at 0.5.
            It matters when you LOWER rel: rel = 0.5 means a cosmic ray on
            the core of a quasar has to exceed half the quasar to be caught,
            so if you want those you reduce rel, and fsys is then the only
            brightness-aware guard left. Raising fsys above rel / nsig makes
            it active; do not do that without re-checking that real cosmic
            rays are still found.
        fgrad (float): the same for the spatial gradient of the estimate,
            which is what a sub-pixel registration error turns into: a frame
            offset by d pixels differs by about d * |grad I|, so fgrad is that
            offset in pixels. Unlike fsys this is not bounded by L -- on a
            steep flank the per-pixel gradient is a large fraction of the peak
            -- so it does bind, in about 28% of bright voxels, and removes
            roughly a quarter of the remaining false positives.
        grow_rel (float): relative threshold used while growing.
        smax (float): the per-spaxel scale is clipped to [1/smax, smax].

    Returns:
        bool ndarray [N, nz, ny, nx]: outlier flags.
    """
    N = I.shape[0]
    Ig = np.where(good, I, np.nan)
    Vg = np.where(good, V, np.nan)
    flags = np.zeros(I.shape, dtype=bool)
    s_xy = np.zeros((3, 3, 3), dtype=bool)
    s_xy[1] = True
    s_z = np.zeros((3, 3, 3), dtype=bool)
    s_z[:, 1, 1] = True
    if loo is None:
        loo = (N < 5)
    Esm_all = None
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        if not loo:
            # median of all cubes: robust to one outlier, computed once
            n_all = np.sum(good, axis=0)
            est_all = np.nanmedian(Ig, axis=0)
            var_all = (np.where(n_all >= 3, np.pi / 2., 1.) *
                       np.nanmean(Vg, axis=0) / np.maximum(n_all, 1))
        for k in range(N):
            if not loo:
                est, var_est, n_oth = est_all, var_all, n_all - good[k]
            else:
                oth = [j for j in range(N) if j != k]
                n_oth = np.sum(good[oth], axis=0)
                if len(oth) == 1:
                    est, var_est = Ig[oth[0]], Vg[oth[0]]
                else:
                    est = np.nanmedian(Ig[oth], axis=0)
                    corr = np.where(n_oth >= 3, np.pi / 2., 1.)
                    var_est = corr * np.nanmean(Vg[oth], axis=0) / np.maximum(n_oth, 1)
            testable = good[k] & (n_oth > 0) & np.isfinite(est) & np.isfinite(var_est)
            e0 = np.nan_to_num(est)
            # Scale match before differencing. Seeing and transparency make one
            # frame a smoothly wavelength-dependent multiple of the others at
            # every spaxel, and high-passing that DIFFERENCE turns every narrow
            # spectral feature into a positive spike: where an absorption line
            # cuts the flux, the deficit of a fainter frame shrinks, and the
            # running median leaves a cosmic-ray-shaped residual. Dividing the
            # mismatch out first removes it at the source. The scale is the
            # ratio of the two cubes smoothed over zratio channels (and 3x3
            # spatially, so the faint halo is matched too), which cosmic rays,
            # being a few channels wide, cannot move.
            if zratio and zratio > 1:
                nsm = int(zratio) * 9
                Ism = ndimage.uniform_filter(
                    _crr_zsmooth(np.nan_to_num(Ig[k]), zratio, zstep),
                    size=(1, 3, 3), mode='nearest')
                if loo:
                    Esm = ndimage.uniform_filter(
                        _crr_zsmooth(e0, zratio, zstep), size=(1, 3, 3), mode='nearest')
                else:
                    # est is the same cube for every k, so smooth it once
                    if Esm_all is None:
                        Esm_all = ndimage.uniform_filter(
                            _crr_zsmooth(e0, zratio, zstep), size=(1, 3, 3), mode='nearest')
                    Esm = Esm_all
                ok_s = Esm > 3. * np.sqrt(np.nan_to_num(var_est) / nsm)
                scale = np.where(ok_s & (Esm != 0), Ism / np.where(Esm == 0, 1., Esm), 1.)
                scale = np.clip(np.nan_to_num(scale, nan=1.), 1. / smax, smax)
            else:
                scale = np.ones_like(e0)
            se = scale * e0
            D = np.where(testable, I[k] - se, 0.)
            if zmed and zmed > 1:
                D = D - ndimage.median_filter(D, size=(int(zmed), 1, 1), mode='nearest')
            # Noise budget. The photon term alone scales as sqrt(I), so on a
            # bright quasar any residual percent-level mismatch is formally
            # significant; fgrad and fsys add the systematic floor the
            # frame-to-frame scatter actually has. fgrad is the one that does
            # the work -- fsys is dominated by the relative test below unless
            # nsig * fsys > rel, see the docstring.
            var = V[k] + scale ** 2 * np.nan_to_num(var_est)
            if fsys:
                var = var + (fsys * np.abs(se)) ** 2
            if fgrad:
                var = var + (fgrad * np.maximum(np.abs(np.gradient(se, axis=1)),
                                                np.abs(np.gradient(se, axis=2)))) ** 2
            sig = np.sqrt(np.where(testable, var, np.inf))
            L = ndimage.maximum_filter(np.abs(se),
                                       size=(int(zwindow), int(window), int(window)),
                                       mode='nearest')
            f = testable & (D > nsig * sig) & (D > rel * L)
            if f.any():
                cand = testable & (D > grow_nsig * sig)
                if grow_rel:
                    # the relative test has to hold while growing as well, or a
                    # cosmic ray beside a bright source grows into it, where a
                    # 2-sigma condition is met by the mismatch alone
                    cand &= (D > grow_rel * L)
                for _ in range(int(grow_z)):
                    new = ndimage.binary_dilation(f, structure=s_z) & cand & ~f
                    if not new.any():
                        break
                    f |= new
                for _ in range(int(grow_xy)):
                    new = ndimage.binary_dilation(f, structure=s_xy) & cand & ~f
                    if not new.any():
                        break
                    f |= new
                if pad_xy > 0:
                    f = ndimage.binary_dilation(f, structure=s_xy, iterations=int(pad_xy)) & testable
            flags[k] = f
    return flags


def _crr_frame_reject(sub, subv, ok, flag_g, est, vest, frame_nsig=3., fsys=0.05):
    """
    Decide which frame(s) of a group cause the outliers flagged in that group
    (used when crr_split groups several frames together, e.g. split-half).

    A frame j is rejected at a flagged voxel if
        (I_j - est_other)                                   > frame_nsig * sigma   and
        (I_j - median of the other frames of the same group) > frame_nsig * sigma,
    where the second test is skipped when the group has only one valid frame
    there. If no frame passes (e.g. in the grown wings of a cosmic ray), the
    frame with the largest significance is rejected, so every flagged voxel is
    cleaned.

    Args:
        sub, subv (ndarray [m, ...]): intensity and variance of the group's frames.
        ok (bool ndarray [m, ...]): valid voxels of those frames.
        flag_g (bool ndarray [...]): voxels flagged for this group.
        est, vest (ndarray [...]): comparison estimate from the other group(s).

    Returns:
        bool ndarray [m, ...]: rejection per frame of the group.
    """
    rej = np.zeros(sub.shape, dtype=bool)
    base = flag_g & np.isfinite(est) & np.isfinite(vest)
    if not base.any():
        return rej
    m = sub.shape[0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        # the same fractional floor as in _crr_flag_outliers: without it the
        # frame with the best seeing is the one rejected at a bright source
        vsys = (fsys * np.abs(np.nan_to_num(est))) ** 2 if fsys else 0.
        z = np.where(ok, (sub - est) / np.sqrt(subv + vest + vsys), -np.inf)
        subg = np.where(ok, sub, np.nan)
        subvg = np.where(ok, subv, np.nan)
        for j in range(m):
            oth = [k for k in range(m) if k != j]
            if len(oth) == 0:
                rej[j] = base & ok[j] & (z[j] > frame_nsig)
                continue
            n_o = np.sum(ok[oth], axis=0)
            med = np.nanmedian(subg[oth], axis=0)
            corr = np.where(n_o >= 3, np.pi / 2., 1.)
            vmed = corr * np.nanmean(subvg[oth], axis=0) / np.maximum(n_o, 1)
            vsys2 = (fsys * np.abs(np.nan_to_num(med))) ** 2 if fsys else 0.
            z2 = (sub[j] - med) / np.sqrt(subv[j] + vmed + vsys2)
            loo = (n_o == 0) | (np.nan_to_num(z2, nan=-np.inf) > frame_nsig)
            rej[j] = base & ok[j] & (z[j] > frame_nsig) & loo
        # fallback: nothing passed, so reject the most deviant valid frame
        left = base & ~rej.any(axis=0) & ok.any(axis=0)
        if left.any():
            jmax = np.argmax(z, axis=0)
            for j in range(m):
                rej[j] |= left & (jmax == j)
    return rej



def _crr_mask_fn(cubefn, crr_mask_dir='', kind='crmsk'):
    """Where the per-frame CR mask (or pre-mask) for one input cube lives.

    <dir>/<frame base>.crmsk.fits, with <dir> defaulting to the directory the
    input cube sits in.
    """
    base = os.path.basename(cubefn)
    for suf in ('_icubes.fits', '_icubed.fits', '_icube.fits', '.fits'):
        if base.endswith(suf):
            base = base[:-len(suf)]
            break
    d = crr_mask_dir if crr_mask_dir else os.path.dirname(cubefn)
    if d and not path.isdir(d):
        os.makedirs(d, exist_ok=True)
    return os.path.join(d, base + '.' + kind + '.fits')


# What moves a FRAME relative to the grid. These are the trap: kcwi_stack
# applies the shift list, the pre-shift and the astrometry correction to each
# frame's own CRPIX/CRVAL and never touches hdr0, so the grid signature is
# identical whether or not the data underneath it moved. A mask reused across
# a change of shift list would land on the wrong voxels and say nothing.
_CRR_PLACE_KEYS = ('CRMXSH', 'CRMYSH', 'CRMPRA', 'CRMPDEC', 'CRMASRA', 'CRMASDE')

_CRR_GRID_KEYS = (('CRMNX', None), ('CRMNY', None), ('CRMNZ', None),
                  ('CRMCRV1', 'CRVAL1'), ('CRMCRV2', 'CRVAL2'), ('CRMCRV3', 'CRVAL3'),
                  ('CRMCRP1', 'CRPIX1'), ('CRMCRP2', 'CRPIX2'), ('CRMCRP3', 'CRPIX3'),
                  ('CRMCD11', 'CD1_1'), ('CRMCD12', 'CD1_2'),
                  ('CRMCD21', 'CD2_1'), ('CRMCD22', 'CD2_2'), ('CRMCD33', 'CD3_3'))


def _crr_grid(hdr0, dimension):
    """The stacked grid a CR mask is tied to, key by key.

    A CR mask is a set of voxels on the COMMON grid, so it only means anything
    on the grid it was made on: change the pixel scale, the dimensions, the
    orientation or the wavelength sampling and voxel (i,j,k) no longer refers
    to the same piece of sky.

    Worth knowing when you plan to reuse a mask: kcwi_stack does NOT build the
    grid from the union of the frames, it takes the tangent point from the
    FIRST frame in the list and the rest from the .par. So a stack of a subset
    lands on the same grid as the full stack only if the first frame and the
    .par are the same. Pin the grid explicitly -- stack_ad in the .par, plus
    wave_ref / nwave / dwave -- and it stops depending on the frame list at
    all, which is what makes 'detect on all frames, combine a subset' work.
    """
    g = {}
    for kw, src in _CRR_GRID_KEYS:
        if src is None:
            g[kw] = {'CRMNX': int(dimension[0]), 'CRMNY': int(dimension[1]),
                     'CRMNZ': int(hdr0['NAXIS3'])}[kw]
        else:
            g[kw] = round(float(hdr0[src]), 10) if src in hdr0 else None
    return g


def _crr_grid_id(hdr0, dimension):
    """One hash over _crr_grid, written into the mask and checked on reuse."""
    import hashlib
    g = _crr_grid(hdr0, dimension)
    return hashlib.md5(repr([g[k] for k, _ in _CRR_GRID_KEYS]).encode()).hexdigest()[:16]


def _crr_read_mask(fnmask, gridid, shape, grid=None, place=None):
    """Load a per-frame mask, refusing it if it was not made on this grid.

    Returns (mask, '') or (None, reason). When the grid differs the reason
    names the keys that moved, because 'different grid' on its own leaves you
    guessing which knob to pin.
    """
    if not path.isfile(fnmask):
        return None, 'no file'
    try:
        hdul = fits.open(fnmask)
        hdu = next((h for h in hdul if h.data is not None), hdul[0])
        got = hdu.header.get('CRGRIDID', '')
        if got != gridid:
            diff = []
            if grid is not None:
                for kw, _ in _CRR_GRID_KEYS:
                    was, now = hdu.header.get(kw, None), grid[kw]
                    if was is None or now is None:
                        continue
                    if abs(float(was) - float(now)) > 1e-9:
                        diff.append('{} {:g}->{:g}'.format(kw, float(was), float(now)))
            return None, ('different output grid' +
                          (' [' + ', '.join(diff) + ']' if diff else
                           ' ({} != {})'.format(got, gridid)))
        if tuple(hdu.data.shape) != tuple(shape):
            return None, 'shape {} != {}'.format(tuple(hdu.data.shape), tuple(shape))
        if place is not None:
            moved = []
            for kw in _CRR_PLACE_KEYS:
                was, now = hdu.header.get(kw, None), place.get(kw, None)
                if was is None or now is None:
                    continue
                if abs(float(was) - float(now)) > 1e-6:
                    moved.append('{} {:g}->{:g}'.format(kw, float(was), float(now)))
            if moved:
                return None, 'the frame moved on the grid [' + ', '.join(moved) + ']'
            src = hdu.header.get('CRMSRC', None)
            if src is not None and place.get('CRMSRC') and src != place['CRMSRC']:
                return None, 'belongs to {}, not {}'.format(src, place['CRMSRC'])
        return hdu.data.astype(bool), ''
    except Exception as e:
        return None, 'unreadable ({})'.format(e)


def _crr_write_mask(fnmask, rej, gridid, dimension, hdr0, frames, params, place=None):
    """Write one frame's CR mask, with everything needed to validate reuse.

    Tile-compressed: the mask is one byte per voxel and almost entirely zero,
    so a 100x120x2596 grid is 31 MB raw and 0.05 MB here. Over a 34-frame
    stack that is the difference between ~1 GB of masks and ~2 MB.
    astropy reads it back like any other image.
    """
    data = rej.astype(np.uint8)
    # Tile along wavelength in slabs. The astropy default is one tile per ROW,
    # which for a 100-pixel row is almost all per-tile overhead: the same mask
    # is 10 MB that way and 0.05 MB in slabs of 64 channels. The keyword is
    # tile_shape on astropy >= 5.3 and tile_size (reversed axes) before that,
    # and an uncompressed image is still a correct mask if neither works.
    tile = (min(64, data.shape[0]), data.shape[1], data.shape[2])
    try:
        hdu = fits.CompImageHDU(data, compression_type='GZIP_1', tile_shape=tile)
    except TypeError:
        try:
            hdu = fits.CompImageHDU(data, compression_type='GZIP_1',
                                    tile_size=tuple(reversed(tile)))
        except Exception:
            hdu = fits.CompImageHDU(data, compression_type='GZIP_1')
    h = hdu.header
    h['CRGRIDID'] = (gridid, 'grid signature; valid only on this grid')
    for kw, val in _crr_grid(hdr0, dimension).items():
        if val is not None:
            h[kw] = val
    for kw, val in (place or {}).items():
        h[kw] = val
    h['CRMNREJ'] = (int(rej.sum()), 'voxels rejected in this frame')
    h['CRMNFRM'] = (len(frames), 'frames compared to make this mask')
    for j, f in enumerate(frames[:99]):
        h['CRMFR%02d' % j] = os.path.basename(f)[:60]
    for k, v in params.items():
        h[k] = v
    h['COMMENT'] = 'Per-frame CR mask from kcwi_stack(crr_split=True).'
    h['COMMENT'] = '1 = voxel rejected in this frame.'
    h['COMMENT'] = 'Reuse with crr_reuse_mask=True on the same grid.'
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(fnmask, overwrite=True)


def _iss6_scale(hdu, dtype=np.int32):
    """
    Work-around for KcwiKit Issue #6 / Montage Issue #80: the MontagePy
    build of mProjectCube reads |values| < 1 as 0, so the trimmed cubes handed
    to it are scaled to a large integer range (only used when montagepy=True). The factor is stored as ISS6FIX and
    divided out again by _iss6_factor() when the drizzled cubes are read.
    (Restored from Yuguang's 37b6dd9, lost in 95de13b.)
    """
    oflow_limit = np.iinfo(dtype).max
    peak = np.nanmax(np.abs(hdu.data))
    fac = np.float64(oflow_limit) / peak / 100. if np.isfinite(peak) and peak > 0 else 1.0
    hdu.data = np.nan_to_num(hdu.data * fac).astype(dtype)
    hdu.header['ISS6FIX'] = (fac, 'scale applied before mProjectCube; divide to undo')
    return hdu


def _iss6_factor(trimfile):
    """ISS6FIX factor of a trimmed cube, 1.0 if it was not scaled."""
    return float(fits.getheader(trimfile).get('ISS6FIX', 1.0))


def kcwi_stack(fnlist,shiftlist='',preshiftfn='',fluxfn='',pixscale_x=0.,pixscale_y=0.,
               dimension=[0,0],orientation=-1000.,cubed=False,drizzle=0,weights=[],
               wave_ref=[0, 0], dwave=0, nwave=0, wave_interp_method='cubic',
               overwrite=False,keep_trim=False,keep_mont=False,method='drizzle',use_astrom=False,
               use_regmask=True, low_mem=False, montagepy=flag_montagepy, crr=False, crr_save_files=False,
               crrthresh=100, medcube=False, nsigma_clip=10, npix_trim = 3,
               crr_split=False, crr_split_mode='auto', crr_nsig=5., crr_rel=0.5,
               crr_window=5, crr_zwindow=3, crr_zmed=15, crr_grow_xy=2, crr_grow_z=2,
               crr_grow_nsig=2., crr_pad_xy=0, crr_frame_nsig=3., crr_zblock=None,
               crr_maxmem=2.0,
               crr_zratio=51, crr_zstep=8, crr_fsys=0.05, crr_fgrad=0.5,
               crr_grow_rel=0.1, crr_pergroup=7,
               crr_mask_dir='', crr_save_mask=False, crr_reuse_mask=False,
               crr_premask=False):
    """
    Stacking the individual data cubes.

    Args:
        fnlist (str): file name of the list that contains individual files.
        shiftlist (str): set if the location of the shift file is non-default.
        preshiftfn (str): set if there is a preshift file.
        fluxfn (str): set if the normalizing file is non-default.
        pixscale_x (float): override the parfile pixel size in the x-direction
            in arcsec.
        pixscale_y (float): override the parfile pixel size in the y-direction
            in arcsec.
        dimension ([float, float]): numbers of pixels in x- and y-directions.
            (Set to override the par file.)
        orientation (float): position angle of the up direction in degrees.
            (Set to override the par file.)
        cubed (bool): using cubed?
        drizzle (float): drizzle factor.
        weights (list-like): weights of the individual cubes if using non-default.
        overwrite (bool): overwrite the cached files that was generated in
            previous runs?
        wave_ref ([float, float]): [CRPIX3, CRVAL3] of the fianl wavelength grid.
            Will override the parfile.
        dwave (float): CD3_3 of the final wavelength grid. Will override par file.
        nawave (int): NAXIS3 of the final wavelength grid. Will override par file.
        wave_interp_method (str): interpolation method for wavelength direction.
            Only applies when spatial method is 'drizzle'.
        keep_trim (bool): cache the trimmed data cubes?
        keep_mont (bool): cache the resampled data cubes?
        method (str): 'drizzle' (default), 'nearest-neighbor', 'bilinear',
            'biquadratic', or 'bicubic'.
        use_astrom (bool): set to use the list generated by 'kcwi_astrometry' to
            correct to absolute astometry.
        use_regmask (bool): set if certain pixels in individual files need to be
            masked by Region files.
        low_mem (bool): turn on low-memory mode. Siginificantly reduce memory
            usage but increase computational time. Useful when working with the
            small slicer.
        montagepy (bool): use MontagePy for drizzling? Otherwise, use the command
            line Montage installation. Both require proper installation.
        #------------------------- Legacy CRR parameters (deprecated) -------------------------
        crr (bool): perform final CRR when stacking based on flux outliers?
        crr_save_files (bool): save trimmed cubes with CRs masked?
        crrthresh (float): Default value = 100. Sets the threshold level above
            which pixels are flagged as CRs. Not super reliable yet - better to
            flag them in 2D images at the beginning of KCWI_DRP.
        nsigma_clip (float): only used for red cameras for sigma clipping to
            remove residual cosmic rays. Default=1.5.
        npix_trim (int): number of pixels to trim from the edges of a cube. Default = 3.
            May want npix_trim = 1 for Large slicer.
        #-------------------------- THE ACTUAL CRR PARAMETERS --------------------------
        crr_split (bool): cosmic-ray rejection after the frames have been shifted
            and projected onto the common grid. Each frame (or each group of
            frames, see crr_split_mode) is compared with the median of the
            others. A voxel is flagged only if its excess, after a running
            median along wavelength is removed, is both significant compared
            with the noise (crr_nsig) and large compared with the local source
            level (crr_rel); see _crr_flag_outliers. Only positive outliers are
            rejected. The flagged voxels are masked in that frame alone before
            the final stack. Can be combined with crr (which runs first).
        crr_split_mode (str, int or array): how the frames are compared.
            An integer N splits the frames into N groups of consecutive frames
            (34 into 5 gives 7,7,7,7,6); 'auto' (default) picks N so that each
            group holds about crr_pergroup frames. The run time is linear in the
            number of groups, so 'frame' on a 34-frame stack costs about seven
            times what five groups do; the frame responsible for each flagged
            voxel is still identified individually inside its group by
            _crr_frame_reject, so grouping costs accuracy only where two
            cosmic rays from different frames of the same group land on the
            same voxel.
            'frame': every frame is compared with the median of all
                the other frames. The most sensitive option; needs >= 3 frames
                to be robust (>= 4-5 is comfortable).
            'half': first half of the list against the second half, each
                co-added with the stacking weights. Less sensitive (the cosmic
                ray is diluted by the frames of its own half) but usable with
                few frames or noisy individual frames; the responsible frame(s)
                of the flagged half are then identified (see _crr_frame_reject).
            'alternate': even frames against odd frames.
            array: one group label per frame, for custom groupings.
        crr_nsig (float): absolute threshold on the excess [sigma]. In 'frame'
            mode the comparison is per frame, so raise it (6-7) if the noise
            of the individual frames is underestimated.
        crr_rel (float): relative threshold, excess > crr_rel * local source level.
        crr_window (int): spatial size (pixels) of the local-level box, about the
            seeing FWHM or larger.
        crr_zwindow (int): spectral size (pixels) of the local-level box.
        crr_zmed (int): width (pixels) of the running median along wavelength
            subtracted from the excess. 0 disables it.
        crr_grow_xy (int): maximum spatial growth (pixels, 8-connected) into
            connected neighbors with excess > crr_grow_nsig * sigma.
        crr_grow_z (int): maximum spectral growth (pixels), same rule.
        crr_grow_nsig (float): threshold for growing.
        crr_pad_xy (int): unconditional spatial margin (pixels). 0 = off.
        crr_mask_dir (str): directory holding the per-frame CR masks written by
            crr_save_mask and read back by crr_reuse_mask / crr_premask. The
            default ('') puts them next to the input cubes. A mask is named
            <frame base>.crmsk.fits and a pre-mask <frame base>.premsk.fits.
        crr_save_mask (bool): write each frame's rejected voxels to its own
            file after the comparison. The mask records the signature of the
            grid it was made on and the frames that took part, so a later run
            can tell whether it still applies.
        crr_reuse_mask (bool): load the saved masks instead of comparing the
            frames again. This is all-or-nothing on purpose: if any frame has
            no valid mask the comparison is redone for every frame, because
            the rejection in one frame depends on which other frames it was
            compared with, so a stack with a different frame list deserves
            fresh statistics rather than a silent mixture.

            Its real use is deliberately decoupling the frames you detect
            cosmic rays with from the frames you combine. Run once over every
            frame of a field with crr_save_mask=True, then stack whatever
            subset you want with crr_reuse_mask=True: the subset inherits the
            rejection quality of the full set. On a field where saturation or
            bad seeing leaves only two or three usable frames, that is the
            difference between a median of two -- which cannot tell which of
            the two is the outlier -- and a median of everything you took.
        crr_premask (bool): before the statistics, read <base>.premsk.fits for
            each frame and treat those voxels as bad. They are left out of the
            median, so a known-bad value cannot drag the estimate and make the
            good frames look like outliers, and they are also rejected in the
            stack. This is where a saturation mask belongs: a clipped QSO core
            sits BELOW the true flux, so leaving it in the comparison biases
            the median down exactly where the source is brightest. Pre-masks
            are applied fresh on every run and are never folded into the saved
            .crmsk.fits, so editing one takes effect without re-running the
            comparison.
        crr_pergroup (int): target frames per group for crr_split_mode='auto'.
        crr_zstep (int): channels between evaluations of the per-spaxel scale,
            which is smooth along wavelength; 1 evaluates it at every channel
            (slow). 8 is about eight times faster and agrees to ~1%, which
            then largely cancels in the ratio.
        crr_zratio (int): channels over which each frame is scaled to the
            comparison estimate per spaxel before differencing, so that
            seeing and transparency differences do not leave cosmic-ray-shaped
            residuals at the narrow features of a bright source. 0 = off.
        crr_fsys (float): fractional frame-to-frame systematic added to the
            noise. It has no effect while crr_nsig * crr_fsys <= crr_rel
            (at the defaults, 0.25 against 0.5), because the relative test is
            already a stronger fractional floor -- it is here for when you
            lower crr_rel to catch cosmic rays on a bright source. See
            _crr_flag_outliers.
        crr_fgrad (float): the same for the spatial gradient, i.e. the error a
            sub-pixel registration offset produces; this one does bind.
        crr_grow_rel (float): relative threshold while growing, so a cosmic ray
            beside a bright source cannot grow into it. Default 0.1 (the same
            default is used by _crr_flag_outliers and kcwi_combine_stacks).
        crr_frame_nsig (float): threshold for picking the responsible frame(s)
            at a flagged voxel. Only used when a group holds several frames.
        crr_zblock (int): number of wavelength channels processed at a time
            when flagging, or None (the default) to size it automatically from
            crr_maxmem. Each block also reads `pad` channels on either side so
            the wavelength filters see a full window, and pad is now
            max(crr_zwindow, crr_zmed, crr_zratio)//2 + crr_grow_z + 1, i.e. 28
            with the defaults. A small block therefore spends most of its time
            on overlap: at 100 channels, 56% of the work is redundant, at 300
            it is 19%. Memory scales as
                (n_frames * 13 + n_groups * 9) * nx * ny * (crr_zblock + 2*pad)
            bytes, on top of a fixed n_frames * nz * nx * ny bits for the
            rejection cube itself.
        crr_maxmem (float): target working set per block [GB] when crr_zblock
            is None. Raising it makes the blocks bigger, which is strictly
            better for accuracy -- a single block over the whole cube is the
            exact answer and every extra boundary is a place where the
            wavelength filters see a window assembled from a margin rather
            than real data. The only costs of a big block are memory and a
            coarser progress bar. The chosen value, the number of blocks, the
            extra work from the margins and the peak memory are printed.

    Returns:
        None
        (Stacked cube with be generated as FITS files.)

    """



    if cubed:
        suffix="cubed"
    else:
        suffix="cubes"

    if medcube:
        suffix='cube'

    if method.lower()!='drizzle':
        if method.lower()=='nearest-neighbor':
            method_flag='nei'
        elif method.lower()=='bilinear':
            method_flag='lin'
        elif method.lower()=='biquadratic':
            method_flag='qua'
        elif method.lower()=='bicubic':
            method_flag='cub'
        else:
            print('Error: Method not found.')
            return 0

    parfn=fnlist.replace(".list",".par")
    par=kcwi_stack_readpar(parfn)

    if shiftlist=="":
        shiftlist=fnlist.replace(".list",".shift.list")

    if pixscale_x==0:
        pixscale_x=par["stack_xpix"]
        if pixscale_x==-1:
            pixscale_x=0.3
    pixscale_x=pixscale_x/3600.

    if pixscale_y==0:
        pixscale_y=par["stack_ypix"]
        if pixscale_y==-1:
            pixscale_y=0.3
    pixscale_y=pixscale_y/3600.


    if dimension[0]==0:
        dimension=par["stack_dimension"]
        if dimension[0]==-1:
            dimension=[100,100]

    if drizzle==0:
        drizzle=par["drizzle"]
        if drizzle==0:
            drizzle=0.7

    if wave_ref[1]==0:
        wave_ref = par['wave_ref']

    if nwave==0:
        nwave = par['nwave']

    if dwave==0:
        dwave = par['dwave']


    # make tmp directory
    if not os.path.exists('kcwi_stack'):
        os.makedirs('kcwi_stack')

    fnhdr='kcwi_stack/'+fnlist.replace('.list','.hdr')

    trim_tab=ascii.read(fnlist,format="no_header")
    fn=trim_tab["col1"]
    trim=np.array([trim_tab["col2"],trim_tab["col3"]])

    shift_tab=ascii.read(shiftlist,format="no_header")
    xshift=shift_tab["col2"]
    yshift=shift_tab["col3"]

    if path.isfile(fn[0]+"_i"+suffix+".fits") == False:
        suffix="cubed"


    vfn=[i+'_v'+suffix+'.fits' for i in fn]
    mfn=[i+'_m'+suffix+'.fits' for i in fn]
    efn=[i+'_e'+suffix+'.fits' for i in fn]
    fn=[i+'_i'+suffix+'.fits' for i in fn]


    if preshiftfn=='':
        preshiftfn=fnlist.replace('.list','.preshift.list')
        if path.isfile(preshiftfn)==False:
            preshiftfn=fnlist.replace('.list','.pre.list')
            if path.isfile(preshiftfn)==False:
                preshiftfn=''
    if preshiftfn!='':
        pre_tab=ascii.read(preshiftfn,format='no_header')
        prefn=[i+'_i'+suffix+'.fits' for i in pre_tab['col1']]
        prera=pre_tab['col2']
        predec=pre_tab['col3']

    # Absolute astrometry offsets, one per frame (restored from Yuguang's
    # df92496/7d43e0d, lost in 95de13b). The .astrom.list file is either
    #   2 columns: dRA dDec            -> one offset for every frame, or
    #   3 columns: frame dRA dDec      -> per-frame offsets (mosaics).
    astrom_rashift = np.zeros(len(fn))
    astrom_decshift = np.zeros(len(fn))
    if use_astrom:
        astromfn = fnlist.replace('.list','.astrom.list')
        astrom_tab=ascii.read(astromfn)
        if len(astrom_tab.colnames)==2:
            astrom_rashift[:] = astrom_tab['col1'][0]
            astrom_decshift[:] = astrom_tab['col2'][0]
        else:
            astrom_basenames = np.array([os.path.basename(str(row['col1'])) for row in astrom_tab])
            for i in range(len(fn)):
                base = os.path.basename(fn[i])
                index = (astrom_basenames == base) | (astrom_basenames + '_i' + suffix + '.fits' == base)
                if np.sum(index) != 1:
                    raise ValueError('%s: expected exactly one entry for %s, found %d'
                                     % (astromfn, base, int(np.sum(index))))
                astrom_rashift[i] = astrom_tab['col2'][index][0]
                astrom_decshift[i] = astrom_tab['col3'][index][0]
        overwrite=True

    # flux weight
    if fluxfn=='':
        fluxfn = fnlist.replace('.list', '.flx.list')
    if os.path.isfile(fluxfn):
        tmp = ascii.read(fluxfn)
        fluxnorm = tmp['col2']
    else:
        fluxnorm = np.ones(len(fn))


    # construct wcs
    hdulist=fits.open(fn[0])
    hdrtmp=hdulist[0].header.copy()
    hdulist.close()
    wcstmp=wcs.WCS(hdrtmp).copy()
    center=wcstmp.wcs_pix2world((wcstmp.pixel_shape[0]-1)/2.,(wcstmp.pixel_shape[1]-1)/2.,0,0,ra_dec_order=True)

    # camera
    camera = hdrtmp['CAMERA']
    if camera != 'RED' and camera != 'BLUE':
        raise ValueError('Unknown camera type - {}'.format(camera))
        return

    if par['stack_ad'][0]!=-1:
        center=par['stack_ad']

    hdr0=hdrtmp.copy()
    hdr0['NAXIS1']=dimension[0]
    hdr0['NAXIS2']=dimension[1]
    hdr0['CRPIX1']=(dimension[0]+1)/2.
    hdr0['CRPIX2']=(dimension[1]+1)/2.
    hdr0['CRVAL1']=float(center[0])
    hdr0['CRVAL2']=float(center[1])
    old_cd11=hdr0['CD1_1']
    old_cd12=hdr0['CD1_2']
    old_cd21=hdr0['CD2_1']
    old_cd22=hdr0['CD2_2']
    hdr0['CD1_1']=-pixscale_x
    hdr0['CD2_2']=pixscale_y
    hdr0['CD1_2']=0
    hdr0['CD2_1']=0
    #hdr0['CTYPE3']='WAVE'
    #hdr0['BUNIT']='10^(-16)erg/s/cm2/Angstrom'
    hdr0['BUNIT']='1e-16 erg/s/cm2/Angstrom/arcsec2'
    if suffix!='cubes':
        #hdr0['BUNIT']='adu/s'
        hdr0['BUNIT']='count/s/arcsec2'

    # orientation
    if orientation==-1000:
        orientation=par['stack_orientation']
        if orientation==-1000:
            orientation=np.rad2deg(np.arctan(old_cd21/(-old_cd11)))
    hdr0['CD1_1']=-pixscale_x*np.cos(np.deg2rad(orientation))
    hdr0['CD2_1']=pixscale_x*np.sin(np.deg2rad(orientation))
    hdr0['CD1_2']=pixscale_y*np.sin(np.deg2rad(orientation))
    hdr0['CD2_2']=pixscale_y*np.cos(np.deg2rad(orientation))

    # wavelength
    if wave_ref[1]!=0:
        hdr0['CRPIX3'] = wave_ref[0]
        hdr0['CRVAL3'] = wave_ref[1]
    if nwave!=0:
        hdr0['NAXIS3'] = int(nwave)
    if dwave!=0:
        hdr0['CD3_3'] = dwave

    hdr0.totextfile(fnhdr,overwrite=1)

    # project
    #void=mProjectCube(fn[0],outfn[0],'kcwi_stack/tmp.hdr',drizzle=0.7,energyMode=True)
    start=ostime.time()
    print('Projecting...')

    # preprocessing
    trimfn=['kcwi_stack/'+os.path.basename(i).replace('.fits','.trim.fits') for i in fn]
    trimvfn=['kcwi_stack/'+os.path.basename(i).replace('.fits','.trim.fits') for i in vfn]
    trimmfn=['kcwi_stack/'+os.path.basename(i).replace('.fits','.trim.fits') for i in mfn]
    trimefn=['kcwi_stack/'+os.path.basename(i).replace('.fits','.trim.fits') for i in efn]
    montfns=[]
    montvfns=[]
    montmfns=[]
    montefns=[]
    etime=np.zeros(len(fn))
    airmass_list = []
    expt_list = []
    for i in range(len(fn)):
        print(os.path.basename(fn[i]))

        # check availability
        if (not os.path.isfile(trimfn[i])) or overwrite==True:
            # science cube
            hdulist=fits.open(fn[i])
            hdr=hdulist[0].header.copy()
            print("     CRVAL (original) = "+str(hdr['CRVAL3']))
            # IDL or Python?

            if len(hdulist) == 1:
                reduxflag = 'idl'
            elif len(hdulist) >= 4:
                reduxflag = 'py'
            else:
                raise ValueError('Reduction pipeline not recognized.')
                return

            if medcube == True:
                reduxflag = 'medcube'

            if reduxflag == 'idl':
                print("idl reduced data")
                hdu_i,vcorr=kcwi_vachelio(hdulist[0],hdr_ref=hdr0)
                print('     Vcorr = '+str(vcorr))
                hdu_i.header['VCORR'] = (vcorr, 'Heliocentric Velocity Correction')
                hdulist.close()

                # variance cube -> sigma cube
                hdulist=fits.open(vfn[i])
                hdu_v,vcorr=kcwi_vachelio(hdulist[0],hdr_ref=hdr0)
                hdulist.close()

                # mask cube
                hdulist=fits.open(mfn[i])
                hdu_m,vcorr=kcwi_vachelio(hdulist[0],hdr_ref=hdr0,mask=True)
                hdulist.close()

            elif reduxflag == 'py':

                hdulist, vcorr = kcwi_vachelio(hdulist, hdr_ref=hdr0)
                print('     Vcorr = '+str(vcorr))
                hdu_i = hdulist[0]
                hdu_i.header['VCORR'] = (vcorr, 'Radial Velocity Correction')

                hdu_v = fits.PrimaryHDU(hdulist['UNCERT'].data**2, hdu_i.header)
                hdu_m = fits.PrimaryHDU(hdulist['FLAGS'].data, hdu_i.header)
                #hdu_m.writeto('testfile/'+os.path.basename(fn[i]).replace('.fits','_vachelio_FLAGS.fits'),overwrite=True)
                #hdu_i.writeto('testfile/'+os.path.basename(fn[i]).replace('.fits','_vachelio_icubes.fits'),overwrite=True)
            elif reduxflag == 'medcube':
                #reduxflag == 'py'
                hdulist, vcorr = kcwi_vachelio(hdulist, hdr_ref=hdr0)
                print('     Vcorr = '+str(vcorr))
                hdu_i = fits.open(fn[i].replace(suffix,'cube.med'))[0] #hdulist[0]

                hdu_v = fits.PrimaryHDU(hdulist['UNCERT'].data**2, hdu_i.header)
                hdu_m = fits.PrimaryHDU(hdulist['FLAGS'].data, hdu_i.header)

            else:
                raise ValueError('reduxflag not assigned')
                return

            # region masks
            regfn = fn[i].replace('.fits','.thum.reg')
            if os.path.isfile(regfn) and use_regmask==True:
                hdr2d=hdu_i.header.copy()
                del hdr2d['CD3_3']
                del hdr2d['CRVAL3']
                del hdr2d['CRPIX3']
                del hdr2d['NAXIS3']
                del hdr2d['CTYPE3']
                del hdr2d['CNAME3']
                del hdr2d['CUNIT3']
                hdr2d['NAXIS']=2

                mask_reg = regions_to_mask(regfn, header=hdr2d, shape=hdu_i.data.shape[1:])

                hdu_i.data[:,mask_reg] = np.nan
                hdu_v.data[:,mask_reg] = np.nan
                hdu_m.data[:,mask_reg] = 128


            # Infinity check
            if medcube == False:
                q=((hdu_i.data==0) | (~np.isfinite(hdu_i.data)) | (hdu_v.data==0) | (~np.isfinite(hdu_v.data)) )
                hdu_i.data[q]=np.nan
                hdu_v.data[q]=np.nan
                hdu_m.data[q]=128


            # check EXPTIME
            hdu_i=kcwi_checkexptime(hdu_i)
            exptime=hdu_i.header['XPOSURE']
            print('     EXPTIME = '+str(exptime))
            etime[i]=exptime
            edata=hdu_i.data*0.+exptime
            q=(hdu_m.data > 0.5) # set to 0.5 to avoid interpolation artifacts, where mask could be a small finite value
            edata[q]=0
            if medcube == False:
                hdu_i.data[q] = np.nan
            hdu_e=fits.PrimaryHDU(edata,header=hdu_i.header)
            hdu_e.header['BUNIT']='s'
            #hdu_e.writeto('testfile/'+os.path.basename(fn[i]).replace('.fits','_vachelio_ecube.fits'),overwrite=True)

            # have to use surface brightness for now, mProjectCube has bug with brightness units combined with drizzle scale
            dx=np.sqrt(hdu_i.header['CD1_1']**2+hdu_i.header['CD2_1']**2)*3600.
            dy=np.sqrt(hdu_i.header['CD1_2']**2+hdu_i.header['CD2_2']**2)*3600.
            area=dx*dy
            if suffix!='cubes':
                hdu_i.data=hdu_i.data/exptime/area
                hdu_v.data=hdu_v.data/exptime**2/area**2
                hdu_i.header['BUNIT']='count/s/arcsec2'
                hdu_v.header['BUNIT']='count2/s2/arcsec4'
            else:
                hdu_i.data=hdu_i.data/area
                hdu_v.data=hdu_v.data/area**2
                hdu_i.header['BUNIT']='1e-16 erg/s/cm2/Angstrom/arcsec2'
                hdu_v.header['BUNIT']='1e-32 erg2/s2/cm4/Angstrom2/arcsec4'


            # preshift
            if preshiftfn!='':
                index=np.where(np.array(prefn)==os.path.basename(fn[i]))
                index=index[0]
                if len(index)>0:
                    index=index[0]
                    hdu_i.header['CRVAL1']=hdu_i.header['CRVAL1']+prera[index]/3600.
                    hdu_i.header['CRVAL2']=hdu_i.header['CRVAL2']+predec[index]/3600.

                    hdu_v.header['CRVAL1']=hdu_v.header['CRVAL1']+prera[index]/3600.
                    hdu_v.header['CRVAL2']=hdu_v.header['CRVAL2']+predec[index]/3600.

                    hdu_m.header['CRVAL1']=hdu_m.header['CRVAL1']+prera[index]/3600.
                    hdu_m.header['CRVAL2']=hdu_m.header['CRVAL2']+predec[index]/3600.

                    hdu_e.header['CRVAL1']=hdu_e.header['CRVAL1']+prera[index]/3600.
                    hdu_e.header['CRVAL2']=hdu_e.header['CRVAL2']+predec[index]/3600.

            # astrometry correction
            if use_astrom:
                hdu_i.header['CRVAL1']=hdu_i.header['CRVAL1']+astrom_rashift[i]/3600.
                hdu_i.header['CRVAL2']=hdu_i.header['CRVAL2']+astrom_decshift[i]/3600.

                hdu_v.header['CRVAL1']=hdu_v.header['CRVAL1']+astrom_rashift[i]/3600.
                hdu_v.header['CRVAL2']=hdu_v.header['CRVAL2']+astrom_decshift[i]/3600.

                hdu_m.header['CRVAL1']=hdu_m.header['CRVAL1']+astrom_rashift[i]/3600.
                hdu_m.header['CRVAL2']=hdu_m.header['CRVAL2']+astrom_decshift[i]/3600.

                hdu_e.header['CRVAL1']=hdu_e.header['CRVAL1']+astrom_rashift[i]/3600.
                hdu_e.header['CRVAL2']=hdu_e.header['CRVAL2']+astrom_decshift[i]/3600.

            # shift
            hdu_i.header['CRPIX1']=hdu_i.header['CRPIX1']+xshift[i]
            hdu_i.header['CRPIX2']=hdu_i.header['CRPIX2']+yshift[i]

            hdu_v.header['CRPIX1']=hdu_v.header['CRPIX1']+xshift[i]
            hdu_v.header['CRPIX2']=hdu_v.header['CRPIX2']+yshift[i]

            hdu_m.header['CRPIX1']=hdu_m.header['CRPIX1']+xshift[i]
            hdu_m.header['CRPIX2']=hdu_m.header['CRPIX2']+yshift[i]

            hdu_e.header['CRPIX1']=hdu_e.header['CRPIX1']+xshift[i]
            hdu_e.header['CRPIX2']=hdu_e.header['CRPIX2']+yshift[i]


            # trim
            for kk in range(hdu_i.header['NAXIS3']):
                img=hdu_i.data[kk,:,:]
                var=hdu_v.data[kk,:,:]
                mask=hdu_m.data[kk,:,:]
                expimg=hdu_e.data[kk,:,:]

                # remove PyDRP edge problem
                if reduxflag=='py':
                    n_pix = npix_trim #3 by default
                    flag_dim = mask.shape
                    mask[0:n_pix,:] = 128
                    mask[:,0:n_pix] = 128
                    mask[flag_dim[0]-n_pix:flag_dim[0],:] = 128
                    mask[:,flag_dim[1]-n_pix:flag_dim[1]] = 128

                index_y,index_x=np.where(mask<0.5) # set to 0.5 to avoid interpolation artifact after wavelength correction
                if len(index_y)==0:
                    continue
                xrange=[index_x.min(),index_x.max()]
                yrange=[index_y.min(),index_y.max()]
                if yrange[0]+trim[0,i] >= yrange[1]-trim[1,i]:
                    continue

                img[yrange[1]-trim[1,i]+1:,:]=np.nan
                img[:yrange[0]+trim[0,i],:]=np.nan

                var[yrange[1]-trim[1,i]+1:,:]=np.nan
                var[:yrange[0]+trim[0,i],:]=np.nan

                mask[yrange[1]-trim[1,i]+1:,:]=128
                mask[:yrange[0]+trim[0,i],:]=128

                expimg[yrange[1]-trim[1,i]+1:,:]=0
                expimg[:yrange[0]+trim[0,i],:]=0

                # Additional trimming due to PyDRP cubic interp
                if reduxflag=='py':
                    img[:,xrange[1]+1:]=np.nan # need to deal with vertical stripes along side
                    img[:,:xrange[0]]=np.nan

                    var[:,xrange[1]+1:]=np.nan
                    var[:,:xrange[0]]=np.nan

                    mask[:,xrange[1]+1:]=128
                    mask[:,:xrange[0]]=128

                    expimg[:,xrange[1]+1:]=0
                    expimg[:,:xrange[0]]=0

                hdu_i.data[kk,:,:]=img
                hdu_v.data[kk,:,:]=var
                hdu_m.data[kk,:,:]=mask
                hdu_e.data[kk,:,:]=expimg

            # flux correction
            hdu_i.data = hdu_i.data * fluxnorm[i]
            hdu_v.data = hdu_v.data * fluxnorm[i]**2

            # wavelength interpolation
            if method=='drizzle':
                # "reproject" package can handle wavelength interpolation, only
                # Montage need this step.
                if not kcwi_check_samewave(hdu_i.header, hdr0):
                    print('  Different wavelength grid: interpolating...')

                    hdu_i = kcwi_resample_wave(hdu_i, hdr0, method=wave_interp_method,plot=False)
                    hdu_v = kcwi_resample_wave(hdu_v, hdr0, method=wave_interp_method)
                    hdu_m = kcwi_resample_wave(hdu_m, hdr0, method='mask')
                    hdu_e = kcwi_resample_wave(hdu_e, hdr0, method='linear')

                # Issue #6 hot fix, MontagePy only: the MontagePy wheel reads
                # |values| < 1 as 0 (Montage issue #80), so integer-scale the
                # cubes handed to it. The command-line mProjectCube does not
                # need it, and the int32 quantisation would cost precision on
                # faint pixels next to a bright source, so it is skipped there.
                if montagepy:
                    hdu_i = _iss6_scale(hdu_i)
                    hdu_v = _iss6_scale(hdu_v)
                    hdu_m = _iss6_scale(hdu_m)
                    hdu_e = _iss6_scale(hdu_e)

            # write
            hdu_i.writeto(trimfn[i],overwrite=True)
            hdu_v.writeto(trimvfn[i],overwrite=True)
            hdu_m.writeto(trimmfn[i],overwrite=True)
            hdu_e.writeto(trimefn[i],overwrite=True)

        # airmass / exposure-time bookkeeping, read from the trimmed cube so it
        # is also filled when trimmed cubes from a previous run are reused
        hdr_trim = fits.getheader(trimfn[i])
        airmass_list.append(hdr_trim['AIRMASS'])
        expt_list.append(hdr_trim['XPOSURE'])

        # Montage
        if method.lower()=='drizzle':
            montfn=trimfn[i].replace('.trim.fits','.mont.fits')
            montvfn=trimvfn[i].replace('.trim.fits','.mont.fits')
            montmfn=trimmfn[i].replace('.trim.fits','.mont.fits')
            montefn=trimefn[i].replace('.trim.fits','.mont.fits')
        else:
            montfn=trimfn[i].replace('.trim.fits','.'+method_flag+'.fits')
            montvfn=trimvfn[i].replace('.trim.fits','.'+method_flag+'.fits')
            montmfn=trimmfn[i].replace('.trim.fits','.'+method_flag+'.fits')
            montefn=trimefn[i].replace('.trim.fits','.'+method_flag+'.fits')

        montfns.append(montfn)
        montvfns.append(montvfn)
        montmfns.append(montmfn)
        montefns.append(montefn)


        if (not os.path.isfile(montfn)) or overwrite==True:
        #if True:
            if method.lower()!='drizzle':
                newhdr=fits.Header.fromtextfile(fnhdr)

                hdut=(fits.open(trimfn[i]))[0]
                hdut.data[np.isfinite(hdut.data)==False]=0.
                newi,newa=reproject_interp(hdut,newhdr,order=method)
                hdui=fits.PrimaryHDU(newi,newhdr)
                hdua=fits.PrimaryHDU(newa,newhdr)
                hdui.writeto(montfn,overwrite=True)
                hdua.writeto(montfn.replace('.'+method_flag,'.'+method_flag+'_area'),overwrite=True)

                hdutv=(fits.open(trimvfn[i]))[0]
                hdutv.data[np.isfinite(hdutv.data)==False]=0.
                newv,newa=reproject_interp(hdutv,newhdr,order=method)
                hduv=fits.PrimaryHDU(newv,newhdr)
                hdua=fits.PrimaryHDU(newa,newhdr)
                hduv.writeto(montvfn,overwrite=True)
                hdua.writeto(montvfn.replace('.'+method_flag,'.'+method_flag+'_area'),overwrite=True)

                hdutm=(fits.open(trimmfn[i]))[0]
                hdutm.data[np.isfinite(hdutm.data)==False]=0.
                newm,newa=reproject_interp(hdutm,newhdr,order='bilinear')
                hdum=fits.PrimaryHDU(newm,newhdr)
                hdua=fits.PrimaryHDU(newa,newhdr)
                hdum.writeto(montmfn,overwrite=True)
                hdua.writeto(montmfn.replace('.'+method_flag,'.'+method_flag+'_area'),overwrite=True)

                hdute=(fits.open(trimefn[i]))[0]
                hdute.data[np.isfinite(hdute.data)==False]=0.
                newe,newa=reproject_interp(hdute,newhdr,order='bilinear')
                hdue=fits.PrimaryHDU(newe,newhdr)
                hdua=fits.PrimaryHDU(newa,newhdr)
                hdue.writeto(montefn,overwrite=True)
                hdua.writeto(montefn.replace('.'+method_flag,'.'+method_flag+'_area'),overwrite=True)

            else:
                if not montagepy:
                    # Command line version
                    exe="mProjectCube -z "+str(drizzle)+" -f "+trimfn[i]+" "+montfn+" "+fnhdr
                    void=os.system(exe)
                    exev="mProjectCube -z "+str(drizzle)+" -f  "+trimvfn[i]+" "+montvfn+" "+fnhdr
                    voidv=os.system(exev)
                    exem="mProjectCube -z "+str(drizzle)+" -f  "+trimmfn[i]+" "+montmfn+" "+fnhdr
                    voidm=os.system(exem)
                    exee="mProjectCube -z "+str(drizzle)+" -f  "+trimefn[i]+" "+montefn+" "+fnhdr
                    voide=os.system(exee)

                else:
                    # MontagePy
                    void=mProjectCube(trimfn[i],montfn,fnhdr,drizzle=drizzle,energyMode=False,fullRegion=True)
                    voidv=mProjectCube(trimvfn[i],montvfn,fnhdr,drizzle=drizzle,energyMode=False,fullRegion=True)
                    voidm=mProjectCube(trimmfn[i],montmfn,fnhdr,drizzle=drizzle,energyMode=False,fullRegion=True)
                    voide=mProjectCube(trimefn[i],montefn,fnhdr,drizzle=drizzle,energyMode=False,fullRegion=True)


    if low_mem==False:
        # cache all cubes
        data0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.float64).T
        vdata0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.float64).T
        mdata0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.int16).T+128
        edata0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.float64).T
        for i in range(len(fn)):
            newcube=fits.open(montfns[i])[0].data / _iss6_factor(trimfn[i])
            newcube[~np.isfinite(newcube)]=0.
            newcubev=fits.open(montvfns[i])[0].data / _iss6_factor(trimvfn[i])
            newcubev[~np.isfinite(newcubev)]=0.
            newcubem=np.ceil(fits.open(montmfns[i])[0].data / _iss6_factor(trimmfn[i]))
            newcubem[~np.isfinite(newcubem)]=128
            newcubee=fits.open(montefns[i])[0].data / _iss6_factor(trimefn[i])
            newcubee[~np.isfinite(newcubee)]=0.
            data0[i,:,:,:]=newcube
            vdata0[i,:,:,:]=newcubev
            mdata0[i,:,:,:]=newcubem
            edata0[i,:,:,:]=newcubee
            #data0.append(newcube)
            #vdata0.append(newcubev)
            #mdata0.append(newcubem)
            #edata0.append(newcubee)
            #print(newcube.shape)
            #print(newcubev.shape)
            #print(newcubem.shape)
            #print(newcubee.shape)

    # Legacy CR Rejection, if needed
    if crr:
        crrfn=trimfn[i].replace('.trim.fits','.crr.fits')
        crrvfn=trimvfn[i].replace('.trim.fits','.crr.fits')
        crrmfn=trimmfn[i].replace('.trim.fits','.crr.fits')
        crrefn=trimefn[i].replace('.trim.fits','.crr.fits')

        fits.PrimaryHDU(data0).writeto(crrfn, overwrite=True)
        fits.PrimaryHDU(vdata0).writeto(crrvfn, overwrite=True)
        fits.PrimaryHDU(mdata0).writeto(crrmfn, overwrite=True)
        fits.PrimaryHDU(edata0).writeto(crrefn, overwrite=True)

        # running in this manner to untilize multiprocessing
        os.system(f'python {pathlib.Path(__file__).parent.resolve()}/kcwi_crr.py {crrfn} {crrvfn} {crrmfn} {crrefn} {crrthresh}')

        data0 = fits.open(crrfn)[0].data
        vdata0 = fits.open(crrvfn)[0].data
        mdata0 = fits.open(crrmfn)[0].data
        edata0 = fits.open(crrefn)[0].data

        if crr_save_files==False:
            os.remove(crrfn)
            os.remove(crrvfn)
            os.remove(crrmfn)
            os.remove(crrefn)



    # ------------------------------------------------------------------
    # crr_split: cosmic-ray rejection by comparing frames (or groups of
    # frames) with each other on the common grid, after shifting/projection.
    # ------------------------------------------------------------------
    if crr_split:
        nfn = len(fn)
        if isinstance(crr_split_mode, str):
            if crr_split_mode == 'frame':
                crr_groups = np.arange(nfn)                     # each frame on its own
            elif crr_split_mode == 'half':
                crr_groups = (np.arange(nfn) >= nfn // 2).astype(int)
            elif crr_split_mode == 'alternate':
                crr_groups = (np.arange(nfn) % 2).astype(int)
            elif crr_split_mode == 'auto':
                # keep the run time sane on a long stack: the cost is linear in
                # the number of groups, so aim at crr_pergroup frames each
                ng = max(2, min(nfn, int(round(nfn / float(crr_pergroup)))))
                crr_groups = (np.arange(nfn) * ng) // nfn
            else:
                raise ValueError("crr_split_mode must be 'frame', 'half', 'alternate', "
                                 "'auto', an integer number of groups, or an array")
        elif np.isscalar(crr_split_mode):
            # an integer: that many groups of consecutive frames, as equal as
            # they divide. 34 frames into 5 gives 7,7,7,7,6.
            ng = int(crr_split_mode)
            if ng < 2 or ng > nfn:
                raise ValueError('crr_split_mode as a number of groups must be '
                                 'between 2 and the number of frames ({})'.format(nfn))
            crr_groups = (np.arange(nfn) * ng) // nfn
        else:
            crr_groups = np.asarray(crr_split_mode).astype(int)
            if len(crr_groups) != nfn:
                raise ValueError('crr_split_mode array must have one label per frame')
            crr_groups = np.unique(crr_groups, return_inverse=True)[1]
        crr_ngroup = int(crr_groups.max()) + 1
        if nfn < 2 or crr_ngroup < 2:
            print('crr_split: need at least two groups; skipping.')
            crr_split = False

    if crr_split:
        crr_members = [np.where(crr_groups == g)[0] for g in range(crr_ngroup)]
        print('CR rejection (crr_split, mode={}): {} group(s): {}'.format(
            crr_split_mode if isinstance(crr_split_mode, str) else 'custom',
            crr_ngroup, [list(m) for m in crr_members]))
        nz3 = hdr0['NAXIS3']
        ny3 = dimension[1]
        nx3 = dimension[0]
        crr_rej = np.zeros((nfn, nz3, ny3, nx3), dtype=bool)

        # ---- per-frame masks: pre-masks in, CR masks out -----------------
        crr_gridid = _crr_grid_id(hdr0, dimension)
        crr_gridkeys = _crr_grid(hdr0, dimension)
        crr_shape = (nz3, ny3, nx3)

        def _crr_place(i):
            """Everything that moves frame i relative to the common grid."""
            pr = pd = 0.
            if preshiftfn != '':
                ix = np.where(np.array(prefn) == os.path.basename(fn[i]))[0]
                if len(ix):
                    pr, pd = float(prera[ix[0]]), float(predec[ix[0]])
            return {'CRMXSH': float(xshift[i]), 'CRMYSH': float(yshift[i]),
                    'CRMPRA': pr, 'CRMPDEC': pd,
                    'CRMASRA': float(np.mean(astrom_rashift)),
                    'CRMASDE': float(np.mean(astrom_decshift)),
                    'CRMSRC': os.path.basename(fn[i])[:60]}
        crr_mskfn = [_crr_mask_fn(f, crr_mask_dir, 'crmsk') for f in fn]

        crr_pre = None
        if crr_premask:
            crr_pre = np.zeros((nfn, nz3, ny3, nx3), dtype=bool)
            npre = 0
            for i in range(nfn):
                m, why = _crr_read_mask(_crr_mask_fn(fn[i], crr_mask_dir, 'premsk'),
                                        crr_gridid, crr_shape, crr_gridkeys, _crr_place(i))
                if m is not None:
                    crr_pre[i] = m
                    npre += 1
                elif why != 'no file':
                    print('  [WARNING] pre-mask for {} ignored: {}'.format(
                        os.path.basename(fn[i]), why))
            if npre:
                print('  pre-masks: {}/{} frames, {} voxels excluded from the '
                      'statistics and from the stack'.format(npre, nfn, int(crr_pre.sum())))
            else:
                print('  crr_premask is set but no usable .premsk.fits was found.')
                crr_pre = None

        crr_reused = False
        if crr_reuse_mask:
            loaded, why_not = [], []
            for i in range(nfn):
                m, why = _crr_read_mask(crr_mskfn[i], crr_gridid, crr_shape,
                                        crr_gridkeys, _crr_place(i))
                loaded.append(m)
                if m is None:
                    why_not.append('{}: {}'.format(os.path.basename(fn[i]), why))
            if not why_not:
                for i in range(nfn):
                    crr_rej[i] = loaded[i]
                crr_reused = True
                # the mask is tile-compressed, so its keywords live in ext 1
                try:
                    nmade = fits.getheader(crr_mskfn[0], 1).get('CRMNFRM', nfn)
                except Exception:
                    nmade = fits.getheader(crr_mskfn[0]).get('CRMNFRM', nfn)
                print('  reusing saved CR masks for all {} frames (made by comparing {} '
                      'frames); skipping the comparison'.format(nfn, nmade))
            else:
                print('  cannot reuse the saved CR masks, recomputing for every frame:')
                for w in why_not[:10]:
                    print('    ' + w)
                if len(why_not) > 10:
                    print('    ... and {} more'.format(len(why_not) - 10))
                if any('the frame moved on the grid' in w for w in why_not):
                    print('    [hint] the shift list, the pre-shift and use_astrom move '
                          'each frame on the grid without changing the grid itself, so a '
                          'mask made before they changed points at the wrong voxels. '
                          'Use the same .shift.list for both runs.')
                if any('different output grid' in w for w in why_not):
                    print('    [hint] the grid follows the FIRST frame in the list plus '
                          'the .par. Pin stack_ad in the .par and pass wave_ref / nwave / '
                          'dwave to make it independent of which frames you stack.')

        # Margin for the widest filter applied along wavelength, which is the
        # zratio scale window (51 channels), not just zmed. The effect is
        # small -- the truncation at a block edge biases the numerator and the
        # denominator of the scale the same way and largely cancels, and a
        # test with marginal cosmic rays piled on a block boundary moved 1
        # flag in 70 between pad 10 and pad 28, none in realistic cases -- but
        # with the block sized by crr_maxmem the margin is cheap and the
        # blocked result then matches the unblocked one exactly.
        pad = max(int(crr_zwindow), int(crr_zmed), int(crr_zratio)) // 2 \
            + int(crr_grow_z) + 1
        per_chan = (ny3 * nx3) * (nfn * 13 + crr_ngroup * 9)
        if crr_zblock is None:
            # the largest block that keeps the working set under crr_maxmem,
            # then evened out: ceil(nz/blk) blocks of ceil(nz/nblocks) channels
            # each. Taking blk straight from the memory budget leaves a long
            # first block and a short last one, which costs the same total
            # work but sets the peak memory by the long one -- 2596 channels
            # at blk=2216 is 2244 + 408, where two blocks of 1298 do the same
            # work with 59% of the peak.
            blk = int(crr_maxmem * 1e9 / max(per_chan, 1)) - 2 * pad
            blk = int(min(max(blk, 4 * pad), nz3))
            nblk = int(np.ceil(nz3 / float(blk)))
            crr_zblock = int(np.ceil(nz3 / float(nblk)))
        edges = list(range(0, nz3, crr_zblock))
        nproc = sum(min(z + crr_zblock + pad, nz3) - max(z - pad, 0) for z in edges)
        if not crr_reused:
            print('  crr_zblock = {} ({} block(s), {:.1f}% extra work from the {}-channel '
              'margins, {:.1f} GB peak)'.format(
                  crr_zblock, len(edges), 100. * (nproc / float(nz3) - 1.), pad,
                  per_chan * min(crr_zblock + 2 * pad, nz3) / 1e9))
        n_test = 0
        for z0 in (tqdm(range(0, nz3, crr_zblock)) if not crr_reused else []):
            z1 = min(z0 + crr_zblock, nz3)
            p0, p1 = max(z0 - pad, 0), min(z1 + pad, nz3)
            nzb = p1 - p0
            f_i = np.zeros((nfn, nzb, ny3, nx3), dtype=np.float32)
            f_v = np.zeros_like(f_i)
            f_w = np.zeros_like(f_i)
            f_ok = np.zeros(f_i.shape, dtype=bool)
            for i in range(nfn):
                if low_mem == False:
                    c_i, c_v = data0[i, p0:p1], vdata0[i, p0:p1] #reading the montage drizzled cubes already in memory
                    c_m, c_e = mdata0[i, p0:p1], edata0[i, p0:p1]
                else:
                    c_i = np.nan_to_num(fits.open(montfns[i], memmap=True)[0].data[p0:p1] / _iss6_factor(trimfn[i]))
                    c_v = np.nan_to_num(fits.open(montvfns[i], memmap=True)[0].data[p0:p1] / _iss6_factor(trimvfn[i]))
                    c_m = np.ceil(fits.open(montmfns[i], memmap=True)[0].data[p0:p1] / _iss6_factor(trimmfn[i]))
                    c_m[~np.isfinite(c_m)] = 128
                    c_e = np.nan_to_num(fits.open(montefns[i], memmap=True)[0].data[p0:p1] / _iss6_factor(trimefn[i]))
                ok = (c_m == 0) & np.isfinite(c_i) & np.isfinite(c_v) & (c_v != 0) & (c_e > 0)
                if crr_pre is not None:
                    ok &= ~crr_pre[i, p0:p1]
                f_i[i], f_v[i], f_ok[i] = c_i, c_v, ok
                if medcube:
                    f_w[i] = ok
                elif len(weights) == 0:
                    f_w[i] = np.where(ok, c_e / fluxnorm[i]**2, 0.)
                else:
                    f_w[i] = np.where(ok, float(weights[i]), 0.)

            # group cubes (a group of one frame is just that frame)
            if crr_ngroup == nfn:
                g_i, g_v, g_good = f_i, f_v, f_ok
            else:
                g_i = np.zeros((crr_ngroup, nzb, ny3, nx3), dtype=np.float32)
                g_v = np.zeros_like(g_i)
                g_good = np.zeros(g_i.shape, dtype=bool)
                for g in range(crr_ngroup):
                    mem = crr_members[g]
                    w = f_w[mem]
                    wsum = w.sum(axis=0)
                    with np.errstate(invalid='ignore', divide='ignore'):
                        g_i[g] = np.where(wsum > 0, (w * np.where(f_ok[mem], f_i[mem], 0.)).sum(axis=0) / wsum, np.nan)
                        g_v[g] = np.where(wsum > 0, (w**2 * np.where(f_ok[mem], f_v[mem], 0.)).sum(axis=0) / wsum**2, np.nan)
                    g_good[g] = (wsum > 0) & np.isfinite(g_i[g]) & np.isfinite(g_v[g]) & (g_v[g] > 0)

            flags = _crr_flag_outliers(g_i, g_v, g_good, nsig=crr_nsig, rel=crr_rel,
                                       window=crr_window, zwindow=crr_zwindow, zmed=crr_zmed,
                                       grow_xy=crr_grow_xy, grow_z=crr_grow_z,
                                       grow_nsig=crr_grow_nsig, pad_xy=crr_pad_xy,
                                       zratio=crr_zratio, zstep=crr_zstep,
                                       fsys=crr_fsys, fgrad=crr_fgrad,
                                       grow_rel=crr_grow_rel)
            c = slice(z0 - p0, z0 - p0 + (z1 - z0))
            n_test += int((g_good.sum(axis=0) > 1)[c].sum())
            for g in range(crr_ngroup):
                mem = crr_members[g]
                if not flags[g].any():
                    continue
                if len(mem) == 1:
                    crr_rej[mem[0], z0:z1] = (flags[g] & f_ok[mem[0]])[c]
                    continue
                # a group of several frames: find the responsible frame(s)
                oth = [h for h in range(crr_ngroup) if h != g]
                if len(oth) == 1:
                    est, vest = g_i[oth[0]], g_v[oth[0]]
                else:
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', RuntimeWarning)
                        n_o = g_good[oth].sum(axis=0)
                        est = np.nanmedian(np.where(g_good[oth], g_i[oth], np.nan), axis=0)
                        vest = (np.where(n_o >= 3, np.pi / 2., 1.) *
                                np.nanmean(np.where(g_good[oth], g_v[oth], np.nan), axis=0) /
                                np.maximum(n_o, 1))
                rej = _crr_frame_reject(f_i[mem], f_v[mem], f_ok[mem], flags[g], est, vest,
                                        frame_nsig=crr_frame_nsig, fsys=crr_fsys)
                for jj, fr in enumerate(mem):
                    crr_rej[fr, z0:z1] = rej[jj][c]
            del f_i, f_v, f_w, f_ok, g_i, g_v, g_good, flags

        # Save the cosmic-ray flags BEFORE folding in the pre-mask: the
        # pre-mask is the user's own input, is re-read on every run and may be
        # edited between runs, so baking it into the saved .crmsk.fits would
        # make a stale copy of it impossible to get rid of.
        if crr_save_mask and not crr_reused:
            par = dict(CRMNSIG=(float(crr_nsig), 'crr_nsig'),
                       CRMREL=(float(crr_rel), 'crr_rel'),
                       CRMZRAT=(int(crr_zratio), 'crr_zratio'),
                       CRMFSYS=(float(crr_fsys), 'crr_fsys'),
                       CRMFGRAD=(float(crr_fgrad), 'crr_fgrad'),
                       CRMGREL=(float(crr_grow_rel), 'crr_grow_rel'),
                       CRMNGRP=(int(crr_ngroup), 'number of comparison groups'),
                       CRMPRE=(bool(crr_pre is not None), 'a pre-mask was applied'))
            for i in range(nfn):
                _crr_write_mask(crr_mskfn[i], crr_rej[i], crr_gridid, dimension,
                                hdr0, list(fn), par, _crr_place(i))
            print('  CR masks written to {}'.format(
                os.path.dirname(crr_mskfn[0]) or '.'))

        if crr_pre is not None:
            crr_rej |= crr_pre

        crr_nrej_frame = crr_rej.reshape(nfn, -1).sum(axis=1).astype(np.int64)
        for i in range(nfn):
            if crr_reused:
                print('  {}: {} voxels masked (reused)'.format(
                    os.path.basename(fn[i]), crr_nrej_frame[i]))
                continue
            frac = crr_nrej_frame[i] / max(n_test, 1)
            print('  {}: {} voxels rejected ({:.3g}% of testable voxels)'.format(
                os.path.basename(fn[i]), crr_nrej_frame[i], 100. * frac))
            if frac > 0.01:
                print('  [WARNING] >1% rejected. Check the alignment and flux scaling, '
                      'or raise crr_rel / crr_nsig.')
        crr_nrej_3d = np.zeros((nx3, ny3, nz3), dtype=np.int16)

    # Stacking!!!
    print('Stacking...')
    data_3d=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3']),dtype=np.float64)
    vdata_3d=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3']),dtype=np.float64)
    edata_3d=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3']),dtype=np.float64)
    mdata_3d=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3']),dtype=np.int16)+1
    for ii in tqdm(range(dimension[0])):

        if low_mem==False:
            img=data0[:,:,:,ii]
            var=vdata0[:,:,:,ii]
            mask=mdata0[:,:,:,ii]
            exp=edata0[:,:,:,ii]
        else:
            #cache columns
            img = np.zeros((len(fn), hdr0['NAXIS3'], dimension[1]))
            var = np.zeros((len(fn), hdr0['NAXIS3'], dimension[1]))
            mask = np.zeros((len(fn), hdr0['NAXIS3'], dimension[1]))
            exp = np.zeros((len(fn), hdr0['NAXIS3'], dimension[1]))
            for i in range(len(fn)):
                newcube=fits.open(montfns[i])[0].data / _iss6_factor(trimfn[i])
                newcube[~np.isfinite(newcube)]=0.
                newcubev=fits.open(montvfns[i])[0].data / _iss6_factor(trimvfn[i])
                newcubev[~np.isfinite(newcubev)]=0.
                newcubem=np.ceil(fits.open(montmfns[i])[0].data / _iss6_factor(trimmfn[i]))
                newcubem[~np.isfinite(newcubem)]=128
                newcubee=fits.open(montefns[i])[0].data / _iss6_factor(trimefn[i])
                newcubee[~np.isfinite(newcubee)]=0.
                img[i,:,:]=newcube[:, :, ii]
                var[i,:,:]=newcubev[:, :, ii]
                mask[i,:,:]=newcubem[:, :, ii]
                exp[i,:,:]=newcubee[:, :, ii]


        mask[~np.isfinite(img)]=1
        mask[~np.isfinite(var)]=1
        mask[var==0]=1

        if camera=='RED' and nsigma_clip > 0:
            # additional sigma clipping
            clip_img = img.copy()
            clip_img[mask!=0] = np.nan
            clip_std = np.nanstd(clip_img, axis=0)
            clip_std = np.repeat(clip_std[np.newaxis, :, :], len(fn), axis=0)
            clip_med = np.nanmedian(clip_img, axis=0)
            clip_med = np.repeat(clip_med[np.newaxis, :, :], len(fn), axis=0)

            clip_mask = (np.abs(clip_img - clip_med) > clip_std * nsigma_clip)
            mask[clip_mask]=1

        # crr_split: mask the voxels rejected in each frame
        if crr_split:
            rej_col = crr_rej[:, :, :, ii]
            for fr in range(len(fn)):
                if rej_col[fr].any():
                    mask[fr][rej_col[fr]] = 1
            crr_nrej_3d[ii, :, :] = np.transpose(rej_col.sum(axis=0)).astype(np.int16)

        q=(mask==0)
        if np.sum(q)==0:
            continue

        weight=np.zeros(var.shape)
        #weight[var!=0]=1/np.abs(var[var!=0])
        fluxweight = 1 / np.repeat(np.repeat(np.array(fluxnorm**2)[:,np.newaxis],
                             hdr0['NAXIS3'],axis=1)[:,:,np.newaxis],dimension[1],axis=2)
        if len(weights)==0:
            weight = exp.copy() * fluxweight
        else:
            weight=np.repeat(np.repeat(np.array(weights)[:,np.newaxis],
                             hdr0['NAXIS3'],axis=1)[:,:,np.newaxis],dimension[1],axis=2).astype(float)
        weight[~q]=np.nan

        if medcube == True:
            weight = np.ones(var.shape)

        #weight[~np.isfinite(weight)]=0

        #q2=stats.sigma_clip(img[q],sigma=5,masked=True)
        #weight[q][q2.mask]=0

        data_3d[ii,:,:]=np.transpose(np.nansum(img*weight,axis=0)/np.nansum(weight,axis=0))
        vdata_3d[ii,:,:]=np.transpose(np.nansum(weight**2*var,axis=0)/np.nansum(weight,axis=0)**2)
        if len(weights)==0:
            edata_3d[ii,:,:] = np.transpose(np.sum(exp * fluxweight * np.isfinite(weight),axis=0))
        else:
            edata_3d[ii,:,:] = np.transpose(np.sum(exp * weight * np.isfinite(weight), axis=0))
        mdata_3d[ii,:,:]=(edata_3d[ii,:,:]==0).astype(int)



    # remove temp files
    for i in range(len(fn)):
        if keep_mont==False:
            os.remove(montfns[i])
            os.remove(montvfns[i])
            os.remove(montmfns[i])
            os.remove(montefns[i])
            if method.lower()=='drizzle':
                os.remove(montfns[i].replace('mont','mont_area'))
                os.remove(montvfns[i].replace('mont','mont_area'))
                os.remove(montmfns[i].replace('mont','mont_area'))
                os.remove(montefns[i].replace('mont','mont_area'))
            else:
                os.remove(montfns[i].replace('.'+method_flag,'.'+method_flag+'_area'))
                os.remove(montvfns[i].replace('.'+method_flag,'.'+method_flag+'_area'))
                os.remove(montmfns[i].replace('.'+method_flag,'.'+method_flag+'_area'))
                os.remove(montefns[i].replace('.'+method_flag,'.'+method_flag+'_area'))

        if keep_trim==False:
            os.remove(trimfn[i])
            os.remove(trimvfn[i])
            os.remove(trimmfn[i])
            os.remove(trimefn[i])

    # exposure-time weighted airmass
    hdr0['AIRMASS']=np.sum(np.array(airmass_list)*np.array(expt_list)/np.sum(expt_list))
    # write
    vhdr0=hdr0.copy()
    if suffix=='cubes':
        #vhdr0['BUNIT']='10^(-32)erg2/s2/cm4/Angstrom2'
        vhdr0['BUNIT']='1e-32 erg2/s2/cm4/Angstrom2/arcsec4'
    else:
        #vhdr0['BUNIT']='adu2/s2'
        vhdr0['BUNIT']='count2/s2/arcsec4'

    mhdr0=hdr0.copy()
    del mhdr0['BUNIT']
    mhdr0['BITPIX']=16

    ehdr0=hdr0.copy()
    ehdr0['BUNIT']='s'

    if method.lower()!='drizzle':
        suffix_all=suffix+'_'+method_flag[0]
    else:
        suffix_all=suffix

    data_3d=np.nan_to_num(data_3d)
    hdu_i=fits.PrimaryHDU(data_3d.T,header=hdr0)
    if use_astrom:
        icubefn=fnlist.replace('.list','_i'+suffix_all+'_wcs.fits')
    else:
        icubefn=fnlist.replace('.list','_i'+suffix_all+'.fits')
    hdu_i.writeto(icubefn,overwrite=True)
    vdata_3d=np.nan_to_num(vdata_3d)
    hdu_v=fits.PrimaryHDU(vdata_3d.T,header=vhdr0)
    hdu_v.writeto(fnlist.replace('.list','_v'+suffix_all+'.fits'),overwrite=True)
    hdu_m=fits.PrimaryHDU(mdata_3d.T,header=mhdr0)
    hdu_m.writeto(fnlist.replace('.list','_m'+suffix_all+'.fits'),overwrite=True)
    hdu_e=fits.PrimaryHDU(edata_3d.T,header=ehdr0)
    hdu_e.writeto(fnlist.replace('.list','_e'+suffix_all+'.fits'),overwrite=True)

    if crr_split:
        chdr = mhdr0.copy()
        chdr['COMMENT'] = 'number of frames rejected per voxel (crr_split)'
        chdr['CRRMODE'] = (str(crr_split_mode) if isinstance(crr_split_mode, str) else 'custom',
                           'crr_split comparison mode')
        chdr['CRRNSIG'] = (crr_nsig, 'absolute threshold [sigma]')
        chdr['CRRREL'] = (crr_rel, 'relative threshold [local level]')
        chdr['CRRWIN'] = (crr_window, 'local-level window, spatial [pix]')
        chdr['CRRZWIN'] = (crr_zwindow, 'local-level window, spectral [pix]')
        chdr['CRRZMED'] = (crr_zmed, 'spectral median high-pass [pix]')
        chdr['CRRGRXY'] = (crr_grow_xy, 'max spatial growth [pix]')
        chdr['CRRGRZ'] = (crr_grow_z, 'max spectral growth [pix]')
        chdr['CRRGRSIG'] = (crr_grow_nsig, 'growth threshold [sigma]')
        chdr['CRRPADXY'] = (crr_pad_xy, 'unconditional spatial pad [pix]')
        chdr['CRRZRAT'] = (crr_zratio, 'channels for the per-spaxel scale match')
        chdr['CRRZSTEP'] = (crr_zstep, 'channels between scale evaluations')
        chdr['CRRNGRP'] = (crr_ngroup, 'number of crr_split groups')
        chdr['CRRFSYS'] = (crr_fsys, 'fractional systematic in the noise')
        chdr['CRRFGRAD'] = (crr_fgrad, 'gradient term in the noise')
        chdr['CRRGREL'] = (crr_grow_rel, 'relative threshold while growing')
        hdu_c0 = fits.PrimaryHDU(crr_nrej_3d.T, header=chdr)
        tab_c = fits.BinTableHDU(table.Table([[os.path.basename(f) for f in fn], crr_groups,
                                              crr_nrej_frame],
                                             names=['file', 'group', 'n_rejected']), name='FRAMES')
        hdul_c = fits.HDUList([hdu_c0, tab_c])
        if crr_save_files:
            mhdr = mhdr0.copy()
            mhdr['EXTNAME'] = 'FRAMEREJ'
            mhdr['COMMENT'] = 'per-frame rejection mask, axis 4 = frame (see FRAMES table)'
            hdul_c.append(fits.ImageHDU(crr_rej.astype(np.uint8), header=mhdr))
        hdul_c.writeto('kcwi_stack/' + fnlist.replace('.list', '.crrsplit.fits'), overwrite=True)

    if use_astrom:
        # wavelength range
        wavebin=par['wavebin']
        if wavebin[0]==-1:
            wavebin=[4000.,5000.]

        cube=hdu_i.data.T
        sz=cube.shape
        wcs_cube=wcs.WCS(hdu_i.header)
        wave=wcs_cube.all_pix2world(np.zeros(sz[2]),np.zeros(sz[2]),np.arange(sz[2]),0)
        wave=wave[2]*1e10

        # collapsing
        qwave=(wave>wavebin[0]) & (wave<wavebin[1])

        hdr_img=hdu_i.header.copy()
        del hdr_img['CD3_3']
        del hdr_img['CRVAL3']
        del hdr_img['CRPIX3']
        del hdr_img['NAXIS3']
        del hdr_img['CTYPE3']
        del hdr_img['CNAME3']
        del hdr_img['CUNIT3']
        hdr_img['NAXIS']=2

        img=np.zeros((sz[0],sz[1]))
        for ii in range(sz[0]):
            for jj in range(sz[1]):
                q=(cube[ii,jj,qwave]!=0) & (np.isfinite(cube[ii,jj,qwave])==1)
                if np.sum(q)>0:
                    img[ii,jj]=np.mean(cube[ii,jj,qwave][q])

        hdu_best=fits.PrimaryHDU(img.T,header=hdr_img)
        hdu_best.writeto('kcwi_astrom/'+fnlist.replace('.list','_i'+suffix_all+'.thum.fits'),overwrite=True)

    end=ostime.time()
    #print(end-start)

    return




def _align_boxes(box):
    """Normalise an alignment box specification to an (n, 4) integer array.

    Accepts the historical [x0, x1, y0, y1], a list of such boxes, or an
    (n, 4) array. Several boxes let the cross-correlation use more than one
    source in the field, which matters when no single object is bright enough
    on its own, when the brightest one is extended or sits on a gradient, or
    when a cosmic ray or a detector artefact can out-weigh the real source and
    lock the solution onto a spurious peak.

    Args:
        box: [x0, x1, y0, y1], or a sequence of them.

    Returns:
        ndarray [n, 4]: one row per box.
    """
    b = np.atleast_2d(np.asarray(box))
    if b.ndim != 2 or b.shape[1] != 4:
        raise ValueError('box must be [x0, x1, y0, y1] or a list of them; '
                         'got shape {}'.format(np.shape(box)))
    return b.astype(int)


def _align_crls(img0, img, boxes, dx, dy, background_subtraction=False,
                background_level=-1000., scale=1, combine='norm'):
    """Cross-correlation surface of two thumbnails over one or more boxes.

    Each box gives its own surface from the statistic kcwi_align has always
    used -- the mean of the product of the two cutouts over the pixels where
    that product is non-zero -- so with a single box the result is unchanged.

    With several boxes the surfaces are combined. Adding them raw would let the
    brightest source decide the shift by itself, which is the thing several
    boxes are meant to avoid, so each surface is first put on a common footing:
    its median is removed and it is divided by its own robust scatter, which
    makes the height of a peak a significance rather than a flux. A faint
    source then counts as much as a bright one, and a box whose peak is
    ambiguous contributes a flat surface that neither helps nor hurts.
    combine='sum' keeps the flux weighting instead.

    Args:
        img0, img (ndarray): reference and current thumbnail, NaNs replaced.
        boxes (ndarray [n, 4]): [x0, x1, y0, y1] per box, unscaled.
        dx, dy (ndarray): integer shift grids, same shape as the output.
        background_subtraction (bool), background_level (float): as in
            kcwi_align.
        scale (int): multiply the box coordinates by this -- the upsampling
            factor in the second pass.
        combine (str): 'norm' (above), 'median' (the median of the
            standardised surfaces instead of their mean -- immune to one bad
            box, needs at least three), or 'sum' (raw, flux-weighted).

    Returns:
        (ndarray, list of ndarray): the combined surface and the per-box ones.
    """
    surfs = []
    for (bx0, bx1, by0, by1) in boxes:
        bx0, bx1 = int(bx0) * int(scale), int(bx1) * int(scale)
        by0, by1 = int(by0) * int(scale), int(by1) * int(scale)
        cut0 = img0[bx0:bx1, by0:by1]
        c = np.zeros(dx.shape)
        for ii in range(c.shape[0]):
            for jj in range(c.shape[1]):
                cut = img[bx0 - dx[ii, jj]:bx1 - dx[ii, jj],
                          by0 - dy[ii, jj]:by1 - dy[ii, jj]]
                if cut.shape != cut0.shape or cut.size == 0:
                    continue          # the shifted box ran off the thumbnail
                cut = cut.copy()
                cut0_ = cut0.copy()
                if background_subtraction:
                    if background_level == -1000:
                        back_val = np.median(cut)
                        back_val0 = np.median(cut0_)
                    else:
                        back_val = float(background_level)
                        back_val0 = back_val
                    cut = cut - back_val
                    cut0_ = cut0_ - back_val0
                else:
                    if background_level != -1000:
                        cut[cut < background_level] = 0
                        cut0_[cut0_ < background_level] = 0
                cut[cut < 0] = 1e-10
                cut0_[cut0_ < 0] = 1e-10
                mult = cut0_ * cut
                n = np.nansum(mult != 0)
                if n > 0:
                    c[ii, jj] = np.nansum(mult) / n
        surfs.append(c)

    if len(surfs) == 1:
        return surfs[0], surfs
    if combine == 'sum':
        return np.sum(surfs, axis=0), surfs
    # standardise each surface: subtract its median and divide by its own
    # robust scatter, so a peak's height becomes a significance rather than a
    # flux and a faint source counts as much as a bright one
    z = []
    for c in surfs:
        f = c[np.isfinite(c) & (c != 0)]
        if f.size < 4:
            continue
        med = np.median(f)
        sig = 1.4826 * np.median(np.abs(f - med))
        if not np.isfinite(sig) or sig <= 0:
            sig = np.std(f)
        if not np.isfinite(sig) or sig <= 0:
            continue
        z.append((c - med) / sig)
    if len(z) == 0:
        return np.sum(surfs, axis=0), surfs
    if combine == 'median' and len(z) >= 3:
        # one box can still carry a very sharp spurious peak -- a cosmic ray,
        # an artefact, a neighbour that moved -- and a sum lets its enormous
        # significance outvote the others. The median across boxes cannot be
        # moved by one of them, so it is the choice when a box might be bad.
        return np.median(np.stack(z, axis=0), axis=0), surfs
    return np.mean(np.stack(z, axis=0), axis=0), surfs


def kcwi_align(fnlist,wavebin=[-1.,-1.],box=[-1,-1,-1,-1],pixscale_x=-1.,pixscale_y=-1.,
    orientation=-1000.,dimension=[-1.,-1.],preshiftfn='',trim=[-1,-1],cubed=False,
    noalign=False,display=True,search_size=-1000,conv_filter=-1000,upfactor=-1000.,
    background_subtraction=False,background_level=-1000.,method='interp',
    use_regmask=True, npix_trim = 4, box_combine='norm'):

    """
    Align individual data cubes and correct their relative astrometry using
        cross-correlation.

    Args:
        fnlist (str): file name of the list that contains individual files.
        wavebin ([float, float]): wavelength range for white-light image. Set to
            override par file.
        box ([int x 4] or [[int x 4], ...]): alignment box as
            [x0, x1, y0, y1] in thumbnail pixels, or a list of such boxes to
            use several sources at once. Overrides the par file, where one
            align_box line per box does the same thing. Several boxes help
            when no single object is bright enough on its own, when the
            brightest one is extended or sits on a gradient, and when an
            artefact could otherwise out-weigh the real source; see
            _align_crls for how the boxes are combined.
        pixscale_x (float): override the parfile pixel size in the x-direction
            in arcsec.
        pixscale_y (float): override the parfile pixel size in the y-direction
            in arcsec.
        dimension ([float, float]): numbers of pixels in x- and y-directions.
            (Set to override the par file.)
        orientation (float): position angle of the up direction in degrees.
            (Set to override the par file.)
        preshiftfn (str): set if there is a preshift file.
        trim ([float, float]): number of pixels to be trimmed ([lower, upper]).
            Override the list file.
        cubed (bool): using cubed?
        noalign(bool): only generate the thum files without actually conducting
            the alignment.
        display (bool): generate diagnositc plots?
        search_size (float): size of the search box. Override the par file.
        conv_filter (float): size of the convolution filter when looking for the
            local maximum. Override the par file.
        upfactor (float): oversampling factor for the second pass. Override the
            par file.
        background_subtraction (bool): set for background subtraction before
            the cross correlation.
        background_level (float): set to use a different background level other
            than the medium.
        method (str): 'interp' (bilinear interpolation) or 'exact' (drizzling).
        use_regmask (bool): set if certain pixels in individual files need to be
            masked by Region files.
        npix_trim (int): number of pixels to trim cubes by. Default = 4. May want
            npix_trim = 1 for Large slicer.
        box_combine (str): how several boxes are combined. 'norm' (default)
            standardises each box's surface so every source counts equally;
            'median' takes the median of those instead of the mean, which one
            bad box cannot move (needs three or more); 'sum' adds the raw
            surfaces, i.e. flux-weighted. Ignored for a single box.
            On synthetic fields with three sources, recovering the right
            integer shift over 60 noise realisations:
                                     1 box  norm  median  sum
              faint field, S/N ~ 6      31    42      35   43
              same + a cosmic ray in
              the first box              0     8      25    0
            So 'norm' for an ordinary field and 'median' when a box might be
            contaminated; 'median' costs a little in the clean case because
            taking the median of three surfaces throws away signal.


    Returns:
        None
        (A table of relative pixel shifts are generated as '*.shift.list'.)

    """

    # support for direct putting in FITS
    if fnlist.endswith('.fits'):
        fnlist=fnlist.replace('.fits','.list')

    # use cubes or cubed?
    if cubed==False:
        suffix='cubes'
    else:
        suffix='cubed'

    parfn=fnlist.replace('.list','.par')
    par=kcwi_stack_readpar(parfn)

    # wavelength range
    if wavebin[0]==-1:
        wavebin=par['wavebin']

    # define alignment box(es)
    if np.ndim(box) == 1 and box[0] == -1:
        box = par['align_box']
        if np.ndim(box) == 1 and box[0] == -1:
            box = [25, 50, 25, 40]
    boxes = _align_boxes(box)
    if len(boxes) > 1:
        print('Alignment: {} boxes, combined with \'{}\''.format(len(boxes), box_combine))
        for ib, bb in enumerate(boxes):
            print('   box {}: x {}-{}, y {}-{}'.format(ib, bb[0], bb[1], bb[2], bb[3]))

    # size of pixels post projection
    if pixscale_x==-1:
        pixscale_x=par['align_xpix']
        if pixscale_x==-1:
            pixscale_x=0.3
    pixscale_x=pixscale_x/3600.

    if pixscale_y==-1:
        pixscale_y=par['align_ypix']
        if pixscale_y==-1:
            pixscale_y=0.3
    pixscale_y=pixscale_y/3600.

    # post-projection image size
    if dimension[0]==-1:
        dimension=par['align_dimension']
        if dimension[0]==-1:
            dimension=[100,100]

    # size of search steps in x and y directions in pixel units after projection
    if search_size==-1000:
        search_size=par['align_search_size']
        if search_size==-1000:
            search_size=10
    search_size=int(search_size)

    # size of the convolution kernel for searching local maxima
    if conv_filter==-1000:
        conv_filter=par['align_conv_filter']
        if conv_filter==-1000:
            conv_filter=2
    conv_filter=int(conv_filter)

    # upsample factor for fine-alignment
    if upfactor==-1000.:
        upfactor=par['align_upfactor']
        if upfactor==-1000.:
            upfactor=10.
    upfactor=np.ceil(upfactor).astype(int)

    # background subtraction in the alignment cut?
    if background_subtraction==False:
        background_subtraction=par['background_subtraction']
    background_subtraction=bool(background_subtraction)

    if background_level==-1000:
        background_level=par['background_level']

    # make tmp directory
    if not os.path.exists('kcwi_align'):
        os.makedirs('kcwi_align')

    # read fnlist
    trimtab=ascii.read(fnlist,format="no_header")
    fn=trimtab['col1']
    fn=[i+'_i'+suffix+'.fits' for i in fn]
    trim=np.array([trimtab['col2'],trimtab['col3']])

    if preshiftfn=='':
        preshiftfn=fnlist.replace('.list','.preshift.list')
        if path.isfile(preshiftfn)==False:
            preshiftfn=fnlist.replace('.list','.pre.list')
            if path.isfile(preshiftfn)==False:
                preshiftfn=''
    if preshiftfn!='':
        pretab=ascii.read(preshiftfn,format='no_header')
        prefn=[i+'_i'+suffix+'.fits' for i in pretab['col1']]
        prera=pretab['col2']
        predec=pretab['col3']

    if display==False:
        oldbackend=matplotlib.get_backend()
        matplotlib.use('Agg')

    # construct WCS
    hdulist=fits.open(fn[0])
    hdrtmp=hdulist[0].header.copy()
    hdulist.close()
    wcstmp=wcs.WCS(hdrtmp).copy()
    center=wcstmp.wcs_pix2world((wcstmp.pixel_shape[0]-1)/2.,(wcstmp.pixel_shape[1]-1)/2.,0,0,ra_dec_order=True)

    # red or blue?
    camera = hdrtmp['CAMERA']
    if camera != 'RED' and camera !='BLUE':
        raise ValueError('Unknown camera type - {}'.format(camera))
        return

    # wavebin
    if wavebin[0]==-1:
        wavebin=[hdrtmp['WAVGOOD0'],hdrtmp['WAVGOOD1']]

    if par['align_ad'][0]!=-1:
        center=par['align_ad']

    hdr0=hdrtmp.copy()
    hdr0['NAXIS1']=dimension[0]
    hdr0['NAXIS2']=dimension[1]
    hdr0['CRPIX1']=(dimension[0]+1)/2.
    hdr0['CRPIX2']=(dimension[1]+1)/2.
    hdr0['CRVAL1']=float(center[0])
    hdr0['CRVAL2']=float(center[1])
    old_cd11=hdr0['CD1_1']
    old_cd12=hdr0['CD1_2']
    old_cd21=hdr0['CD2_1']
    old_cd22=hdr0['CD2_2']
    hdr0['CD1_1']=-pixscale_x
    hdr0['CD2_2']=pixscale_y
    hdr0['CD1_2']=0.
    hdr0['CD2_1']=0.
    hdr0['NAXIS']=2
    del hdr0['CD3_3']
    del hdr0['CRVAL3']
    del hdr0['CRPIX3']
    del hdr0['NAXIS3']
    del hdr0['CTYPE3']
    del hdr0['CNAME3']
    del hdr0['CUNIT3']

    # orientation
    if orientation==-1000:
        orientation=par['align_orientation']
        if orientation==-1000:
            orientation=np.rad2deg(np.arctan(old_cd21/(-old_cd11)))
    hdr0['CD1_1']=-pixscale_x*np.cos(np.deg2rad(orientation))
    hdr0['CD2_1']=pixscale_x*np.sin(np.deg2rad(orientation))
    hdr0['CD1_2']=pixscale_y*np.sin(np.deg2rad(orientation))
    hdr0['CD2_2']=pixscale_y*np.cos(np.deg2rad(orientation))

    # align
    data_thum=np.zeros((dimension[0],dimension[1],len(fn)))
    data0_thum=np.zeros((dimension[0],dimension[1],len(fn)))
    xshift=np.zeros(len(fn))
    yshift=np.zeros(len(fn))
    xshift_xy=np.zeros(len(fn))
    yshift_xy=np.zeros(len(fn))
    pngfn=[]
    for i in range(len(fn)):
        print(os.path.basename(fn[i]))

        hdulist = fits.open(fn[i])

        # Intensity cube
        hdu=hdulist[0]
        img=hdu.data.T.copy()
        sz=img.shape
        wcs_i=wcs.WCS(hdu.header)
        wave=wcs_i.wcs_pix2world(np.zeros(sz[2]),np.zeros(sz[2]),np.arange(sz[2]),0)
        wave=wave[2]*1e10
        qwave=(wave > wavebin[0]) & (wave < wavebin[1])

        hdr=hdu.header.copy()
        del hdr['CD3_3']
        del hdr['CRVAL3']
        del hdr['CRPIX3']
        del hdr['NAXIS3']
        del hdr['CTYPE3']
        del hdr['CNAME3']
        del hdr['CUNIT3']

        hdr['NAXIS']=2

        # mask?
        regfn1 = fn[i].replace('.fits','_rm.reg')
        regfn2 = fn[i].replace('.fits','_rm.reg')
        if "zap" in fn[i]:
            regfn = regfn1
        else:
            regfn = regfn2
        print("Checking for region file: {}".format(regfn))
        if os.path.isfile(regfn) and use_regmask==True and ".reg" in regfn:
            print('Using region mask: {}'.format(regfn))
            # img is (x, y, wave) here, the mask is (ny, nx)
            mask_reg = regions_to_mask(regfn, header=hdr, shape=(img.shape[1], img.shape[0]))
            img[mask_reg.T]=np.nan

        thum=np.zeros((sz[0],sz[1]))
        for ii in range(sz[0]):
            for jj in range(sz[1]):
                q=(img[ii,jj,qwave]!=0) & (np.isfinite(img[ii,jj,qwave])==1)
                if np.sum(q)>0:
                    if camera=='BLUE':
                        thum[ii,jj]=np.mean(img[ii,jj,qwave][q])
                    elif camera=='RED':
                        thum[ii,jj]=np.median(img[ii,jj,qwave][q])

        # additional trimming for PyDRP
        if len(hdulist)>=4:
            mask = hdulist['FLAGS']
            mask_dim = mask.shape
            # set edge pixels = 128
            n_pix = npix_trim # nominally 4 pix, really only need 1
            mask.data[:,0:n_pix,:] = 128
            mask.data[:,:,0:n_pix] = 128
            mask.data[:,mask_dim[1]-n_pix:mask_dim[1],:] = 128
            mask.data[:,:,mask_dim[2]-n_pix:mask_dim[2]] = 128

            twod_mask = np.nanmedian(mask.data[qwave,:,:], axis = 0).T # don't really need [qwave,:,:]
            thum[twod_mask==128] = 0
        elif len(hdulist)!=1:
            raise ValueError('Unknown reduction type')

        # trim
        index_x,index_y=np.where(thum!=0)
        xrange=[index_x.min(),index_x.max()]
        yrange=[index_y.min(),index_y.max()]
        thum[:,yrange[1]-trim[1,i]+1:]=np.nan
        thum[:,:yrange[0]+trim[0,i]]=np.nan
        thum[:xrange[0],:]=np.nan
        thum[xrange[1]:,:]=np.nan
        #thum=np.nan_to_num(thum)
        thum[thum == 0] = np.nan

        # preshift
        if preshiftfn!='':
            index=np.where(np.array(prefn)==os.path.basename(fn[i]))
            index=index[0]
            if len(index)>0:
                index=index[0]
                hdr['CRVAL1']=hdr['CRVAL1']+prera[index]/3600.
                hdr['CRVAL2']=hdr['CRVAL2']+predec[index]/3600.

        # initial projection
        if method=='interp':
            newthum,coverage=reproject_interp((thum.T,hdr),hdr0,order='bilinear')
        elif method=='exact':
            newthum,coverage=reproject_exact((thum.T,hdr),hdr0)
        newthum=newthum.T
        newthum[np.isfinite(newthum)==0]=np.nan
        #hdutmp=fits.PrimaryHDU(newthum)
        #hdutmp.writeto('test.fits',overwrite=True)
        data0_thum[:,:,i]=newthum
        hdr_preshift=hdr.copy()
        wcs_preshift=wcs.WCS(hdr_preshift)

        if i==0:
            thum_1=thum.copy()
            hdr_1=hdr.copy()
            data_thum[:,:,i]=newthum
        else:
            if noalign==False:
                #img0=data0_thum[box[0]:box[1],box[2]:box[3],0]
                #img=newthum[box[0]:box[1],box[2]:box[3]]
                img0=np.nan_to_num(data0_thum[:,:,0])
                img=np.nan_to_num(newthum)

                #dx,dy,ex,ey=chi2_shift(img,img0,return_error='True',verbose=True)
                #img_test=shift.shiftnd(img0,(-dy,-dx))
                #hdu=fits.PrimaryHDU(img.T)
                #hdu.writeto('test.fits',overwrite=True)
                #hdu=fits.PrimaryHDU(newthum.T)
                #hdu.writeto('test2.fits',overwrite=True)
                # turns out image_registration does not work well in this case, will adopt the key algorithm.

                # +/- 10 pixels
                crls_size=search_size+conv_filter
                #print(search_size,crls_size)
                xx=np.linspace(-crls_size,crls_size,2*crls_size+1)
                yy=np.linspace(-crls_size,crls_size,2*crls_size+1)
                dy,dx=np.meshgrid(yy,xx)
                dx=dx.astype(int)
                dy=dy.astype(int)
                crls, crls_box = _align_crls(
                    img0, img, boxes, dx, dy,
                    background_subtraction=background_subtraction,
                    background_level=background_level, combine=box_combine)

                from astropy.visualization import ImageNormalize, ZScaleInterval, LinearStretch, AsinhStretch

                finite = crls[np.isfinite(crls)]
                vmin, vmax = ZScaleInterval(contrast=0.5).get_limits(finite)   # compute limits on finite values only

                norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
                fig,ax=plt.subplots(figsize=(4,4))
                xplot=np.append(xx,xx[1]-xx[0]+xx[-1])-0.5
                yplot=np.append(yy,yy[1]-yy[0]+yy[-1])-0.5
                
                im = ax.imshow(crls.T, origin='lower', norm=norm, cmap='Greys',
                    extent=[xplot[0], xplot[-1], yplot[0], yplot[-1]],
                    interpolation='nearest', aspect='equal')                #fig.colorbar(im, ax=ax)

                # find closest local maximum
                max_conv=ndimage.filters.maximum_filter(crls,2*conv_filter+1)
                maxima=(crls==max_conv) & (crls != 0)
                labeled, num_objects=ndimage.label(maxima)
                slices=ndimage.find_objects(labeled)
                xindex,yindex=[],[]
                for dx,dy in slices:
                    x_center=(dx.start+dx.stop-1)/2
                    xindex.append(x_center)
                    y_center=(dy.start+dy.stop-1)/2
                    yindex.append(y_center)
                #    print(x_center,y_center)
                #print(xindex,yindex)
                xindex=np.array(xindex).astype(int)
                yindex=np.array(yindex).astype(int)
                index=((xindex>=conv_filter) & (xindex<2*crls_size-conv_filter) &
                        (yindex>=conv_filter) & (yindex<2*crls_size-conv_filter))
                xindex=xindex[index]
                yindex=yindex[index]
                # filter out weak ones
                #print(index)
                
                #print(2*crls_size-conv_filter)
                #print(max_conv[xindex,yindex])
                max=np.nanmax(max_conv[xindex,yindex])
                med=np.nanmedian(crls)
                index=np.where(max_conv[xindex,yindex] > 0.3*(max-med)+med)
                xindex=xindex[index]
                yindex=yindex[index]
                r=(xx[xindex]**2+yy[yindex]**2)
                index=r.argmin()
                xshift[i]=xx[xindex[index]]
                yshift[i]=yy[yindex[index]]

                #tmp=np.unravel_index(crls.argmax(),crls.shape)
                #xshift[i]=xx[tmp[0]]
                #yshift[i]=yy[tmp[1]]


                # upsample
                hdr0_up=hdr0.copy()
                hdr0_up['NAXIS1']=hdr0_up['NAXIS1']*upfactor
                hdr0_up['NAXIS2']=hdr0_up['NAXIS2']*upfactor
                hdr0_up['CRPIX1']=(hdr0_up['CRPIX1']-0.5)*upfactor+0.5
                hdr0_up['CRPIX2']=(hdr0_up['CRPIX2']-0.5)*upfactor+0.5
                hdr0_up['CD1_1']=hdr0_up['CD1_1']/upfactor
                hdr0_up['CD2_1']=hdr0_up['CD2_1']/upfactor
                hdr0_up['CD1_2']=hdr0_up['CD1_2']/upfactor
                hdr0_up['CD2_2']=hdr0_up['CD2_2']/upfactor
                if method=='interp':
                    newthum1,coverage=reproject_interp((thum_1.T,hdr_1),hdr0_up,order='bilinear')
                elif method=='exact':
                    newthum1,coverage=reproject_exact((thum_1.T,hdr_1),hdr0_up)
                newthum1=newthum1.T

                # do the shift from last iteration
                wcs_hdr0=wcs.WCS(hdr0)
                tmp=wcs_hdr0.all_pix2world(hdr0['CRPIX1']+xshift[i],hdr0['CRPIX2']+yshift[i],1)
                hdr['CRVAL1']=hdr['CRVAL1']+(float(tmp[0])-hdr0['CRVAL1'])
                hdr['CRVAL2']=hdr['CRVAL2']+(float(tmp[1])-hdr0['CRVAL2'])
                if method=='interp':
                    newthum2,coverage=reproject_interp((thum.T,hdr),hdr0_up,order='bilinear')
                elif method=='exact':
                    newthum2,coverage=reproject_exact((thum.T,hdr),hdr0_up)
                newthum2=newthum2.T

                img0=np.nan_to_num(newthum1)
                img=np.nan_to_num(newthum2)

                # +/-1 pix
                ncrl=np.ceil(upfactor).astype(int)
                xx=np.linspace(-ncrl,ncrl,2*ncrl+1)
                yy=np.linspace(-ncrl,ncrl,2*ncrl+1)
                dy,dx=np.meshgrid(yy,xx)
                dx=dx.astype(int)
                dy=dy.astype(int)
                crls, crls_box = _align_crls(
                    img0, img, boxes, dx, dy,
                    background_subtraction=background_subtraction,
                    background_level=background_level, scale=upfactor,
                    combine=box_combine)

                #plt.figure(1)
                #plt.clf()
                xplot=(np.append(xx,xx[1]-xx[0]+xx[-1])-0.5)/upfactor+xshift[i]
                yplot=(np.append(yy,yy[1]-yy[0]+yy[-1])-0.5)/upfactor+yshift[i]
                #ax.pcolormesh(xplot,yplot,crls.T,cmap='plasma')

                ax.imshow(crls.T, origin='lower', norm=norm, cmap='plasma',
                    extent=[xplot[0], xplot[-1], yplot[0], yplot[-1]],
                    interpolation='nearest', aspect='equal')                #fig.colorbar(im, ax=ax)


                tmp=np.unravel_index(crls.argmax(),crls.shape)
                xshift[i]+=xx[tmp[0]]/upfactor
                yshift[i]+=yy[tmp[1]]/upfactor
                ax.plot(xshift[i],yshift[i],'+',color='r')
                plt.title('{}{}'.format(os.path.basename(fn[i]),
                          '' if len(boxes)==1 else '  ({} boxes)'.format(len(boxes))))


                # get pixel shift
                tmp=wcs_hdr0.all_pix2world(hdr0['CRPIX1']+xshift[i],hdr0['CRPIX2']+yshift[i],1)
                ashift=float(tmp[0])-hdr0['CRVAL1']
                dshift=float(tmp[1])-hdr0['CRVAL2']
                tmp=wcs_preshift.all_world2pix(hdr_preshift['CRVAL1']-ashift,hdr_preshift['CRVAL2']-dshift,1)
                xshift_xy[i]=tmp[0]-hdr_preshift['CRPIX1']
                yshift_xy[i]=tmp[1]-hdr_preshift['CRPIX2']
                print(xshift_xy[i],yshift_xy[i])

                # make shifted thumnail
                hdr_shift=hdr_preshift.copy()
                hdr_shift['CRPIX1']=hdr_shift['CRPIX1']+xshift_xy[i]
                hdr_shift['CRPIX2']=hdr_shift['CRPIX2']+yshift_xy[i]
                thum_shift,coverage=reproject_interp((thum.T,hdr_shift),hdr0)
                thum_shift=thum_shift.T
                thum_shift[np.isfinite(thum_shift)==0]=np.nan
                data_thum[:,:,i]=thum_shift

                #hdr_shift=hdr_preshift.copy()
                #hdr_shift['CRVAL1']=hdr_shift['CRVAL1']+ashift
                #hdr_shift['CRVAL2']=hdr_shift['CRVAL2']+dshift
                #thumshift,coverage=reproject_interp((thum.T,hdr),hdr0,order='bilinear')
                #thumshift=thumshift.T
                #data_thum[:,:,i]=thum_shift

                plt.show()
                fig.tight_layout()

                pngfn.append('kcwi_align/'+os.path.basename(fn[i]).replace('.fits','_align.png'))
                fig.savefig(pngfn[-1])



    if noalign==False:
        hdu=fits.PrimaryHDU(data_thum.T)
        hdu.writeto('kcwi_align/'+fnlist.replace('.list','.thum.fits'),overwrite=True)

        pdf=FPDF()
        for i in pngfn:
            pdf.add_page()
            pdf.image(i,w=180,h=180)
        pdf.output('kcwi_align/'+fnlist.replace('.list','.align.pdf'))


    hdr0['NAXIS']=3
    hdr0['NAXIS3']=len(fn)
    hdu=fits.PrimaryHDU(data0_thum.T,header=hdr0)
    hdu.writeto('kcwi_align/'+fnlist.replace('.list','.thum0.fits'),overwrite=True)

    writefn=[i.replace('.fits','') for i in fn]
    xytable=table.Table([np.array(writefn),xshift_xy,yshift_xy])
    ascii.write(xytable,fnlist.replace('.list','.shift.list'),overwrite=True,format='no_header')

    if display==False:
        matplotlib.use(oldbackend)
    return


#def kcwi_medfilter(fnlist,med_x=-1.,med_y=-1.,med_z=-1,regfn=''):

#    print(os.path.basename(fnlist))

#    suffix='cubes'
#    cubefn=fnlist.replace('.list','_i'+suffix+'.fits')
#    if path.isfile(cubefn)==False:
#        suffix='cubed'
#        cubefn=fnlist.replace('.list','_i'+suffix+'.fits')

#    parfn=fnlist.replace('.list','.par')
#    par=kcwi_stack_readpar(parfn)




def kcwi_astrometry(fnlist,imgfn='',wavebin=[-1.,-1.],display=True,search_size=-1000,
    conv_filter=-1000,upfactor=-1000,box=[-1.,-1.,-1.,-1.],nocrl=0,method='drizzle',
    save_shift=False,interp_order='bilinear',background_subtraction=False,
    background_kcwi=0., background_ref=0.,nhdu=0):

    """
    Conduct astrometry correction of the stacked cube by cross-correlating the
        cube with existing imaging data.

    Args:
        fnlist (str): file name of the list that contains individual files.
        imgfn (str): location of the imaging data to be aligned to.
        wavebin ([float, float]): wavelength range for white-light image. Set to
            override par file.
        display (bool): generate diagnositc plots?
        search_size (float): size of the search box. Override the par file.
        conv_filter (float): size of the convolution filter when looking for the
            local maximum. Override the par file.
        upfactor (float): oversampling factor for the second pass. Override the
            par file.
        box ([float x 4]): location of the alignment box ([x-lower-left,
            y-lower-left, x-upper-right, y-upper-right]). Override the par file.

    Returns:
        None
        (Astrometry-corrected data cubes are generated as '*_wcs.fits')


    """

    print(os.path.basename(fnlist))

    if method.lower()!='drizzle':
        if method.lower()=='nearest-neighbor':
            method_flag='nei'
        elif method.lower()=='bilinear':
            method_flag='lin'
        elif method.lower()=='biquadratic':
            method_flag='qua'
        elif method.lower()=='bicubic':
            method_flag='cub'
        else:
            print('Error: Method not found.')
            return 0

    suffix='cubes'
    if method.lower()!='drizzle':
        suffix=suffix+'_'+method_flag[0]
    cubefn=fnlist.replace('.list','_i'+suffix+'.fits')
    if path.isfile(cubefn)==False:
        suffix='cubed'
        if method.lower()!='drizzle':
            suffix=suffix+'_'+method_flag[0]
        cubefn=fnlist.replace('.list','_i'+suffix+'.fits')

    parfn=fnlist.replace('.list','.par')
    par=kcwi_stack_readpar(parfn)

    # imgfn
    if imgfn=='':
        imgfn=par['ref_fn']
        if imgfn=='':
            print('[ERROR] kcwi_astrometry: Specify alignment image.')
            return

    if par['ref_xy'][0]==-1:
        print('[ERROR] kcwi_astrometry: Set reference x-y coordinate')
        return
    if par['ref_ad'][0]==-1:
        print('[ERROR] kcwi_astrometry: Set reference RA-DEC coordinate')
        return

    # wavelength range
    if wavebin[0]==-1:
        wavebin=par['wavebin']

    # search size
    if search_size==-1000:
        search_size=par['ref_search_size']
        if search_size==-1000:
            search_size=10
    search_size=int(search_size)

    # convoution kernel
    if conv_filter==-1000:
        conv_filter=par['ref_conv_filter']
        if conv_filter==-1000:
            conv_filter=2
    conv_filter=int(conv_filter)

    # upsample factor
    if upfactor==-1000:
        upfactor=par['ref_upfactor']
        if upfactor==-1000:
            upfactor=10
    upfactor=int(upfactor)

    # alignment box
    # default is none, unless specifically defined. This is needed if there is a strong contamination source nearby
    if box[0]==-1:
        box=par['ref_box']

    # nocrl - for QSOs, just put in a number w/o doing cross-correlation
    if nocrl==0:
        nocrl=par['ref_nocrl']


    # make tmp directory
    if not os.path.exists('kcwi_astrom'):
        os.makedirs('kcwi_astrom')


    if display==False:
        oldbackend=matplotlib.get_backend()
        matplotlib.use('Agg')


    hdu_cube=fits.open(cubefn)[0]
    cube=hdu_cube.data.T
    sz=cube.shape
    wcs_cube=wcs.WCS(hdu_cube.header)
    wave=wcs_cube.all_pix2world(np.zeros(sz[2]),np.zeros(sz[2]),np.arange(sz[2]),0)
    wave=wave[2]*1e10

    # initial position
    # record old solution
    oref_ra,oref_dec,_=wcs_cube.all_pix2world(par['ref_xy'][0],par['ref_xy'][1],1,1)
    dref_ra=par['ref_ad'][0]-oref_ra
    dref_dec=par['ref_ad'][1]-oref_dec

    hdu_cube.header['CRPIX1']=par['ref_xy'][0]
    hdu_cube.header['CRPIX2']=par['ref_xy'][1]
    hdu_cube.header['CRVAL1']=par['ref_ad'][0]
    hdu_cube.header['CRVAL2']=par['ref_ad'][1]


    # collapsing
    if wavebin[0]==-1:
        wavebin=[hdu_cube.header['WAVGOOD0'],hdu_cube.header['WAVGOOD1']]
    qwave=(wave>wavebin[0]) & (wave<wavebin[1])

    hdr_img=hdu_cube.header.copy()
    del hdr_img['CD3_3']
    del hdr_img['CRVAL3']
    del hdr_img['CRPIX3']
    del hdr_img['NAXIS3']
    del hdr_img['CTYPE3']
    del hdr_img['CNAME3']
    del hdr_img['CUNIT3']
    hdr_img['NAXIS']=2

    img=np.zeros((sz[0],sz[1]))
    for ii in range(sz[0]):
        for jj in range(sz[1]):
            q=(cube[ii,jj,qwave]!=0) & (np.isfinite(cube[ii,jj,qwave])==1)
            if np.sum(q)>0:
                img[ii,jj]=np.mean(cube[ii,jj,qwave][q])

    img_withbkg = img.copy()
    if background_subtraction:
        img = img - background_kcwi
        img[img < 0] = 0

    # Reference image: use the requested extension, otherwise the first one
    # that actually holds a 2-D image (handles MEF and compressed files).
    hdl_img0 = fits.open(imgfn)
    if nhdu is not None and nhdu < len(hdl_img0) and hdl_img0[nhdu].data is not None \
            and hdl_img0[nhdu].data.ndim == 2:
        hdu_img0 = hdl_img0[nhdu]
    else:
        cands = [h for h in hdl_img0 if h.data is not None and h.data.ndim == 2]
        if len(cands) == 0:
            raise ValueError('%s: no 2-D image extension found' % imgfn)
        if nhdu is not None and nhdu != 0:
            print('  [WARNING] extension %d has no 2-D image, using the first one that does' % nhdu)
        hdu_img0 = cands[0]

    img0=hdu_img0.data.T
    hdr0=hdu_img0.header

    if background_subtraction:
        img0 = img0 - background_ref
        img0[img0 < 0] = 0

    print("reading ref image: "+imgfn)
    print("search_size,cov_filter,upfactor,box:",search_size,conv_filter,upfactor,box)
    #plt.imshow(img,origin='lower',cmap='gray')
    #plt.show()
    if nocrl==0:

        hdr_shift=hdr_img.copy()
        # entry 1:
        print('Iter: 1/2')
        crls_size=search_size+conv_filter
        xx=np.linspace(-crls_size,crls_size,2*crls_size+1)
        yy=np.linspace(-crls_size,crls_size,2*crls_size+1)
        dy,dx=np.meshgrid(yy,xx)
        dx=dx.astype(int)
        dy=dy.astype(int)
        crls=np.zeros(dx.shape)
        for ii in range(crls.shape[0]):
            for jj in range(crls.shape[1]):
                hdr_shift['CRPIX1']=hdr_img['CRPIX1']+dx[ii,jj]
                hdr_shift['CRPIX2']=hdr_img['CRPIX2']+dy[ii,jj]

                img0_shift,coverage=reproject_interp((img0.T,hdr0),hdr_shift,order=interp_order)
                img0_shift=img0_shift.T
                img0_shift=np.nan_to_num(img0_shift)

                if box[0]==-1:
                    mult=img0_shift*img
                else:
                    mult=img0_shift[box[0]:box[1],box[2]:box[3]]*img[box[0]:box[1],box[2]:box[3]]
                if np.sum(mult!=0)>0:
                    crls[ii,jj]=np.sum(mult)/np.sum(mult!=0)
                else:
                    plt.imshow(img0_shift,origin='lower',cmap='gray')
                    plt.show()
                    print('Warning: no overlap region found!')

        fig=plt.figure(1)
        plt.clf()
        xplot=np.append(xx,xx[1]-xx[0]+xx[-1])-0.5
        yplot=np.append(yy,yy[1]-yy[0]+yy[-1])-0.5
        plt.pcolormesh(xplot,yplot,crls.T)

        max_conv=ndimage.filters.maximum_filter(crls,2*conv_filter+1)
        maxima=(crls==max_conv) & (crls != 0)
        labeled, num_objects=ndimage.label(maxima)
        slices=ndimage.find_objects(labeled)
        xindex,yindex=[],[]
        for dx,dy in slices:
            x_center=(dx.start+dx.stop-1)/2
            xindex.append(x_center)
            y_center=(dy.start+dy.stop-1)/2
            yindex.append(y_center)
        xindex=np.array(xindex).astype(int)
        yindex=np.array(yindex).astype(int)
        index=((xindex>=conv_filter) & (xindex<2*crls_size-conv_filter) &
            (yindex>=conv_filter) & (yindex<2*crls_size-conv_filter))
        xindex=xindex[index]
        yindex=yindex[index]
        max=np.max(max_conv[xindex,yindex])
        med=np.median(crls)
        index=np.where(max_conv[xindex,yindex] > 0.3*(max-med)+med)
        xindex=xindex[index]
        yindex=yindex[index]
        r=(xx[xindex]**2+yy[yindex]**2)
        index=r.argmin()
        xmax=xx[xindex[index]]
        ymax=yy[yindex[index]]

        # upscale
        print('Iter: 2/2')
        ncrl=np.ceil(upfactor).astype(int)
        xx=np.linspace(-1,1,2*ncrl+1)+xmax
        yy=np.linspace(-1,1,2*ncrl+1)+ymax
        dy,dx=np.meshgrid(yy,xx)
        crls=np.zeros(dx.shape)
        for ii in range(crls.shape[0]):
            for jj in range(crls.shape[1]):
                hdr_shift['CRPIX1']=hdr_img['CRPIX1']+dx[ii,jj]
                hdr_shift['CRPIX2']=hdr_img['CRPIX2']+dy[ii,jj]

                img0_shift,coverage=reproject_interp((img0.T,hdr0),hdr_shift,order=interp_order)
                img0_shift=img0_shift.T
                img0_shift=np.nan_to_num(img0_shift)

                if box[0]==-1:
                    mult=img0_shift*img
                else:
                    mult=img0_shift[box[0]:box[1],box[2]:box[3]]*img[box[0]:box[1],box[2]:box[3]]
                if np.sum(mult!=0)>0:
                    crls[ii,jj]=np.sum(mult)/np.sum(mult!=0)

        xplot=np.append(xx,xx[1]-xx[0]+xx[-1])-0.5*(xx[1]-xx[0])
        yplot=np.append(yy,yy[1]-yy[0]+yy[-1])-0.5*(xx[1]-xx[0])
        plt.pcolormesh(xplot,yplot,crls.T,cmap='plasma')

        tmp=np.unravel_index(crls.argmax(),crls.shape)
        xmax=xx[tmp[0]]
        ymax=yy[tmp[1]]
        plt.plot(xmax,ymax,'+',color='r')
        plt.title(os.path.basename(cubefn))
        print(xmax,ymax)

        # test
        #hdr_shift['CRPIX1']=hdr_img['CRPIX1']+xmax
        #hdr_shift['CRPIX2']=hdr_img['CRPIX2']+ymax
        #img0_shift,coverage=reproject_interp((img0.T,hdr0),hdr_shift,order='bilinear')
        #hdu_test=fits.PrimaryHDU(img0_shift,header=hdr_shift)
        #hdu_test.writeto('test.fits',overwrite=True)


        # write plot
        if display:
            fig.tight_layout()
            plt.show()

        fig.savefig('kcwi_astrom/'+cubefn.replace('.fits','.astrom.pdf'))
    else:
        xmax=0
        ymax=0

    # best fit
    hdr_best=hdr_img.copy()
    hdr_best['CRPIX1']=hdr_img['CRPIX1']+xmax
    hdr_best['CRPIX2']=hdr_img['CRPIX2']+ymax
    hdu_best=fits.PrimaryHDU(img_withbkg.T, header=hdr_best)
    hdu_best.writeto('kcwi_astrom/'+cubefn.replace('.fits','.thum.fits'),overwrite=True)

    hdu_cube.header['CRPIX1']+=xmax
    hdu_cube.header['CRPIX2']+=ymax
    hdu_cube.writeto(cubefn.replace('.fits','_wcs.fits'),overwrite=True)

    if save_shift:
        wcs_best=wcs.WCS(hdr_best)
        ra,dec=wcs_best.all_pix2world([hdr_best['CRPIX1'],hdr_best['CRPIX1']+xmax],
                                      [hdr_best['CRPIX2'],hdr_best['CRPIX2']+ymax],1)
        rashift=(dref_ra+ra[0]-ra[1])*3600
        decshift=(dref_dec+dec[0]-dec[1])*3600
        shifttab=table.Table([[rashift],[decshift]])
        ascii.write(shifttab,fnlist.replace('.list','.astrom.list'),overwrite=True,format='no_header')


    if display==False:
        matplotlib.use(oldbackend)
    return
