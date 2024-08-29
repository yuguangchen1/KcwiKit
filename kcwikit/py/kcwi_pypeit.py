import numpy as np
import os.path as path
import os
import shutil
import glob
from astropy.io import ascii
from astropy.io import fits
from astropy import wcs
from astropy import coordinates
from astropy import units
from astropy import time
from astropy import stats
from astropy import table
import pyregion
from reproject import reproject_interp
from reproject import reproject_exact
#from MontagePy.main import mProjectCube
from PyAstronomy import pyasl
from scipy import interpolate
from scipy import signal
from scipy import ndimage
#from image_registration import chi2_shift
#from image_registration.fft_tools import shift
import matplotlib
import matplotlib.pyplot as plt
import pyregion
from fpdf import FPDF
from tqdm import tqdm
import pdb
import time as ostime


# Read parameter files for alignment, stacking, and astrometry
def kcwi_stack_readpar(parname='q0100-bx172.par'):

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
        "background_level":-1000}

    with open(parname,'r') as file:
        lins=file.readlines()
        keys=[]
        for i in range(len(lins)):
            if not lins[i].isspace():
                keys.append(lins[i].split()[0].casefold())


    # Alignment keywords
    q=np.where(np.array(keys)=="align_box")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        tmp=np.array(ele[1:5]).astype(np.float)
        par['align_box']=np.array([tmp[0],tmp[2],tmp[1],tmp[3]]).astype(int)

    q=np.where(np.array(keys)=="align_dimension")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["align_dimension"]=np.array(ele[1:3]).astype(np.int)

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
        par["align_ad"]=np.array(ele[1:3]).astype(np.float)

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
        par["stack_dimension"]=np.array(ele[1:3]).astype(np.int)

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
        par["stack_ad"]=np.array(ele[1:3]).astype(np.float)


    # Global keywords
    q=np.where(np.array(keys)=="dimension")[0]
    if len(q) > 0:
        q=q[0]
        ele=lins[q].split()
        par["stack_dimension"]=np.array(ele[1:3]).astype(np.int)
        par["align_dimension"]=np.array(ele[1:3]).astype(np.int)

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
        par["stack_ad"]=np.array(ele[1:3]).astype(np.float)
        par["align_ad"]=np.array(ele[1:3]).astype(np.float)

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
        par["wavebin"]=np.array(ele[1:3]).astype(np.float)

    # astrometry
    q=np.where(np.array(keys)=="ref_xy")[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par["ref_xy"]=np.array(ele[1:3]).astype(np.float)

    q=np.where(np.array(keys)=="ref_ad")[0]
    if len(q) >0:
        q=q[0]
        ele=lins[q].split()
        par["ref_ad"]=np.array(ele[1:3]).astype(np.float)

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
        tmp=np.array(ele[1:5]).astype(np.float)
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


def kcwi_wave_zeropoint(fnlist, correct_wave = 3226.5, cubed=False):
    tab = ascii.read(fnlist)
    fn = tab['col1']

    if cubed:
        suffix="cubed"
    else:
        suffix="cubes"
    for i in range(len(fn)):
        hdu = fits.open(fn[i]+'_i'+suffix+'.fits')
        hdr = hdu[0].header
        print(f"CRVAL3 (old, new): {hdr['CRVAL3']}, {correct_wave}")
        hdr.set('CRVAL3', correct_wave)
        hdu.writeto(fn[i]+'_i'+suffix+'.fits', overwrite=True)
        hdu.close()


# Perform air-to-vac, and barycentric correction in wavelength
def kcwi_vachelio(hdu,hdr_ref='',mask=False,uncert=False,method='heliocentric'):


    if hdr_ref=='':
        hdr_new=hdu['FLUX'].header.copy()
    else:
        hdr_new=hdu['FLUX'].header.copy()
        hdr_new['NAXIS3']=hdr_ref['NAXIS3']
        hdr_new['CTYPE3']=hdr_ref['CTYPE3']
        hdr_new['CUNIT3']=hdr_ref['CUNIT3']
        hdr_new['CNAME3']=hdr_ref['CNAME3']
        hdr_new['CRVAL3']=hdr_ref['CRVAL3']
        hdr_new['CRPIX3']=hdr_ref['CRPIX3']
        hdr_new['PC3_3']=hdr_ref['PC3_3']


    hdr_new['CTYPE3']='WAVE'

    hdr_old=hdu['FLUX'].header.copy()
    # print(hdr_old)
    if mask==True:
        hdu = hdu['MASK']
        # print('mask=True')
    elif uncert==True:
        hdu = hdu['VARIANCE']
        # print('uncert=True')
    else:
        # print('science frame')
        hdu = hdu['FLUX']

    cube_old=hdu.data.copy()
    cube_old=np.nan_to_num(cube_old)
    shape_old=[hdr_old['NAXIS1'],hdr_old['NAXIS2'],hdr_old['NAXIS3']]

    wcs_old=wcs.WCS(hdr_old)
    wave_old=wcs_old.wcs_pix2world(np.zeros(shape_old[2]),np.zeros(shape_old[2]),np.arange(shape_old[2]),0)
    wave_old=wave_old[2]*1e10

    shape_new=[hdr_new['NAXIS1'],hdr_new['NAXIS2'],hdr_new['NAXIS3']]
    wcs_new=wcs.WCS(hdr_new)
    wave_new=wcs_new.wcs_pix2world(np.zeros(shape_new[2]),np.zeros(shape_new[2]),np.arange(shape_new[2]),0)
    wave_new=wave_new[2]*1e10

    if 'VCORR' in hdr_new:
        if mask==False and uncert==False:
            print('Skipping Wavelength Correction') # python version of DRP already performs correction
        hdu_new=fits.PrimaryHDU(cube_old,header=hdr_new)
        vcorr = hdr_new['VCORR']
        return (hdu_new,vcorr)

    # if skip_wave == True:
    #     hdu_new=fits.PrimaryHDU(cube_old,header=hdr_new)
    #     return hdu_new

    # print(wave_old)

    wave_vac=pyasl.airtovac2(wave_old)

    targ=coordinates.SkyCoord(hdr_old['TARGRA'],hdr_old['TARGDEC'],unit='deg',obstime=hdr_old['DATE-BEG'])
    keck=coordinates.EarthLocation.of_site('Keck Observatory')
    vcorr=targ.radial_velocity_correction(kind=method,location=keck)

    hdr_new['VCORR'] = (vcorr.to('km/s').value, 'Heliocentric Velocity Correction')

    wave_hel=wave_vac*(1+vcorr.value/2.99792458e8)


    #resample
    cube_new=np.zeros(shape_new)
    for i in range(shape_new[0]):
        for j in range(shape_new[1]):
            spec=cube_old[:,j,i]
            if mask==False:
                f_cubic=interpolate.interp1d(wave_hel,spec,kind='cubic',fill_value='extrapolate')
                spec_new=f_cubic(wave_new)
                # testing
                if (spec_new.shape[0]-np.sum(np.isfinite(spec_new)))>0:
                    pdb.set_trace()
            else:
                f_pre=interpolate.interp1d(wave_hel,spec,kind='previous',bounds_error=False,fill_value=128)
                spec_pre=f_pre(wave_new)
                f_nex=interpolate.interp1d(wave_hel,spec,kind='next',bounds_error=False,fill_value=128)
                spec_nex=f_nex(wave_new)

                spec_new=np.zeros(shape_new[2])
                for k in range(shape_new[2]):
                    spec_new[k]=max(spec_pre[k],spec_nex[k])
            cube_new[i,j,:]=spec_new

    hdu_new=fits.PrimaryHDU(cube_new.T,header=hdr_new)

#   plt.clf()
#   plt.plot(wave_old,cube_old[:,45,15],drawstyle='steps-mid')
#   plt.plot(wave_new,cube_new.T[:,45,15],drawstyle='steps-mid')
#   plt.xlim(3300,3500)
#   plt.ylim(-0.01,0.01)
    return (hdu_new,vcorr.to('km/s').value)



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
    # Calculate the relative flux of overlapping sources.

    tab = ascii.read(fnlist)
    fn = tab['col1']

    if thumfn is None:
        thumfn = 'kcwi_align/'+fnlist.replace('.list','.thum.fits')
    hdu_thum = fits.open(thumfn)[0]

    # convert to SB units
    if cubed==False:
        for i in range(hdu_thum.shape[0]):
            hdu_i = fits.open(fn[i]+'_icubes.fits')['FLUX']
            # print(hdu_i.header)
            dx=np.sqrt(hdu_i.header['PC1_1']**2+hdu_i.header['PC2_1']**2)*3600.
            dy=np.sqrt(hdu_i.header['PC1_2']**2+hdu_i.header['PC2_2']**2)*3600.
            area=dx*dy
            #hdu_thum.data[i,:,:] = hdu_thum.data[i,:,:] / area

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

    return


def kcwi_norm_flux(fnlist, frame=[], thumfn=None, nsig=1.5, cubed=False):
    # Generate correction factor.

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
        hdu_i = fits.open(fn[i]+'_i'+suffix+'.fits')['FLUX']
        dx=np.sqrt(hdu_i.header['PC1_1']**2+hdu_i.header['PC2_1']**2)*3600.
        dy=np.sqrt(hdu_i.header['PC1_2']**2+hdu_i.header['PC2_2']**2)*3600.
        area=dx*dy
        #hdu_thum.data[i,:,:] = hdu_thum.data[i,:,:] / area


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




def kcwi_stack(fnlist,shiftlist='',preshiftfn='',fluxfn='',pixscale_x=0.,pixscale_y=0.,
               dimension=[0,0],orientation=-1000.,cubed=False,stepsig=0,drizzle=0,weights=[],
               overwrite=False,keep_trim=True,keep_mont=False,method='drizzle',use_astrom=False,
               use_regmask=True, low_mem=False):
#   fnlist="q0100-bx172.list"
#   shiftlist=""
#   preshiftfn=""
#   pixscale_x=0
#   pixscale_y=0
#   dimension=[0,0]
#   orientation=-1000.
#   cubed=False
#   stepsig=0
#   overwrite=False


    if cubed:
        suffix="cubed"
    else:
        suffix="cubes"

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

    if stepsig==0:
        stepsig=par["stepsig"]

    if drizzle==0:
        drizzle=par["drizzle"]
        if drizzle==0:
            drizzle=0.7

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

    if use_astrom:
        astrom_tab=ascii.read(fnlist.replace('.list','.astrom.list'))
        astrom_rashift=astrom_tab['col1'][0]
        astrom_decshift=astrom_tab['col2'][0]
        overwrite=True
    else:
        astrom_rashift=0.
        astrom_decshift=0.

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
    hdrtmp=hdulist['FLUX'].header.copy()
    hdulist.close()
    wcstmp=wcs.WCS(hdrtmp).copy()
    center=wcstmp.wcs_pix2world((wcstmp.pixel_shape[0]-1)/2.,(wcstmp.pixel_shape[1]-1)/2.,0,0,ra_dec_order=True)

    if par['stack_ad'][0]!=-1:
        center=par['stack_ad']

    hdr0=hdrtmp.copy()
    hdr0['NAXIS1']=dimension[0]
    hdr0['NAXIS2']=dimension[1]
    hdr0['CRPIX1']=(dimension[0]+1)/2.
    hdr0['CRPIX2']=(dimension[1]+1)/2.
    hdr0['CRVAL1']=float(center[0])
    hdr0['CRVAL2']=float(center[1])
    old_cd11=hdr0['PC1_1']
    old_cd12=hdr0['PC1_2']
    old_cd21=hdr0['PC2_1']
    old_cd22=hdr0['PC2_2']
    hdr0['PC1_1']=-pixscale_x
    hdr0['PC2_2']=pixscale_y
    hdr0['PC1_2']=0
    hdr0['PC2_1']=0
    hdr0['CTYPE3']='WAVE'
    #hdr0['BUNIT']='10^(-16)erg/s/cm2/Angstrom'
    hdr0['BUNIT']='10^(-8)erg/s/cm3/arcsec2'
    if suffix!='cubes':
        #hdr0['BUNIT']='adu/s'
        hdr0['BUNIT']='count/s/arcsec2'

    # orientation
    if orientation==-1000:
        orientation=par['stack_orientation']
        if orientation==-1000:
            orientation=np.rad2deg(np.arctan(old_cd21/(-old_cd11)))
    hdr0['PC1_1']=-pixscale_x*np.cos(np.deg2rad(orientation))
    hdr0['PC2_1']=pixscale_x*np.sin(np.deg2rad(orientation))
    hdr0['PC1_2']=pixscale_y*np.sin(np.deg2rad(orientation))
    hdr0['PC2_2']=pixscale_y*np.cos(np.deg2rad(orientation))
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
    # montmfns=[]
    montefns=[]
    etime=np.zeros(len(fn))
    for i in range(len(fn)):
        print(os.path.basename(fn[i]))

        # check availability
        if (not os.path.isfile(trimfn[i])) or overwrite==True:
            # if wave_corr == False:
            #     print('Skipping Heliocentric Wavelength Correction')
            #     # science cube
            #     hdulist=fits.open(fn[i])
            #     hdu_i=kcwi_vachelio(hdulist,hdr_ref=hdr0, skip_wave=True)
            #     hdulist.close()
            #
            #     # variance cube -> sigma cube
            #     hdulist=fits.open(fn[i]) #fits.open(vfn[i])
            #     hdu_v=kcwi_vachelio(hdulist,hdr_ref=hdr0, uncert=True, skip_wave=True)
            #     hdulist.close()
            #
            #     # mask cube
            #     hdulist=fits.open(fn[i]) #fits.open(mfn[i])
            #     hdu_m=kcwi_vachelio(hdulist,hdr_ref=hdr0,mask=True, skip_wave=True)
            #     hdulist.close()
            # else:
                # print('skip_wave == False')
                # science cube
            hdulist=fits.open(fn[i])
            hdu_i,vcorr=kcwi_vachelio(hdulist,hdr_ref=hdr0)
            print('     Vcorr = '+ str(vcorr)+ ' km/s')
            # hdu_i.header['VCORR'] = (vcorr, 'Heliocentric Velocity Correction')
            hdulist.close()

            # variance cube -> sigma cube
            hdulist=fits.open(fn[i]) #fits.open(vfn[i])
            hdu_v,vcorr=kcwi_vachelio(hdulist,hdr_ref=hdr0, uncert=True)
            hdulist.close()

            # mask cube
            # hdulist=fits.open(fn[i]) #fits.open(mfn[i])
            # hdu_m,vcorr=kcwi_vachelio(hdulist,hdr_ref=hdr0,mask=True)
            # hdulist.close()

            # region masks
            regfn = fn[i].replace('.fits','.thum.reg')
            if os.path.isfile(regfn) and use_regmask==True:
                hdr2d=hdu_i.header.copy()
                del hdr2d['PC3_3']
                del hdr2d['CRVAL3']
                del hdr2d['CRPIX3']
                del hdr2d['NAXIS3']
                del hdr2d['CTYPE3']
                del hdr2d['CNAME3']
                del hdr2d['CUNIT3']
                hdr2d['NAXIS']=2

                region = pyregion.open(regfn).as_imagecoord(hdr2d)
                tmp = np.mean(hdu_i.data,axis=0)
                mask_reg = region.get_mask(hdu=fits.PrimaryHDU(tmp,header=hdr2d))

                hdu_i.data[:,mask_reg] = np.nan
                hdu_v.data[:,mask_reg] = np.Inf
                # hdu_m.data[:,mask_reg] = 128

            # Infinity check
            q=((hdu_i.data==0) | (~np.isfinite(hdu_i.data)) | (hdu_v.data==0) | (~np.isfinite(hdu_v.data)) )
            hdu_i.data[q]=np.nan
            hdu_v.data[q]=np.Inf
            # hdu_m.data[q]=128


            # check EXPTIME
            hdu_i=kcwi_checkexptime(hdu_i)
            exptime=hdu_i.header['XPOSURE']
            print('     EXPTIME = '+str(exptime))
            etime[i]=exptime
            edata=hdu_i.data*0.+exptime
            # q=(hdu_m.data != 0)
            # edata[q]=0
            hdu_e=fits.PrimaryHDU(edata,header=hdu_i.header)
            hdu_e.header['BUNIT']='s'


            cube=hdu_i.data.T
            sz=cube.shape
            wcs_cube=wcs.WCS(hdu_i.header)
            wave=wcs_cube.all_pix2world(np.zeros(sz[2]),np.zeros(sz[2]),np.arange(sz[2]),0)
            wave=wave[2]*1e10

            # have to use surface brightness for now, mProjectCube has bug with brightness units combined with drizzle scale
            dx=np.sqrt(hdu_i.header['PC1_1']**2+hdu_i.header['PC2_1']**2)*3600.
            dy=np.sqrt(hdu_i.header['PC1_2']**2+hdu_i.header['PC2_2']**2)*3600.
            area=dx*dy
            if suffix!='cubes': # i.e. suffix == cubed
                # for lam in range(len(wave)):
                #     hdu_i.data[lam,:,:] *= wave[lam]*(0.69/0.29)**2
                #     # print(wave[lam])
                #     hdu_v.data[lam,:,:] *= wave[lam]*(0.69/0.29)**2
                hdu_i.header['BUNIT']='count/s/arcsec2'
                hdu_v.header['BUNIT']='count2/s2/arcsec4'
            #else:
                #hdu_i.data=hdu_i.data/area
                #hdu_v.data=hdu_v.data/area**2
                #hdu_i.header['BUNIT']='10^(-8)erg/s/cm3/arcsec2'
                #hdu_v.header['BUNIT']='10^(-16)erg2/s2/cm6/arcsec4'


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

                    # hdu_m.header['CRVAL1']=hdu_m.header['CRVAL1']+prera[index]/3600.
                    # hdu_m.header['CRVAL2']=hdu_m.header['CRVAL2']+predec[index]/3600.

                    hdu_e.header['CRVAL1']=hdu_e.header['CRVAL1']+prera[index]/3600.
                    hdu_e.header['CRVAL2']=hdu_e.header['CRVAL2']+predec[index]/3600.

            # astrometry correction
            if use_astrom:
                hdu_i.header['CRVAL1']=hdu_i.header['CRVAL1']+astrom_rashift/3600.
                hdu_i.header['CRVAL2']=hdu_i.header['CRVAL2']+astrom_decshift/3600.

                hdu_v.header['CRVAL1']=hdu_v.header['CRVAL1']+astrom_rashift/3600.
                hdu_v.header['CRVAL2']=hdu_v.header['CRVAL2']+astrom_decshift/3600.

                # hdu_m.header['CRVAL1']=hdu_m.header['CRVAL1']+astrom_rashift/3600.
                # hdu_m.header['CRVAL2']=hdu_m.header['CRVAL2']+astrom_decshift/3600.

                hdu_e.header['CRVAL1']=hdu_e.header['CRVAL1']+astrom_rashift/3600.
                hdu_e.header['CRVAL2']=hdu_e.header['CRVAL2']+astrom_decshift/3600.

            # shift
            hdu_i.header['CRPIX1']=hdu_i.header['CRPIX1']+xshift[i]
            hdu_i.header['CRPIX2']=hdu_i.header['CRPIX2']+yshift[i]

            hdu_v.header['CRPIX1']=hdu_v.header['CRPIX1']+xshift[i]
            hdu_v.header['CRPIX2']=hdu_v.header['CRPIX2']+yshift[i]

            # hdu_m.header['CRPIX1']=hdu_m.header['CRPIX1']+xshift[i]
            # hdu_m.header['CRPIX2']=hdu_m.header['CRPIX2']+yshift[i]

            hdu_e.header['CRPIX1']=hdu_e.header['CRPIX1']+xshift[i]
            hdu_e.header['CRPIX2']=hdu_e.header['CRPIX2']+yshift[i]


            # trim
            for kk in range(hdr0['NAXIS3']):
                img=hdu_i.data[kk,:,:]
                var=hdu_v.data[kk,:,:]
                # mask=hdu_m.data[kk,:,:]
                expimg=hdu_e.data[kk,:,:]

                index_y,index_x=np.where(expimg!=0)
                if len(index_y)==0:
                    continue
                xrange=[index_x.min(),index_x.max()]
                yrange=[index_y.min(),index_y.max()]

                if yrange[0]+trim[0,i] >= yrange[1]-trim[1,i]:
                    continue

                img[yrange[1]-trim[1,i]+1:,:]=np.nan
                img[:yrange[0]+trim[0,i],:]=np.nan

                var[yrange[1]-trim[1,i]+1:,:]=np.Inf
                var[:yrange[0]+trim[0,i],:]=np.Inf

                # mask[yrange[1]-trim[1,i]+1:,:]=128
                # mask[:yrange[0]+trim[0,i],:]=128

                expimg[yrange[1]-trim[1,i]+1:,:]=0
                expimg[:yrange[0]+trim[0,i],:]=0

                hdu_i.data[kk,:,:]=img
                hdu_v.data[kk,:,:]=var
                # hdu_m.data[kk,:,:]=mask
                hdu_e.data[kk,:,:]=expimg

            # flux correction
            hdu_i.data = hdu_i.data * fluxnorm[i]
            hdu_v.data = hdu_v.data * fluxnorm[i]**2


            # write
            hdu_i.writeto(trimfn[i],overwrite=True)
            hdu_v.writeto(trimvfn[i],overwrite=True)
            # hdu_m.writeto(trimmfn[i],overwrite=True)
            hdu_e.writeto(trimefn[i],overwrite=True)

        # Montage
        if method.lower()=='drizzle':
            montfn=trimfn[i].replace('.trim.fits','.mont.fits')
            montvfn=trimvfn[i].replace('.trim.fits','.mont.fits')
            # montmfn=trimmfn[i].replace('.trim.fits','.mont.fits')
            montefn=trimefn[i].replace('.trim.fits','.mont.fits')
        else:
            montfn=trimfn[i].replace('.trim.fits','.'+method_flag+'.fits')
            montvfn=trimvfn[i].replace('.trim.fits','.'+method_flag+'.fits')
            # montmfn=trimmfn[i].replace('.trim.fits','.'+method_flag+'.fits')
            montefn=trimefn[i].replace('.trim.fits','.'+method_flag+'.fits')

        montfns.append(montfn)
        montvfns.append(montvfn)
        # montmfns.append(montmfn)
        montefns.append(montefn)


        if (not os.path.isfile(montfn)) or overwrite==True:
        #if True:
            if method.lower()!='drizzle':
                newhdr=fits.Header.fromtextfile(fnhdr)

                hdut=(fits.open(trimfn[i]))[0]
                hdut.data[np.isfinite(hdut.data)==False]=0.
                newi,newa=reproject_interp(hdut,newhdr,order=method,independent_celestial_slices=True)
                hdui=fits.PrimaryHDU(newi,newhdr)
                hdua=fits.PrimaryHDU(newa,newhdr)
                hdui.writeto(montfn,overwrite=True)
                hdua.writeto(montfn.replace('.'+method_flag,'.'+method_flag+'_area'),overwrite=True)

                hdutv=(fits.open(trimvfn[i]))[0]
                hdutv.data[np.isfinite(hdutv.data)==False]=0.
                newv,newa=reproject_interp(hdutv,newhdr,order=method,independent_celestial_slices=True)
                hduv=fits.PrimaryHDU(newv,newhdr)
                hdua=fits.PrimaryHDU(newa,newhdr)
                hduv.writeto(montvfn,overwrite=True)
                hdua.writeto(montvfn.replace('.'+method_flag,'.'+method_flag+'_area'),overwrite=True)

                # hdutm=(fits.open(trimmfn[i]))[0]
                # hdutm.data[np.isfinite(hdutm.data)==False]=0.
                # newm,newa=reproject_interp(hdutm,newhdr,order='bilinear',independent_celestial_slices=True)
                # hdum=fits.PrimaryHDU(newm,newhdr)
                # hdua=fits.PrimaryHDU(newa,newhdr)
                # hdum.writeto(montmfn,overwrite=True)
                # hdua.writeto(montmfn.replace('.'+method_flag,'.'+method_flag+'_area'),overwrite=True)

                hdute=(fits.open(trimefn[i]))[0]
                hdute.data[np.isfinite(hdute.data)==False]=0.
                newe,newa=reproject_interp(hdute,newhdr,order='bilinear',independent_celestial_slices=True)
                hdue=fits.PrimaryHDU(newe,newhdr)
                hdua=fits.PrimaryHDU(newa,newhdr)
                hdue.writeto(montefn,overwrite=True)
                hdua.writeto(montefn.replace('.'+method_flag,'.'+method_flag+'_area'),overwrite=True)

            else:
                # using shell version for now because of memory leakage of MontagePy
                exe="mProjectCube -z "+str(drizzle)+" -f "+trimfn[i]+" "+montfn+" "+fnhdr
                void=os.system(exe)
                exev="mProjectCube -z "+str(drizzle)+" -f  "+trimvfn[i]+" "+montvfn+" "+fnhdr
                voidv=os.system(exev)
                # exem="mProjectCube -z "+str(drizzle)+" -f  "+trimmfn[i]+" "+montmfn+" "+fnhdr
                # voidm=os.system(exem)
                exee="mProjectCube -z "+str(drizzle)+" -f  "+trimefn[i]+" "+montefn+" "+fnhdr
                voide=os.system(exee)

                #void=mProjectCube(trimfn[i],montfn,fnhdr,drizzle=drizzle,energyMode=True,fullRegion=True)
                #voidv=mProjectCube(trimvfn[i],montvfn,fnhdr,drizzle=drizzle,energyMode=True,fullRegion=True)
                #voidm=mProjectCube(trimmfn[i],montmfn,fnhdr,drizzle=drizzle,energyMode=False,fullRegion=True)
                #voide=mProjectCube(trimefn[i],montefn,fnhdr,drizzle=drizzle,energyMode=False,fullRegion=True)


    if low_mem==False:
        # cache all cubes
        data0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.float64).T
        vdata0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.float64).T
        # mdata0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.int16).T+128
        edata0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.float64).T
        for i in range(len(fn)):
            newcube=fits.open(montfns[i])[0].data
            newcube[~np.isfinite(newcube)]=0.
            newcubev=fits.open(montvfns[i])[0].data
            newcubev[~np.isfinite(newcubev)]=0.
            # newcubem=np.ceil(fits.open(montmfns[i])[0].data)
            # newcubem[~np.isfinite(newcubem)]=128
            newcubee=fits.open(montefns[i])[0].data
            newcubee[~np.isfinite(newcubee)]=0.
            data0[i,:,:,:]=newcube
            vdata0[i,:,:,:]=newcubev
            # mdata0[i,:,:,:]=newcubem
            edata0[i,:,:,:]=newcubee
            #data0.append(newcube)
            #vdata0.append(newcubev)
            #mdata0.append(newcubem)
            #edata0.append(newcubee)
            #print(newcube.shape)
            #print(newcubev.shape)
            #print(newcubem.shape)
            #print(newcubee.shape)


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
            # mask=mdata0[:,:,:,ii]
            exp=edata0[:,:,:,ii]
        else:
            #cache columns
            img = np.zeros((len(fn), hdr0['NAXIS3'], dimension[1]))
            var = np.zeros((len(fn), hdr0['NAXIS3'], dimension[1]))
            # mask = np.zeros((len(fn), hdr0['NAXIS3'], dimension[1]))
            exp = np.zeros((len(fn), hdr0['NAXIS3'], dimension[1]))
            for i in range(len(fn)):
                newcube=fits.open(montfns[i])[0].data
                newcube[~np.isfinite(newcube)]=0.
                newcubev=fits.open(montvfns[i])[0].data
                newcubev[~np.isfinite(newcubev)]=0.
                # newcubem=np.ceil(fits.open(montmfns[i])[0].data)
                # newcubem[~np.isfinite(newcubem)]=128
                newcubee=fits.open(montefns[i])[0].data
                newcubee[~np.isfinite(newcubee)]=0.
                img[i,:,:]=newcube[:, :, ii]
                var[i,:,:]=newcubev[:, :, ii]
                # mask[i,:,:]=newcubem[:, :, ii]
                exp[i,:,:]=newcubee[:, :, ii]


        # mask[~np.isfinite(img)]=1
        # mask[~np.isfinite(var)]=1
        # mask[var==0]=1

        # q=(mask==0)
        # if np.sum(q)==0:
            # continue

        weight=np.zeros(var.shape)
        #weight[var!=0]=1/np.abs(var[var!=0])
        fluxweight = 1 / np.repeat(np.repeat(np.array(fluxnorm**2)[:,np.newaxis],
                             hdr0['NAXIS3'],axis=1)[:,:,np.newaxis],dimension[1],axis=2)
        if len(weights)==0:
            weight = exp.copy() * fluxweight
        else:
            weight=np.repeat(np.repeat(np.array(weights)[:,np.newaxis],
                             hdr0['NAXIS3'],axis=1)[:,:,np.newaxis],dimension[1],axis=2).astype(float)
        # weight[~q]=np.nan


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
            # os.remove(montmfns[i])
            os.remove(montefns[i])
            if method.lower()=='drizzle':
                os.remove(montfns[i].replace('mont','mont_area'))
                os.remove(montvfns[i].replace('mont','mont_area'))
                # os.remove(montmfns[i].replace('mont','mont_area'))
                os.remove(montefns[i].replace('mont','mont_area'))
            else:
                os.remove(montfns[i].replace('.'+method_flag,'.'+method_flag+'_area'))
                os.remove(montvfns[i].replace('.'+method_flag,'.'+method_flag+'_area'))
                # os.remove(montmfns[i].replace('.'+method_flag,'.'+method_flag+'_area'))
                os.remove(montefns[i].replace('.'+method_flag,'.'+method_flag+'_area'))

        if keep_trim==False:
            os.remove(trimfn[i])
            os.remove(trimvfn[i])
            # os.remove(trimmfn[i])
            os.remove(trimefn[i])

    # write
    hdr0['CD1_1'] = hdr0['PC1_1']
    hdr0['CD1_2'] = hdr0['PC1_2']
    hdr0['CD2_1'] = hdr0['PC2_1']
    hdr0['CD2_2'] = hdr0['PC2_2']
    hdr0['CD3_3'] = hdr0['PC3_3']*1e10
    hdr0['CUNIT3'] = 'Angstrom'
    hdr0['CRVAL3'] *= 1e10

    del hdr0['PC1_1']
    del hdr0['PC1_2']
    del hdr0['PC2_1']
    del hdr0['PC2_2']
    del hdr0['PC3_3']


    vhdr0=hdr0.copy()
    if suffix=='cubes':
        #vhdr0['BUNIT']='10^(-32)erg2/s2/cm4/Angstrom2'
        vhdr0['BUNIT']='10^(-16)erg2/s2/cm6/arcsec4'
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
    # hdu_m=fits.PrimaryHDU(mdata_3d.T,header=mhdr0)
    # hdu_m.writeto(fnlist.replace('.list','_m'+suffix_all+'.fits'),overwrite=True)
    hdu_e=fits.PrimaryHDU(edata_3d.T,header=ehdr0)
    hdu_e.writeto(fnlist.replace('.list','_e'+suffix_all+'.fits'),overwrite=True)

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
        del hdr_img['PC3_3']
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




def kcwi_align(fnlist,wavebin=[-1.,-1.],box=[-1,-1,-1,-1],pixscale_x=-1.,pixscale_y=-1.,orientation=-1000.,dimension=[-1.,-1.],preshiftfn='',trim=[-1,-1],cubed=False,noalign=False,display=True,search_size=-1000,conv_filter=-1000,upfactor=-1000.,background_subtraction=False,background_level=-1000.,method='interp',
              use_regmask=True):

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
        if wavebin[0]==-1:
            wavebin=[4000.,5000.]

    # define alignment box
    if box[0]==-1:
        box=par['align_box']
        if box[0]==-1:
            box=[25,50,25,40]

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
    hdrtmp=hdulist['FLUX'].header.copy()
    hdulist.close()
    wcstmp=wcs.WCS(hdrtmp).copy()
    center=wcstmp.wcs_pix2world((wcstmp.pixel_shape[0]-1)/2.,(wcstmp.pixel_shape[1]-1)/2.,0,0,ra_dec_order=True)

    if par['align_ad'][0]!=-1:
        center=par['align_ad']

    hdr0=hdrtmp.copy()
    hdr0['NAXIS1']=dimension[0]
    hdr0['NAXIS2']=dimension[1]
    hdr0['CRPIX1']=(dimension[0]+1)/2.
    hdr0['CRPIX2']=(dimension[1]+1)/2.
    hdr0['CRVAL1']=float(center[0])
    hdr0['CRVAL2']=float(center[1])
    old_cd11=hdr0['PC1_1']
    old_cd12=hdr0['PC1_2']
    old_cd21=hdr0['PC2_1']
    old_cd22=hdr0['PC2_2']
    hdr0['PC1_1']=-pixscale_x
    hdr0['PC2_2']=pixscale_y
    hdr0['PC1_2']=0.
    hdr0['PC2_1']=0.
    hdr0['NAXIS']=2
    del hdr0['PC3_3']
    del hdr0['CRVAL3']
    del hdr0['CRPIX3']
    del hdr0['NAXIS3']
    del hdr0['CTYPE3']
    del hdr0['CNAME3']
    del hdr0['CUNIT3']
    del hdr0['CDELT3']

    hdr0['WCSAXES'] = 2

    # orientation
    if orientation==-1000:
        orientation=par['align_orientation']
        if orientation==-1000:
            orientation=np.rad2deg(np.arctan(old_cd21/(-old_cd11)))
    hdr0['PC1_1']=-pixscale_x*np.cos(np.deg2rad(orientation))
    hdr0['PC2_1']=pixscale_x*np.sin(np.deg2rad(orientation))
    hdr0['PC1_2']=pixscale_y*np.sin(np.deg2rad(orientation))
    hdr0['PC2_2']=pixscale_y*np.cos(np.deg2rad(orientation))

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

        # Intensity cube
        hdu=fits.open(fn[i])['FLUX']
        img=hdu.data.T.copy()
        sz=img.shape
        wcs_i=wcs.WCS(hdu.header)
        wave=wcs_i.wcs_pix2world(np.zeros(sz[2]),np.zeros(sz[2]),np.arange(sz[2]),0)
        wave=wave[2]*1e10
        qwave=(wave > wavebin[0]) & (wave < wavebin[1])

        hdr=hdu.header.copy()
        del hdr['PC3_3']
        del hdr['CRVAL3']
        del hdr['CRPIX3']
        del hdr['NAXIS3']
        del hdr['CTYPE3']
        del hdr['CNAME3']
        del hdr['CUNIT3']
        del hdr['CDELT3']

        hdr['NAXIS']=2
        hdr['WCSAXES'] = 2

        # mask?
        regfn = fn[i].replace('.fits','.thum.reg')
        if os.path.isfile(regfn) and use_regmask==True:
            region = pyregion.open(regfn).as_imagecoord(hdr)
            tmp = np.mean(img,axis=2)
            mask_reg = region.get_mask(hdu=fits.PrimaryHDU(tmp.T,header=hdr))
            img[mask_reg.T!=0]=np.nan

        thum=np.zeros((sz[0],sz[1]))
        for ii in range(sz[0]):
            for jj in range(sz[1]):
                q=(img[ii,jj,qwave]!=0) & (np.isfinite(img[ii,jj,qwave])==1)
                if np.sum(q)>0:
                    thum[ii,jj]=np.mean(img[ii,jj,qwave][q])

        # trim
        index_x,index_y=np.where(thum!=0)
        xrange=[index_x.min(),index_x.max()]
        yrange=[index_y.min(),index_y.max()]
        thum[:,yrange[1]-trim[1,i]+1:]=np.nan
        thum[:,:yrange[0]+trim[0,i]]=np.nan
        thum[:xrange[0],:]=np.nan
        thum[xrange[1]:,:]=np.nan
        #thum=np.nan_to_num(thum)

        # preshift
        if preshiftfn!='':
            index=np.where(np.array(prefn)==os.path.basename(fn[i]))
            index=index[0]
            if len(index)>0:
                index=index[0]
                hdr['CRVAL1']=hdr['CRVAL1']+prera[index]/3600.
                hdr['CRVAL2']=hdr['CRVAL2']+predec[index]/3600.
        # print(hdr)
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
                xx=np.linspace(-crls_size,crls_size,2*crls_size+1)
                yy=np.linspace(-crls_size,crls_size,2*crls_size+1)
                dy,dx=np.meshgrid(yy,xx)
                dx=dx.astype(int)
                dy=dy.astype(int)
                crls=np.zeros(dx.shape)
                for ii in range(crls.shape[0]):
                    for jj in range(crls.shape[1]):
                        cut0=img0[box[0]:box[1],box[2]:box[3]]
                        cut=img[box[0]-dx[ii,jj]:box[1]-dx[ii,jj],box[2]-dy[ii,jj]:box[3]-dy[ii,jj]]
                        if background_subtraction:
                            if background_level==-1000:
                                back_val=np.median(cut)
                                back_val0=np.median(cut0)
                            else:
                                back_val=float(background_level)
                                back_val0=back_val
                            cut=cut-back_val
                            cut0=cut0-back_val0
                        else:
                            if background_level!=-1000:
                                cut[cut<background_level]=0
                                cut0[cut0<background_level]=0
                        cut[cut<0]=0
                        cut0[cut0<0]=0
                        mult=cut0*cut
                        if np.sum(mult!=0)>0:
                            crls[ii,jj]=np.sum(mult)/np.sum(mult!=0)

                fig,ax=plt.subplots(figsize=(4,4))
                xplot=np.append(xx,xx[1]-xx[0]+xx[-1])-0.5
                yplot=np.append(yy,yy[1]-yy[0]+yy[-1])-0.5
                ax.pcolormesh(xplot,yplot,crls.T)

                # find closest local maximum
                max_conv=ndimage.filters.maximum_filter(crls,2*conv_filter+1)
                maxima=(crls==max_conv)
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
                # filter out weak ones
                max=np.max(max_conv[xindex,yindex])
                med=np.median(crls)
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
                hdr0_up['PC1_1']=hdr0_up['PC1_1']/upfactor
                hdr0_up['PC2_1']=hdr0_up['PC2_1']/upfactor
                hdr0_up['PC1_2']=hdr0_up['PC1_2']/upfactor
                hdr0_up['PC2_2']=hdr0_up['PC2_2']/upfactor
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
                crls=np.zeros(dx.shape)
                for ii in range(crls.shape[0]):
                    for jj in range(crls.shape[1]):
                        cut0=img0[box[0]*upfactor:box[1]*upfactor,box[2]*upfactor:box[3]*upfactor]
                        cut=img[box[0]*upfactor-dx[ii,jj]:box[1]*upfactor-dx[ii,jj],box[2]*upfactor-dy[ii,jj]:box[3]*upfactor-dy[ii,jj]]
                        if background_subtraction:
                            if background_level==-1000:
                                back_val=np.median(cut)
                                back_val0=np.median(cut0)
                            else:
                                back_val=float(background_level)
                                back_val0=back_val
                            cut0=cut0-back_val0
                            cut=cut-back_val
                        else:
                            if background_level!=-1000:
                                cut[cut<background_level]=0
                                cut0[cut0<background_level]=0
                        cut[cut<0]=0
                        cut0[cut0<0]=0
                        mult=cut0*cut
                        crls[ii,jj]=np.sum(mult)/np.sum(mult!=0)

                #plt.figure(1)
                #plt.clf()
                xplot=(np.append(xx,xx[1]-xx[0]+xx[-1])-0.5)/upfactor+xshift[i]
                yplot=(np.append(yy,yy[1]-yy[0]+yy[-1])-0.5)/upfactor+yshift[i]
                ax.pcolormesh(xplot,yplot,crls.T,cmap='plasma')

                tmp=np.unravel_index(crls.argmax(),crls.shape)
                xshift[i]+=xx[tmp[0]]/upfactor
                yshift[i]+=yy[tmp[1]]/upfactor
                ax.plot(xshift[i],yshift[i],'+',color='r')
                plt.title(os.path.basename(fn[i]))


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




def kcwi_astrometry(fnlist,imgfn='',wavebin=[-1.,-1.],display=True,search_size=-1000,conv_filter=-1000,upfactor=-1000,box=[-1.,-1.,-1.,-1.],nocrl=0,method='drizzle',save_shift=False):

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
        if wavebin[0]==-1:
            wavebin=[4000.,5000.]

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
    qwave=(wave>wavebin[0]) & (wave<wavebin[1])

    hdr_img=hdu_cube.header.copy()
    del hdr_img['PC3_3']
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


    hdu_img0=fits.open(imgfn)[0]
    img0=hdu_img0.data.T
    hdr0=hdu_img0.header


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

                img0_shift,coverage=reproject_interp((img0.T,hdr0),hdr_shift,order='bilinear')
                img0_shift=img0_shift.T
                img0_shift=np.nan_to_num(img0_shift)

                if box[0]==-1:
                    mult=img0_shift*img
                else:
                    mult=img0_shift[box[0]:box[1],box[2]:box[3]]*img[box[0]:box[1],box[2]:box[3]]
                if np.sum(mult!=0)>0:
                    crls[ii,jj]=np.sum(mult)/np.sum(mult!=0)

        fig=plt.figure(1)
        plt.clf()
        xplot=np.append(xx,xx[1]-xx[0]+xx[-1])-0.5
        yplot=np.append(yy,yy[1]-yy[0]+yy[-1])-0.5
        plt.pcolormesh(xplot,yplot,crls.T)

        max_conv=ndimage.filters.maximum_filter(crls,2*conv_filter+1)
        maxima=(crls==max_conv)
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

                img0_shift,coverage=reproject_interp((img0.T,hdr0),hdr_shift,order='bilinear')
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
        fig.savefig('kcwi_astrom/'+cubefn.replace('.fits','.astrom.pdf'))
    else:
        xmax=0
        ymax=0

    # best fit
    hdr_best=hdr_img.copy()
    hdr_best['CRPIX1']=hdr_img['CRPIX1']+xmax
    hdr_best['CRPIX2']=hdr_img['CRPIX2']+ymax
    hdu_best=fits.PrimaryHDU(img.T,header=hdr_best)
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


def kcwi_package_fn(fnlist):
    # Get the names of the relevant files to be released.

    basename = fnlist.replace('.list', '')

    default_fns = [basename + '_icubes_wcs.fits',
                   basename + '_vcubes.fits',
                   basename + '_mcubes.fits',
                   basename + '_ecubes.fits',
                   basename + '.list',
                   basename + '.par',
                   basename + '.pre.list',
                   basename + '.preshift.list',
                   basename + '.flx.list',
                   basename + '.shift.list',
                   basename + '.ipynb',
                   'kcwi_align/' + basename + '.align.pdf',
                   'kcwi_astrom/' + basename + '_icubes.astrom.pdf']

    final_fns = []
    for fn in default_fns:
        if os.path.isfile(fn):
            final_fns.append(fn)

    return final_fns


def kcwi_upload_google(allfiles, googledirid, jsonfn='../py/pydrive/client_secrets.json'):
    # automatically upload to google

    from pydrive.auth import GoogleAuth
    from pydrive.drive import GoogleDrive

    # Authentication
    gauth = GoogleAuth()
    gauth.DEFAULT_SETTINGS['client_config_file'] = jsonfn
    gauth.LocalWebserverAuth() # Creates local webserver and auto handles authentication.

    drive = GoogleDrive(gauth)

    for fn in allfiles:
        # check existing
        existing_f = drive.ListFile({'q': "'{0:s}' in parents and trashed=false".format(googledirid)}).GetList()

        existing_flag = 0
        for f in existing_f:
            if f['title'] == os.path.basename(fn):
                existing_flag = 1
                break
        if existing_flag == 0:
            f = drive.CreateFile({"parents": [{"kind": "drive#fileLink", "id": googledirid}]})

        f.SetContentFile(fn)
        f['title'] = os.path.basename(fn)
        f.Upload()
        print('title: %s, id: %s' % (f['title'], f['id']))

    return
