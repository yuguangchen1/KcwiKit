import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import cosmology
import reproject
import pyregion
import os
import warnings
import pdb
from matplotlib import pyplot as plt
import kcwi_stats


def subcube(hdu,wave,writefn='',box=[-1,-1,-1,-1],pixel_wave=False,pixel_box=True):
    """
    Trim a big cube down to a smaller one.

    Parameters
    ----------
    hdu: HDU object or str
        The input data cube in header-data unit or a string specifying its
        path.

    wave: array_like with 2 elements
        The lower and higher boundaries in the wavelength direction.

    writefn: str, optional
        The file name of the output HDU.

    box: array_like with 4 elements, optional
        Coordinates of the lower-left and upper-right corners of the sub-cube.

    pixel_wave: bool, optional
        Using pixel coordinate in the wavelength direction? Default: False.

    pixel_box: bool, optional
        Using pixel coordinate in the spatial directions? Default: True.

    Returns
    -------
        newhdu: HDU object
            The extracted sub-cube.

    """

    if type(hdu)==type(''):
        tmp=fits.open(hdu)
        hdu=tmp[0]

    shape=hdu.data.shape
    wcs0=wcs.WCS(hdu.header)
    shape0=hdu.data.shape
    wave0=wcs0.wcs_pix2world(np.zeros(shape0[0]),np.zeros(shape0[0]),np.arange(shape0[0]),0)
    wave0=wave0[2]*1e10

    newhdu=hdu.copy()


    if pixel_wave==False:
        qwave=(wave0 >= wave[0]) & (wave0 < wave[1])
    else:
        qwave=(np.arange(shape0[2]) >= wave[0]) & (np.arange(shape0[2]) < wave[1])

    if np.sum(qwave)==0:
        return -1

    newwave=wave0[qwave]
    newhdu.data=newhdu.data[qwave,:,:]

    if np.sum(newhdu.data)==0:
        return -1

    newhdu.header['NAXIS3']=newhdu.data.shape[0]
    newhdu.header['CRPIX3']=1
    newhdu.header['CRVAL3']=newwave[0]

    if box[0]!=-1:
        ra0=wcs0.wcs_pix2world(np.arange(shape0[2]),np.zeros(shape0[2]),np.zeros(shape0[2]),0)
        ra0=ra0[0]
        dec0=wcs0.wcs_pix2world(np.zeros(shape0[1]),np.arange(shape0[1]),np.zeros(shape0[1]),0)
        dec0=dec0[1]

        if pixel_box==False:
            # real RA DEC from WCS

            qra=(ra0 <= box[0]) & (ra0 > box[2])
            qdec=(dec0 >= box[1]) & (dec0 < box[3])
        else:
            # pixel RA DEC
            qra=(np.arange(shape[2]) >= box[0]) & (np.arange(shape[2]) < box[2])
            qdec=(np.arange(shape[1]) >= box[1]) & (np.arange(shape[1])< box[3])

        newra=ra0[qra]
        newdec=dec0[qdec]
        newhdu.data=newhdu.data[:,qdec,:]
        newhdu.data=newhdu.data[:,:,qra]

        newhdu.header['NAXIS1']=newhdu.data.shape[2]
        newhdu.header['NAXIS2']=newhdu.data.shape[1]
        newhdu.header['CRPIX1']=1
        newhdu.header['CRPIX2']=1
        newhdu.header['CRVAL1']=newra[0]
        newhdu.header['CRVAL2']=newdec[0]



    if writefn!='':
        newhdu.writeto(writefn,overwrite=True)

    return newhdu


def collapse_header(hdr):
    """
    Quick wrapper to collapse a 3-D header into a 2-D one.

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



def collapse(hdu,wavebin=[-1.,-1.],usepix=False,var=False,weight=False,usemean=False,usesum=False,writefn='',ignore_blank=False):
    """
    Collapse the cube into a 2-D image whitelight/narrowband image.

    Parameters
    ----------
    hdu: HDU object or string of the file name

    wavebin: array_like (n*2 elements), optional
        The range of which the cube is collapsed into.
        Can be split in n seperate ranges. The final image will be collapsed into one.

    usepix: bool, optional
        Use pixel indices for wavebin?

    var: bool, optional
        variance cube?

    weight: bool, optional
        Output image in weight, instead of variance.

    usemean: bool, optional
        Using mean instead of median.

    ignore_blank: bool, optional
        Ignore blank images without writing files.
    """

    # default wavebin
    tab_grating=np.array(['BL','BM'])
    tab_wave=np.array([500,300])

    if type(hdu)==type(''):
        ofn=hdu
        tmp=fits.open(hdu)
        hdu=tmp[0]
    else:
        ofn=''


    if type(wavebin)==type([]) or type(wavebin)==type(()):
        wavebin=np.array(wavebin)

    if weight==True:
        var=True

    if len(wavebin.shape)==1:
        wavebin=np.array([wavebin])


    # get cube parameters
    wcs0=wcs.WCS(hdu.header)
    shape0=hdu.data.shape
    wave0=wcs0.wcs_pix2world(np.zeros(shape0[0]),np.zeros(shape0[0]),np.arange(shape0[0]),0)
    wave0=wave0[2]*1e10
    cwave=hdu.header['BCWAVE']

    # get pixel indices
    if wavebin[0,0]==-1:
        grat=hdu.header['BGRATNAM']
        qg=(tab_grating==grat)

        if np.sum(qg)==0:
            wrange=[np.max([(cwave-500),3500]),np.min([(cwave+500),5500])]
        else:
            wrange=[np.max([(cwave-tab_wave[qg][0]),3500]),
                            np.min([(cwave+tab_wave[qg][0]),5500])]
            qwave=(wave0>wrange[0]) & (wave0<wrange[1])
    else:
        if usepix==False:
            qwave=np.zeros(wave0.shape[0],dtype=bool)
            for i in range(wavebin.shape[0]):
                qwave_tmp=(wave0>wavebin[i,0]) & (wave0<wavebin[i,1])
                qwave=np.bitwise_or(qwave_tmp,qwave)

        else:
            qwave=np.zeros(wave0.shape[0],dtype=bool)
            windex=np.arange(wave0.shape[0])
            for i in range(wavebin.shape[0]):
                qwave_tmp=(windex>wavebin[i,0]) & (windex<wavebin[i,1])
                qwave=np.bitwise_or(qwave_tmp,qwave)

    # collapse
    if np.sum(qwave)==0:
        if ignore_blank:
            return -1
        else:
            img=np.zeros(hdu.shape[1:3])

    else:
        cube_0=hdu.data.copy()
        cube_0[cube_0==0]=np.nan
        if var==False:
            if usemean:
                img=np.nanmean(cube_0[qwave,:,:],axis=0)
            elif usesum:
                img=np.nansum(cube_0[qwave,:,:],axis=0)
            else:
                img=np.nanmedian(cube_0[qwave,:,:],axis=0)
        else:
            if usemean==False:
                img=np.nanmean(cube_0[qwave,:,:],axis=0)/np.sum(np.isfinite(cube_0[qwave,:,:]),axis=0)
            elif usesum:
                img=np.nansum(cube_0[qwave,:,:],axis=0)
            else:
                img=np.nanmedian(cube_0[qwave,:,:],axis=0)/np.sum(np.isfinite(cube_0[qwave,:,:]),axis=0)

    # convert to weight
    if weight==True:
        img=np.nan_to_num(1./img)


    hdr_img=collapse_header(hdu.header)
    if weight==True:
        hdr_img['BUNIT']='weight'

    hdu_img=fits.PrimaryHDU(img,header=hdr_img)

    if writefn!='':
        hdu_img.writeto(writefn,overwrite=True)

    return hdu_img


def onedspec(hdu,center=None,radius=None,writefn='',maskfn='',sourcemap='',mcubefn='',c_radec=False,
                r_arcsec=False,source_seg=False):
    """
    Extract 1-D spectrum from data cubes.

    Parameters
    ----------
    hdu: HDU object or str
        The input data cube in header-data unit or a string specifying its
        path.

    center: array_like, optional
        Center of the source in pixel position.

    radius: float, optional
        Pixel radius of the source.

    writefn: string, optional
        Filename of the output spectra.

    maskfn: string, optional
        Continuum mask generated by SExtractor.

    sourcemap: string, optional
        Filename of the source map in FITS.

    mcubefn: string, optional
        Name of the mask cube.

    c_radec: bool, optional
        Using RA Dec in decimal degrees, instead of pixel postion, for the center of the source.

    r_arcsec: bool, optional
        Using arcsec in radius, instead of pixel radius.

    var: bool, optional
        Variance cube?

    source_seg: bool, optional
        Accompanying sourcemap. When True, the file name points to the segmented map generated by
        SExtractor.

    """


    if ((center is None) or (radius is None)) and (sourcemap==''):
        print('[Error] Specify extraction area.')
        return -1

    if type(hdu)==type(''):
        ofn=hdu
        tmp=fits.open(hdu)
        hdu=tmp[0]
    else:
        ofn=''

    # cosmology
    cos=cosmology.LambdaCDM(70.,0.3,0.7)

    # 0 - the original HDU
    hdu0=hdu
    wcs0=wcs.WCS(hdu0.header)
    sz=hdu0.data.shape

    dx=np.sqrt(hdu0.header['CD1_1']**2+hdu0.header['CD2_1']**2)*3600
    dy=np.sqrt(hdu0.header['CD1_2']**2+hdu0.header['CD2_2']**2)*3600

    # ra dec
    if c_radec==True:
        center_ad=center
        tmp=wcs0.wcs_world2pix(center[0],center[1],0,0)
        center_pix=[float(tmp[0]),float(tmp[1])]
    else:
        center_pix=center
        tmp=wcs0.wcs_pix2world(center[0],center[1],0,0)
        center_ad=[float(tmp[0]),float(tmp[1])]

    # masking
    #mask_img=np.zeros(hdu0.shape[1:3],dtype=int)
    fitsmask=None
    if maskfn!='':

        fitsmask_hdu=fits.open(maskfn)[0]
        fitsmask=fitsmask_hdu.data.copy()

        # remove the central souce itself
        center_int=np.round(np.flip(center_pix,axis=-1)).astype(int)
        mask_card=fitsmask[center_int[0],center_int[1]]

        if mask_card!=0:
            fitsmask[fitsmask==mask_card]=0

        fitsmask=fitsmask.astype(bool)

        # expand by 1 pix
        tmpmask=fitsmask.copy()
        xindex,yindex=np.where(tmpmask==True)
        for i in range(xindex.shape[0]):
            if xindex[i]-1>=0:
                fitsmask[xindex[i]-1,yindex[i]]=True
            if xindex[i]+1<fitsmask.shape[0]:
                fitsmask[xindex[i]+1,yindex[i]]=True
            if yindex[i]-1>=0:
                fitsmask[xindex[i],yindex[i]-1]=True
            if yindex[i]+1<fitsmask.shape[1]:
                fitsmask[xindex[i],yindex[i]+1]=True

    if fitsmask is not None:
        mask_3d=np.repeat([fitsmask],sz[0],axis=0)
    else:
        mask_3d=np.zeros_like(hdu0.data)

    if mcubefn!='':
        hdu_mcube=fits.open(mcubefn)[0]
        mask_3d=np.bitwise_or(mask_3d,hdu_mcube.data)

    hdu1=hdu0.copy()
    hdu1.data[(mask_3d==True)]=0


    # source map
    if sourcemap=='':
        xx,yy=np.meshgrid(np.arange(sz[2]),np.arange(sz[1]))
        xx=xx-center_pix[0]
        yy=yy-center_pix[1]

        if r_arcsec==True:
            xx=xx*dx
            yy=yy*dy

        rr=np.sqrt(xx**2+yy**2)
        source=(rr<radius).astype(int)
    else:
        tmphdu=fits.open(sourcemap)[0]
        if source_seg==True:
            center_int=np.round(np.flip(center_pix,axis=-1)).astype(int)
            source_card=tmphdu.data[center_int[0],center_int[1]]
            source=(tmphdu.data==source_card).astype(int)
        else:
            source=(tmphdu.data>0).astype(int)

    # extract spec
    data2=hdu1.data.copy()
    data2[data2==0]=np.nan


    spec=np.nansum(data2*np.repeat([source],sz[0],axis=0),axis=(1,2))
    spec=np.nan_to_num(spec)*dx*dy
    #plt.plot(spec,drawstyle='steps-mid')

    # header
    hdr=hdu1.header.copy()
    hdr['NAXIS']=1
    hdr['NAXIS1']=hdr['NAXIS3']
    hdr['CTYPE1']=hdr['CTYPE3']
    hdr['CNAME1']=hdr['CNAME3']
    hdr['CRVAL1']=hdr['CRVAL3']
    hdr['CRPIX1']=hdr['CRPIX3']
    hdr['CUNIT1']=hdr['CUNIT3']
    hdr['CD1_1']=hdr['CD3_3']
    hdr['CDELT1']=hdr['CD3_3']
    hdr['BUNIT']='10^(-8)erg/s/cm3'
    del hdr['NAXIS2']
    del hdr['NAXIS3']
    del hdr['CTYPE2']
    del hdr['CTYPE3']
    del hdr['CNAME2']
    del hdr['CNAME3']
    del hdr['CRVAL2']
    del hdr['CRVAL3']
    del hdr['CRPIX2']
    del hdr['CRPIX3']
    del hdr['CUNIT2']
    del hdr['CUNIT3']
    del hdr['CD1_2']
    del hdr['CD2_1']
    del hdr['CD2_2']
    del hdr['CD3_3']
    del hdr['WCSDIM']
    if center is None:
        hdr['S1D_SFN']=(sourcemap,'Source map')
    else:
        hdr['S1D_SRA']=(center_ad[0],'Source RA')
        hdr['S1D_SDEC']=(center_ad[1],'Source DEC')
        hdr['S1D_SX']=(center_pix[0],'Source X')
        hdr['S1D_SY']=(center_pix[1],'Source Y')
        if r_arcsec==True:
            hdr['S1D_RAS']=(radius,'Radius in arcsec')
        else:
            hdr['S1D_RPIX']=(radius,'Radius in pix')
    hdr['S1D_OFN']=(ofn,'Previous filename')

    hdu_write=fits.PrimaryHDU(spec,header=hdr)

    if writefn!='':
        hdu_write.writeto(writefn,overwrite=True)

    return hdu_write




def cont_sub(hdu,wrange,writefn='',fit_order=1,w_center=None,w_vel=False,auto_reduce=True):
    """
    Conduct continuum-subtraction to cubes.

    Parameters
    ----------
    hdu: HDU object or str
        The input data cube in header-data unit or a string specifying its
        path.

    wrange: N*2 array
        The range of wavelength in Angstrom for the polynomial fit.

    writefn: str, optional
        The file name of the output HDU.

    fit_order: int, optional
        Order of the polynomial fit.

    w_center: float, optional
        If w_vel==True, this specifies the center of the wavelength.

    w_vel: bool, optional
        Use velocity in km/s for the wrange bins, intead of Angstrom.

    auto_reduce: bool, optional
        Automatically reduce the fitting order to 0, if all valid wavelength bins are on
        one side of the central wavelength.

    """

    if w_vel==True and w_center is None:
        print('[Error] Central wavelength required for w_vel=True.')
        return -1

    if type(hdu)==type(''):
        ofn=hdu
        tmp=fits.open(hdu)
        hdu=tmp[0]
    else:
        ofn=''

    data_new=hdu.data.copy()
    data_new=np.nan_to_num(data_new)

    sz=hdu.shape
    wcs0=wcs.WCS(hdu.header)
    wave0=wcs0.wcs_pix2world(np.zeros(sz[0]),np.zeros(sz[0]),np.arange(sz[0]),0)
    wave0=wave0[2]*1e10
    if not (w_center is None):
        v0=(wave0-w_center)/w_center*3e5


    # Get wavelength bins
    if type(wrange)==type([]):
        wrange=np.array(wrange)

    if len(wrange.shape)==1:
        wrange=np.array([wrange])

    sz_wr=wrange.shape

    index0=np.zeros(sz[0],dtype=bool)
    for i in range(sz_wr[0]):
        if w_vel==True:
            index=((v0>=np.min(wrange[i,:])) & (v0<np.max(wrange[i,:])))
        else:
            index=((wave0>np.min(wrange[i,:])) & (wave0<np.max(wrange[i,:])))

        index0=np.bitwise_or(index0,index)



    # Fitting
    cube_cont=np.zeros(data_new.shape)
    for i in range(sz[2]):
        for j in range(sz[1]):
            fit_order_tmp=fit_order
            if not (w_center is None):
                index_r=((wave0[index0] > w_center) & (data_new[index0,j,i]!=0))
                index_l=((wave0[index0] < w_center) & (data_new[index0,j,i]!=0))

                if (np.sum(index_r)==0) or (np.sum(index_r)==0):
                    if auto_reduce==True:
                        fit_order_tmp=0


            fit_func=kcwi_stats.iter_polyfit(wave0[index0],hdu.data[index0,j,i],fit_order_tmp)
            dp=fit_func(wave0)
            cube_cont[:,j,i]=dp

    # subtract
    index_num=(data_new!=0)
    data_new[index_num]=data_new[index_num]-cube_cont[index_num]

    # write file
    hdr_new=hdu.header.copy()
    hdr_new['CSUB']=1
    hdu_new=fits.PrimaryHDU(data_new,header=hdr_new)

    if writefn!='':
        hdu_new.writeto(writefn,overwrite=True)

    return hdu_new









def cart2cyli(hdu,center,writefn='',maskfn='',ellip=1.,pa=0.,nr=None,npa=None,r_range=None,pa_range=[0,360],dr=None,dpa=None,redshift=-1.,drizzle=1.0,c_radec=False,clean=True,compress=False,automask_basename='',segnum=0,exact_area=False,montage=False):
    """
    Resample a cube in cartesian coordinate to cylindrical coordinate (with ellipticity)
    using drizzling algorithm.
    This function can be used to extract 2D spectrum along radial/angular direction.

    Parameters
    ----------
    hdu: HDU object or str
        The input data cube in header-data unit or a string specifying its
        path.

    center: array_like,
        Specify the central of the cylindrical projection, in pixel coordinate.
        Can be switched to [RA, DEC] with 'c_radec' parameter.

    writefn: str, optional
        The file name of the output HDU.

    maskfn: str, optional
        The file name of the DS9 region file that indicates the position of the
        continuum source.

    ellip: float, optional
        The axis-ratio of the series of ellipses to which the cartesian cube will be
        project to.
        Default: 1. - circle.

    pa: float, optional
        The position angle of the *MAJOR* axis that is east of the north direction.
        Default: 0.

    nr: int, optional, CHANGE TO KPC
        Number of pixels of the post-projection cubes in the radial direction.
        Default: The closest integer that makes the size of individual pixels to
        be 0.3 arcsec in the major axis.

    npa: int, optional
        Number of pixels of the post-projection cubes in the angular direction.
        Default: The closest integer that makes the size of individual pixels to
        be 1 degree.

    r_range: array_like with 2 elements, optional
        The boundaries of the radii to be extracted for the major axis.
        Default: 0 to the minimum radius to include all signals in the input cube.

    pa_range: array_like with 2 elements, optional
        The boundaries of the position angles to be extracted.
        Default: [0,360]

    drizzle: float, optional
        The drizzling factor.
        Default: 1.

    c_radec: bool, optional
        Use [RA, DEC] in decimal degrees to specify the center.

    clean: bool, optional
        Clean up the temp files.
        Default: True

    compress: bool, optional
        Remove the axis with only 1 pixel. This will lose the WCS information on
        the corresponding axis, but is convenient for displaying purposes.
        Default: False

    automask_basefn: str, optional
        Automatically apply default *_seg.fits and *.reg masks. This keyword is
        overridden by maskfn keyword. Directory + base name of the *_seg.fits and
        *.reg masks.

    exact_area: bool, optional
        Get the exact 2D-area map. Default is a 1D collapsed approximate.

    montage: bool, optional
        Use Montage? Defalt = True. If set False, use reproject.interp

    Returns
    -------
        (hdu5, ahdu5)
        hdu5: HDU object
            The resampled data cube.
        ahdu5: HDU object
            The resampled area coverage cube.
    """

    if not os.path.exists('kcwi_tools'):
        os.makedirs('kcwi_tools')

    if type(hdu)==type(''):
        ofn=hdu
        tmp=fits.open(hdu)
        hdu=tmp[0]
    else:
        ofn=''

    if montage==False:
        exact_area=True


    # cosmology
    cos=cosmology.LambdaCDM(70.,0.3,0.7)

    # 0 - the original HDU
    hdu0=hdu
    wcs0=wcs.WCS(hdu0.header)
    sz=hdu0.data.shape


    # masking
    mask_3d=np.zeros(hdu0.shape,dtype=bool)
    if maskfn!='':
        if '.reg' in maskfn:
            # make 2-D header
            hdr_tmp=collapse_header(hdu0.header)

            reg=pyregion.open(maskfn)
            reg_mask=reg.get_mask(header=hdr_tmp,shape=hdu0.shape[1:3])

            for i in range(mask_3d.shape[0]):
                mask_3d[i,:,:]=np.bitwise_or(mask_3d[i,:,:],reg_mask)

        if '.fits' in maskfn:
            fitsmask_hdu=fits.open(maskfn)[0]
            fitsmask=fitsmask_hdu.data.copy()

            # remove the central souce itself
            center_int=np.round(np.flip(center,axis=-1)).astype(int)-1
            if segnum==0:
                mask_card=fitsmask[center_int[0],center_int[1]]
            else:
                mask_card=int(segnum)


            if mask_card!=0:
                fitsmask[fitsmask==mask_card]=0

            fitsmask=fitsmask.astype(bool)

            # expand by 1 pix
            tmpmask=fitsmask.copy()
            xindex,yindex=np.where(tmpmask==True)
            for i in range(xindex.shape[0]):
                if xindex[i]-1>=0:
                    fitsmask[xindex[i]-1,yindex[i]]=True
                if xindex[i]+1<fitsmask.shape[0]:
                    fitsmask[xindex[i]+1,yindex[i]]=True
                if yindex[i]-1>=0:
                    fitsmask[xindex[i],yindex[i]-1]=True
                if yindex[i]+1<fitsmask.shape[1]:
                    fitsmask[xindex[i],yindex[i]+1]=True

            all_mask=fitsmask
            for i in range(mask_3d.shape[0]):
                mask_3d[i,:,:]=np.bitwise_or(mask_3d[i,:,:],all_mask)

    elif automask_basename!='':
        regmask_fn=automask_basename+'.reg'
        fitsmask_fn=automask_basename+'_seg.fits'

        all_mask=np.zeros(hdu0.shape[1:3],dtype=int)

        if os.path.isfile(fitsmask_fn):
            fitsmask_hdu=fits.open(fitsmask_fn)[0]
            fitsmask=fitsmask_hdu.data.copy()

            # remove the central souce itself
            center_int=np.round(np.flip(center,axis=-1)).astype(int)-1
            if segnum==0:
                mask_card=fitsmask[center_int[0],center_int[1]]
            else:
                mask_card=segnum

            if mask_card!=0:
                fitsmask[fitsmask==mask_card]=0

            fitsmask=fitsmask.astype(bool)

            # expand by 1 pix
            tmpmask=fitsmask.copy()
            xindex,yindex=np.where(tmpmask==True)
            for i in range(xindex.shape[0]):
                if xindex[i]-1>=0:
                    fitsmask[xindex[i]-1,yindex[i]]=True
                if xindex[i]+1<fitsmask.shape[0]:
                    fitsmask[xindex[i]+1,yindex[i]]=True
                if yindex[i]-1>=0:
                    fitsmask[xindex[i],yindex[i]-1]=True
                if yindex[i]+1<fitsmask.shape[1]:
                    fitsmask[xindex[i],yindex[i]+1]=True

            all_mask=fitsmask


        if os.path.isfile(regmask_fn):

            hdr_tmp=collapse_header(hdu0.header)

            reg=pyregion.open(regmask_fn)
            reg_mask=reg.get_mask(header=hdr_tmp,shape=hdu0.shape[1:3])

            all_mask=np.bitwise_or(all_mask,reg_mask)


        for i in range(mask_3d.shape[0]):
            mask_3d[i,:,:]=np.bitwise_or(mask_3d[i,:,:],all_mask)

    hdu0_mask=hdu0.copy()
    hdu0_mask.data[mask_3d==True]=np.nan



    # Skewing
    hdu1=hdu0_mask.copy()
    hdr1=hdu1.header
    if ellip!=1:
        CD0=np.array([[hdr1['CD1_1'],hdr1['CD1_2']],
            [hdr1['CD2_1'],hdr1['CD2_2']]])
        rot=np.radians(pa)
        ROT=np.array([[np.cos(rot),-np.sin(rot)],
            [np.sin(rot),np.cos(rot)]])
        CD_rot=np.matmul(ROT,CD0)
        CD_shr=CD_rot
        if ellip<1:
            CD_shr[0,0]=CD_rot[0,0]/ellip
            CD_shr[0,1]=CD_rot[0,1]/ellip
        elif ellip>1:
            CD_shr[0,0]=CD_rot[0,0]*ellip
            CD_shr[0,1]=CD_rot[0,1]*ellip

        ROT=np.array([[np.cos(-rot),-np.sin(-rot)],
            [np.sin(-rot),np.cos(-rot)]])
        CD1=np.matmul(ROT,CD_shr)

        hdr1['CD1_1']=CD1[0,0]
        hdr1['CD1_2']=CD1[0,1]
        hdr1['CD2_1']=CD1[1,0]
        hdr1['CD2_2']=CD1[1,1]

    #cube1fn='kcwi_tools/cart2cyli_cube1.fits'
    #hdu1.writeto(cube1fn,overwrite=True)

    # Shift hdu to the south pole
    wcs1=wcs.WCS(hdr1)
    hdu2=hdu1.copy()
    hdr2=hdu2.header
    if c_radec==True:
        tmp=wcs1.wcs_world2pix(center[0],center[1]+0.3/3600.,0,0)
        ref_pix=[float(tmp[0]),float(tmp[1])]

        center_ad=center
        tmp=wcs0.wcs_world2pix(center[0],center[1],0,0)
        center_pix=[float(tmp[0]),float(tmp[1])]
    else:
        center_pix=center
        tmp=wcs0.all_pix2world(center[0],center[1],0,0)
        center_ad=[float(tmp[0]),float(tmp[1])]
        tmp=wcs1.all_world2pix(tmp[0],tmp[1]+0.3/3600.,0,0)
        ref_pix=[float(tmp[0]),float(tmp[1])]


    hdr2['CRPIX1']=ref_pix[0]+1
    hdr2['CRPIX2']=ref_pix[1]+1
    hdr2['CRVAL1']=0.
    hdr2['CRVAL2']=-90+0.3/3600.
    wcs2=wcs.WCS(hdr2)

    cube2fn='kcwi_tools/cart2cyli_cube2.fits'
    hdu2.writeto(cube2fn,overwrite=True)


    # Calculate the x dimension of the final cube
    if (npa==None) and (dpa==None):
        dx0=1.
        nx0=int(np.round((pa_range[1]-pa_range[0])/dx0))
        pa_range[1]=pa_range[0]+dx0*nx0
    else:
        if npa!=None:
            nx0=int(npa)
            dx0=(pa_range[1]-pa_range[0])/nx0
        if dpa!=None:
            nx0=int(np.round(pa_range[1]-pa_range[0]/dpa))
            dx0=dpa
            pa_range[1]=pa_range[0]+dx0*nx0

    # Split too large pixels
    if dx0>1:
        nx=int(dx0)*nx0
        dx=(pa_range[1]-pa_range[0])/nx
    else:
        nx=nx0
        dx=dx0
    # Montage can't handle DPA>180, splitting. Also need padding to avoid edge effect.
    nx3=np.round(nx/3)
    nx3=np.array([nx3,nx3,nx-2*nx3])
    xr3=np.zeros((2,3))
    xr3[0,0]=pa_range[1]
    xr3[1,0]=pa_range[1]-nx3[0]*dx
    xr3[0,1]=xr3[1,0]
    xr3[1,1]=pa_range[1]-(nx3[0]+nx3[1])*dx
    xr3[0,2]=xr3[1,1]
    xr3[1,2]=pa_range[1]-np.sum(nx3)*dx
    nx3=nx3+10
    xr3[0,:]=xr3[0,:]+5*dx
    xr3[1,:]=xr3[1,:]-5*dx


    # Calculate the y dimension
    if r_range==None:
        p_corner_x=np.array([0,0,sz[2]-1,sz[2]-1])
        p_corner_y=np.array([0,sz[1]-1,sz[1]-1,0])
        p_corner_z=np.array([0,0,0,0])
        tmp=wcs2.wcs_pix2world(p_corner_x,p_corner_y,p_corner_z,0)
        a_corner=tmp[0]
        d_corner=tmp[1]
        r_corner=np.abs(d_corner+90)*3600.
        r_range=[0,np.max(r_corner)]
    if (nr==None) and (dr==None):
        dy=0.3
        ny=int(np.round((r_range[1]-r_range[0])/dy))
        r_range[1]=r_range[0]+dy*ny
    else:
        if nr!=None:
            ny=int(nr)
            dy=(r_range[1]-r_range[0])/ny
        if dr!=None:
            ny=int(np.round((r_range[1]-r_range[0])/dr))
            dy=dr
            r_range[1]=r_range[0]+dy*ny


    # Set up headers
    hdr3_1=hdu2.header.copy()
    hdr3_1['NAXIS1']=int(nx3[0])
    hdr3_1['NAXIS2']=ny
    hdr3_1['CTYPE1']='RA---CAR'
    hdr3_1['CTYPE2']='DEC--CAR'
    hdr3_1['CUNIT1']='deg'
    hdr3_1['CUNIT2']='deg'
    hdr3_1['CRVAL1']=xr3[0,0]
    hdr3_1['CRVAL2']=0.
    hdr3_1['CRPIX1']=0.5
    hdr3_1['CRPIX2']=(90.-r_range[0])/(dy/3600.)+0.5
    hdr3_1['CD1_1']=-dx
    hdr3_1['CD2_1']=0.
    hdr3_1['CD1_2']=0.
    hdr3_1['CD2_2']=dy/3600.
    hdr3_1['LONPOLE']=180.
    hdr3_1['LATPOLE']=0.
    hdr3_1fn='kcwi_tools/cart2cyli_hdr3_1.hdr'
    hdr3_1.totextfile(hdr3_1fn,overwrite=True)

    hdr3_2=hdr3_1.copy()
    hdr3_2['NAXIS1']=int(nx3[1])
    hdr3_2['CRVAL1']=xr3[0,1]
    hdr3_2fn='kcwi_tools/cart2cyli_hdr3_2.hdr'
    hdr3_2.totextfile(hdr3_2fn,overwrite=True)

    hdr3_3=hdr3_1.copy()
    hdr3_3['NAXIS1']=int(nx3[2])
    hdr3_3['CRVAL1']=xr3[0,2]
    hdr3_3fn='kcwi_tools/cart2cyli_hdr3_3.hdr'
    hdr3_3.totextfile(hdr3_3fn,overwrite=True)


    if montage==True:
            # montage
            cube3_1fn='kcwi_tools/cart2cyli_cube3_1.fits'
            cube3_2fn='kcwi_tools/cart2cyli_cube3_2.fits'
            cube3_3fn='kcwi_tools/cart2cyli_cube3_3.fits'
            exe='mProjectCube -z '+str(drizzle)+' -f '+cube2fn+' '+cube3_1fn+' '+hdr3_1fn
            void=os.system(exe)
            exe='mProjectCube -z '+str(drizzle)+' -f '+cube2fn+' '+cube3_2fn+' '+hdr3_2fn
            void=os.system(exe)
            exe='mProjectCube -z '+str(drizzle)+' -f '+cube2fn+' '+cube3_3fn+' '+hdr3_3fn
            void=os.system(exe)

            # Merge
            hdu3_1=fits.open(cube3_1fn)[0]
            hdu3_2=fits.open(cube3_2fn)[0]
            hdu3_3=fits.open(cube3_3fn)[0]

            data4=np.zeros((hdu3_1.shape[0],hdu3_1.shape[1],
                hdu3_1.shape[2]+hdu3_2.shape[2]+hdu3_3.shape[2]-30))
            data4[:,:,0:hdu3_1.shape[2]-10]=hdu3_1.data[:,:,5:hdu3_1.shape[2]-5]
            data4[:,:,hdu3_1.shape[2]-10:
                    hdu3_1.shape[2]+hdu3_2.shape[2]-20]=hdu3_2.data[:,:,5:hdu3_2.shape[2]-5]
            data4[:,:,hdu3_1.shape[2]+hdu3_2.shape[2]-20:
                    hdu3_1.shape[2]+hdu3_2.shape[2]+hdu3_3.shape[2]-30]=hdu3_3.data[:,:,5:hdu3_3.shape[2]-5]
            data4[data4==0]=np.nan
            #hdutmp=fits.PrimaryHDU(data4)
            #hdutmp.writeto('kcwi_tools/tmp.fits',overwrite=True)

            # Area
            area3_1fn=cube3_1fn.replace('.fits','_area.fits')
            area3_2fn=cube3_2fn.replace('.fits','_area.fits')
            area3_3fn=cube3_3fn.replace('.fits','_area.fits')
            ahdu3_1=fits.open(area3_1fn)[0]
            ahdu3_2=fits.open(area3_2fn)[0]
            ahdu3_3=fits.open(area3_3fn)[0]

            area4=np.zeros((hdu3_1.shape[1],
                hdu3_1.shape[2]+hdu3_2.shape[2]+hdu3_3.shape[2]-30))
            area4[:,0:hdu3_1.shape[2]-10]=ahdu3_1.data[:,5:hdu3_1.shape[2]-5]
            area4[:,hdu3_1.shape[2]-10:
                    hdu3_1.shape[2]+hdu3_2.shape[2]-20]=ahdu3_2.data[:,5:hdu3_2.shape[2]-5]
            area4[:,hdu3_1.shape[2]+hdu3_2.shape[2]-20:
                    hdu3_1.shape[2]+hdu3_2.shape[2]+hdu3_3.shape[2]-30]=ahdu3_3.data[:,5:hdu3_3.shape[2]-5]
            if exact_area:
                area4=np.repeat([area4],data4.shape[0],axis=0)
                area4[data4==0]=0
                area4[~np.isfinite(data4)]=0
            #hdutmp=fits.PrimaryHDU(area4)
            #hdutmp.writeto('kcwi_tools/tmp.fits',overwrite=True)
    else:
        cube3_1,area3_1=reproject.reproject_interp(hdu2,hdr3_1)
        cube3_2,area3_2=reproject.reproject_interp(hdu2,hdr3_2)
        cube3_3,area3_3=reproject.reproject_interp(hdu2,hdr3_3)

        data4=np.zeros((cube3_1.shape[0],cube3_1.shape[1],
                cube3_1.shape[2]+cube3_2.shape[2]+cube3_3.shape[2]-30))
        data4[:,:,0:cube3_1.shape[2]-10]=cube3_1[:,:,5:cube3_1.shape[2]-5]
        data4[:,:,cube3_1.shape[2]-10:cube3_1.shape[2]+cube3_2.shape[2]-20]=cube3_2[:,:,5:cube3_2.shape[2]-5]
        data4[:,:,cube3_1.shape[2]+cube3_2.shape[2]-20:cube3_1.shape[2]+cube3_2.shape[2]+cube3_3.shape[2]-30]=cube3_3[:,:,5:cube3_3.shape[2]-5]
        data4[data4==0]=np.nan

        area4=np.zeros((area3_1.shape[0],area3_1.shape[1],
                area3_1.shape[2]+area3_2.shape[2]+area3_3.shape[2]-30))
        area4[:,:,0:area3_1.shape[2]-10]=area3_1[:,:,5:area3_1.shape[2]-5]
        area4[:,:,area3_1.shape[2]-10:area3_1.shape[2]+area3_2.shape[2]-20]=area3_2[:,:,5:area3_2.shape[2]-5]
        area4[:,:,area3_1.shape[2]+area3_2.shape[2]-20:area3_1.shape[2]+area3_2.shape[2]+area3_3.shape[2]-30]=area3_3[:,:,5:area3_3.shape[2]-5]
        area4=np.nan_to_num(area4)
    area4[~np.isfinite(data4)]=0

    # Resample to final grid
    hdr5=hdr3_1.copy()
    hdr5['NAXIS1']=nx0
    hdr5['CTYPE1']='PA'
    hdr5['CTYPE2']='Radius'
    hdr5['CNAME1']='PA'
    hdr5['CNAME2']='Radius'
    hdr5['CRVAL1']=pa_range[1]
    hdr5['CRVAL2']=r_range[0]
    hdr5['CRPIX2']=0.5
    hdr5['CUNIT2']='arcsec'
    hdr5['CD1_1']=-dx0
    hdr5['CD2_2']=dy
    hdr5['C2C_ORA']=(center_ad[0],'RA of origin')
    hdr5['C2C_ODEC']=(center_ad[1],'DEC of origin')
    hdr5['C2C_OX']=(center_pix[0]+1,'X of origin')
    hdr5['C2C_OY']=(center_pix[1]+1,'Y of origin')
    hdr5['C2C_E']=(ellip,'Axis raio')
    hdr5['C2C_EPA']=(pa,'PA of major-axis')
    hdr5['C2C_DRIZ']=(drizzle,'Drizzle factor')
    hdr5['C2C_MTD']=('montage','Algorithm')
    if ofn!='':
        hdr5['C2C_OFN']=(ofn,'Previous filename')
    # 2nd WCS
    hdr5['CTYPE1A']='PA'
    hdr5['CTYPE2A']='Radius'
    hdr5['CTYPE3A']=hdr5['CTYPE3']
    hdr5['CNAME1A']='PA'
    hdr5['CNAME2A']='Radius'
    hdr5['CNAME3A']=hdr5['CNAME3']
    hdr5['CRVAL1A']=pa_range[1]
    hdr5['CRVAL2A']=r_range[0]
    hdr5['CRVAL3A']=hdr5['CRVAL3']
    hdr5['CRPIX1A']=hdr5['CRPIX1']
    hdr5['CRPIX2A']=0.5
    hdr5['CRPIX3A']=hdr5['CRPIX3']
    hdr5['CUNIT1A']=hdr5['CUNIT1']
    hdr5['CUNIT2A']=hdr5['CUNIT2']
    hdr5['CUNIT3A']=hdr5['CUNIT3']
    hdr5['CD1_1A']=-dx0
    hdr5['CD2_2A']=dy
    hdr5['CD3_3A']=hdr5['CD3_3']


    if redshift>0:
        a_dis=(cos.arcsec_per_kpc_proper(redshift)).value

        hdr5['CRVAL2A']=r_range[0]/a_dis
        hdr5['CUNIT2A']='kpc'
        hdr5['CD2_2A']=dy/a_dis

        hdr5['CRVAL3A']=hdr5['CRVAL3']/(1+redshift)
        hdr5['CD3_3A']=hdr5['CD3_3']/(1+redshift)



    ahdr5=hdr5.copy()
    if ~exact_area:
        ahdr5['NAXIS']=2
        del ahdr5['NAXIS3']
        del ahdr5['CTYPE3']
        del ahdr5['CUNIT3']
        del ahdr5['CNAME3']
        del ahdr5['CRVAL3']
        del ahdr5['CRPIX3']
        del ahdr5['CD3_3']
        del ahdr5['CTYPE3A']
        del ahdr5['CUNIT3A']
        del ahdr5['CNAME3A']
        del ahdr5['CRVAL3A']
        del ahdr5['CRPIX3A']
        del ahdr5['CD3_3A']

    data5=np.zeros((data4.shape[0],data4.shape[1],nx0))
    if exact_area:
        area5=np.zeros((data4.shape[0],data4.shape[1],nx0))
    else:
        area5=np.zeros((data4.shape[1],nx0))
    ratio=int(dx0/dx)
    if ratio!=1:
        #warnings.filterwarnings('ignore')
        for i in range(nx0):
            tmp=data4[:,:,i*ratio:(i+1)*ratio]
            data5[:,:,i]=np.nanmean(tmp,axis=2)
            if exact_area:
                tmp=area4[:,:,i*ratio:(i+1)*ratio]
                area5[:,:,i]=np.sum(tmp,axis=2)
            else:
                tmp=area4[:,i*ratio:(i+1)*ratio]
                area5[:,i]=np.sum(tmp,axis=1)
        #warnings.resetwarnings()


    # Compress
    if compress==True:
        if nx0==1:
            tmp=hdr5.copy()
            hdr5['NAXIS']=2
            hdr5['NAXIS1']=tmp['NAXIS3']
            hdr5['NAXIS2']=tmp['NAXIS2']
            del hdr5['NAXIS3']
            hdr5['CTYPE1']=tmp['CTYPE3']
            hdr5['CTYPE2']=tmp['CTYPE2']
            hdr5['CTYPE1A']=tmp['CTYPE3A']
            hdr5['CTYPE2A']=tmp['CTYPE2A']
            del hdr5['CTYPE3']
            del hdr5['CTYPE3A']
            hdr5['CUNIT1']=tmp['CUNIT3']
            hdr5['CUNIT2']=tmp['CUNIT2']
            hdr5['CUNIT1A']=tmp['CUNIT3A']
            hdr5['CUNIT2A']=tmp['CUNIT2A']
            del hdr5['CUNIT3']
            del hdr5['CUNIT3A']
            hdr5['CNAME1']=tmp['CNAME3']
            hdr5['CNAME2']=tmp['CNAME2']
            hdr5['CNAME1A']=tmp['CNAME3A']
            hdr5['CNAME2A']=tmp['CNAME2A']
            del hdr5['CNAME3']
            del hdr5['CNAME3A']
            hdr5['CRVAL1']=tmp['CRVAL3']
            hdr5['CRVAL2']=tmp['CRVAL2']
            hdr5['CRVAL1A']=tmp['CRVAL3A']
            hdr5['CRVAL2A']=tmp['CRVAL2A']
            del hdr5['CRVAL3']
            del hdr5['CRVAL3A']
            hdr5['CRPIX1']=tmp['CRPIX3']
            hdr5['CRPIX2']=tmp['CRPIX2']
            hdr5['CRPIX1A']=tmp['CRPIX3A']
            hdr5['CRPIX2A']=tmp['CRPIX2A']
            del hdr5['CRPIX3']
            del hdr5['CRPIX3A']
            hdr5['CD1_1']=tmp['CD3_3']
            hdr5['CD2_2']=tmp['CD2_2']
            hdr5['CD1_1A']=tmp['CD3_3A']
            hdr5['CD2_2A']=tmp['CD2_2A']
            del hdr5['CD3_3']
            del hdr5['CD3_3A']
            data5=np.transpose(np.squeeze(data5,axis=2))

            if exact_area:
                ahdr5=hdr5.copy()
                area5=np.transpose(np.squeeze(area5,axis=2))
            else:
                tmp=ahdr5.copy()
                ahdr5['NAXIS']=1
                ahdr5['NAXIS1']=tmp['NAXIS2']
                del ahdr5['NAXIS2']
                ahdr5['CTYPE1']=tmp['CTYPE2']
                ahdr5['CTYPE1A']=tmp['CTYPE2A']
                del ahdr5['CTYPE2']
                del ahdr5['CTYPE2A']
                ahdr5['CUNIT1']=tmp['CUNIT2']
                ahdr5['CUNIT1A']=tmp['CUNIT2A']
                del ahdr5['CUNIT2']
                del ahdr5['CUNIT2A']
                ahdr5['CNAME1']=tmp['CNAME2']
                ahdr5['CNAME1A']=tmp['CNAME2A']
                del ahdr5['CNAME2']
                del ahdr5['CNAME2A']
                ahdr5['CRVAL1']=tmp['CRVAL2']
                ahdr5['CRVAL1A']=tmp['CRVAL2A']
                del ahdr5['CRVAL2']
                del ahdr5['CRVAL2A']
                ahdr5['CRPIX1']=tmp['CRPIX2']
                ahdr5['CRPIX1A']=tmp['CRPIX2A']
                del ahdr5['CRPIX2']
                del ahdr5['CRPIX2A']
                ahdr5['CD1_1']=tmp['CD2_2']
                ahdr5['CD1_1A']=tmp['CD2_2A']
                del ahdr5['CD2_2']
                del ahdr5['CD2_2A']
                area5=np.transpose(np.squeeze(area5))

        elif ny==1:
            tmp=hdr5.copy()
            hdr5['NAXIS']=2
            hdr5['NAXIS1']=tmp['NAXIS3']
            hdr5['NAXIS2']=tmp['NAXIS1']
            del hdr5['NAXIS3']
            hdr5['CTYPE1']=tmp['CTYPE3']
            hdr5['CTYPE2']=tmp['CTYPE1']
            hdr5['CTYPE1A']=tmp['CTYPE3A']
            hdr5['CTYPE2A']=tmp['CTYPE1A']
            del hdr5['CTYPE3']
            del hdr5['CTYPE3A']
            hdr5['CUNIT1']=tmp['CUNIT3']
            hdr5['CUNIT2']=tmp['CUNIT1']
            hdr5['CUNIT1A']=tmp['CUNIT3A']
            hdr5['CUNIT2A']=tmp['CUNIT1A']
            del hdr5['CUNIT3']
            del hdr5['CUNIT3A']
            hdr5['CNAME1']=tmp['CNAME3']
            hdr5['CNAME2']=tmp['CNAME1']
            hdr5['CNAME1A']=tmp['CNAME3A']
            hdr5['CNAME2A']=tmp['CNAME1A']
            del hdr5['CNAME3']
            del hdr5['CNAME3A']
            hdr5['CRVAL1']=tmp['CRVAL3']
            hdr5['CRVAL2']=tmp['CRVAL1']
            hdr5['CRVAL1A']=tmp['CRVAL3A']
            hdr5['CRVAL2A']=tmp['CRVAL1A']
            del hdr5['CRVAL3']
            del hdr5['CRVAL3A']
            hdr5['CRPIX1']=tmp['CRPIX3']
            hdr5['CRPIX2']=tmp['CRPIX1']
            hdr5['CRPIX1A']=tmp['CRPIX3A']
            hdr5['CRPIX2A']=tmp['CRPIX1A']
            del hdr5['CRPIX3']
            del hdr5['CRPIX3A']
            hdr5['CD1_1']=tmp['CD3_3']
            hdr5['CD2_2']=tmp['CD1_1']
            hdr5['CD1_1A']=tmp['CD3_3A']
            hdr5['CD2_2A']=tmp['CD1_1A']
            del hdr5['CD3_3']
            del hdr5['CD3_3A']
            data5=np.transpose(np.squeeze(data5,axis=1))


    hdu5=fits.PrimaryHDU(data5,header=hdr5)
    ahdu5=fits.PrimaryHDU(area5,header=ahdr5)



    if writefn!='':
        hdu5.writeto(writefn,overwrite=True)
        ahdu5.writeto(writefn.replace('.fits','_area.fits'),overwrite=True)

    if clean==True:
        os.system('rm -f kcwi_tools/cart2cyli*')

    return (hdu5,ahdu5)



def cart2cyli_pix(hdu,center,writefn='',maskfn='',ellip=1.,pa=0.,nr=None,npa=None,r_range=None,pa_range=[0,360],
                 dr=None,dpa=None,redshift=-1,c_radec=False,compress=False,automask_basename=''):
    """
    Resample a cube in cartesian coordinate to cylindrical coordinate (with ellipticity)
    using drizzling algorithm.
    This function can be used to extract 2D spectrum along radial/angular direction.

    Parameters
    ----------
    hdu: HDU object or str
        The input data cube in header-data unit or a string specifying its
        path.

    center: array_like,
        Specify the central of the cylindrical projection, in pixel coordinate.
        Can be switched to [RA, DEC] with 'c_radec' parameter.

    writefn: str, optional
        The file name of the output HDU.

    maskfn: str, optional
        The file name of the DS9 region file that indicates the position of the
        continuum source.

    ellip: float, optional
        The axis-ratio of the series of ellipses to which the cartesian cube will be
        project to.
        Default: 1. - circle.

    pa: float, optional
        The position angle of the *MAJOR* axis that is east of the north direction.
        Default: 0.

    nr: int, optional
        Number of pixels of the post-projection cubes in the radial direction.
        Default: The closest integer that makes the size of individual pixels to
        be 0.3 arcsec in the major axis.

    npa: int, optional
        Number of pixels of the post-projection cubes in the angular direction.
        Default: The closest integer that makes the size of individual pixels to
        be 1 degree.

    r_range: array_like with 2 elements, optional
        The boundaries of the radii to be extracted for the major axis.
        Default: 0 to the minimum radius to include all signals in the input cube.

    pa_range: array_like with 2 elements, optional
        The boundaries of the position angles to be extracted.
        Default: [0,360]

    c_radec: bool, optional
        Use [RA, DEC] in decimal degrees to specify the center.

    compress: bool, optional
        Remove the axis with only 1 pixel. This will lose the WCS information on
        the corresponding axis, but is convenient for displaying purposes.
        Default: False

    automask_basefn: str, optional
        Automatically apply default *_seg.fits and *.reg masks. This keyword is
        overridden by maskfn keyword. Directory + base name of the *_seg.fits and
        *.reg masks.

    Returns
    -------
        (hdu5, ahdu5)
        hdu5: HDU object
            The resampled data cube.
        ahdu5: HDU object
            The resampled area coverage cube.
    """

    if not os.path.exists('kcwi_tools/'):
        os.makedirs('kcwi_tools')

    if type(hdu)==type(''):
        ofn=hdu
        tmp=fits.open(hdu)
        hdu=tmp[0]
    else:
        ofn=''

    # cosmology
    cos=cosmology.LambdaCDM(70.,0.3,0.7)

    # 0 - the original HDU
    hdu0=hdu
    wcs0=wcs.WCS(hdu0.header)
    sz=hdu0.data.shape

    # masking
    mask_3d=np.zeros(hdu0.shape,dtype=bool)
    if maskfn!='':
        if '.reg' in maskfn:
            # make 2-D header
            hdr_tmp=collapse_header(hdu0.header)

            reg=pyregion.open(maskfn)
            reg_mask=reg.get_mask(header=hdr_tmp,shape=hdu0.shape[1:3])

            for i in range(mask_3d.shape[0]):
                mask_3d[i,:,:]=np.bitwise_or(mask_3d[i,:,:],reg_mask)

        if '.fits' in maskfn:
            hdu_mask=fits.open(maskfn)[0]
            mask=bool(hdu_mask.data)
            for i in range(mask_3d.shape[0]):
                mask_3d[i,:,:]=np.bitwise_or(mask_3d[i,:,:],mask)

    elif automask_basename!='':
        regmask_fn=automask_basename+'.reg'
        fitsmask_fn=automask_basename+'_seg.fits'

        all_mask=np.zeros(hdu0.shape[1:3],dtype=int)

        if os.path.isfile(fitsmask_fn):
            fitsmask_hdu=fits.open(fitsmask_fn)[0]
            fitsmask=fitsmask_hdu.data.copy()

            # remove the central souce itself
            center_int=np.round(np.flip(center,axis=-1)).astype(int)-1
            mask_card=fitsmask[center_int[0],center_int[1]]

            if mask_card!=0:
                fitsmask[fitsmask==mask_card]=0

            fitsmask=fitsmask.astype(bool)

            # expand by 1 pix
            tmpmask=fitsmask.copy()
            xindex,yindex=np.where(tmpmask==True)
            for i in range(xindex.shape[0]):
                if xindex[i]-1>=0:
                    fitsmask[xindex[i]-1,yindex[i]]=True
                if xindex[i]+1<fitsmask.shape[0]:
                    fitsmask[xindex[i]+1,yindex[i]]=True
                if yindex[i]-1>=0:
                    fitsmask[xindex[i],yindex[i]-1]=True
                if yindex[i]+1<fitsmask.shape[1]:
                    fitsmask[xindex[i],yindex[i]+1]=True

            all_mask=fitsmask


        if os.path.isfile(regmask_fn):

            hdr_tmp=collapse_header(hdu0.header)

            reg=pyregion.open(regmask_fn)
            reg_mask=reg.get_mask(header=hdr_tmp,shape=hdu0.shape[1:3])

            all_mask=np.bitwise_or(all_mask,reg_mask)


        for i in range(mask_3d.shape[0]):
            mask_3d[i,:,:]=np.bitwise_or(mask_3d[i,:,:],all_mask)

    hdu0_mask=hdu0.copy()
    hdu0_mask.data[mask_3d==True]=0

    # get center
    if c_radec==True:
        center_ad=center
        tmp=wcs0.all_world2pix(center[0],center[1],0,0)
        center_pix=[float(tmp[0]),float(tmp[1])]
    else:
        center_pix=center
        tmp=wcs0.all_pix2world(center[0],center[1],0,0)
        center_ad=[float(tmp[0]),float(tmp[1])]

    # build location vectors of cube pixels
    xx,yy=np.meshgrid(np.arange(sz[2]),np.arange(sz[1]))
    xx=(xx-center_pix[0])*np.sqrt(hdu0.header['CD1_1']**2+hdu0.header['CD2_1']**2)*3600
    yy=(yy-center_pix[1])*np.sqrt(hdu0.header['CD1_2']**2+hdu0.header['CD2_2']**2)*3600
    if ellip!=1:
        rot=np.radians(pa)

        xx_circ=np.zeros(xx.shape)
        yy_circ=np.zeros(yy.shape)
        for i in range(xx.shape[1]):
            for j in range(yy.shape[0]):
                vec4=np.array([xx[j,i],yy[j,i]])
                ROT1=np.array([[np.cos(rot),np.sin(rot)],
                              [-np.sin(rot),np.cos(rot)]])
                if ellip>1:
                    FRAC=np.array([[ellip,0],
                                   [0,1]])
                else:
                    FRAC=np.array([[1/ellip,0],
                                   [0,1]])
                ROT2=np.array([[np.cos(rot),-np.sin(rot)],
                              [np.sin(rot),np.cos(rot)]])
                vec1=np.matmul(ROT2,np.matmul(FRAC,np.matmul(ROT1,vec4)))
                #vec1=np.matmul(FRAC,np.matmul(ROT1,vec4))
                xx_circ[j,i]=vec1[0]
                yy_circ[j,i]=vec1[1]

    else:
        xx_circ=xx
        yy_circ=yy

    rr_circ=np.sqrt(xx_circ**2+yy_circ**2)
    theta_circ=np.zeros(rr_circ.shape)
    theta_circ[(yy_circ==0) & (xx_circ>0)]=270
    theta_circ[(yy_circ==0) & (xx_circ<0)]=90
    theta_circ[yy_circ!=0]==np.degrees(np.arctan(-xx_circ[yy_circ!=0]/yy_circ[yy_circ!=0]))
    theta_circ[(xx_circ>0) & (yy_circ>0)]=theta_circ[(xx_circ>0) & (yy_circ>0)]+360
    theta_circ[(xx_circ<=0) & (yy_circ<0)]=theta_circ[(xx_circ<=0) & (yy_circ<0)]+180
    theta_circ[(xx_circ>0) & (yy_circ<0)]=theta_circ[(xx_circ>0) & (yy_circ<0)]+180
    #plt.pcolormesh(yy_circ)



    # setup target grid
    # x dimension
    if (npa==None) and (dpa==None):
        dx0=1.
        nx0=int(np.round((pa_range[1]-pa_range[0])/dx0))
        pa_range[1]=pa_range[0]+dx0*nx0
    else:
        if npa!=None:
            nx0=int(npa)
            dx0=(pa_range[1]-pa_range[0])/nx0
        if dpa!=None:
            nx0=int(np.round(pa_range[1]-pa_range[0]/dpa))
            dx0=dpa
            pa_range[1]=pa_range[0]+dx0*nx0

    nx=nx0
    dx=dx0

    # y dimension
    if r_range==None:
        p_corner_x=np.array([0,0,sz[2]-1,sz[2]-1])
        p_corner_y=np.array([0,sz[1]-1,sz[1]-1,0])
        p_corner_z=np.array([0,0,0,0])
        tmp=wcs2.wcs_pix2world(p_corner_x,p_corner_y,p_corner_z,0)
        a_corner=tmp[0]
        d_corner=tmp[1]
        r_corner=np.abs(d_corner+90)*3600.
        r_range=[0,np.max(r_corner)]
    if (nr==None) and (dr==None):
        dy=0.3
        ny=int(np.round((r_range[1]-r_range[0])/dy))
        r_range[1]=r_range[0]+dy*ny
    else:
        if nr!=None:
            ny=int(nr)
            dy=(r_range[1]-r_range[0])/ny
        if dr!=None:
            ny=int(np.round((r_range[1]-r_range[0])/dr))
            dy=dr
            r_range[1]=r_range[0]+dy*ny

    # project to 2D
    x_index=np.linspace(pa_range[0],pa_range[1],nx+1)
    y_index=np.linspace(r_range[0],r_range[1],ny+1)
    data5=np.zeros((sz[0],y_index.shape[0]-1,x_index.shape[0]-1))

    for i in range(x_index.shape[0]-1):
        for j in range(y_index.shape[0]-1):
            qq=(theta_circ >= x_index[i]) & (theta_circ < x_index[i+1]) & (rr_circ >= y_index[j]) & (rr_circ < y_index[j+1])
            if np.sum(qq)==0:
                continue

            element=hdu0_mask.data[:,qq]
            element[element==0]=np.nan

            data5[:,j,i]=np.nanmean(element,axis=1)


    # header
    hdr5=hdu0.header.copy()
    hdr5['NAXIS1']=nx0
    hdr5['CTYPE1']='PA'
    hdr5['CTYPE2']='Radius'
    hdr5['CNAME1']='PA'
    hdr5['CNAME2']='Radius'
    hdr5['CRVAL1']=pa_range[1]
    hdr5['CRVAL2']=r_range[0]
    hdr5['CRPIX2']=0.5
    hdr5['CUNIT2']='arcsec'
    hdr5['CD1_1']=-dx0
    hdr5['CD2_2']=dy
    hdr5['C2C_ORA']=(center_ad[0],'RA of origin')
    hdr5['C2C_ODEC']=(center_ad[1],'DEC of origin')
    hdr5['C2C_OX']=(center_pix[0]+1,'X of origin')
    hdr5['C2C_OY']=(center_pix[1]+1,'Y of origin')
    hdr5['C2C_E']=(ellip,'Axis raio')
    hdr5['C2C_EPA']=(pa,'PA of major-axis')
    hdr5['C2C_MTD']=('pix_mean','Algorithm')
    if ofn!='':
        hdr5['C2C_OFN']=(ofn,'Previous filename')
    # 2nd WCS
    hdr5['CTYPE1A']='PA'
    hdr5['CTYPE2A']='Radius'
    hdr5['CTYPE3A']=hdr5['CTYPE3']
    hdr5['CNAME1A']='PA'
    hdr5['CNAME2A']='Radius'
    hdr5['CNAME3A']=hdr5['CNAME3']
    hdr5['CRVAL1A']=pa_range[1]
    hdr5['CRVAL2A']=r_range[0]
    hdr5['CRVAL3A']=hdr5['CRVAL3']
    hdr5['CRPIX1A']=hdr5['CRPIX1']
    hdr5['CRPIX2A']=0.5
    hdr5['CRPIX3A']=hdr5['CRPIX3']
    hdr5['CUNIT1A']=hdr5['CUNIT1']
    hdr5['CUNIT2A']=hdr5['CUNIT2']
    hdr5['CUNIT3A']=hdr5['CUNIT3']
    hdr5['CD1_1A']=-dx0
    hdr5['CD2_2A']=dy
    hdr5['CD3_3A']=hdr5['CD3_3']


    if redshift>0:
        a_dis=(cos.arcsec_per_kpc_proper(redshift)).value

        hdr5['CRVAL2A']=r_range[0]/a_dis
        hdr5['CUNIT2A']='kpc'
        hdr5['CD2_2A']=dy/a_dis

        hdr5['CRVAL3A']=hdr5['CRVAL3']/(1+redshift)
        hdr5['CD3_3A']=hdr5['CD3_3']/(1+redshift)

    # Compress
    if compress==True:
        if nx0==1:
            tmp=hdr5.copy()
            hdr5['NAXIS']=2
            hdr5['NAXIS1']=tmp['NAXIS3']
            hdr5['NAXIS2']=tmp['NAXIS2']
            del hdr5['NAXIS3']
            hdr5['CTYPE1']=tmp['CTYPE3']
            hdr5['CTYPE2']=tmp['CTYPE2']
            hdr5['CTYPE1A']=tmp['CTYPE3A']
            hdr5['CTYPE2A']=tmp['CTYPE2A']
            del hdr5['CTYPE3']
            del hdr5['CTYPE3A']
            hdr5['CUNIT1']=tmp['CUNIT3']
            hdr5['CUNIT2']=tmp['CUNIT2']
            hdr5['CUNIT1A']=tmp['CUNIT3A']
            hdr5['CUNIT2A']=tmp['CUNIT2A']
            del hdr5['CUNIT3']
            del hdr5['CUNIT3A']
            hdr5['CNAME1']=tmp['CNAME3']
            hdr5['CNAME2']=tmp['CNAME2']
            hdr5['CNAME1A']=tmp['CNAME3A']
            hdr5['CNAME2A']=tmp['CNAME2A']
            del hdr5['CNAME3']
            del hdr5['CNAME3A']
            hdr5['CRVAL1']=tmp['CRVAL3']
            hdr5['CRVAL2']=tmp['CRVAL2']
            hdr5['CRVAL1A']=tmp['CRVAL3A']
            hdr5['CRVAL2A']=tmp['CRVAL2A']
            del hdr5['CRVAL3']
            del hdr5['CRVAL3A']
            hdr5['CRPIX1']=tmp['CRPIX3']
            hdr5['CRPIX2']=tmp['CRPIX2']
            hdr5['CRPIX1A']=tmp['CRPIX3A']
            hdr5['CRPIX2A']=tmp['CRPIX2A']
            del hdr5['CRPIX3']
            del hdr5['CRPIX3A']
            hdr5['CD1_1']=tmp['CD3_3']
            hdr5['CD2_2']=tmp['CD2_2']
            hdr5['CD1_1A']=tmp['CD3_3A']
            hdr5['CD2_2A']=tmp['CD2_2A']
            del hdr5['CD3_3']
            del hdr5['CD3_3A']
            data5=np.transpose(np.squeeze(data5,axis=2))

        elif ny==1:
            tmp=hdr5.copy()
            hdr5['NAXIS']=2
            hdr5['NAXIS1']=tmp['NAXIS3']
            hdr5['NAXIS2']=tmp['NAXIS1']
            del hdr5['NAXIS3']
            hdr5['CTYPE1']=tmp['CTYPE3']
            hdr5['CTYPE2']=tmp['CTYPE1']
            hdr5['CTYPE1A']=tmp['CTYPE3A']
            hdr5['CTYPE2A']=tmp['CTYPE1A']
            del hdr5['CTYPE3']
            del hdr5['CTYPE3A']
            hdr5['CUNIT1']=tmp['CUNIT3']
            hdr5['CUNIT2']=tmp['CUNIT1']
            hdr5['CUNIT1A']=tmp['CUNIT3A']
            hdr5['CUNIT2A']=tmp['CUNIT1A']
            del hdr5['CUNIT3']
            del hdr5['CUNIT3A']
            hdr5['CNAME1']=tmp['CNAME3']
            hdr5['CNAME2']=tmp['CNAME1']
            hdr5['CNAME1A']=tmp['CNAME3A']
            hdr5['CNAME2A']=tmp['CNAME1A']
            del hdr5['CNAME3']
            del hdr5['CNAME3A']
            hdr5['CRVAL1']=tmp['CRVAL3']
            hdr5['CRVAL2']=tmp['CRVAL1']
            hdr5['CRVAL1A']=tmp['CRVAL3A']
            hdr5['CRVAL2A']=tmp['CRVAL1A']
            del hdr5['CRVAL3']
            del hdr5['CRVAL3A']
            hdr5['CRPIX1']=tmp['CRPIX3']
            hdr5['CRPIX2']=tmp['CRPIX1']
            hdr5['CRPIX1A']=tmp['CRPIX3A']
            hdr5['CRPIX2A']=tmp['CRPIX1A']
            del hdr5['CRPIX3']
            del hdr5['CRPIX3A']
            hdr5['CD1_1']=tmp['CD3_3']
            hdr5['CD2_2']=tmp['CD1_1']
            hdr5['CD1_1A']=tmp['CD3_3A']
            hdr5['CD2_2A']=tmp['CD1_1A']
            del hdr5['CD3_3']
            del hdr5['CD3_3A']
            data5=np.transpose(np.squeeze(data5,axis=1))

    hdu5=fits.PrimaryHDU(data5,header=hdr5)

    if writefn!='':
        hdu5.writeto(writefn,overwrite=True)

    return hdu5
