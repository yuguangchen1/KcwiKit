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
from reproject import reproject_interp
from MontagePy.main import mProjectCube
from PyAstronomy import pyasl
from scipy import interpolate
from scipy import signal
from scipy import ndimage
#from image_registration import chi2_shift
#from image_registration.fft_tools import shift
import matplotlib
import matplotlib.pyplot as plt
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
		"drizzle":0.}

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


	return par


# Perform air-to-vac, and barycentric correction in wavelength
def kcwi_vachelio(hdu,hdr_ref='',mask=False):

	if hdr_ref=='':
		hdr_new=hdu.header.copy()
	else:
		hdr_new=hdu.header.copy()
		hdr_new['NAXIS3']=hdr_ref['NAXIS3']
		hdr_new['CTYPE3']=hdr_ref['CTYPE3']
		hdr_new['CUNIt3']=hdr_ref['CUNIT3']
		hdr_new['CNAME3']=hdr_ref['CNAME3']
		hdr_new['CRVAL3']=hdr_ref['CRVAL3']
		hdr_new['CRPIX3']=hdr_ref['CRPIX3']
		hdr_new['CD3_3']=hdr_ref['CD3_3']
		

	hdr_new['CTYPE3']='WAVE'

	hdr_old=hdu.header.copy()
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


	wave_vac=pyasl.airtovac2(wave_old)
	
	targ=coordinates.SkyCoord(hdr_old['TARGRA'],hdr_old['TARGDEC'],unit='deg',obstime=hdr_old['DATE-BEG'])
	keck=coordinates.EarthLocation.of_site('Keck Observatory')
	vcorr=targ.radial_velocity_correction(kind='barycentric',location=keck)
	
	wave_hel=wave_vac*(1-vcorr.value/2.99792458e8)

	
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
	
#	plt.clf()
#	plt.plot(wave_old,cube_old[:,45,15],drawstyle='steps-mid')
#	plt.plot(wave_new,cube_new.T[:,45,15],drawstyle='steps-mid')
#	plt.xlim(3300,3500)
#	plt.ylim(-0.01,0.01)
	return (hdu_new,vcorr)



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





def kcwi_stack(fnlist,shiftlist='',preshiftfn='',pixscale_x=0.,pixscale_y=0.,dimension=[0,0],orientation=-1000.,cubed=False,stepsig=0,drizzle=0,overwrite=False,keepmont=False):
#	fnlist="q0100-bx172.list"
#	shiftlist=""
#	preshiftfn=""
#	pixscale_x=0
#	pixscale_y=0
#	dimension=[0,0]
#	orientation=-1000.
#	cubed=False
#	stepsig=0
#	overwrite=False


	if cubed:
		suffix="cubed"
	else:
		suffix="cubes"

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




	# construct wcs
	hdulist=fits.open(fn[0])
	hdrtmp=hdulist[0].header.copy()
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
	old_cd11=hdr0['CD1_1']
	old_cd12=hdr0['CD1_2']
	old_cd21=hdr0['CD2_1']
	old_cd22=hdr0['CD2_2']
	hdr0['CD1_1']=-pixscale_x
	hdr0['CD2_2']=pixscale_y
	hdr0['CD1_2']=0
	hdr0['CD2_1']=0
	hdr0['CTYPE3']='WAVE'
	#hdr0['BUNIT']='10^(-16)erg/s/cm2/Angstrom'
	hdr0['BUNIT']='10^(-8)erg/s/cm3/arcsec2'
	if suffix!='cubes':
		#hdr0['BUNIT']='adu/s'
		hdr0['BUNIT']='adu/s/arcsec2'

	# orientation
	if orientation==-1000:
		orientation=par['stack_orientation']
		if orientation==-1000:
			orientation=np.rad2deg(np.arctan(old_cd21/(-old_cd11)))
	hdr0['CD1_1']=-pixscale_x*np.cos(np.deg2rad(orientation))
	hdr0['CD2_1']=pixscale_x*np.sin(np.deg2rad(orientation))
	hdr0['CD1_2']=pixscale_y*np.sin(np.deg2rad(orientation))
	hdr0['CD2_2']=pixscale_y*np.cos(np.deg2rad(orientation))
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
	data0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.float64).T
	vdata0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.float64).T
	mdata0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.int16).T+128
	edata0=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3'],len(fn)),dtype=np.float64).T
	#data0=[]
	#vdata0=[]
	#mdata0=[]
	#edata0=[]
	etime=np.zeros(len(fn))
	for i in range(len(fn)):
		print(os.path.basename(fn[i]))

		# check availability
		if (not os.path.isfile(trimfn[i])) or overwrite==True:
			# science cube
			hdulist=fits.open(fn[i])
			hdu_i,vcorr=kcwi_vachelio(hdulist[0],hdr_ref=hdr0)
			print('     Vcorr = '+str(vcorr.to('km/s')))
			hdulist.close()

			# variance cube -> sigma cube
			hdulist=fits.open(vfn[i])
			hdu_v,vcorr=kcwi_vachelio(hdulist[0],hdr_ref=hdr0)
			hdulist.close()

			# mask cube
			hdulist=fits.open(mfn[i])
			hdu_m,vcorr=kcwi_vachelio(hdulist[0],hdr_ref=hdr0,mask=True)
			hdulist.close()

			# Infinity check
			q=((hdu_i.data==0) | (~np.isfinite(hdu_i.data)) | (hdu_v.data==0) | (~np.isfinite(hdu_v.data)) )
			hdu_i.data[q]=np.nan
			hdu_v.data[q]=np.Inf
			hdu_m.data[q]=128
			

			# check EXPTIME
			hdu_i=kcwi_checkexptime(hdu_i)
			exptime=hdu_i.header['XPOSURE']
			print('     EXPTIME = '+str(exptime))
			etime[i]=exptime
			edata=hdu_i.data*0.+exptime
			q=(hdu_m.data != 0)
			edata[q]=0
			hdu_e=fits.PrimaryHDU(edata,header=hdu_i.header)
			hdu_e.header['BUNIT']='s'


			# have to use surface brightness for now, mProjectCube has bug with brightness units combined with drizzle scale
			dx=np.sqrt(hdu_i.header['CD1_1']**2+hdu_i.header['CD2_1']**2)*3600.
			dy=np.sqrt(hdu_i.header['CD1_2']**2+hdu_i.header['CD2_2']**2)*3600.
			area=dx*dy
			if suffix!='cubes':
				hdu_i.data=hdu_i.data/exptime/area
				hdu_v.data=hdu_v.data/exptime**2/area**2
				hdu_i.header['BUNIT']='adu/s/arcsec2'
				hdu_v.header['BUNIT']='adu2/s2/arcsec4'
			else:
				hdu_i.data=hdu_i.data/area
				hdu_v.data=hdu_v.data/area**2
				hdu_i.header['BUNIT']='10^(-8)erg/s/cm3/arcsec2'
				hdu_v.header['BUNIT']='10^(-16)erg2/s2/cm6/arcsec4'


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
			for kk in range(hdr0['NAXIS3']):
				img=hdu_i.data[kk,:,:]
				var=hdu_v.data[kk,:,:]
				mask=hdu_m.data[kk,:,:]
				expimg=hdu_e.data[kk,:,:]

				index_y,index_x=np.where(mask==0)
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

				mask[yrange[1]-trim[1,i]+1:,:]=128
				mask[:yrange[0]+trim[0,i],:]=128

				expimg[yrange[1]-trim[1,i]+1:,:]=0
				expimg[:yrange[0]+trim[0,i],:]=0

				hdu_i.data[kk,:,:]=img
				hdu_v.data[kk,:,:]=var
				hdu_m.data[kk,:,:]=mask
				hdu_e.data[kk,:,:]=expimg


			# write
			hdu_i.writeto(trimfn[i],overwrite=True)
			hdu_v.writeto(trimvfn[i],overwrite=True)
			hdu_m.writeto(trimmfn[i],overwrite=True)
			hdu_e.writeto(trimefn[i],overwrite=True)

		# Montage
		montfn=trimfn[i].replace('.trim.fits','.mont.fits')
		montvfn=trimvfn[i].replace('.trim.fits','.mont.fits')
		montmfn=trimmfn[i].replace('.trim.fits','.mont.fits')
		montefn=trimefn[i].replace('.trim.fits','.mont.fits')
		if (not os.path.isfile(montfn)) or overwrite==True:
			# using shell version for now, see if energy leakage is still a problem
			exe="mProjectCube -z "+str(drizzle)+" -f "+trimfn[i]+" "+montfn+" "+fnhdr
			void=os.system(exe)
			exev="mProjectCube -z "+str(drizzle)+" -f  "+trimvfn[i]+" "+montvfn+" "+fnhdr
			voidv=os.system(exev)
			exem="mProjectCube -z "+str(drizzle)+" -f  "+trimmfn[i]+" "+montmfn+" "+fnhdr
			voidm=os.system(exem)
			exee="mProjectCube -z "+str(drizzle)+" -f  "+trimefn[i]+" "+montefn+" "+fnhdr
			voide=os.system(exee)

			#void=mProjectCube(trimfn[i],montfn,fnhdr,drizzle=drizzle,energyMode=True,fullRegion=True)
			#voidv=mProjectCube(trimvfn[i],montvfn,fnhdr,drizzle=drizzle,energyMode=True,fullRegion=True)
			#voidm=mProjectCube(trimmfn[i],montmfn,fnhdr,drizzle=drizzle,energyMode=False,fullRegion=True)
			#voide=mProjectCube(trimefn[i],montefn,fnhdr,drizzle=drizzle,energyMode=False,fullRegion=True)

		# cach cubes
		newcube=fits.open(montfn)[0].data
		newcube[~np.isfinite(newcube)]=0.
		newcubev=fits.open(montvfn)[0].data
		newcubev[~np.isfinite(newcubev)]=0.
		newcubem=np.ceil(fits.open(montmfn)[0].data)
		newcubem[~np.isfinite(newcubem)]=128
		newcubee=fits.open(montefn)[0].data
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

		if keepmont==False:
			os.remove(montfn)
			os.remove(montvfn)
			os.remove(montmfn)
			os.remove(montefn)
			os.remove(montfn.replace('mont','mont_area'))
			os.remove(montvfn.replace('mont','mont_area'))
			os.remove(montmfn.replace('mont','mont_area'))
			os.remove(montefn.replace('mont','mont_area'))


	# Stacking!!!
	print('Stacking...')
	data_3d=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3']),dtype=np.float64)
	vdata_3d=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3']),dtype=np.float64)
	edata_3d=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3']),dtype=np.float64)
	mdata_3d=np.zeros((dimension[0],dimension[1],hdr0['NAXIS3']),dtype=np.int16)+1
	for ii in tqdm(range(dimension[0])):
		#print(ii)
		for jj in range(dimension[1]):
			for kk in range(hdr0['NAXIS3']):
				#img=np.array([cube[kk,jj,ii] for cube in data0])
				#var=np.array([cube[kk,jj,ii] for cube in vdata0])
				#mask=np.array([cube[kk,jj,ii] for cube in mdata0])
				#exp=np.array([cube[kk,jj,ii] for cube in edata0])
				img=data0[:,kk,jj,ii]
				var=vdata0[:,kk,jj,ii]
				mask=mdata0[:,kk,jj,ii]
				exp=edata0[:,kk,jj,ii]

				mask[~np.isfinite(img)]=1
				mask[~np.isfinite(var)]=1
				mask[var==0]=1

				q=(mask==0)
				if np.sum(q)==0:
					continue
				
				weight=np.zeros(var.shape[0])
				weight[var!=0]=1/np.abs(var[var!=0])
				#weight=exp
				weight[~q]=0

				weight[~np.isfinite(weight)]=0
				
				#q2=stats.sigma_clip(img[q],sigma=5,masked=True)
				#weight[q][q2.mask]=0

				data_3d[ii,jj,kk]=np.sum(img*weight)/np.sum(weight)
				vdata_3d[ii,jj,kk]=np.sum(weight**2*var)/np.sum(weight)**2
				edata_3d[ii,jj,kk]=np.sum(exp[weight!=0])
				mdata_3d[ii,jj,kk]=0


	# write
	vhdr0=hdr0.copy()
	if suffix=='cubes':
		#vhdr0['BUNIT']='10^(-32)erg2/s2/cm4/Angstrom2'
		vhdr0['BUNIT']='10^(-8)erg/s/cm3/arcsec2'
	else:
		#vhdr0['BUNIT']='adu2/s2'
		vhdr0['BUNIT']='adu2/s2/arcsec4'

	mhdr0=hdr0.copy()
	del mhdr0['BUNIT']
	mhdr0['BITPIX']=16

	ehdr0=hdr0.copy()
	ehdr0['BUNIT']='s'


	hdu_i=fits.PrimaryHDU(data_3d.T,header=hdr0)
	hdu_i.writeto(fnlist.replace('.list','_i'+suffix+'.fits'),overwrite=True)
	hdu_v=fits.PrimaryHDU(vdata_3d.T,header=vhdr0)
	hdu_v.writeto(fnlist.replace('.list','_v'+suffix+'.fits'),overwrite=True)
	hdu_m=fits.PrimaryHDU(mdata_3d.T,header=mhdr0)
	hdu_m.writeto(fnlist.replace('.list','_m'+suffix+'.fits'),overwrite=True)
	hdu_e=fits.PrimaryHDU(edata_3d.T,header=ehdr0)
	hdu_e.writeto(fnlist.replace('.list','_e'+suffix+'.fits'),overwrite=True)

	end=ostime.time()
	#print(end-start)
	
	return 
	



def kcwi_align(fnlist,wavebin=[-1.,-1.],box=[-1,-1,-1,-1],pixscale_x=-1.,pixscale_y=-1.,orientation=-1000.,dimension=[-1.,-1.],preshiftfn='',trim=[-1,-1],cubed=False,noalign=False,display=True,search_size=-1000,conv_filter=-1000,upfactor=-1000.):

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
			orientation=np.ra2deg(np.arctan(old_cd21/(-old_cd11)))
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
		
		hdu=fits.open(fn[i])[0]
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

		# preshift
		if preshiftfn!='':
			index=np.where(np.array(prefn)==os.path.basename(fn[i]))
			index=index[0]
			if len(index)>0:
				index=index[0]
				hdr['CRVAL1']=hdr['CRVAL1']+prera[index]/3600.
				hdr['CRVAL2']=hdr['CRVAL2']+predec[index]/3600.
		
		# initial projection
		newthum,coverage=reproject_interp((thum.T,hdr),hdr0,order='bilinear')
		newthum=newthum.T
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
						mult=cut0*cut
						if np.sum(mult!=0)>0:
							crls[ii,jj]=np.sum(mult)/np.sum(mult!=0)
				
				fig=plt.figure(1)
				plt.clf()
				xplot=np.append(xx,xx[1]-xx[0]+xx[-1])-0.5
				yplot=np.append(yy,yy[1]-yy[0]+yy[-1])-0.5
				plt.pcolormesh(xplot,yplot,crls.T)

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
				hdr0_up['CD1_1']=hdr0_up['CD1_1']/upfactor
				hdr0_up['CD2_1']=hdr0_up['CD2_1']/upfactor
				hdr0_up['CD1_2']=hdr0_up['CD1_2']/upfactor
				hdr0_up['CD2_2']=hdr0_up['CD2_2']/upfactor
				newthum1,coverage=reproject_interp((thum_1.T,hdr_1),hdr0_up,order='bilinear')
				newthum1=newthum1.T

				# do the shift from last iteration
				wcs_hdr0=wcs.WCS(hdr0)
				tmp=wcs_hdr0.all_pix2world(hdr0['CRPIX1']+xshift[i],hdr0['CRPIX2']+yshift[i],1)
				hdr['CRVAL1']=hdr['CRVAL1']+(float(tmp[0])-hdr0['CRVAL1'])
				hdr['CRVAL2']=hdr['CRVAL2']+(float(tmp[1])-hdr0['CRVAL2'])
				newthum2,coverage=reproject_interp((thum.T,hdr),hdr0_up,order='bilinear')
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
						mult=cut*cut0
						crls[ii,jj]=np.sum(mult)/np.sum(mult!=0)
				
				#plt.figure(1)
				#plt.clf()
				xplot=(np.append(xx,xx[1]-xx[0]+xx[-1])-0.5)/upfactor+xshift[i]
				yplot=(np.append(yy,yy[1]-yy[0]+yy[-1])-0.5)/upfactor+yshift[i]
				plt.pcolormesh(xplot,yplot,crls.T,cmap='plasma')

				tmp=np.unravel_index(crls.argmax(),crls.shape)
				xshift[i]+=xx[tmp[0]]/upfactor
				yshift[i]+=yy[tmp[1]]/upfactor
				plt.plot(xshift[i],yshift[i],'+',color='r')
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
				data_thum[:,:,i]=thum_shift
				
				#hdr_shift=hdr_preshift.copy()
				#hdr_shift['CRVAL1']=hdr_shift['CRVAL1']+ashift
				#hdr_shift['CRVAL2']=hdr_shift['CRVAL2']+dshift
				#thumshift,coverage=reproject_interp((thum.T,hdr),hdr0,order='bilinear')
				#thumshift=thumshift.T
				#data_thum[:,:,i]=thum_shift

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


	hdu=fits.PrimaryHDU(data0_thum.T)
	hdu.writeto('kcwi_align/'+fnlist.replace('.list','.thum0.fits'),overwrite=True)

	writefn=[i.replace('.fits','') for i in fn]
	xytable=table.Table([np.array(writefn),xshift_xy,yshift_xy])
	ascii.write(xytable,fnlist.replace('.list','.shift.list'),overwrite=True,format='no_header')

	if display==False:
		matplotlib.use(oldbackend)
	return



def kcwi_astrometry(fnlist,imgfn='',wavebin=[-1.,-1.],display=True,search_size=-1000,conv_filter=-1000,upfactor=-1000,box=[-1.,-1.,-1.,-1.],nocrl=0):
	
	print(os.path.basename(fnlist))

	suffix='cubes'
	cubefn=fnlist.replace('.list','_i'+suffix+'.fits')
	if path.isfile(cubefn)==False:
		suffix='cubed'
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
	hdu_cube.header['CRPIX1']=par['ref_xy'][0]
	hdu_cube.header['CRPIX2']=par['ref_xy'][1]
	hdu_cube.header['CRVAL1']=par['ref_ad'][0]
	hdu_cube.header['CRVAL2']=par['ref_ad'][1]


	# collapsing
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
	

	if display==False:
		matplotlib.use(oldbackend)
	return





