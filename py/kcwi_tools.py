import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
import os
import warnings
import pdb


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

	newwave=wave0[qwave]
	newhdu.data=newhdu.data[qwave,:,:]
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



def cart2cyli(hdu,center,writefn='',ellip=1.,pa=0.,nr=None,npa=None,r_range=None,pa_range=[0,360],dr=None,dpa=None,drizzle=1.0,c_radec=False,clean=True,compress=False):
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
	
	# 0 - the original HDU
	hdu0=hdu
	wcs0=wcs.WCS(hdu0.header)
	sz=hdu0.data.shape

		
	# Skewing
	hdu1=hdu0.copy()
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

	cube1fn='kcwi_tools/cart2cyli_cube1.fits'
	hdu1.writeto(cube1fn,overwrite=True)

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
	hdr3_1['NAXIS1']=nx3[0]
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
	hdr3_2['NAXIS1']=nx3[1]
	hdr3_2['CRVAL1']=xr3[0,1]
	hdr3_2fn='kcwi_tools/cart2cyli_hdr3_2.hdr'
	hdr3_2.totextfile(hdr3_2fn,overwrite=True)

	hdr3_3=hdr3_1.copy()
	hdr3_3['NAXIS1']=nx3[2]
	hdr3_3['CRVAL1']=xr3[0,2]
	hdr3_3fn='kcwi_tools/cart2cyli_hdr3_3.hdr'
	hdr3_3.totextfile(hdr3_3fn,overwrite=True)


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
	hdutmp=fits.PrimaryHDU(area4)
	hdutmp.writeto('kcwi_tools/tmp.fits',overwrite=True)

	
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
	hdr5['CD1_1']=-dx0
	hdr5['CD2_2']=dy
	hdr5['C2C_ORA']=(center_ad[0],'RA of origin')
	hdr5['C2C_ODEC']=(center_ad[1],'DEC of origin')
	hdr5['C2C_OX']=(center_pix[0]+1,'X of origin')
	hdr5['C2C_OY']=(center_pix[1]+1,'Y of origin')
	hdr5['C2C_E']=(ellip,'Axis raio')
	hdr5['C2C_EPA']=(pa,'PA of major-axis')
	hdr5['C2C_DRIZ']=(drizzle,'Drizzle factor')
	if ofn!='':
		hdr5['C2C_OFN']=(ofn,'Previous filename')

	ahdr5=hdr5.copy()
	ahdr5['NAXIS']=2
	del ahdr5['NAXIS3']
	del ahdr5['CTYPE3']
	del ahdr5['CUNIT3']
	del ahdr5['CNAME3']
	del ahdr5['CRVAL3']
	del ahdr5['CRPIX3']
	del ahdr5['CD3_3']

	data5=np.zeros((data4.shape[0],data4.shape[1],nx0))
	area5=np.zeros((data4.shape[1],nx0))
	ratio=int(dx0/dx)
	if ratio!=1:
		#warnings.filterwarnings('ignore')
		for i in range(nx0):
			tmp=data4[:,:,i*ratio:(i+1)*ratio]
			data5[:,:,i]=np.nanmean(tmp,axis=2)
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
			del hdr5['CTYPE3']
			hdr5['CUNIT1']=tmp['CUNIT3']
			hdr5['CUNIT2']=tmp['CUNIT2']
			del hdr5['CUNIT3']
			hdr5['CNAME1']=tmp['CNAME3']
			hdr5['CNAME2']=tmp['CNAME2']
			del hdr5['CNAME3']
			hdr5['CRVAL1']=tmp['CRVAL3']
			hdr5['CRVAL2']=tmp['CRVAL2']
			del hdr5['CRVAL3']
			hdr5['CRPIX1']=tmp['CRPIX3']
			hdr5['CRPIX2']=tmp['CRPIX2']
			del hdr5['CRPIX3']
			hdr5['CD1_1']=tmp['CD3_3']
			hdr5['CD2_2']=tmp['CD2_2']
			del hdr5['CD3_3']
			data5=np.transpose(np.squeeze(data5,axis=2))
			
			tmp=ahdr5.copy()
			ahdr5['NAXIS']=1
			ahdr5['NAXIS1']=tmp['NAXIS2']
			del ahdr5['NAXIS2']
			ahdr5['CTYPE1']=tmp['CTYPE2']
			del ahdr5['CTYPE2']
			ahdr5['CUNIT1']=tmp['CUNIT2']
			del ahdr5['CUNIT2']
			ahdr5['CNAME1']=tmp['CNAME2']
			del ahdr5['CNAME2']
			ahdr5['CRVAL1']=tmp['CRVAL2']
			del ahdr5['CRVAL2']
			ahdr5['CRPIX1']=tmp['CRPIX2']
			del ahdr5['CRPIX2']
			ahdr5['CD1_1']=tmp['CD2_2']
			del ahdr5['CD2_2']
			area5=np.transpose(np.squeeze(area5))

		elif ny==1:
			tmp=hdr5.copy()
			hdr5['NAXIS']=2
			hdr5['NAXIS1']=tmp['NAXIS3']
			hdr5['NAXIS2']=tmp['NAXIS1']
			del hdr5['NAXIS3']
			hdr5['CTYPE1']=tmp['CTYPE3']
			hdr5['CTYPE2']=tmp['CTYPE1']
			del hdr5['CTYPE3']
			hdr5['CUNIT1']=tmp['CUNIT3']
			hdr5['CUNIT2']=tmp['CUNIT1']
			del hdr5['CUNIT3']
			hdr5['CNAME1']=tmp['CNAME3']
			hdr5['CNAME2']=tmp['CNAME1']
			del hdr5['CNAME3']
			hdr5['CRVAL1']=tmp['CRVAL3']
			hdr5['CRVAL2']=tmp['CRVAL1']
			del hdr5['CRVAL3']
			hdr5['CRPIX1']=tmp['CRPIX3']
			hdr5['CRPIX2']=tmp['CRPIX1']
			del hdr5['CRPIX3']
			hdr5['CD1_1']=tmp['CD3_3']
			hdr5['CD2_2']=tmp['CD1_1']
			del hdr5['CD3_3']
			data5=np.transpose(np.squeeze(data5,axis=1))

			
	hdu5=fits.PrimaryHDU(data5,header=hdr5)
	ahdu5=fits.PrimaryHDU(area5,header=ahdr5)

			

	if writefn!='':
		hdu5.writeto(writefn,overwrite=True)
		ahdu5.writeto(writefn.replace('.fits','_area.fits'),overwrite=True)

	if clean==True:
		os.system('rm -f kcwi_tools/cart2cyli*')
		
	return (hdu5,ahdu5)
	

