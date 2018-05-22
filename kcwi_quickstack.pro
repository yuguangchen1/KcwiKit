
; quick stack in 3D space
pro kcwi_quickstack,fnlist,shiftlist,preshiftfn=preshiftfn,$
	pixscale=pixscale,dimension=dimension,orientaion=orientaion,$
	cubed=cubed

if ~keyword_set(cubed) then begin
	suffix='cubes'
endif else begin
	suffix='cubed'
endelse

parfn=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.par'
par=kcwi_stack_readpar(parfn)

if ~keyword_set(shiftlist) then begin
	shiftlist=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.shift.list'
endif

if ~keyword_set(pixscale) then begin
	pixscale=par.pixscale
	if pixscale eq -1 then begin
		pixscale=0.3
	endif 
endif
pixscale=pixscale/3600.

if ~keyword_set(dimension) then begin
	dimension=par.dimension
	if dimension[0] eq -1 then begin
		dimension=[100,100]
	endif
endif

readcol,fnlist,fn,trim1,trim2,format='a,f,f',comment='#'
trim=fltarr(2,n_elements(fn))
trim[0,*]=trim1
trim[1,*]=trim2

readcol,shiftlist,fn,xshift,yshift,format='a,f,f',comment='#'
if file_test(fn[0]+'_i'+suffix+'.fits') eq 0 then begin
	suffix='cubed'
endif
vfn=fn+'_v'+suffix+'.fits'
mfn=fn+'_m'+suffix+'.fits'
fn=fn+'_i'+suffix+'.fits'

if ~keyword_set(preshiftfn) then begin
	preshiftfn=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+$
		'.preshift.list'
	if file_test(preshiftfn) eq 0 then begin
		preshiftfn=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+$
			'.pre.list'
		if file_test(preshiftfn) eq 0 then begin
			preshiftfn=''
		endif
	endif
endif

if preshiftfn ne '' then begin
	readcol,preshiftfn,prefn,prera,predec,format='a,f,f',comment='#'
	prefn=prefn+'_i'+suffix+'.fits'
endif

; construct WCS
hdr0=headfits(fn[0])
sxaddpar,hdr0,'naxis1',dimension[0]
sxaddpar,hdr0,'naxis2',dimension[1]
sxaddpar,hdr0,'crpix1',dimension[0]/2.
sxaddpar,hdr0,'crpix2',dimension[1]/2.
old_cd11=sxpar(hdr0,'cd1_1')
old_cd12=sxpar(hdr0,'cd1_2')
old_cd21=sxpar(hdr0,'cd2_1')
old_cd22=sxpar(hdr0,'cd2_2')
sxaddpar,hdr0,'cd1_1',-pixscale
sxaddpar,hdr0,'cd2_2',pixscale
sxaddpar,hdr0,'cd1_2',0
sxaddpar,hdr0,'cd2_1',0
sxaddpar,hdr0,'naxis',2
old_cd33=sxpar(hdr0,'cd3_3')
old_crval3=sxpar(hdr0,'crval3')
old_crpix3=sxpar(hdr0,'crpix3')
old_naxis3=sxpar(hdr0,'naxis3')
sxdelpar,hdr0,'cd3_3'
sxdelpar,hdr0,'crval3'
sxdelpar,hdr0,'crpix3'
sxdelpar,hdr0,'naxis3'


; orientation
if ~isa(orientation) then begin
	orientation=par.orientation
	if orientation eq -1000 then begin
		orientation=atan(old_cd21/(-old_cd11))*!radeg
		;orientation=atan(old_cd12/old_cd22)*!radeg
	endif
endif
sxaddpar,hdr0,'cd1_1',-pixscale*cos(orientation/!radeg)
sxaddpar,hdr0,'cd2_1',pixscale*sin(orientation/!radeg)
sxaddpar,hdr0,'cd1_2',pixscale*sin(orientation/!radeg)
sxaddpar,hdr0,'cd2_2',pixscale*cos(orientation/!radeg)


; distribute flux
old_pixscale=fltarr(2)
old_pixscale[0]=sqrt(old_cd11^2+old_cd21^2)
old_pixscale[1]=sqrt(old_cd12^2+old_cd22^2)
fluxscale=pixscale^2/(old_pixscale[0]*old_pixscale[1])


print,'Projecting...'
data0=dblarr(dimension[0],dimension[1],old_naxis3,n_elements(fn))
vdata0=data0
mdata0=data0+1
edata0=data0
etime=fltarr(n_elements(fn))
for i=0,n_elements(fn)-1 do begin
	print,file_basename(fn[i])

	data=mrdfits(fn[i],0,hdr,/silent)
	data=double(data)
	sz=size(data,/dim)
	sxdelpar,hdr,'cd3_3'
	sxdelpar,hdr,'crval3'
	sxdelpar,hdr,'crpix3'
	sxdelpar,hdr,'naxis3'
	sxaddpar,hdr,'naxis',2
	
	vdata=mrdfits(vfn[i],0,vhdr,/silent)
	vdata=double(vdata)
	sxdelpar,vhdr,'cd3_3'
	sxdelpar,vhdr,'crval3'
	sxdelpar,vhdr,'crpix3'
	sxdelpar,vhdr,'naxis3'
	sxaddpar,vhdr,'naxis',2

	mdata=mrdfits(mfn[i],0,mhdr,/silent)
	sxdelpar,mhdr,'cd3_3'
	sxdelpar,mhdr,'crval3'
	sxdelpar,mhdr,'crpix3'
	sxdelpar,mhdr,'naxis3'
	sxaddpar,mhdr,'naxis',2

	; infinity check
	q=where(finite(data) eq 0 or finite(vdata) eq 0)
	if q[0] ne -1 then begin
		data[q]=0.
		vdata[q]=0.
		mdata[q]=1
	endif


	; exptime
	exptime=sxpar(hdr,'xposure')
	expend=sxpar(hdr,'date-end')
	expend=date_conv(expend,'MODIFIED')
	rdend=sxpar(hdr,'daterend')
	rdend=date_conv(rdend,'MODIFIED')
	if rdend le expend then begin
		print,'Warning XPOSURE incorrect...'
		expbeg=date_conv(sxpar(hdr,'date-beg'),'MODIFIED')
		expend=rdend-53.64/3600./24.
		exptime=(expend-expbeg)*3600.*24.
		
		;expend=date_conv(expend,'FITS')
		;sxaddpar,hdr,'date-end',expend
		sxaddpar,hdr,'xposure',exptime
		print,'   Setting to'+string(exptime)+'.' 
	endif
	print,exptime
	etime[i]=exptime
	edata=data*0.+exptime
	q=where(data eq 0 or mdata ne 0)
	if q[0] ne -1 then edata[q]=0.

	; count/s
	if suffix ne 'cubes' then begin
		data=data/exptime
		vdata=vdata/exptime^2
	endif
	
	; flux
	data=data*fluxscale
	vdata=vdata/fluxscale^2
	

	; pre-shift
	if keyword_set(preshiftfn) then begin
		q=where(prefn eq fn[i])
		if q[0] ne -1 then begin
			sxaddpar,hdr,'crval1',sxpar(hdr,'crval1')+prera[q[0]]/3600.
			sxaddpar,hdr,'crval2',sxpar(hdr,'crval2')+predec[q[0]]/3600.

			sxaddpar,vhdr,'crval1',sxpar(vhdr,'crval1')+prera[q[0]]/3600.
			sxaddpar,vhdr,'crval2',sxpar(vhdr,'crval2')+predec[q[0]]/3600.

			sxaddpar,mhdr,'crval1',sxpar(mhdr,'crval1')+prera[q[0]]/3600.
			sxaddpar,mhdr,'crval2',sxpar(mhdr,'crval2')+predec[q[0]]/3600.
		endif
	endif


	for kk=0,sz[2]-1 do begin
		hdrshift=hdr
		vhdrshift=vhdr
		mhdrshift=mhdr

		sxaddpar,hdrshift,'crpix1',sxpar(hdr,'crpix1')+xshift[i]
		sxaddpar,hdrshift,'crpix2',sxpar(hdr,'crpix2')+yshift[i]
		sxaddpar,vhdrshift,'crpix1',sxpar(vhdr,'crpix1')+xshift[i]
		sxaddpar,vhdrshift,'crpix2',sxpar(vhdr,'crpix2')+yshift[i]
		sxaddpar,mhdrshift,'crpix1',sxpar(mhdr,'crpix1')+xshift[i]
		sxaddpar,mhdrshift,'crpix2',sxpar(mhdr,'crpix2')+yshift[i]

		img=data[*,*,kk]
		var=vdata[*,*,kk]
		mask=mdata[*,*,kk]
		expimg=edata[*,*,kk]

		;weight=sig*0.
		;q=where(sig ne 0)
		;weight[q]=1./sig[q]^2

		;q=where(mask ne 0)
		;if q[0] ne -1 then begin
		;	sig[q]=0.
		;	expimg[q]=0.
		;endif

		q=where(img eq 0)
		if q[0] ne -1 then begin
			var[q]=0.
			expimg[q]=0.
			mask[q]=1
		endif

		ehdrshift=mhdrshift

		; trim
		q=where(img ne 0)
		if q[0] eq -1 then continue
		qx=q mod sz[0]
		qy=q/sz[0]
		xrange=[min(qx),max(qx)]
		yrange=[min(qy),max(qy)]
		if yrange[0]+trim[0,i] ge yrange[1]-trim[1,i] then continue
		hextract,img,hdrshift,xrange[0],xrange[1],$
			yrange[0]+trim[0,i],yrange[1]-trim[1,i],/silent
		hextract,var,vhdrshift,xrange[0],xrange[1],$
			yrange[0]+trim[0,i],yrange[1]-trim[1,i],/silent
		hextract,mask,mhdrshift,xrange[0],xrange[1],$
			yrange[0]+trim[0,i],yrange[1]-trim[1,i],/silent
		hextract,expimg,ehdrshift,xrange[0],xrange[1],$
			yrange[0]+trim[0,i],yrange[1]-trim[1,i],/silent
		

		hastrom,img,hdrshift,newimg,newhdr,hdr0,missing=0,interp=2,cubic=-0.5
		if (size(img))[0] ne 2 then stop
		data0[*,*,kk,i]=newimg
		
		hastrom,var,vhdrshift,newvar,newvhdr,hdr0,missing=0,interp=2,cubic=-0.5
		vdata0[*,*,kk,i]=newvar

		hastrom,mask,mhdrshift,newmask,newmhdr,hdr0,missing=1,interp=1,cubic=-0.5
		mdata0[*,*,kk,i]=newmask

		hastrom,expimg,ehdrshift,newexpimg,newehdr,hdr0,missing=0,interp=1,cubic=-0.5
		edata0[*,*,kk,i]=newexpimg
		
	endfor

endfor

print,'Stacking...'
data_3d=dblarr(dimension[0],dimension[1],old_naxis3)
vdata_3d=data_3d
edata_3d=data_3d
mdata_3d=long(data_3d)+1
for ii=0,dimension[0]-1 do begin
	for jj=0,dimension[1]-1 do begin
		for kk=0,old_naxis3-1 do begin
			q=where(vdata0[ii,jj,kk,*] ne 0 and mdata0[ii,jj,kk,*] eq 0)
			if q[0] ne -1 then begin
				weight=1./abs(vdata0[ii,jj,kk,*])
				index=where(vdata0[ii,jj,kk,*] eq 0 or mdata0[ii,jj,kk,*] ne 0)
				if index[0] ne -1 then weight[index]=0.
				;weight=edata0[ii,jj,kk,*]

				data_3d[ii,jj,kk]=total(data0[ii,jj,kk,*]*weight)/$
					total(weight)
				vdata_3d[ii,jj,kk]=total((weight*$
					vdata0[ii,jj,kk,*])^2)/total(weight)^2
				edata_3d[ii,jj,kk]=total(edata0[ii,jj,kk,q])
				mdata_3d[ii,jj,kk]=0

				;if finite(data_3d[ii,jj,kk]) eq 0 then stop
			endif
		endfor
	endfor
endfor

sxaddpar,hdr0,'cd3_3',old_cd33
sxaddpar,hdr0,'crval3',old_crval3
sxaddpar,hdr0,'crpix3',old_crpix3
sxaddpar,hdr0,'naxis3',old_naxis3
sxaddpar,hdr0,'bitpix',-64
sxaddpar,hdr0,'naxis',3

writefits,file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'_i'+suffix+'.fits',data_3d,hdr0
writefits,file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'_v'+suffix+'.fits',vdata_3d,hdr0
writefits,file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'_e'+suffix+'.fits',edata_3d,hdr0

sxaddpar,hdr0,'bitpix',32
writefits,file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'_m'+suffix+'.fits',mdata_3d,hdr0

end
