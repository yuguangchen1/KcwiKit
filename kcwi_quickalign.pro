
; quick alignment
pro kcwi_quickalign,fnlist,wavebin=wavebin,box=box,prefix=prefix,num=num,$
	pixscale_x=pixscale_x,pixscale_y=pixscale_y,orientation=orientation,$
	dimension=dimension,preshiftfn=preshiftfn,$
	trim=trim,cubed=cubed,noalign=noalign
; Users can use prefix+num method instead of providing fnlist.
; Unit in preshiftfn is RA, DEC in arcsec.

if ~keyword_set(cubed) then begin
	suffix='cubes'
endif else begin
	suffix='cubed'
endelse

parfn=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.par'
par=kcwi_stack_readpar(parfn)

if ~keyword_set(wavebin) then begin
	wavebin=par.wavebin
	if wavebin[0] eq -1 then begin
		wavebin=[4000.,5000.]
	endif
endif


if ~keyword_set(box) then begin
	box=par.align_box
	if box[0] eq -1 then begin
		box=[25,50,25,40]
	endif
endif

if ~keyword_set(pixscale_x) then begin
	pixscale_x=par.align_xpix
	if pixscale_x eq -1 then begin
		pixscale_x=0.3
	endif
endif
pixscale_x=pixscale_x/3600.

if ~keyword_set(pixscale_y) then begin
	pixscale_y=par.align_ypix
	if pixscale_y eq -1 then begin
		pixscale_y=0.3
	endif
endif
pixscale_y=pixscale_y/3600.

if ~keyword_set(dimension) then begin
	dimension=par.align_dimension
	if dimension[0] eq -1 then begin
		dimension=[100,100]
	endif
endif

if ~keyword_set(fnlist) then begin
;dir='/scr/yuguangchen/obs/kcwi/kcwi_sep17/2017sep22/redux_don_yc/redux/'
	if n_elements(prefix) ne n_elements(num) or n_elements(prefix) ne 1 then begin
		print,'% QUICKASTROM Error: Incorrect number of dir.'
	      	return
	endif
	fn0=prefix+num
	fn=fn0+'_i'+suffix+'.fits'
	;fn=dir+'kb170922_'+num+'_icuber_3d.fits'
endif else begin
	;readcol,'bx418.list',num,format='a'
	readcol,fnlist,fn0,trim1,trim2,format='a',comment='#'
	fn=fn0+'_i'+suffix+'.fits'
	trim=fltarr(2,n_elements(fn))
	trim[0,*]=trim1
	trim[1,*]=trim2
endelse

if ~keyword_set(preshiftfn) then begin
	tmp=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.preshift.list'
	if file_test(tmp) then begin
		preshiftfn=tmp
	endif else begin
		tmp=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.pre.list'
		if file_test(tmp) then begin
			preshiftfn=tmp
		endif else preshiftfn=''
	endelse
endif
if preshiftfn ne '' then begin
	readcol,preshiftfn,prefn,prera,predec,format='a,f,f',comment='#'
	prefn=prefn+'_i'+suffix+'.fits'
endif

starteps,file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.align.ps',xs=6,ys=6
loadct,13

; construct WCS
hdrtmp=headfits(fn[0])
xyad,hdrtmp,float(sxpar(hdrtmp,'naxis1'))/2.-1,float(sxpar(hdrtmp,'naxis2'))/2.-1,ratmp,dectmp
center=[ratmp,dectmp]

hdr0=hdrtmp
sxaddpar,hdr0,'naxis1',dimension[0]
sxaddpar,hdr0,'naxis2',dimension[1]
sxaddpar,hdr0,'crpix1',float(dimension[0])/2.
sxaddpar,hdr0,'crpix2',float(dimension[1])/2.
sxaddpar,hdr0,'crval1',center[0]
sxaddpar,hdr0,'crval2',center[1]
old_cd11=sxpar(hdr0,'cd1_1')
old_cd12=sxpar(hdr0,'cd1_2')
old_cd21=sxpar(hdr0,'cd2_1')
old_cd22=sxpar(hdr0,'cd2_2')
sxaddpar,hdr0,'cd1_1',-pixscale_x
sxaddpar,hdr0,'cd2_2',pixscale_y
sxaddpar,hdr0,'cd1_2',0
sxaddpar,hdr0,'cd2_1',0
sxaddpar,hdr0,'naxis',2
sxdelpar,hdr0,'cd3_3'
sxdelpar,hdr0,'crval3'
sxdelpar,hdr0,'crpix3'
sxdelpar,hdr0,'naxis3'


; orientation
if ~isa(orientation) then begin
	orientation=par.align_orientation
	if orientation eq -1000 then begin
		orientation=atan(old_cd21/(-old_cd11))*!radeg
	endif
endif
sxaddpar,hdr0,'cd1_1',-pixscale_x*cos(orientation/!radeg)
sxaddpar,hdr0,'cd2_1',pixscale_x*sin(orientation/!radeg)
sxaddpar,hdr0,'cd1_2',pixscale_y*sin(orientation/!radeg)
sxaddpar,hdr0,'cd2_2',pixscale_y*cos(orientation/!radeg)


data_thum=fltarr(dimension[0],dimension[1],n_elements(fn))
data0_thum=data_thum
xshift=fltarr(n_elements(fn))
yshift=fltarr(n_elements(fn))
for i=0,n_elements(fn)-1 do begin
	print,file_basename(fn[i])

	data=mrdfits(fn[i],0,hdr,/silent)
	sz=size(data,/dim)

	wave=(findgen(sz[2])-sxpar(hdr,'crpix3')+1)*sxpar(hdr,'cd3_3')+sxpar(hdr,'crval3')
	qwave=where(wave gt wavebin[0] and wave le wavebin[1])
	
	sxdelpar,hdr,'cd3_3'
	sxdelpar,hdr,'crval3'
	sxdelpar,hdr,'crpix3'
	sxdelpar,hdr,'naxis3'
	sxaddpar,hdr,'naxis',2

	thum=fltarr(sz[0],sz[1])
	for ii=0,sz[0]-1 do begin
		for jj=0,sz[1]-1 do begin
			q=where(data[ii,jj,qwave] ne 0 and finite(data[ii,jj,qwave]) eq 1)
			thum[ii,jj]=mean(data[ii,jj,qwave[q]])
		endfor
	endfor

	q=where(thum ne 0)
	qx=q mod sz[0]
	qy=q/sz[0]
	xrange=[min(qx),max(qx)]
	yrange=[min(qy),max(qy)]
	hextract,thum,hdr,xrange[0],xrange[1],yrange[0]+trim[0,i],yrange[1]-trim[1,i],/silent
	;cgimage,thum<40,/scale,/axes

	; apply pre-shifts
	if preshiftfn ne '' then begin
		q=where(fn[i] eq prefn)
		if q[0] ne -1 then begin
			sxaddpar,hdr,'crval1',sxpar(hdr,'crval1')+prera[q[0]]/3600.
			sxaddpar,hdr,'crval2',sxpar(hdr,'crval2')+predec[q[0]]/3600.

			;sxaddpar,hdr,'crval1',sxpar(hdr,'crval1')+2.1/3600.
			;sxaddpar,hdr,'crval2',sxpar(hdr,'crval2')-9.6/3600.
	
			;sxaddpar,hdr,'crpix1',sxpar(hdr,'crpix1')+4
			;sxaddpar,hdr,'crpix2',sxpar(hdr,'crpix2')+4
		endif
	endif

	; initial projection
	hastrom,thum,hdr,newthum,newhdr,hdr0,missing=0,interp=2,cubic=-0.5
	data0_thum[*,*,i]=newthum


	if i eq 0 then begin
		hastrom,thum,hdr,newthum,newhdr,hdr0,missing=0,interp=2,cubic=-0.5
		;cgimage,newthum<40,/scale,/axes
		data_thum[*,*,i]=newthum
	endif else begin
	
		if ~keyword_set(noalign) then begin
			; iteration 1
			; +/- 4 pixels
			dy=((findgen(17)-8)/2.) ## (fltarr(17)+1)
			dx=(fltarr(17)+1) ## ((findgen(17)-8)/2.)
			rmss=fltarr(17,17)
			crls=fltarr(17,17)

			for ii=0,16 do begin
				for jj=0,16 do begin
					hdrshift=hdr
					sxaddpar,hdrshift,'crpix1',sxpar(hdr,'crpix1')+dx[ii,jj]
					sxaddpar,hdrshift,'crpix2',sxpar(hdr,'crpix2')+dy[ii,jj]
	
					hastrom,thum,hdrshift,newthum,newhdr,hdr0,missing=0,$
						interp=2,cubic=-0.5
					
					crls[ii,jj]=correlate(data0_thum[box[0]:box[1],$
						box[2]:box[3],0],$
						newthum[box[0]:box[1],box[2]:box[3]])
					;residual=data0_thum[*,*,0]-newthum
					;res=residual[box[0]:box[1],box[2]:box[3]]
					;rmss[ii,jj]=rms(res)
				endfor
			endfor
			
			cgimage,crls,/scale,/axes,title=file_basename(fn[i])+' 1'
	
			q=where(crls eq max(crls))
			ii1=q[0] mod 17
			jj1=q[0]/17
			xshift[i]=dx[ii1,jj1]
			yshift[i]=dy[ii1,jj1]

	
			; iteration 2
			; +/- 1 pixel
			dy=((findgen(17)-8)/8.+yshift[i]) ## (fltarr(17)+1)
			dx=(fltarr(17)+1) ## ((findgen(17)-8)/8.+xshift[i])
			crls=fltarr(17,17)
			
			for ii=0,16 do begin
				for jj=0,16 do begin
					hdrshift=hdr
					sxaddpar,hdrshift,'crpix1',sxpar(hdr,'crpix1')+dx[ii,jj]
					sxaddpar,hdrshift,'crpix2',sxpar(hdr,'crpix2')+dy[ii,jj]
	
					hastrom,thum,hdrshift,newthum,newhdr,hdr0,missing=0,$
						interp=2,cubic=-0.5
	
					crls[ii,jj]=correlate(data0_thum[box[0]:box[1],$
						box[2]:box[3],0],$
						newthum[box[0]:box[1],box[2]:box[3]])
				endfor
			endfor
	
			cgimage,crls,/scale,/axes,title=file_basename(fn[i])+' 2'
	
			q=where(crls eq max(crls))
			ii1=q[0] mod 17
			jj1=q[0]/17
			xshift[i]=dx[ii1,jj1]
			yshift[i]=dy[ii1,jj1]


			hdrshift=hdr
			sxaddpar,hdrshift,'crpix1',sxpar(hdr,'crpix1')+xshift[i]
			sxaddpar,hdrshift,'crpix2',sxpar(hdr,'crpix2')+yshift[i]

			hastrom,thum,hdrshift,newthum,newhdr,hdr0,missing=0,interp=2
			data_thum[*,*,i]=newthum
		endif
	endelse

	print,xshift[i],yshift[i]

endfor

writecol,file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.shift.list',fn0,$
	replicate(string(9b),n_elements(fn)),xshift,yshift


sxaddpar,hdr0,'naxis',3
sxaddpar,hdr0,'naxis3',n_elements(fn)
writefits,file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.thum.fits',data_thum,hdr0
writefits,file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.thum0.fits',data0_thum,hdr0

endeps
loadct,0


end
