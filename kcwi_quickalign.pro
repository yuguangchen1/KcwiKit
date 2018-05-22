
; quick alignment
pro kcwi_quickalign,fnlist,zrange=zrange,box=box,prefix=prefix,num=num,$
	pixscale=pixscale,dimension=dimension,preshiftfn=preshiftfn,$
	trim=trim
; Users can use prefix+num method instead of providing fnlist.
; Unit in preshiftfn is RA, DEC in arcsec.

parfn=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.par'
par=kcwi_stack_readpar(parfn)

if ~keyword_set(zrange) then zrange=[774,1300]
; TODO: substitute zrange to wavebin

if ~keyword_set(box) then begin
	box=par.box
	if box[0] eq -1 then begin
		box=[25,50,25,40]
	endif
endif
if ~keyword_set(pixscale) then pixscale=0.3
pixscale=pixscale/3600.
if ~keyword_set(dimension) then begin
	dimension=par.dimension
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
	fn=fn0+'_icubes.fits'
	;fn=dir+'kb170922_'+num+'_icuber_3d.fits'
endif else begin
	;readcol,'bx418.list',num,format='a'
	readcol,fnlist,fn0,trim1,trim2,format='a',comment='#'
	fn=fn0+'_icubes.fits'
	trim=fltarr(2,n_elements(fn))
	trim[0,*]=trim1
	trim[1,*]=trim2
endelse

if ~keyword_set(preshiftfn) then begin
	tmp=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.preshift.list'
	if file_test(tmp) then begin
		preshiftfn=tmp
	else begin
		tmp=file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.preshift.list'
		if file_test(tmp) then begin
			preshiftfn=tmp
		endif 
	endelse
endif
if keyword_set(preshiftfn) then begin
	readcol,preshiftfn,prefn,prera,predec,format='a,f,f',comment='#'
	prefn=prefn+'_icubes.fits'
endif

starteps,file_dirname(fnlist)+'/'+file_basename(fnlist,'.list')+'.astrom.ps',xs=6,ys=6

; construct WCS
hdr0=headfits(fn[0])
sxaddpar,hdr0,'naxis1',dimension[0]
sxaddpar,hdr0,'naxis2',dimension[1]
sxaddpar,hdr0,'crpix1',dimension[0]/2.
sxaddpar,hdr0,'crpix2',dimension[1]/2.
sxaddpar,hdr0,'cd1_1',pixscale
sxaddpar,hdr0,'cd2_2',pixscale
sxaddpar,hdr0,'cd1_2',0
sxaddpar,hdr0,'cd2_1',0
sxaddpar,hdr0,'naxis',2
sxdelpar,hdr0,'cd3_3'
sxdelpar,hdr0,'crval3'
sxdelpar,hdr0,'crpix3'
sxdelpar,hdr0,'naxis3'


data_thum=fltarr(dimension[0],dimension[1],n_elements(fn))
data0_thum=data_thum
xshift=fltarr(n_elements(fn))
yshift=fltarr(n_elements(fn))
for i=0,n_elements(fn)-1 do begin
	print,file_basename(fn[i])

	data=mrdfits(fn[i],0,hdr,/silent)
	sz=size(data,/dim)
	
	sxdelpar,hdr,'cd3_3'
	sxdelpar,hdr,'crval3'
	sxdelpar,hdr,'crpix3'
	sxdelpar,hdr,'naxis3'
	sxaddpar,hdr,'naxis',2

	thum=fltarr(sz[0],sz[1])
	for ii=0,sz[0]-1 do begin
		for jj=0,sz[1]-1 do begin
			thum[ii,jj]=mean(data[ii,jj,zrange[0]:zrange[1]])
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
	if keyword_set(preshiftfn) then begin
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
	hastrom,thum,hdr,newthum,newhdr,hdr0,missing=0,interp=2
	data0_thum[*,*,i]=newthum

	; +/- 4 pixels
	dy=((findgen(17)-8)/2.) ## (fltarr(17)+1)
	dx=(fltarr(17)+1) ## ((findgen(17)-8)/2.)
	rmss=fltarr(17,17)
	if i eq 0 then begin
		hastrom,thum,hdr,newthum,newhdr,hdr0,missing=0,interp=2
		;cgimage,newthum<40,/scale,/axes
		data_thum[*,*,i]=newthum
	endif else begin
	
		for ii=0,16 do begin
			for jj=0,16 do begin
				hdrshift=hdr
				sxaddpar,hdrshift,'crpix1',sxpar(hdr,'crpix1')+dx[ii,jj]
				sxaddpar,hdrshift,'crpix2',sxpar(hdr,'crpix2')+dy[ii,jj]

				hastrom,thum,hdrshift,newthum,newhdr,hdr0,missing=0,interp=2
				
				residual=data0_thum[*,*,0]-newthum
				res=residual[box[0]:box[1],box[2]:box[3]]
				rmss[ii,jj]=rms(res)
			endfor
		endfor
		
		cgimage,rmss,/scale,/axes,title=file_basename(fn[i])

		q=where(rmss eq min(rmss))
		ii1=q[0] mod 17
		jj1=q[0]/17
		xshift[i]=dx[ii1,jj1]
		yshift[i]=dy[ii1,jj1]

		hdrshift=hdr
		sxaddpar,hdrshift,'crpix1',sxpar(hdr,'crpix1')+xshift[i]
		sxaddpar,hdrshift,'crpix2',sxpar(hdr,'crpix2')+yshift[i]

		hastrom,thum,hdrshift,newthum,newhdr,hdr0,missing=0,interp=2
		data_thum[*,*,i]=newthum
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


end
