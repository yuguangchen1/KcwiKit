pro kcwi_combinestd,listfn,noplot=noplot,extension=ext

if ~keyword_set(noplot) then noplot=0
if ~keyword_set(ext) then ext='' else begin
	if ext[0] ne '_' then ext='_'+ext
endelse

readcol,listfn,invsensfn,format='a',comment='#',/silent


for i=0,n_elements(invsensfn)-1 do begin
	data=mrdfits(invsensfn[i],0,hdr,/silent)
	if i eq 0 then begin
		data0=data
		hdr0=hdr
	endif
	
	; Check instrument setup
	ifunam=strtrim(sxpar(hdr,'IFUNAM'))
	gratnam=strtrim(sxpar(hdr,'BGRATNAM'))
	cwave=round(sxpar(hdr,'BCWAVE'))
	pwave=round(sxpar(hdr,'BPWAVE'))
	if i eq 0 then begin
		ifunam0=ifunam
		gratnam0=gratnam
		cwave0=cwave
		pwave0=pwave
	endif
	if ifunam ne ifunam0 or gratnam ne gratnam0 or $
		cwave ne cwave0 or pwave ne pwave0 then begin
		print,'Warning: Inconsistent instrument setup.'
		print,ifunam0,gratnam0,cwave0,pwave0
		print,file_basename(invsensfn[i]),ifunam,gratnam,cwave,pwave
	endif


	; Check wavelength
	wave=(findgen(sxpar(hdr,'NAXIS1'))-sxpar(hdr,'CRPIX1')+1)*$
		sxpar(hdr,'CDELT1')+sxpar(hdr,'CRVAL1')
	dw=wave[1]-wave[0]
	if i eq 0 then begin
		wave0=wave
		dw0=dw
	endif
	
	; Resample
	if wave0[0] ne wave[0] or dw0 ne dw or $
		n_elements(wave0) ne n_elements(wave) then begin
		linterp,wave,data[*,0],wave0,invsens_0
		linterp,wave,data[*,1],wave0,invsens_1
	endif else begin
		invsens_0=data[*,0]
		invsens_1=data[*,1]
	endelse

	; Master array
	if i eq 0 then begin
		invsens_0_all=fltarr(n_elements(wave0),n_elements(invsensfn))
		invsens_1_all=fltarr(n_elements(wave0),n_elements(invsensfn))
	endif
	invsens_0_all[*,i]=invsens_0
	invsens_1_all[*,i]=invsens_1
	
endfor

if n_elements(invsensfn) ne 1 then begin
	data0[*,0]=mean(invsens_0_all,dimension=2)
	data0[*,1]=mean(invsens_1_all,dimension=2)
endif else begin
	data0[*,0]=invsens_0_all
	data0[*,1]=invsens_0_all
endelse

writefits,file_basename(listfn,'.list')+'_invsens'+ext+'.fits',data0,hdr0

; Plot effective aperture
if not noplot then begin

	device,get_decomposed=color_flag
	device,decomposed=0
	loadct,13

	lg_titles=strarr(n_elements(invsensfn))
	lg_colors=intarr(n_elements(invsensfn))
	for i=0,n_elements(invsensfn)-1 do begin
		data=mrdfits(repstr(invsensfn[i],'_invsens','_ea'),0,hdr,/silent)
		wave=(findgen(sxpar(hdr,'NAXIS1'))-sxpar(hdr,'CRPIX1')+1)*$
			sxpar(hdr,'CDELT1')+sxpar(hdr,'CRVAL1')
		if i eq 0 then begin
			cgplot,wave,data[*,0],yr=[0,2.5e5],psym=10,$
				color=long(float(i)/float(n_elements(invsensfn))*255)
		endif else begin
			cgoplot,wave,data[*,0],psym=10,$
				color=long(float(i)/float(n_elements(invsensfn))*255)
		endelse
		cgoplot,wave,data[*,1],color=long(float(i)/float(n_elements(invsensfn))*255)

		qgood=where(wave gt sxpar(hdr0,'WAVGOOD0') and $
			wave lt sxpar(hdr0,'WAVGOOD1'))

		;lg_titles[i]=file_basename(invsensfn[i],'_invsens.fits')+' '+$
			;str(max(data[qgood,1]))
		;lg_colors[i]=long(float(i)/float(n_elements(invsensfn))*255)

	endfor

	cgoplot,replicate(sxpar(hdr0,'WAVGOOD0'),2),[0,1e6],linestyle=2
	cgoplot,replicate(sxpar(hdr0,'WAVGOOD1'),2),[0,1e6],linestyle=2

	;cglegend,titles=lg_titles,colors=lg_colors,alignment=0

	img=tvrd(true=1)
	write_png,file_basename(listfn,'.list')+ext+'.png',img

	device,decomposed=color_flag
endif




end
