; Collapse data cubes to pseudo-broad band images, using default  wavelength range

pro kcwi_collapse,cubefn,wavebin=wavebin,usepix=usepix,suffix=suffix,var=var


; default wavelength range
tab_grating=['BL','BM']
tab_wave=[500,300]


if ~keyword_set(suffix) then begin
	suffix='icube'
endif

if ~keyword_set(cubefn) then begin
	dirflag=1
	dir='redux/'
	if file_test(dir) eq 0 then begin
		dir='./'
	endif


	cubefn=file_search(dir+'/*_'+suffix+'.fits')
endif

if keyword_set(wavebin) then begin
	if n_elements(wavebin) eq 1 then begin
		wavebin=replicate(wavebin,n_elements(cubefn))
	endif
endif


if strmatch(suffix,'vcube*') then var=1


; Loop the cubes
for i=0,n_elements(cubefn)-1 do begin
	print,file_basename(cubefn[i])

	cube=mrdfits(cubefn[i],0,hdr,/silent)
	sz=size(cube,/dimension)
	wave=(findgen(sz[2])-sxpar(hdr,'crpix3')+1)*sxpar(hdr,'cd3_3')+sxpar(hdr,'crval3')
	cwave=round(sxpar(hdr,'bcwave'))
	
	; get wave pixels
	if keyword_set(wavebin) then begin
		if ~keyword_set(usepix) then begin
			wrange=[cwave-wavebin[i],(cwave+wavebin[i]<5500)]
			qwave=where(wave gt wrange[0] and wave lt wrange[1])
		endif else begin
			tmp=min(abs(wave-cwave),qtmp)
			qwave=where(abs(findgen(sz[2])-qtmp) lt wavebin[i])
		endelse
	endif else begin
		grat=strtrim(sxpar(hdr,'bgratnam'),2)
		q=where(tab_grating eq grat)
		if q[0] eq -1 then begin
			wrange=[sxpar(hdr,'WAVGOOD0'),sxpar(hdr,'WAVGOOD1')<5500]
		endif else begin
			wrange=[cwave-tab_wave[q[0]],(cwave+tab_wave[q[0]])<5500]
		endelse
		qwave=where(wave gt wrange[0] and wave lt wrange[1])
	endelse


	; collapse
	if ~keyword_set(var) then begin
		img=median(cube[*,*,qwave],dimension=3)
	endif else begin
		img=mean(cube[*,*,qwave],dimension=3)/n_elements(qwave)
	endelse

	; header
	sxaddpar,hdr,'NAXIS',2
	sxdelpar,hdr,'NAXIS3'
	sxdelpar,hdr,'CD3_3'
	sxdelpar,hdr,'CTYPE3'
	sxdelpar,hdr,'CUNIT3'
	sxdelpar,hdr,'CNAME3'
	sxdelpar,hdr,'CRVAL3'
	sxdelpar,hdr,'CRPIX3'

	
	writefn=file_dirname(cubefn[i])+'/'+file_basename(cubefn[i],'.fits')+'.thum.fits'
	writefits,writefn,img,hdr
	
endfor




end
