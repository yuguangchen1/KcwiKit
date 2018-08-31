function kcwi_vachelio,cube,hdr_old,hdr_new,mask=mask,vcorr=vcorr
; Convert from air-wavelength to vacuum-heliocentric wavelengt for KCWI cubes. 
; TODO: Better handling for mask cubes. 
; 07/30/2018 - Yuguang Chen

if ~keyword_set(hdr_new) then begin
	hdr_new=hdr_old
endif

sz=[sxpar(hdr_new,'naxis1'),sxpar(hdr_new,'naxis2'),sxpar(hdr_new,'naxis3')]

naxis3=sxpar(hdr_old,'naxis3')
crval3=sxpar(hdr_old,'crval3')
crpix3=sxpar(hdr_old,'crpix3')
cd33=sxpar(hdr_old,'cd3_3')
wave_old=(findgen(naxis3)-crpix3+1)*cd33+crval3



naxis3=sxpar(hdr_new,'naxis3')
crval3=sxpar(hdr_new,'crval3')
crpix3=sxpar(hdr_new,'crpix3')
cd33=sxpar(hdr_new,'cd3_3')
wave_new=(findgen(naxis3)-crpix3+1)*cd33+crval3


; Vacuum wavelength
airtovac,wave_old,wave_vac

; Heliocentric wavelength
observatory,'keck',obs

targra=sxpar(hdr_old,'targra')
targdec=sxpar(hdr_old,'targdec')
jd=sxpar(hdr_old,'mjd')+2400000.5
vcorr=heliocentric(targra,targdec,jd=jd,longitude=obs.longitude,latitude=obs.latitude,$
	altitude=obs.altitude)
wave_hel=(1-vcorr/2.99792458e5)*wave_vac


; Resample
cube_new=dblarr(sz[0],sz[1],sz[2])
if keyword_set(mask) then cube_new=long(cube_new)
for i=0,sz[0]-1 do begin
	for j=0,sz[1]-1 do begin
		spec=cube[i,j,*]
		if ~keyword_set(mask) then begin
			y2deriv=spl_init(wave_hel,spec,/double)
			spec_new=spl_interp(wave_hel,spec,y2deriv,wave_new,/double)
		endif else begin
			spec_new=interpol(spec,wave_hel,wave_new)
			spec_new=ceil(spec_new)
		endelse
		cube_new[i,j,*]=spec_new
	endfor
endfor



return,cube_new


end
