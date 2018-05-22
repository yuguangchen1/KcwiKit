function kcwi_quickastrom_project,img,hdr,img0,hdr0,dx,dy,newhdr=hdrshift

hdrshift=hdr

x0=sxpar(hdr,'crpix1')
y0=sxpar(hdr,'crpix2')
sxaddpar,hdrshift,'crpix1',x0+dx
sxaddpar,hdrshift,'crpix2',y0+dy

hastrom,img0,hdr0,newimg,newhdr,hdrshift,missing=0,interp=2

q=where(img eq 0)
if q[0] ne -1 then begin
	newimg[q]=0
endif
return,newimg

end


function kcwi_quickastrom_crl,img,hdr2d,img0,hdr0,dx,dy,crl=crl,$
	bestimg=bestimg,besthdr=besthdr

crl=fltarr(n_elements(dx),n_elements(dy))
for i=0,n_elements(dx)-1 do begin
	for j=0,n_elements(dy)-1 do begin
		newimg0=kcwi_quickastrom_project(img,hdr2d,img0,hdr0,dx[i],dy[j])
		crl[i,j]=correlate(img,newimg0)
	endfor
endfor
q=where(crl eq max(crl))
qx=q mod n_elements(dx)
qy=q/n_elements(dx)
xmax=dx[qx]
ymax=dy[qy]

bestimg=kcwi_quickastrom_project(img,hdr2d,img0,hdr0,mean(xmax),mean(ymax),newhdr=besthdr)

return,[mean(xmax),mean(ymax)]

end



; ===================================
; Update WCS to ground baesd images.
; ===================================
pro kcwi_quickastrom,cubefn,imgfn,wavebin=wavebin

suffix='_icubes.fits'
if strmatch(cubefn,'*_icubes.fits') eq 0 then begin
	suffix='_icubed.fits'
endif


parfn=file_dirname(cubefn)+'/'+file_basename(cubefn,suffix)+'.par'
par=kcwi_stack_readpar(parfn)

if par.ref_xy[0] eq -1 then begin
	print,'%KCWI_QUICKASTROM: Set reference x-y coordinate.'
	return
endif
if par.ref_ad[0] eq -1 then begin
	print,'%KCWI_QUICKASTROM: Set reference RA-DEC coordinate.'
	return
endif

if ~keyword_set(wavebin) then begin
	wavebin=par.wavebin
	if par.wavebin[0] eq -1 then begin
		wavebin=[4000.,5000.]
	endif
endif



cube=mrdfits(cubefn,0,hdr,/silent)
sz=size(cube,/dimension)
wave=(findgen(sxpar(hdr,'naxis3'))-sxpar(hdr,'crpix3')+1)*sxpar(hdr,'cd3_3')+$
	sxpar(hdr,'crval3')


; collapsing
qwave=where(wave gt wavebin[0] and wave le wavebin[1])
img=fltarr(sz[0],sz[1])
for i=0,sz[0]-1 do begin
	for j=0,sz[1]-1 do begin
		q=where(cube[i,j,qwave] ne 0)
		if q[0] ne -1 then begin
			img[i,j]=mean(cube[i,j,qwave[q]])
		endif
	endfor
endfor

; header
hdr2d=hdr
sxdelpar,hdr2d,'cd3_3'
sxdelpar,hdr2d,'crval3'
sxdelpar,hdr2d,'crpix3'
sxdelpar,hdr2d,'naxis3'
sxaddpar,hdr2d,'naxis',2


img0=mrdfits(imgfn,0,hdr0,/silent)
if sxpar(hdr0,'equinox') eq 0 then begin
	sxaddpar,hdr0,'equinox',2000
endif


; initial guess
sxaddpar,hdr2d,'crpix1',par.ref_xy[0]
sxaddpar,hdr2d,'crpix2',par.ref_xy[1]
sxaddpar,hdr2d,'crval1',par.ref_ad[0]
sxaddpar,hdr2d,'crval2',par.ref_ad[1]


; entry 1: 
dx=(findgen(11)-5)*2
dy=(findgen(11)-5)*2
tmp=kcwi_quickastrom_crl(img,hdr2d,img0,hdr0,dx,dy,crl=crl,bestimg=bestimg)
xmax=tmp[0]
ymax=tmp[1]

; entry 2:
dx=(findgen(7)-3)*0.5+xmax
dy=(findgen(7)-3)*0.5+ymax
tmp=kcwi_quickastrom_crl(img,hdr2d,img0,hdr0,dx,dy,crl=crl,bestimg=bestimg)
xmax=tmp[0]
ymax=tmp[1]

; entry 3:
dx=(findgen(7)-3)*0.1+xmax
dy=(findgen(7)-3)*0.1+ymax
tmp=kcwi_quickastrom_crl(img,hdr2d,img0,hdr0,dx,dy,crl=crl,bestimg=bestimg,$
	besthdr=besthdr)
xmax=tmp[0]
ymax=tmp[1]


; write thumbnail
thumbfn=file_dirname(cubefn)+'/'+file_basename(cubefn,'.fits')+'_wcs.thum.fits'
sxaddpar,besthdr,'bitpix',-32
writefits,thumbfn,img,besthdr

; write cube
wcsfn=file_dirname(cubefn)+'/'+file_basename(cubefn,'.fits')+'_wcs.fits'
hdrwcs=hdr
sxaddpar,hdrwcs,'crval1',sxpar(besthdr,'crval1')
sxaddpar,hdrwcs,'crval2',sxpar(besthdr,'crval2')
sxaddpar,hdrwcs,'crpix1',sxpar(besthdr,'crpix1')
sxaddpar,hdrwcs,'crpix2',sxpar(besthdr,'crpix2')
writefits,wcsfn,cube,hdrwcs



end
