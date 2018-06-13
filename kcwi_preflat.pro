pro kcwi_preflat,twinum,overwrite=overwrite,dir=dir
; This program hacks the KCWI-KDERP pipeline to normalize twiflat. 

if ~keyword_set(dir) then begin
	dir='redux/'
endif

if ~keyword_set(twinum) then begin
	print,'%KCWI_PREFLAT: Must specify image numbers.'
	return
endif
;twinum=35

; read config
ppfname=dir+'/kcwi.ppar'
procfname=dir+'/kcwi.proc'

ppar=kcwi_read_ppar(ppfname)
kpars=kcwi_read_proc(ppar,procfname,imgnum,count=nproc)

for i=0,n_elements(twinum)-1 do begin
	; hacking config
	index=where(imgnum eq twinum[i])
	obfil=kcwi_get_imname(kpars[index],imgnum[index],'_intd',/reduced)
	vfil=repstr(obfil,'_int','_var')

	kcfg=kcwi_read_cfg(obfil)

	gfile=repstr(strtrim(kpars[index].geomcbar,2),'_int','_geom')

	kpars[index].mastersky=repstr(obfil,'_intd','_sky')

	; read
	nocopy=0
	img=mrdfits(obfil,0,hdr,/fscale,/silent)
	var=mrdfits(vfil,0,vhdr,/fscale,/silent)
	if sxpar(hdr,'fnorm') eq 1 then begin
		if ~keyword_set(overwrite) then begin
			print,'Image '+string(twinum[i])+' has already been normalized...'
			continue
		endif else begin
			img=mrdfits(file_dirname(obfil)+'/old/'+file_basename(obfil),0,hdr,/silent)
			var=mrdfits(file_dirname(obfil)+'/old/'+file_basename(obfil),0,hdr,/silent)
			nocopy=1
		endelse
	endif

	; constructing flat model
	kcwi_make_sky,kpars[index],img,hdr,gfile,flat

	; normalize
	iflat=img
	vflat=var

	q=where(flat ne 0)
	iflat[q]=img[q]/flat[q]*5000.
	vflat[q]=var[q]/flat[q]^2*5000.^2
	
	sxaddpar,hdr,'fnorm',1
	sxaddpar,vhdr,'fnorm',1

	if nocopy eq 0 then begin
		if ~file_test(file_dirname(obfil)+'/old') then begin
			spawn,'mkdir '+file_dirname(obfil)+'/old'
		endif
		spawn,'cp '+obfil+' '+file_dirname(obfil)+'/old/'
		spawn,'cp '+vfil+' '+file_dirname(obfil)+'/old/'
	endif
	spawn,'mv '+repstr(obfil,'_intd','_sky')+' '+repstr(obfil,'_intd','_twiflat')

	writefits,obfil,iflat,hdr
	writefits,vfil,vflat,hdr

endfor

end
