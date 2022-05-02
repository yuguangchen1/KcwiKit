function kcwi_medfilter_readpar,parfn,parflag=parflag,fnonly=fnonly
; read parameter files if present

fnonly=-1
result={fnnum:0l,zbin:0l,ybin:0l,ytrim:[0l,0l],skipflag:0l}

if file_test(parfn) eq 0 then begin
	parflag=0
	return,result
endif else begin
	parflag=1
endelse

; read par file into string
openr,lun,parfn,/get_lun
tmp=''
str=''
while ~eof(lun) do begin
	readf,lun,tmp
	tmp=strtrim(tmp,2)
	if tmp ne '' then begin
		str=[str,tmp]
	endif
endwhile
str=str[1:n_elements(str)-1]
close,lun
free_lun,lun


; convert string to structure
ele=result
flag=0
for i=0,n_elements(str)-1 do begin
	tmp=strsplit(str[i],'=',/extract)
	key=strtrim(tmp[0],2)
	value=strtrim(tmp[1],2)

	case key of
	'FNONLY':begin
		fnonly=long(strsplit(value,',',/extract))
	end

	'FNNUM':begin
		ele.fnnum=long(value)
		ele.ybin=0l
		ele.zbin=0l
		ele.ytrim=[0l,0l]
		result=[result,ele]
		flag++
	end
	'ZBIN':begin
		; filter size in the wavelength direction
		result[flag].zbin=long(value)
	end
	'YBIN':begin
		; filter size in the y direction
		result[flag].ybin=long(value)
	end
	'YTRIM':begin
		; number of pixs to be trimmed on the edge of slicers
		trimstr=strsplit(value,',',/extract)
		for j=0,1 do begin
			result[flag].ytrim[j]=long(trimstr[j])
		endfor
	end
	'SKIPFLAG':begin
		; Ignore? 1-Skipping; 2-Subtract single sky value
		result[flag].skipflag=long(value)
	end
	endcase
endfor
if n_elements(result) gt 1 then begin
	result=result[1:n_elements(result)-1]
endif

return,result

end


function kcwi_medfilter_linspl,qv,v
; Determine the interpolation method.

if qv[1] ne qv[0]+1 then return,0

if qv[n_elements(qv)-2] ne qv[n_elements(qv)-1] then return,0

return,1


end



pro kcwi_medfilter2_py,dir=dir,parfn=parfn,zbin=zbin,ybin=ybin,ytrim=ytrim,suffix=suffix,$
	overwrite=overwrite
; Run median filter on cubes

if ~keyword_set(dir) then begin
	dir='redux/'
	if file_test(dir) eq 0 then begin
		dir='./'
	endif
endif
if ~keyword_set(parfn) then begin
	parfn='medfilter.ppar'
endif
if ~keyword_set(suffix) then begin
	suffix='icubes'
endif

medpar=kcwi_medfilter_readpar(dir+'/'+parfn,parflag=parflag,fnonly=fnonly)
; print,medpar

; Read Kderp files
; ppfname=dir+'/kcwi.ppar'
; procname=dir+'/kcwi.proc'
; ppar=kcwi_read_ppar(ppfname)
; kpars=kcwi_read_proc(ppar,procfname,imgnum,count=nproc)
;print,kpars
;for i=75,76 do begin
; print,'n_elements',n_elements(imgnum)

; get files in the directory (assume we're in redux)
files = FILE_SEARCH('*_icubes.fits',/FULLY_QUALIFY_PATH)

for i=0,n_elements(files)-1 do begin
	; print,'this is i ',i
	; print,'now the file name ',files[i]
	if fnonly[0] ne -1 then begin
		q=where(fnonly eq i)
		if q[0] eq -1 then continue
	endif

	nocopy=0

	; get file names
	cubefn = files[i] ;kcwi_get_imname(kpars[i],imgnum[i],'_'+suffix,/reduced)
	; print,'this is the cubefn',cubefn
	; print,file_basename(cubefn)
	if file_test(cubefn) eq 0 then continue
	; mcubefn=repstr(cubefn,'_icube','_mcube')
	mimgfn1=file_dirname(cubefn)+'/'+file_basename(cubefn,'.fits')+'.stackmask.fits'
	print,cubefn
	print,mimgfn1

	; mimgfn2=file_dirname(cubefn)+'/'+file_basename(cubefn,'.fits')+'_2d.mask.fits'
	; print,mcubefn,mimgfn1,mimgfn2
	; read files
	cube=mrdfits(cubefn,0,hdr,/silent)
	if strtrim(sxpar(hdr,'imtype'),2) ne 'OBJECT' then begin
		continue
	endif

	print,file_basename(cubefn)

	; check header
	; if sxpar(hdr,'medfilt') eq 1 then begin
	; 	if keyword_set(overwrite) then begin
	; 		cube=mrdfits(file_dirname(cubefn)+'/old/'+file_basename(cubefn),$
	; 			0,hdr,/silent)
	; 		nocopy=1
	; 		print,' Overwrite=True, Using original icube file'
	; 	endif else begin
	; 		print,' Already processed. Skipping...'
	; 		continue
	; 	endelse
	; endif


	cube0=cube
	sz=size(cube,/dimension)
	cube0=cube
	mcube=mrdfits(cubefn,'mask',mhdr,/silent)

	vcube=mrdfits(cubefn,'uncert',vhdr,/silent) ; rest of the extensions for later
	fcube=mrdfits(cubefn,'flags',fhdr,/silent) ; "

	if file_test(mimgfn1) then begin
		mimg1=mrdfits(mimgfn1,0,/silent)
	endif else begin
		mimg1=fltarr(sz[0],sz[1])
	endelse
	; if file_test(mimgfn2) then begin
	; 	mimg2=mrdfits(mimgfn2,0,/silent)
	; endif else begin
	mimg2=0.
	; endelse


	wave=(findgen(sz[2])-sxpar(hdr,'crpix3')+1)*sxpar(hdr,'cd3_3')+sxpar(hdr,'crval3')
	cwave=round(sxpar(hdr,'bcwave'))

	; get defaults
	ytrim1=[4,4]
	ybin1=16
	zbin1=100

	q=where(medpar.fnnum eq i)
	if q[0] ne -1 then begin
		if medpar[q[0]].skipflag eq 1 then begin
			print,' SKIPFLAG == 1. Skipping...'
			continue
		endif
		if medpar[q[0]].skipflag eq 2 then begin
			print,' SKIPFLAG == 2. Subtracting single value.'
			cube1=cube0
			medcube=fltarr(sz[0],sz[1],sz[2])
;			goto,skipflag2
		endif
		if medpar[q[0]].ytrim[0] ne 0 then begin
			ytrim1=medpar[q[0]].ytrim
		endif
		if medpar[q[0]].ybin ne 0 then begin
			ybin1=medpar[q[0]].ybin
		endif
		if medpar[q[0]].zbin ne 0 then begin
			zbin1=medpar[q[0]].zbin
		endif
	endif
	if keyword_set(ytrim) then begin
		ytrim1=ytrim
	endif
	if keyword_set(ybin) then begin
		ybin1=ybin
	endif
	if keyword_set(zbin) then begin
		zbin1=zbin
	endif



	; masking
	q=where(mcube ne 0)
	if q[0] ne -1 then begin
		cube[q]=0.
	endif

	q=where(~finite(mimg1), /null) ; where there are NaNs
	print,q[0]
	if q[0] ne -1 then begin
		cube[q]=0.
	endif

	;save 0 values
	q0=where(cube eq 0)

	if n_elements(mimg2) ne 1 then begin
		; resemble cube - copied from kcwi_flatten_cube
		mcube2=fltarr(sz[0],sz[1],sz[2])
		trm=0
		for ii=0,sz[0]-1 do begin
			ix0 = 0 + trm
			ix1 = sz[1] - (trm + 1)
			ox0 = ii*sz[1] + trm
			ox1 = (ii+1)*sz[1] - (trm + 1)
			mcube2[ii,ix0:ix1,*] = mimg2[ox0:ox1,*]
		endfor
		q=where(mcube2 eq 1)
		cube[q]=-100
	endif


	q=where(mimg1 eq 1)
	if q[0] ne -1 then begin
		for kk=0,sz[2]-1 do begin
			tmp=cube[*,*,kk]
			tmp[q]=-100.	; flag for interpolation
			cube[*,*,kk]=tmp
		endfor
	endif
	cube[*,0:ytrim1[0]-1,*]=-200   ; flag for nearest neighbor
	cube[*,sz[1]-ytrim1[1]:sz[1]-1,*]=-200

	cube[q0]=0


	; looping in cube
	medcube=fltarr(sz[0],sz[1],sz[2])
	for ii=0,sz[0]-1 do begin
		for jj=0,sz[1]-1 do begin
			for kk=0,sz[2]-1 do begin
				if cube[ii,jj,kk] eq 0 then continue
				if cube[ii,jj,kk] eq -100 then begin
					medcube[ii,jj,kk]=-100
					continue
				endif
				if cube[ii,jj,kk] eq -200 then begin
					medcube[ii,jj,kk]=-200
					continue
				endif

				; median filtering
				yrange=[(ytrim1[0] > (jj-ybin1/2)),$
					((sz[1]-ytrim1[1]-1) < (jj+ybin1/2))]
				zrange=[(0 > (kk-zbin1/2)),$
					((sz[2]-1) < (kk+zbin1/2))]
				tmpcube=cube[ii,yrange[0]:yrange[1],zrange[0]:zrange[1]]
				q=where(tmpcube ne 0 and tmpcube ne -100 and tmpcube ne -200)
				medcube[ii,jj,kk]=median(tmpcube[q])
			endfor
		endfor
	endfor
	;writefits,'tmp.fits',medcube

	; interpolation
	for ii=0,sz[0]-1 do begin
		for kk=0,sz[2]-1 do begin
			line0=medcube[ii,*,kk]
			line=line0

			qv=where(line ne -100 and line ne -200 and line ne 0)
			if qv[0] eq -1 then begin
				q100=where(line eq -100 or line eq -200)
				if q100[0] ne -1 then begin
					line[q100]=-300 ; linear interp in wavelength
				endif
				medcube[ii,*,kk]=line
			endif else begin
				qm=[min(qv),max(qv)]

				; cubic interpolate
				q100=where(line eq -100)
				if q100[0] ne -1 then begin
					q=where(q100 lt qm[0])
					if q[0] ne -1 then begin
						if n_elements(q) eq n_elements(q100) then begin
							line[q100]=0
							goto,escape
						endif else begin
							line[q100[q]]=-200
							remove,q,q100
						endelse
					endif
					q=where(q100 gt qm[1])
					if q[0] ne -1 then begin
						if n_elements(q) eq n_elements(q100) then begin
							line[q100]=0
							goto,escape
						endif else begin
							line[q100[q]]=-200
							remove,q,q100
						endelse
					endif

					if n_elements(qv) gt 3 then begin
						; linear or spline?
						flag=kcwi_medfilter_linspl(qv,line[qv])
						if flag eq 0 then begin
							tmp=interpol(line[qv],qv,q100)
						endif else begin
							tmp=spline(qv,line[qv],q100)
						endelse
						;if ii eq 13 and kk eq 1642 then stop
						line[q100]=tmp
					endif else begin
						line[q100]=-200
					endelse
				endif
				escape:

				;if ii eq 13 and kk eq 1642 then stop

				; nearest neighbor
				q200=where(line eq -200)
				if q200[0] ne -1 then begin
					if qv[0] ne -1 then begin
						for jj=0,n_elements(q200)-1 do begin
							tmp=min(abs(qv-q200[jj]),q)
							line[q200[jj]]=line[qv[q]]
						endfor
					endif
				endif

				medcube[ii,*,kk]=line
			endelse

		endfor
	endfor


	; linear in wavelength
	index3=where(medcube eq -300)
	if index3[0] ne -1 then begin
		for ii=0,sz[0]-1 do begin
			for jj=0,sz[1]-1 do begin
				tmpmed0=medcube[ii,jj,*]
				tmpmed=tmpmed0
				qlin3=where(tmpmed eq -300)
				if qlin3[0] eq -1 then continue

				qlin0=where(tmpmed ne -300 and tmpmed ne -200 and $
					tmpmed ne -100 and tmpmed ne 0)
				qq=where(qlin3 lt min(qlin0))
				if qq[0] ne -1 then begin
					tmpmed[qlin3[qq]]=tmpmed[min(qlin0)]
				endif
				qq=where(qlin3 gt max(qlin0))
				if qq[0] ne -1 then begin
					tmpmed[qlin3[qq]]=tmpmed[max(qlin0)]
				endif

				qlin3=where(tmpmed eq -300)
				if qlin3[0] ne -1 then begin
					tmpinter=interpol(tmpmed[qlin0],qlin0,qlin3)
					tmpmed[qlin3]=tmpinter
				endif

				medcube[ii,jj,*]=tmpmed

			endfor
		endfor
	endif
	;writefits,'tmp.fits',medcube


	; subtract
	cube1=cube0-medcube
	;writefits,'tmp.fits',cube1


	; further offset background
	skipflag2:
	tab_grating=['BL','BM']
	tab_wave=[500,300]

	grat=strtrim(sxpar(hdr,'bgratnam'),2)
	q=where(tab_grating eq grat)
	if q[0] eq -1 then begin
		wrange=[(cwave-500)>3500,(cwave+500)<5500]
	endif else begin
		wrange=[cwave-tab_wave[q[0]],(cwave+tab_wave[q[0]])<5500]
	endelse
	qwave=where(wave gt wrange[0] and wave lt wrange[1])
	img1=mean(cube1[*,*,qwave],dimension=3)
	img1[*,0:ytrim1[0]-1]=0
	img1[*,sz[1]-ytrim1[1]:sz[1]-1]=0
	q=where(mimg1 eq 1)
	img1[q]=0
	q=where(img1 ne 0)

	;print,q -- something wrong with _icube.mask.fits
	;print,img1[q]

	bkg=median(img1[q])
	;print,bkg
	q=where(cube1 ne 0)
	cube1[q]=cube1[q]-bkg
	medcube[q]=medcube[q]+bkg


	; write
	sxaddpar,hdr,'medfilt',2
	if nocopy eq 0 then begin
		if ~file_test(file_dirname(cubefn)+'/old') then begin
			spawn,'mkdir '+file_dirname(cubefn)+'/old'
		endif
		spawn,'cp '+cubefn+' '+file_dirname(cubefn)+'/old'
	endif

	medfn=file_dirname(cubefn)+'/'+file_basename(cubefn,'.fits')+'.med.fits'
	writefits,medfn,medcube,hdr ; .med.fits
	writefits,cubefn,cube1,hdr ; _icube.fits
	writefits,cubefn,mcube,mhdr,/APPEND ; probably a better way to do this in a loop...
	writefits,cubefn,vcube,vhdr,/APPEND
	writefits,cubefn,fcube,fhdr,/APPEND




endfor




end
