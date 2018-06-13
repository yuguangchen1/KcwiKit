function kcwi_stack_readpar,parname

if file_test(parname) then begin
	openr,1,parname
	str=''
	tmp=''
	while ~eof(1) do begin
		readf,1,tmp
		str=[str,tmp]
	endwhile
	str=str[1:n_elements(str)-1]
	close,1
endif else begin
	str=' '
endelse

par={align_box:[-1l,-1l,-1l,-1l],align_dimension:[-1l,-1l],align_xpix:-1.,align_ypix:-1.,$
	align_orientation:-1000.,stack_dimension:[-1l,-1l],stack_xpix:-1.,stack_ypix:-1.,$
	stack_orientation:-1000.,wavebin:[-1.,-1.],ref_xy:[-1.,-1.],ref_ad:[-1d,-1d]}


; alignment keywords
q=where(strmatch(str,'align_box*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.align_box[0]=long(tmp[1])
	par.align_box[2]=long(tmp[2])
	par.align_box[1]=long(tmp[3])
	par.align_box[3]=long(tmp[4])
endif

q=where(strmatch(str,'align_dimension*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.align_dimension[0]=long(tmp[1])
	par.align_dimension[1]=long(tmp[2])
endif

q=where(strmatch(str,'align_xpix*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.align_xpix=float(tmp[1])
endif

q=where(strmatch(str,'align_ypix*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.align_ypix=float(tmp[1])
endif

q=where(strmatch(str,'align_orientation*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.align_orientation=float(tmp[1])
endif


; stacking keywords
q=where(strmatch(str,'stack_dimension*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.stack_dimension[0]=long(tmp[1])
	par.stack_dimension[1]=long(tmp[2])
endif

q=where(strmatch(str,'stack_xpix*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.stack_xpix=float(tmp[1])
endif

q=where(strmatch(str,'stack_ypix*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.stack_ypix=float(tmp[1])
endif

q=where(strmatch(str,'stack_orientation*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.stack_orientation=float(tmp[1])
endif


; overwrite both stack_* and align_* keywords
q=where(strmatch(str,'dimension*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.align_dimension[0]=long(tmp[1])
	par.align_dimension[1]=long(tmp[2])
	par.stack_dimension[0]=long(tmp[1])
	par.stack_dimension[1]=long(tmp[2])
endif

q=where(strmatch(str,'xpix*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.align_xpix=float(tmp[1])
	par.stack_xpix=float(tmp[1])
endif

q=where(strmatch(str,'ypix*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.align_ypix=float(tmp[1])
	par.stack_ypix=float(tmp[1])
endif

q=where(strmatch(str,'orientation*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.align_orientation=float(tmp[1])
	par.stack_orientation=float(tmp[1])
endif


; wavebin for alignment and astrometry
q=where(strmatch(str,'*wavebin*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.wavebin[0]=float(tmp[1])
	par.wavebin[1]=float(tmp[2])
endif


; astrometry
q=where(strmatch(str,'*ref_xy*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.ref_xy[0]=float(tmp[1])
	par.ref_xy[1]=float(tmp[2])
endif

q=where(strmatch(str,'*ref_ad*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.ref_ad[0]=double(tmp[1])
	par.ref_ad[1]=double(tmp[2])
endif


return,par


end
