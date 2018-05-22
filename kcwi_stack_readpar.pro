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

par={box:[-1l,-1l,-1l,-1l],dimension:[-1l,-1l],pixscale:-1.,orientation:-1000.,$
	wavebin:[-1.,-1.],ref_xy:[-1.,-1.],ref_ad:[-1d,-1d]}

q=where(strmatch(str,'*box*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.box[0]=long(tmp[1])
	par.box[2]=long(tmp[2])
	par.box[1]=long(tmp[3])
	par.box[3]=long(tmp[4])
endif

q=where(strmatch(str,'*dimension*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.dimension[0]=long(tmp[1])
	par.dimension[1]=long(tmp[2])
endif

q=where(strmatch(str,'*pixscale*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.pixscale=float(tmp[1])
endif

q=where(strmatch(str,'*orientation*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.orientation=float(tmp[1])
endif

q=where(strmatch(str,'*wavebin*',/fold_case) eq 1)
if q[0] ne -1 then begin
	tmp=strsplit(str[q[0]],/extract)
	par.wavebin[0]=float(tmp[1])
	par.wavebin[1]=float(tmp[2])
endif


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
