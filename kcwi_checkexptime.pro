pro kcwi_checkexptime,dir,redux=redux
; 


print,'Checking exposure time..'

if ~keyword_set(dir) then begin
	dir='./'
endif

if keyword_set(redux) then begin
	dir='./redux/'
endif

fn=file_search(dir+'/kb*.fits')
for i=0,n_elements(fn)-1 do begin
	hdr=headfits(fn[i])
	
	exptime=sxpar(hdr,'xposure')
	expend=sxpar(hdr,'date-end')
	expend=date_conv(expend,'MODIFIED')
	rdend=sxpar(hdr,'daterend')
	rdend=date_conv(rdend,'MODIFIED')
	if rdend le expend then begin
		print,file_basename(fn[i])
		expbeg=date_conv(sxpar(hdr,'date-beg'),'MODIFIED')
		expend=rdend-53.64/3600./24.
		exptime=(expend-expbeg)*3600.*24.
		
		sxaddpar,hdr,'xposure',exptime
		sxaddpar,hdr,'telapse',exptime+0.005
		sxaddpar,hdr,'ttime',exptime
		print,'   Setting to'+string(exptime)+'.'

		data=mrdfits(fn[i],0,/silent)
		
		if ~keyword_set(redux) then begin
			if ~file_test(dir+'/old') then begin
				spawn,'mkdir '+dir+'/old'
			endif
			spawn,'cp '+fn[i]+' '+dir+'/old'
		endif
		writefits,fn[i],data,hdr
	endif

endfor


end
