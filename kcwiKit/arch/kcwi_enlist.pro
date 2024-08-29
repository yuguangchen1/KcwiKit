pro kcwi_enlist,list,pattern=pattern,fold_case=fold_case
if ~keyword_set(fold_case) then fold_case=0

openr,1,list
tmp=''
str=''
while ~eof(1) do begin
	readf,1,tmp
	str=[str,tmp]
endwhile
str=str[1:n_elements(str)-1]
close,1

obj=strarr(n_elements(str))
ra=fltarr(n_elements(str))
dec=ra
for i=0,n_elements(str)-1 do begin
	tmp=strsplit(str[i],/extract)
	obj[i]=tmp[0]
	ra[i]=float(tmp[1])+float(tmp[2])/60.+float(tmp[3])/3600.
	dec[i]=float(tmp[4])+float(tmp[5])/60.+float(tmp[6])/3600.
endfor

index=lindgen(n_elements(ra))
for i=1,n_elements(ra)-1 do begin
	q=where(ra[0:i-1] eq ra[i] and dec[0:i-1] eq dec[i])
	if q[0] ne -1 then begin
		qq=where(index ne i)
		index=index[qq]
	endif
endfor
obj=obj[index]
ra=ra[index]
dec=dec[index]

if keyword_set(pattern) then begin
	index=0
	for i=0,n_elements(pattern)-1 do begin
		q=where(strmatch(obj,pattern[i],fold_case=fold_case) eq 1)
		if q[0] ne -1 then begin
			index=[index,q]
		endif
	endfor
	index=index[1:n_elements(index)-1]

	index0=lindgen(n_elements(index))
	for i=1,n_elements(index)-1 do begin
		q=where(index[0:i-1] eq index[i])
		if q[0] ne -1 then begin
			qq=where(index0 ne i)
			index=index[qq]
		endif
	endfor
	obj=obj[index]
	ra=ra[index]
	dec=dec[index]
endif


xr=[max(ra)+(max(ra)-min(ra))/10.,min(ra)-(max(ra)-min(ra))/10.*2]
yr=[min(dec)-(max(dec)-min(dec))/10.,max(dec)+(max(dec)-min(dec))/10.]
cgplot,ra,dec,psym=16,color='red',xr=xr,yr=yr
for i=0,n_elements(obj)-1 do begin
	cgtext,ra[i],dec[i],obj[i],color='red'
endfor


npair=n_elements(obj)*(n_elements(obj)-1)/2.
dra=fltarr(npair)
ddec=fltarr(npair)
dis=fltarr(npair)
count=0
ra0=dblarr(npair)
ra1=dblarr(npair)
dec0=dblarr(npair)
dec1=dblarr(npair)
for j=1,n_elements(obj)-1 do begin
	for i=0,j-1 do begin
		gcirc,1,ra[i],dec[i],ra[j],dec[j],tmp
		dis[count]=tmp
	
		gcirc,1,ra[i],dec[i],ra[i],dec[j],tmp
		if dec[i] gt dec[j] then tmp=-tmp
		ddec[count]=tmp

		gcirc,1,ra[i],dec[i],ra[j],dec[i],tmp
		if ra[i] gt ra[j] then tmp=-tmp
		dra[count]=tmp

		ra0[count]=ra[i]
		ra1[count]=ra[j]
		dec0[count]=dec[i]
		dec1[count]=dec[j]

		cgarrow,ra0[count],dec0[count],ra1[count],dec1[count],color='black',/data,$
			hthick=3
		cgtext,(ra0[count]+ra1[count])/2.,(dec0[count]+dec1[count])/2.,$
			string(count,format='(i0)'),color='black'

		print,string(count,format='(i0)')+':'+'   en'+$
			string(dra[count])+string(ddec[count])
			
	
		count++
	endfor

endfor


end
