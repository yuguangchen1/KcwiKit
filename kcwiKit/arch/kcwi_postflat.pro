function kcwi_postflat_readpar,parfn
; Read twiflat.ppar

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


result={interflat:0l,twifn:'',twiflat:lonarr(100)}
ele=result
flag=0
for i=0,n_elements(str)-1 do begin
	tmp=strsplit(str[i],'=',/extract)
	key=strtrim(tmp[0],2)
	value=strtrim(tmp[1],2)

	case key of
	'INTERFLAT':begin
		ele.interflat=long(value)
		ele.twifn=''
		ele.twiflat=lonarr(100)
		result=[result,ele]
		flag++
	end
	'TWIFN':begin
		result[flag].twifn=value
	end
	'TWIFLAT':begin
		twinum=strsplit(value,',',/extract)
		for j=0,n_elements(twinum)-1 do begin
			result[flag].twiflat[j]=long(twinum[j])
		endfor
	end
	endcase
endfor
result=result[1:n_elements(result)-1]

return,result

end




pro kcwi_postflat,dir=dir,parfn=parfn,overwrite=overwrite

; Initial setup
if ~keyword_set(dir) then begin
	dir='redux/'
	if file_test(dir) eq 0 then begin
		dir='./'
	endif
endif
if ~keyword_set(parfn) then begin
	parfn='twiflat.ppar'
endif

twistr=kcwi_postflat_readpar(dir+'/'+parfn)

; Read Kderp files
ppfname=dir+'/kcwi.ppar'
procname=dir+'/kcwi.proc'

ppar=kcwi_read_ppar(ppfname)
kpars=kcwi_read_proc(ppar,procfname,imgnum,count=nproc)


for iconfig=0,n_elements(twistr)-1 do begin

	;Get file names 
	flatnum=twistr[iconfig].interflat
	index=where(imgnum eq flatnum)
	interflatfn=kcwi_get_imname(kpars[index],imgnum[index],'_mflat',/reduced)

	
	if twistr[iconfig].twifn eq ''  then begin
		numlist=twistr[iconfig].twiflat
		q=where(numlist ne 0)
		if q[0] eq -1 then begin
			print,'[Error] Please indicate twiflat files in twiflat.ppar.'
			continue
		endif
		numlist=numlist[q]

		index0=where(imgnum eq numlist[0])
		fn0=kcwi_get_imname(kpars[index0],imgnum[index0],'_intf',/reduced)
		tfifn=repstr(fn0,'intf','tfimg')
		tfvfn=repstr(fn0,'intf','tfvar')
		twiskyfn=repstr(fn0,'intf','twisky')
		tflatfn=repstr(fn0,'intf','tflat')
	endif else begin
		tflatfn=twistr[iconfig].twifn
		tflatfn=file_dirname(interflatfn)+'/'+tflatfn
		tfifn=repstr(tflatfn,'tflat','tfimg')
		tfvfn=repstr(tflatfn,'tflat','tfvar')
		twiskyfn=repstr(tflatfn,'tflat','twisky')
	endelse


	;Check file existance
	if file_test(tflatfn) then begin
		iflat=mrdfits(tflatfn,/fscale,/silent)
	endif else begin
		; Make master flat
		for i=0,n_elements(numlist)-1 do begin
			index=where(imgnum eq numlist[i])
			fn=kcwi_get_imname(kpars[index],imgnum[index],'_intf',/reduced)
			vfn=repstr(fn,'_int','_var')
			mfn=repstr(fn,'_int','_msk')
	
			img=mrdfits(fn,0,hdr,/fscale,/silent)
			var=mrdfits(vfn,0,vhdr,/fscale,/silent)
			msk=mrdfits(mfn,0,mhdr,/fscale,/silent)

			if sxpar(hdr,'twiflat') eq 1 then begin
				img=mrdfits(file_dirname(fn)+'/old/'+file_basename(fn),$
					0,hdr,/fscale,/silent)
				var=mrdfits(file_dirname(vfn)+'/old/'+file_basename(vfn),$
					0,vhdr,/fscale,/silent)
			endif

			if i eq 0 then begin
				sz=size(img,/dim)
				icube=fltarr(sz[0],sz[1],n_elements(numlist))
				vcube=icube
				mcube=icube
			endif

			icube[*,*,i]=img
			vcube[*,*,i]=var
			mcube[*,*,i]=msk

		endfor

		istack=fltarr(sz[0],sz[1])
		vstack=fltarr(sz[0],sz[1])
		for i=0,sz[0]-1 do begin
			for j=0,sz[1]-1 do begin
				weight=abs(1./vcube[i,j,*])
				q=where(weight eq 0 or finite(weight) eq 0 or $
					(mcube[i,j,*] and 240) eq 0 and (mcube[i,j,*] ne 0))
				if q[0] ne -1 then begin
					weight[q]=0
				endif
			
				istack[i,j]=total(icube[i,j,*]*weight)/total(weight)
				vstack[i,j]=total(weight^2*vcube[i,j,*])/total(weight)^2
			endfor
		endfor

		; Write stack
		writefits,tfifn,istack,hdr0
		writefits,tfvfn,vstack,hdr0


		; Make sky model
		gfn=repstr(strtrim(kpars[index0].geomcbar,2),'_int','_geom')
		kpars[index0].mastersky=repstr(fn0,'_intf','_sky')
		kcwi_make_sky,kpars[index0],istack,hdr,gfn,flat

		; Make master flat
		iflat=fltarr(sz[0],sz[1])+1
		q=where(flat ne 0 and istack ne 0)
		if q[0] ne -1 then begin
			iflat[q]=flat[q]/istack[q]
		endif
		q=where(finite(iflat) eq 0)
		if q[0] ne -1 then begin
			iflat[q]=1
		endif


		; Write sky model and master flat
		spawn,'mv '+repstr(fn0,'_intf','_sky')+' '+twiskyfn
		writefits,tflatfn,iflat,hdr

	endelse


	; Update intf files
	index=where(kpars.masterflat eq interflatfn)
	for i=0,n_elements(index)-1 do begin
		nocopy=0

		fn=kcwi_get_imname(kpars[index[i]],imgnum[index[i]],'_intf',/reduced)
		vfn=repstr(fn,'_int','_var')
	
		if file_test(fn) eq 0 then begin
			print,'Warning: File '+file_basename(fn)+' not found.'
			continue
		endif

		print,file_basename(fn)+':'

		img=mrdfits(fn,0,hdr,/fscale,/silent)
		var=mrdfits(vfn,0,vhdr,/fscale,/silent)
		
		if sxpar(hdr,'twiflat') eq 1 then begin
			if keyword_set(overwrite) then begin
				img=mrdfits(file_dirname(fn)+'/old/'+file_basename(fn),$
					0,hdr,/fscale,/silent)
				var=mrdfits(file_dirname(vfn)+'/old/'+file_basename(fn),$
					0,vhdr,/fscale,/silent)
				nocopy=1
			endif else begin
				print,' Already processed. Skipping...'
				continue
			endelse
		endif


		q=where(iflat ne 0)
		img[q]=img[q]*iflat[q]
		var[q]=var[q]*iflat[q]^2

		q=where(finite(img) eq 0)
		if q[0] ne -1 then begin
			img[q]=0
		endif

		sxaddpar,hdr,'twiflat',1
		sxaddpar,hdr,'twiflat',1
	
		if nocopy eq 0 then begin
			if ~file_test(file_dirname(fn)+'/old') then begin
				spawn,'mkdir '+file_dirname(fn)+'/old'
			endif
			spawn,'cp '+fn+' '+file_dirname(fn)+'/old/'
			spawn,'cp '+vfn+' '+file_dirname(vfn)+'/old/'

		endif
	
		writefits,fn,img,hdr
		writefits,vfn,var,vhdr
	endfor

	
endfor

end
