pro kcwi_makemask,dir

if ~keyword_set(dir) then begin
	dir='./'
endif

reduxdir=dir+'/redux/'

fn=file_search(reduxdir+'/kb*.reg')

for i=0,n_elements(fn)-1 do begin
	fitsfn=file_dirname(fn[i])+'/'+file_basename(fn[i],'.reg')+'_intf.fits'
	spawn,'python /scr/yuguangchen/Soft/kderp/kderp/devel/kcwi_masksky_ds9.py '+$
		fitsfn+' '+fn[i]

endfor


end
