pro kcwi_gcirc,coord0,coord1

tmp0=strsplit(coord0,' ',/extract)
tmp1=strsplit(coord1,' ',/extract)

ad0=[float(tmp0[0])+float(tmp0[1])/60.+float(tmp0[2])/3600.,$
	float(tmp0[3])+float(tmp0[4])/60.+float(tmp0[5])/3600.]
ad1=[float(tmp1[0])+float(tmp1[1])/60.+float(tmp1[2])/3600.,$
	float(tmp1[3])+float(tmp1[4])/60.+float(tmp1[5])/3600.]

gcirc,1,ad0[0],ad0[1],ad1[0],ad1[1],dis
print,'Tot:'+string(dis)

gcirc,1,ad0[0],ad0[1],ad0[0],ad1[1],dis_dec
if ad0[1] gt ad1[1] then dis_dec=-dis_dec

gcirc,1,ad0[0],ad0[1],ad1[0],ad0[1],dis_ra
if ad0[0] gt ad1[0] then dis_ra=-dis_ra

print,'en'+string(dis_ra)+string(dis_dec)


end
