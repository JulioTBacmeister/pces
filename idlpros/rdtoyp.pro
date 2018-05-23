        

nx=1000
nz=100
;**********************************************************
  close,1
  openr,1,/f77,'fort.211'
  b=fltarr(nx,nz  )
  p=fltarr(nx,nz+1)
  readu,1,b
  readu,1,p
 
 close,1

end
