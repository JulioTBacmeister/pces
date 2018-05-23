 nx=0l
 ny=0l
 nz=0l

close,1
openr,1,/f77_u,'../fort.611'

conread:

readu,1,nx,ny,nz
xx=fltarr(nx)
yy=fltarr(ny)
zz=fltarr(nz)
readu,1,xx,yy,zz

ur2=fltarr(nx,nz)
zax=fltarr(nx,nz)
KE=fltarr(nx,nz)
uczr=fltarr(nx,nz)
uczz=fltarr(nx,nz)
ducz=fltarr(nx,nz)
dke=fltarr(nx,nz)
dkeog=fltarr(nx,nz)

readu,1,ur2,zax,ke
readu,1,uczr,uczz
readu,1,ducz,dke,dkeog

ur3=fltarr(nx,ny,nz)
readu,1,ur3



STOP
 goto,conread

end


