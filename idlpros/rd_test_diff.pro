close,1
 openr,1,'../fort.211',/f77_u
 
 nz=0l
 readu,1,nz
xx=fltarr(nz,2000)
x=fltarr(nz)
rho=fltarr(nz)
readu,1,rho
readu,1,x
xx(*,0)=x
i=1
while not eof(1) do begin

    readu,1,x
    xx(*,i)=x
    i=i+1

 endwhile

end
