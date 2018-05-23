
close,1
openr,1,'../wojtek.dat',/f77_u

nx=0L&ny=0L&nz=0L
time=0.
readu,1,time,nx,ny,nz
rh=fltarr(nz)
gac=fltarr(nz)
th0=fltarr(nz)
xx=fltarr(nx)
yy=fltarr(ny)
zz=fltarr(nz)

    readu,1, xx,yy,zz,rh,th0,gac


close,1
openr,1,'../fort.18',/f77_u


nx=100&ny=100&nz=100

time=0.
u3=fltarr(nx,ny,nz)
v3=fltarr(nx,ny,nz)
w3=fltarr(nx,ny,nz)
th3=fltarr(nx,ny,nz)
p3=fltarr(nx,ny,nz)
pb=fltarr(nx,ny,nz+1)
pd=fltarr(nx,ny,nz+1)
;pb=fltarr(nx,ny,nz)
;pd=fltarr(nx,ny,nz)
phyd=fltarr(nx,ny,nz+1)
buoy=fltarr(nx,ny,nz)
buoy2=fltarr(nx,ny,nz)


readi:

readu,1,time,p3,pb,pd,th3,buoy,phyd
readu,1,u3,v3,w3

p3hv=avg( avg(p3,0),0)
pbhv=avg( avg(pb,0),0)

p3=p3-p3(nx/2,ny/2,nz-1)
pb=pb-pb(nx/2,ny/2,nz-1)
pd=pd-pd(nx/2,ny/2,nz-1)

pt=pb+pd
pd2=p3-pb

th3p=th3*0

for l=0,nz-1 do begin
    th3p(*,*,L)=rh(L)*( th3(*,*,L)-th0(L) )/th0(L)
endfor
for l=1,nz-2 do begin
    buoy2(*,*,L)=(1./rh(L) )* ( th3p(*,*,L+1)-th3p(*,*,L-1) )/(zz(L+1)-zz(L-1))
endfor
;for l=nz-2,1,-1 do begin
;    phyd(*,*,L-1)= phyd(*,*,L+1) - grav * (th3p(*,*,L)/rh(L)) * (zz(L+1)-zz(L-1))
;endfor

if time gt -4.5 then STOP

goto,readi

end
