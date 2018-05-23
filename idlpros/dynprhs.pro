pro dynprhs,u=u3,v=v3,w=w3,x=xx,y=yy,z=zz

   s=size(u3)
   nx=s(1)&ny=s(2)&nz=s(3)
   ztx3=u3*0
   zty3=u3*0
   ztz3=u3*0
   ztmg=u3*0


   ux = u3*0
   uy = u3*0
   uz = u3*0

   vx = u3*0
   vy = u3*0
   vz = u3*0

   wx = u3*0
   wy = u3*0
   wz = u3*0

   e2ij = u3*0



   for k=1,nz-2 do begin
   for j=1,ny-2 do begin
   for i=1,nx-2 do begin
      ztx3(i,j,k) =   ( w3(i,j+1,k) - w3(i,j-1,k) )/(yy(j+1)-yy(j-1)) $
                    - ( v3(i,j,k+1) - v3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      zty3(i,j,k) = - ( w3(i+1,j,k) - w3(i-1,j,k) )/(xx(i+1)-xx(i-1)) $
                    + ( u3(i,j,k+1) - u3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      ztz3(i,j,k) =   ( v3(i+1,j,k) - v3(i-1,j,k) )/(xx(i+1)-xx(i-1)) $ 
                    - ( u3(i,j+1,k) - u3(i,j-1,k) )/(yy(j+1)-yy(j-1)) 
   endfor
   endfor
   endfor

   for k=1,nz-2 do begin
   for j=1,ny-2 do begin
   for i=1,nx-2 do begin
      ux(i,j,k) =     ( u3(i+1,j,k) - u3(i-1,j,k) )/(xx(i+1)-xx(i-1))
      vx(i,j,k) =     ( v3(i+1,j,k) - v3(i-1,j,k) )/(xx(i+1)-xx(i-1))
      wx(i,j,k) =     ( w3(i+1,j,k) - w3(i-1,j,k) )/(xx(i+1)-xx(i-1))
   endfor
   endfor
   endfor

   for k=1,nz-2 do begin
   for j=1,ny-2 do begin
   for i=1,nx-2 do begin
      uy(i,j,k) =      ( u3(i,j+1,k) - u3(i,j-1,k) )/(yy(j+1)-yy(j-1))
      vy(i,j,k) =      ( v3(i,j+1,k) - v3(i,j-1,k) )/(yy(j+1)-yy(j-1))
      wy(i,j,k) =      ( w3(i,j+1,k) - w3(i,j-1,k) )/(yy(j+1)-yy(j-1))
   endfor
   endfor
   endfor

   for k=1,nz-2 do begin
   for j=1,ny-2 do begin
   for i=1,nx-2 do begin
      uz(i,j,k) =      ( u3(i,j,k+1) - u3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      vz(i,j,k) =      ( v3(i,j,k+1) - v3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      wz(i,j,k) =      ( w3(i,j,k+1) - w3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
   endfor
   endfor
   endfor

   for k=1,nz-2 do begin
   for j=1,ny-2 do begin
   for i=1,nx-2 do begin
      e2ij(i,j,k) =      4*(ux(i,j,k)^2)  +  2*((uy(i,j,k)+vx(i,j,k))^2)  +  2*((uz(i,j,k)+wx(i,j,k))^2)  $
                     +   4*(vy(i,j,k)^2)  +  2*((vz(i,j,k)+wy(i,j,k))^2)  $
                     +   4*(wz(i,j,k)^2) 
   endfor
   endfor
   endfor

   e2ij=(1./4.)*e2ij


   ztmg = ztx3^2 + zty3^2+ztz3^2


   for n=0,4 do begin
   tmp=e2ij
   for k=1,nz-2 do begin
   for j=1,ny-2 do begin
   for i=1,nx-2 do begin
      e2ij(i,j,k) = (1./7.)* (tmp(i,j,k) + tmp(i+1,j,k)  + tmp(i-1,j,k)  + tmp(i,j+1,k)  + tmp(i,j-1,k)   + tmp(i,j,k+1)  + tmp(i,j,k-1) )
   endfor
   endfor
   endfor
   endfor
   
   for n=0,4 do begin
   tmp=ztmg
   for k=1,nz-2 do begin
   for j=1,ny-2 do begin
   for i=1,nx-2 do begin
      ztmg(i,j,k) = (1./7.)* (tmp(i,j,k) + tmp(i+1,j,k)  + tmp(i-1,j,k)  + tmp(i,j+1,k)  + tmp(i,j-1,k)   + tmp(i,j,k+1)  + tmp(i,j,k-1) )
   endfor
   endfor
   endfor
   endfor
   

STOP

return
end
