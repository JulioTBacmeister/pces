#define SMOOTHDYNSRC
#define EULAGDYN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program wojtek_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real, allocatable :: xx(:),yy(:),zz(:),rho(:),th0(:),th_e(:),A(:),B(:),C(:),rhzn(:)
real :: time,C1,C2,DX,DZ,DY
integer :: nx,ny,nz,nn,itm,iocon,ipu

integer :: xperod,yperod,zperod,ierror

real, allocatable :: u3(:,:,:),v3(:,:,:),w3(:,:,:),th3(:,:,:),p3(:,:,:),buoysc(:,:,:),dynsc(:,:,:),pdyn(:,:,:),pbyc(:,:,:),phyd(:,:,:)

!   real :: A(0:NZ),B(0:NZ),C(0:NZ),C1,C2   


    OPEN (unit = 17, file= "wojtek.dat" ,form="UNFORMATTED" )

    read(17) time,nx,ny,nz
         print*, time,nx,ny,nz

    allocate( xx(nx) )
    allocate( yy(ny) )
    allocate( zz(nz) )
    allocate( rho(nz) )
    allocate( rhzn(0:nz) )
    allocate( TH0(nz) )
    allocate( TH_e(nz) )
 
    allocate( p3(nx,ny,nz) )
    allocate( pdyn(nx,ny,0:nz) )
    allocate( pbyc(nx,ny,0:nz) )
    allocate( phyd(nx,ny,0:nz) )
 !   allocate( pdyn(nx,ny,1:nz) )
 !   allocate( pbyc(nx,ny,1:nz) )



 !   allocate( A(1:NZ),B(1:NZ),C(1:NZ) )
    allocate( A(0:NZ),B(0:NZ),C(0:NZ) )


    allocate( u3(nx,ny,nz) )
    allocate( v3(nx,ny,nz) )
    allocate( w3(nx,ny,nz) )
    allocate( th3(nx,ny,nz) )
    allocate( dynsc(nx,ny,nz) )
    allocate( buoysc(nx,ny,nz) )

    read(17) xx,yy,zz,rho,th0
    read(17) u3,v3,w3,th_e,th3,p3

    DX = xx(2)-xx(1)
    DY = yy(2)-yy(1)
    DZ = zz(2)-zz(1)


    do k=2,nz-1 
       rhzn(k) = (1./rho(k)) * ( rho(k-1)-rho(k+1) ) / (2*dz)
    end do
    k=0
    rhzn(k)=0.
    k=1
    rhzn(k)=(1./rho(k)) * ( rho(k)-rho(k+1) ) / (1*dz)
    k=nz
    rhzn(k)=(1./rho(k)) * ( rho(k-1)-rho(k) ) / (1*dz)

    !rhzn=-rhzn

! Fishpack soln
!
!                        C1*(X(I-1,J,K) -2.*X(I,J,K) +X(I+1,J,K)) +
!                        C2*(X(I,J-1,K) -2.*X(I,J,K) +X(I,J+1,K)) +
!                        A(K)*X(I,J,K-1) +B(K)*X(I,J,K)+ C(K)*X(I,J,K+1)
!                        = F(I,J,K)
!

  C1 = 1./dx**2
  C2 = 1./dy**2

  A  = 1./dz**2  + rhzn/(2.*dz)  
  B  =-2./dz**2
  C  = 1./dz**2  - rhzn/(2.*dz)  


  xperod=1 ! zero BCs at sides
  yperod=1 !

!  xperod=0 ! Periodic BCs at sides
!  yperod=0 !

!--------------------------------------
! Second order deriv w/ ghost point at NZ+1
! AND j=0 is zero
!--------------------------------------
  A(0)=0.
  C(0) =2./(DZ**2)
!  A(1)=0.
!  C(1) =2./(DZ**2)
  C(NZ)=0.
  A(NZ)=2./(DZ**2)
  zperod=1
!---------------------------------------


! end of initial read
! loop over output times
  
   iocon=0
   !do while ( iocon==0 )
   do ipu=1,30
       read(17) time
       read(17,iostat=iocon) u3,v3,w3,th_e,th3,p3
         print*, time
         print*,"Max TH0-TH_Eq",maxval(th0-th_e)


         call rhs4p( xx,yy,zz,rho,th_e,u3,v3,w3,th3,p3,nx,ny,nz,buoysc,dynsc,phyd )

         pbyc=0.
         pdyn=0.

         pbyc(:,:,1:nz)=buoysc(:,:,:)     
         !write(*,*) " going in t pois3d from fishpack "
           call POIS3D(XPEROD, NX, C1, YPEROD, NY, C2, ZPEROD, NZ+1, A, B, C, NX, NY, pbyc , IERROR)
         write(*,*) " got out pois3d from fishpack, exit condition: ", IERROR


         pdyn(:,:,1:nz)=-dynsc(:,:,:)     
         !write(*,*) " going in t pois3d from fishpack "
           call POIS3D(XPEROD, NX, C1, YPEROD, NY, C2, ZPEROD, NZ+1, A, B, C, NX, NY, pdyn , IERROR)
         write(*,*) "      got out pois3d from fishpack, exit condition: ", IERROR

         !write(*,*) DX,DY,DZ
         !write(*,*) maxval(pbyc),maxval(p3)
         !write(*,*) minval(pbyc),minval(p3)
 
         write(18) time , p3 , pbyc, pdyn, th3, buoysc, phyd
         write(18) u3 , v3, w3

   end do


   close(unit=18)


end program wojtek_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rhs4p( xx,yy,zz,rho,theq,u3,v3,w3,th3,p3,nx,ny,nz,buoysc,dynsc,ph )

integer,  intent(in) :: nx,ny,nz
real, intent(in) :: xx(nx),yy(ny),zz(nz),rho(nz),theq(nz)

real, intent(in)  :: u3(nx,ny,nz),v3(nx,ny,nz),w3(nx,ny,nz),th3(nx,ny,nz),p3(nx,ny,nz)
real, intent(out) :: dynsc(nx,ny,nz),buoysc(nx,ny,nz),ph(nx,ny,0:nz)

   real :: e2ij(nx,ny,nz),ztmg(nx,ny,nz),pref(nz)
   real :: ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),uten(nx,ny,nz)
   real :: vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz),vten(nx,ny,nz)
   real :: wx(nx,ny,nz),wy(nx,ny,nz),wz(nx,ny,nz),wten(nx,ny,nz)
   real :: ztx3(nx,ny,nz),zty3(nx,ny,nz),ztz3(nx,ny,nz),tmp(nx,ny,nz)
   real :: th3p(nx,ny,nz) !,th3pz(nx,ny,nz)
   real ::  d2ph(nx,ny,0:nz)
   real :: dx,dy, dy2

   integer :: i,j,k,n,l,i2
 
   e2ij=0.
   ztmg=0.
   dynsc=0.
   wx=0.
    wy=0.
     wz=0.
   vx=0.
    vy=0.
     vz=0.
   ux=0.
    uy=0.
     uz=0.
 

   do k=2,nz-1
   do j=2,ny-1
   do i=2,nx-1
      ztx3(i,j,k) =   ( w3(i,j+1,k) - w3(i,j-1,k) )/(yy(j+1)-yy(j-1)) & 
                    - ( v3(i,j,k+1) - v3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      zty3(i,j,k) = - ( w3(i+1,j,k) - w3(i-1,j,k) )/(xx(i+1)-xx(i-1)) & 
                    + ( u3(i,j,k+1) - u3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      ztz3(i,j,k) =   ( v3(i+1,j,k) - v3(i-1,j,k) )/(xx(i+1)-xx(i-1)) & 
                    - ( u3(i,j+1,k) - u3(i,j-1,k) )/(yy(j+1)-yy(j-1)) 
   end do
   end do
   end do


   do k=1,nz 
   do j=1,ny 
   do i=2,nx-1 
      ux(i,j,k) =     ( u3(i+1,j,k) - u3(i-1,j,k) )/(xx(i+1)-xx(i-1))
      vx(i,j,k) =     ( v3(i+1,j,k) - v3(i-1,j,k) )/(xx(i+1)-xx(i-1))
      wx(i,j,k) =     ( w3(i+1,j,k) - w3(i-1,j,k) )/(xx(i+1)-xx(i-1))
   end do
   end do
   end do

   do k=1,nz 
   do j=1,ny 
      i=1
      ux(i,j,k) =     ( u3(i+1,j,k) - u3(nx,j,k) )/(xx(i+2)-xx(i))
      vx(i,j,k) =     ( v3(i+1,j,k) - v3(nx,j,k) )/(xx(i+2)-xx(i))
      wx(i,j,k) =     ( w3(i+1,j,k) - w3(nx,j,k) )/(xx(i+2)-xx(i))
      i=nx
      ux(i,j,k) =     ( u3(1,j,k) - u3(i-1,j,k) )/(xx(i)-xx(i-2))
      vx(i,j,k) =     ( v3(1,j,k) - v3(i-1,j,k) )/(xx(i)-xx(i-2))
      wx(i,j,k) =     ( w3(1,j,k) - w3(i-1,j,k) )/(xx(i)-xx(i-2))
   end do
   end do





   do k=1,nz
   do j=2,ny-1 
   do i=1,nx
      uy(i,j,k) =      ( u3(i,j+1,k) - u3(i,j-1,k) )/(yy(j+1)-yy(j-1))
      vy(i,j,k) =      ( v3(i,j+1,k) - v3(i,j-1,k) )/(yy(j+1)-yy(j-1))
      wy(i,j,k) =      ( w3(i,j+1,k) - w3(i,j-1,k) )/(yy(j+1)-yy(j-1))
   end do
   end do
   end do

   dx2=xx(3)-xx(1)
   dy2=yy(3)-yy(1)
   do k=1,nz
   do i=1,nx
      j=1
      uy(i,j,k) =      ( u3(i,j+1,k) - u3(i,ny,k) )/dy2
      vy(i,j,k) =      ( v3(i,j+1,k) - v3(i,ny,k) )/dy2
      wy(i,j,k) =      ( w3(i,j+1,k) - w3(i,ny,k) )/dy2
      j=ny
      uy(i,j,k) =      ( u3(i,1,k) - u3(i,j-1,k) )/dy2
      vy(i,j,k) =      ( v3(i,1,k) - v3(i,j-1,k) )/dy2
      wy(i,j,k) =      ( w3(i,1,k) - w3(i,j-1,k) )/dy2
   end do
   end do

   do k=2,nz-1 
   do j=1,ny 
   do i=1,nx 
      uz(i,j,k) =      ( u3(i,j,k+1) - u3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      vz(i,j,k) =      ( v3(i,j,k+1) - v3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      wz(i,j,k) =      ( w3(i,j,k+1) - w3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
   end do
   end do
   end do

   
   do j=1,ny 
   do i=1,nx 
      k=1
      uz(i,j,k) =      ( u3(i,j,k+1) - u3(i,j,k) )/(zz(k+1)-zz(k)) 
      vz(i,j,k) =      ( v3(i,j,k+1) - v3(i,j,k) )/(zz(k+1)-zz(k)) 
      wz(i,j,k) =      ( w3(i,j,k+1) - w3(i,j,k) )/(zz(k+1)-zz(k)) 
      k=nz
      uz(i,j,k) =      ( u3(i,j,k) - u3(i,j,k-1) )/(zz(k)-zz(k-1)) 
      vz(i,j,k) =      ( v3(i,j,k) - v3(i,j,k-1) )/(zz(k)-zz(k-1)) 
      wz(i,j,k) =      ( w3(i,j,k) - w3(i,j,k-1) )/(zz(k)-zz(k-1)) 
   end do
   end do

#ifdef EULAGDYN

      write(*,*) "  Define EuLag Dynam press "

   uten = u3*ux + v3*uy + w3*uz
   vten = u3*vx + v3*vy + w3*vz
   wten = u3*wx + v3*wy + w3*wz

   do k=2,nz-1 
   do j=2,ny-1 
   do i=2,nx-1 
      !dynsc(i,j,k) =     & 
      !     rho(k) * ( uten(i+1,j,k) - uten(i-1,j,k) )/(xx(i+1)-xx(i-1))  &
      !   + rho(k) * ( vten(i,j+1,k) - vten(i,j-1,k) )/(yy(j+1)-yy(j-1))  &
      !   +    ( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
      dynsc(i,j,k) =     & 
           ( uten(i+1,j,k) - uten(i-1,j,k) )/(xx(i+1)-xx(i-1))  &
         + ( vten(i,j+1,k) - vten(i,j-1,k) )/(yy(j+1)-yy(j-1))  &
         +    (1./rho(k))*( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
   end do
   end do
   end do

   do k=2,nz-1 
   do j=2,ny-1  
     i=1 
     dynsc(i,j,k) =     & 
           ( uten(i+1,j,k) - uten(nx,j,k) )/dx2  &
         + ( vten(i,j+1,k) - vten(i,j-1,k) )/dy2  &
         +    (1./rho(k))*( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
     i=nx
     dynsc(i,j,k) =     & 
           ( uten(1,j,k) - uten(i-1,j,k) )/dx2  &
         + ( vten(i,j+1,k) - vten(i,j-1,k) )/dy2  &
         +    (1./rho(k))*( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
   end do
   end do

   do k=2,nz-1 
   do i=2,nx-1  
     j=1 
     dynsc(i,j,k) =     & 
           ( uten(i+1,j,k) - uten(i-1,j,k) )/dx2  &
         + ( vten(i,j+1,k) - vten(i,ny,k) )/dy2  &
         +    (1./rho(k))*( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
     j=ny
     dynsc(i,j,k) =     & 
           ( uten(i+1,j,k) - uten(i-1,j,k) )/dx2  &
         + ( vten(i,1,k) - vten(i,j-1,k) )/dy2  &
         +    (1./rho(k))*( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
   end do
   end do

   do k=2,nz-1 
     j=1 
     i=1
     dynsc(i,j,k) =     & 
           ( uten(i+1,j,k) - uten(nx,j,k) )/dx2  &
         + ( vten(i,j+1,k) - vten(i,ny,k) )/dy2  &
         +    (1./rho(k))*( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
     i=nx
     dynsc(i,j,k) =     & 
           ( uten(1,j,k) - uten(i-1,j,k) )/dx2  &
         + ( vten(i,j+1,k) - vten(i,ny,k) )/dy2  &
         +    (1./rho(k))*( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
    j=ny
    i=1
     dynsc(i,j,k) =     & 
           ( uten(i+1,j,k) - uten(nx,j,k) )/dx2  &
         + ( vten(i,1,k) - vten(i,j-1,k) )/dy2  &
         +    (1./rho(k))*( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
    i=nx
     dynsc(i,j,k) =     & 
           ( uten(1,j,k) - uten(i-1,j,k) )/dx2  &
         + ( vten(i,1,k) - vten(i,j-1,k) )/dy2  &
         +    (1./rho(k))*( rho(k+1)*wten(i,j,k+1) - rho(k-1)*wten(i,j,k-1) )/(zz(k+1)-zz(k-1))  
   end do

   do j=1,ny 
   do i=1,nx  
      dynsc(i,j,1)  = dynsc(i,j,2)    
      dynsc(i,j,nz) = dynsc(i,j,nz-1)    
   end do
   end do
   

     tmp = dynsc
     do k=2,nz-1 
     do j=1,ny
     do i=1,nx 
        dynsc(i,j,k) =  0.25* ( 2.*tmp(i,j,k) + tmp(i,j,k+1)  + tmp(i,j,k-1) )
     end do
     end do
     end do

#else

   do k=2,nz-1 
   do j=2,ny-1 
   do i=2,nx-1 
      e2ij(i,j,k) =      4*(ux(i,j,k)**2)  +  2*((uy(i,j,k)+vx(i,j,k))**2)  +  2*((uz(i,j,k)+wx(i,j,k))**2)  &
                     +   4*(vy(i,j,k)**2)  +  2*((vz(i,j,k)+wy(i,j,k))**2)  &
                     +   4*(wz(i,j,k)**2) 
   end do
   end do
   end do

   e2ij=(1./4.)*e2ij


   ztmg = ztx3**2 + zty3**2 + ztz3**2

#ifdef SMOOTHDYNSRC
! Smooth 
   do n=1,5
     tmp = e2ij
     do k=2,nz-1 
     do j=2,ny-1 
     do i=2,nx-1 
        e2ij(i,j,k) =  (1./7.)* (tmp(i,j,k) + tmp(i+1,j,k)  + tmp(i-1,j,k)  + tmp(i,j+1,k)  + tmp(i,j-1,k)   + tmp(i,j,k+1)  + tmp(i,j,k-1) )
     end do
     end do
     end do
   end do

   do n=1,5
     tmp = ztmg
     do k=2,nz-1 
     do j=2,ny-1 
     do i=2,nx-1 
        ztmg(i,j,k) =  (1./7.)* (tmp(i,j,k) + tmp(i+1,j,k)  + tmp(i-1,j,k)  + tmp(i,j+1,k)  + tmp(i,j-1,k)   + tmp(i,j,k+1)  + tmp(i,j,k-1) )
     end do
     end do
     end do
   end do
#endif

     write(*,*) "  Define wierd Dynam press "

     dynsc = e2ij - 0.5 * ztmg

     do k=1,nz
        dynsc(:,:,k) = dynsc(:,:,k) ! * Rho(k)  
     end do
     do k=1,nz
        pref(k) = 100000. * Rho(k)/Rho(nz)  
     end do

#endif



   do k=1,nz
   do j=1,ny 
   do i=1,nx 
      !!th3p(i,j,k) =     9.81 * rho(k)* ( (th3(i,j,k) - theq(k))/theq(k) -  p3(i,j,k)/pref(k) )
      th3p(i,j,k) =       9.81 * rho(k)* ( (th3(i,j,k) - theq(k))/theq(k) )
   end do
   end do
   end do

   buoysc=0.
   ph=0.

!#ifndef DAVIESJONES
   do k=1,nz-1 
   do j=1,ny 
   do i=1,nx 
      !buoysc(i,j,k) =      ( th3p(i,j,k+1) - th3p(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      !buoysc(i,j,k) =      ( 1./rho(k) )*( th3p(i,j,k+1) - th3p(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      buoysc(i,j,k) =      ( 1./rho(k) )*( th3p(i,j,k+1) - th3p(i,j,k) )/(zz(k+1)-zz(k)) 
      !buoysc(i,j,k) =      ( 1./rho(k) )*( th3p(i,j,k+1) - th3p(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
   end do
   end do
   end do

     tmp = buoysc
     do k=2,nz-1 
     do j=1,ny
     do i=1,nx 
        buoysc(i,j,k) =  0.25* ( 2.*tmp(i,j,k) + tmp(i,j,k+1)  + tmp(i,j,k-1) )
     end do
     end do
     end do

!#else

  dx=xx(2)-xx(1)
  dy=xx(2)-xx(1)
  do k=nz,1,-1
  do j=1,ny
  do i=1,nx
    ph(i,j,k-1) =    ph(i,j,k  ) - ( th3p(i,j,k)/rho(k) )*( zz(k)-zz(k-1) ) 
  end do
  end do
  end do


  do k=0,nz
  do j=2,ny-1
  do i=2,nx-1
     d2ph(i,j,k) =     ( ph(i+1,j  ,k  ) + ph(i-1,j  ,k  ) - 2*ph(i,j,k)  ) / (dx**2)   &
                    + ( ph(i,  j+1,k  ) + ph(i,  j-1,k  ) - 2*ph(i,j,k)  ) / (dy**2)  
  end do
  end do
  end do
  

  d2ph(nx,:,: )     = 0. !d2ph(2,   :,:)
  d2ph(1, :,: )     = 0. !d2ph(nx-1,:,:)
  d2ph(:,ny,: )     = 0. !d2ph(:,2,   :)
  d2ph(:,1, : )     = 0. !d2ph(:,ny-1,:)

  !buoysc(:,:,1:nz)=-d2ph(:,:,1:nz)
 
!#endif

end subroutine rhs4p


