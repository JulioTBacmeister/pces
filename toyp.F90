program toyp

implicit none

integer, parameter :: NX=100
integer, parameter :: NZ=100

real, parameter :: dx =100.
real, parameter :: dz =100.

real, dimension(NX,0:NZ) :: p,bz,pm1,pp1,res
real, dimension(NX,1:NZ) :: b
real, dimension(NX) :: p00

real :: err,pii
integer :: i,j,n

pii = 2.*asin(1.0)

b=0.

!b( :, 1:int(0.5*NZ) ) =0.1
!b(int(0.15*NX):int(0.75*NX), 1:int(0.5*NZ) ) =0.1
b(int(0.45*NX):int(0.55*NX), 1:int(0.5*NZ) ) =0.1


do j=1,nz
   b(:,j) = b(:,j) - sum( b(:,j) , 1 )/NX
end do

p=0.



bz(:,1:NZ-1) = ( b(:,2:NZ) - b(:, 1:NZ-1) ) / dz
bz(:,0)      = 0.
bz(:,NZ)     = 0.

res=0.
p= 0.
p00(1*nx/4:3*nx/4 )=1.

do i=1,nx-1

   p00(i)=sin( 2.*pii *(i-1) / (nx-1) )

enddo

p00(nx)=p00(1)

pp1=p
do n=1,5000

  pm1=p
 
! $Id: toyp3.F90,v 1.16 2007/02/16 21:23:03 bacmj Exp $
program toyp

implicit none

integer, parameter :: NX=100 ! 100 ! 200
integer, parameter :: NY=100 ! 100 ! 200
integer, parameter :: NZ=100 ! 100
      
REAL, PARAMETER :: PI=3.1415927

real, parameter :: dx =200.
real, parameter :: dy =200.
real, parameter :: dz =200.

real, dimension(NX,NY,0:NZ) :: p,d2ph,pp1,ph,res
real, dimension(NX,NY,1:NZ) :: b

real,dimension(nx,ny) :: p00

real, dimension(NX)   :: x
real, dimension(NY)   :: y
real, dimension(NZ  ) :: z , rho, thbck
real, dimension(0:NZ) :: ze


real  ::  omega_sor ,dx1,dy1,dz1

real :: err,hrad,hscale
integer :: i,j,k,n

!----------------------------------------------------------------------------
! p      dimension(NX,NY,0:NZ)     Nonhydrostatic pressure
! ph     dimension(NX,NY,0:NZ)     Hydrostatic pressure 
!                                  (g*z-intgeral of overlying perturbation mass field)
! b      dimension(NX,NY,NZ)       buoyancy = grav * density perturbation/mean density
!                                    (m s-2)
! d2ph   dimension(NX,NY,0:NZ)     Horizontal Laplacian of hydrostatic pressure

b=0.


do i=1,nx
   x(i) = (i-nx/2)*dx
enddo
do i=1,ny
   y(i) = (i-ny/2)*dy
enddo
do i=0,nz
   ze(i) =  i * dz
enddo
z=(ze(0:nz-1)+ze(1:nz))/2.

rho = 1.2*exp( -z / 7000. )
thbck = 300.*exp( z /30000. )

!#if 0

do k=1,nz
   do j=1,ny
      do i=1,nx
         !hscale    = 4000.  !*exp( - 0.5*(z(k)/2500.)**2 )
         hscale    = 4000.  !*exp( - ((z(k)-200.)/3000.)**2 ) + 0.01
         hrad      = sqrt( ( x(i)-0.5*z(k))**2 + y(j)**2 ) / hscale
         hscale    = 4000.*exp( - ((z(k)-200.)/12000.)**2 )
         if ( (hrad < 1.0 ) .and. (hscale > 200.) ) then
            b(i,j,k)  = 0.5* (hscale/4000.) * (1.0 + cos( PI*  hrad ) ) 
         else 
            b(i,j,k)  = 0.0
         end if
      end do
   end do
end do 

where ( b < 0.)
  b = 0.
endwhere

b=5.0*b


do k=1,nz
   do j=1,ny
      do i=1,nx
        b(i,j,k) = - 9.8 * rho(k) * b(i,j,k) / thbck(k)
      end do
   end do
end do

!#endif

p00=0.
p00(1*nx/4:3*nx/4 , 1*ny/4:3*ny/4)=1.


! implementing Davies-Jones 2000 p_nh / p_h decomposition

ph(:,:,NZ) = 0.
do k=nz-1,0,-1
   ph(:,:,k) = ph(:,:,k+1) + b(:,:,k+1)*DZ
end do

! One source term for elliptic p_nh eq is HORIZONTAL del2 of p_h
!--------------------------------------------------------------- 
  do k=0,nz
  do j=2,ny-1
  do i=2,nx-1
       

    d2ph(i,j,k) =     ( ph(i+1,j  ,k  ) + ph(i-1,j  ,k  ) - 2*ph(i,j,k)  ) / (dx**2)   &
                    + ( ph(i,  j+1,k  ) + ph(i,  j-1,k  ) - 2*ph(i,j,k)  ) / (dy**2)  

  end do
  end do
  end do
  
  !d2ph(:,:,0)       = d2ph(:,:,1   ) 
  !d2ph(:,:,nz)      = d2ph(:,:,nz-1)

  d2ph(nx,:,: )     = d2ph(2,   :,:)
  d2ph(1, :,: )     = d2ph(nx-1,:,:)
  d2ph(:,ny,: )     = d2ph(:,2,   :)
  d2ph(:,1, : )     = d2ph(:,ny-1,:)




pp1= 0.0 ! d2ph / ( 2./(dx**2) + 2./(dy**2)  + 2./(dz**2)  )
pp1(:,:,0)      = p00(:,:)

dx1=1.
dy1=1.
dz1=1.

p=d2ph

do n=1,1000

  ! Periodic BCs  50x50x51 Dirichlet Top/Bottom
  omega_sor = 1.88
  omega_sor = 1.85

  do k=1,nz-1
  do j=2,ny-1
  do i=2,nx-1


!    pp1(i,j,k) =        p(i,j,k)  +    omega_sor *                               &
!                    ( ( p(i+1,j  ,k  ) + p(i-1,j  ,k  ) ) / (dx**2)               &
!                    + ( p(i,  j+1,k  ) + p(i,  j-1,k  ) ) / (dy**2)               &
!                    + ( p(i,  j  ,k+1) + p(i,  j  ,k-1) ) / (dz**2)               &
!                    -   p(i,  j,  k  )*( 2./(dx**2) + 2./(dy**2) + 2./(dz**2) )   &
!                    + d2ph(i,j,k)   ) / ( 2./(dx**2) + 2./(dy**2) + 2./(dz**2) )

!     pp1 (i,j,k) =        p(i,j,k)  +    omega_sor *                               &
!                    ( ( p(i+1,j  ,k  ) + p(i-1,j  ,k  ) )               &
!                    + ( p(i,  j+1,k  ) + p(i,  j-1,k  ) )               &
!                    + ( p(i,  j  ,k+1) + p(i,  j  ,k-1) )               &
!                    -   p(i,  j,  k  )*( 6. )                           &
!                    + d2ph(i,j,k)*(dx**2)   ) / (6. )

     res (i,j,k) =      p(i+1,j  ,k  ) + p(i-1,j  ,k  )                 &
                    +   p(i,  j+1,k  ) + p(i,  j-1,k  )                 &
                    +   p(i,  j  ,k+1) + p(i,  j  ,k-1)                 &
                    -   p(i,  j,  k  )* 6.                              &
                    +   d2ph(i,j,k)*(dx**2)  

     p (i,j,k) =  p(i,j,k) + omega_sor * res(i,j,k) / 6.0

  end do
  end do
  end do
  
#if 0 
  p(:,:,0)      = p00(:,:) ! pp1(:,:,1   ) !- b(:,:,1 )*dz
  p(:,:,nz)     = 0. ! pp1(:,:,nz-1) !+ b(:,:,nz)*dz
#endif

  p(:,:,0)      = p(:,:,1   ) !- b(:,:,1 )*dz
  p(:,:,nz)     = p(:,:,nz-1) !+ b(:,:,nz)*dz

    ! periodic BCs really slow you down

#if 1
  p(nx,:,: )     = 0. ! p(2,   :,:)
  p(1, :,: )     = 0. ! p(nx-1,:,:)
  p(:,ny,: )     = 0. ! p(:,2,   :)
  p(:,1, : )     = 0. ! p(:,ny-1,:)
#endif

#if 0
  p(nx,:,: )     = p(2,   :,:)
  p(1, :,: )     = p(nx-1,:,:)
  p(:,ny,: )     = p(:,2,   :)
  p(:,1, : )     = p(:,ny-1,:)
#endif

  !err = sum( (pp1-p)**2 ) / sum(p**2)

  err = maxval( abs(res) )/maxval( abs(p) )

  write(*,*) n, maxval( abs(p) ), err, omega_sor

  if ( (err <= 1.0e-6) .and. (n > 5 ) ) exit

end do

write(211) NX,NY,NZ
write(211) X,  Y, Z, ZE
write(211) b
write(211) p
write(211) ph
write(211) d2ph

end program toyp
  do j=1,nz-1
  do i=2,nx-1


    ! pp1(i,j)  =  ( p(i+1,j  ) + p(i-1,j  )   &
    !            + p(i,  j+1) + p(i,  j-1) ) / 4.0

     res(i,j)  =  ( p(i+1,j  ) + p(i-1,j  )   &
                + p(i,  j+1) + p(i,  j-1) - 4.0*p(i,j) ) 
                                         
     p(i,j)  = p(i,j) + 1.0*res(i,j)/4.0

  end do
  end do
  
  p(:,0)      = p00
  p(:,nz)     = 0.0 !p(:,nz-1) + b(:,nz)*dz


  p(nx,: )     = p(2,   :)
  p(1, : )     = p(nx-1,:)


  err=  maxval(abs(res)) / maxval(abs(p))


  write(*,*) n,err,maxval(abs(res)),maxval(abs(p))
  if ( (err < 1.e-6) .and. (n > 5) )  exit

end do

write(211) b
write(211) p


end program toyp
