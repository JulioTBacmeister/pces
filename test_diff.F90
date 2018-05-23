program test_diff

use CE1_TYPES
use CE1_UTILS
use CE1_VDIFF

implicit none

type (T_ONE_PCE), allocatable, save, dimension(:) ::  PCES

integer, parameter :: NP=1
integer, parameter :: NZ=100 ! 100
      
REAL, PARAMETER :: PI=3.1415927
REAL, PARAMETER :: DZ=300.


real  :: ZE(NZ+1),ZL(NZ),KZZ(NZ+1),RHO(NZ),DT,X(NZ),ZCEN

integer :: L,I,J,KP,N


allocate( PCES(NP) )
   KP=1
   do N=1,NP
      call CREATE_CE1 ( NZ , KP , PCES(N) )  ! in ce1_utils
   end do



DT=30.

do L=1,NZ+1
   ZE(L) = (NZ+1-L)*DZ
end do
do L=1,NZ
   ZL(L) = 0.5*(ZE(L+1)+ZE(L))
end do



#if 1
PCES(1)%ZGE=ZE
PCES(1)%PPE=1000.*EXP( -ZE/7000. )
ZCEN=5000.

do L=1,NZ
   PCES(1)%A(L) = 2000.*EXP( -( (ZL(L)-ZCEN)/2000.)**2    )
   PCES(1)%W(L) = 10.*EXP( -( (ZL(L)-ZCEN)/2000.)**2    )
   PCES(1)%Q(L) = 0.02*EXP( -( (ZL(L)-ZCEN)/300.)**2    )
end do

call VERTICAL_DIFFUSION_TEST( PCES, DT, NZ )

#endif



#if 1
KZZ = 0.00001

KZZ(40:80) = 10000.

X   = 0.
!X(50:70)=10.
X(55:60)=10.

RHO = 1.2*EXP( -ZL/70000000. )
!rho(1:40)=.00001
!rho(71:NZ) = 0.00001


write(211) NZ
write(211) RHO
write(211) X


do i=1,1000
   call vdiff(  X, KZZ, ZE, ZL, RHO , NZ, DT )
   write(211) X
end do   
#endif

end program test_diff
