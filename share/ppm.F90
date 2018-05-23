!  $Id: ppm.F90,v 1.2 2006/07/25 22:06:29 bacmj Exp $
!
!-----------------------------
module ppm

IMPLICIT NONE
PRIVATE

PUBLIC FXPPM
PUBLIC FXPPM2
PUBLIC UPST1
PUBLIC UPSTAD3
PUBLIC LMTPPM

CONTAINS


SUBROUTINE fxppm( U , Q , I1 , IN  & 
                , UNCONSTRAINED )

IMPLICIT NONE

INTEGER, INTENT(IN)         :: I1, IN

#include "ppms.h"

REAL, INTENT(IN)            :: U(BOUNDS)
REAL, INTENT(INOUT)         :: Q(BOUNDS)
INTEGER, INTENT(IN),OPTIONAL   :: UNCONSTRAINED


REAL              :: UU(PADDED_BOUNDS)
REAL              :: FF(PADDED_BOUNDS)
REAL              :: QQ(PADDED_BOUNDS)
REAL              :: D4(PADDED_BOUNDS)
REAL              :: DC(PADDED_BOUNDS)
REAL              :: PMAX(PADDED_BOUNDS)
REAL              :: PMIN(PADDED_BOUNDS)
REAL              :: TEST(PADDED_BOUNDS)

REAL              :: AR(PADDED_BOUNDS)
REAL              :: AL(PADDED_BOUNDS)
REAL              :: A6(PADDED_BOUNDS)
REAL              :: UPOS(BOUNDS)

INTEGER           :: I

REAL              :: R23=2.0/3.0

UU(I1:IN)=U
QQ(I1:IN)=Q

UU(I1-5:I1-1)=U(IN-4:IN)
UU(IN+1:IN+5)=U(I1:I1+4)

QQ(I1-5:I1-1)=0. !Q(IN-4:IN)
QQ(IN+1:IN+5)=0. !Q(I1:I1+4)

!! PRINT *, ' UU ========= '
!! PRINT *, UU 


!! 4th ORDER DERIVATIVE OF Q AT Q POINTS

D4(BOUNDS) = (  8.0 * (  QQ(X_1_RT_OF_X) - QQ(X_1_LF_OF_X) )    &
         -    (  QQ(X_2_RT_OF_X) - QQ(X_2_LF_OF_X) )  )/24.0


DO I=I1,IN
    PMAX(I) = AMAX1( QQ(I-1),QQ(I),QQ(I+1) ) - QQ(I)
    PMIN(I) = QQ(I) - AMIN1( QQ(I-1),QQ(I),QQ(I+1) ) 
    DC(I)   = SIGN( AMIN1( ABS(D4(I)),PMIN(I),PMAX(I) ), D4(I) )
END DO

DC(I1-5:I1-1)=DC(IN-4:IN)
DC(IN+1:IN+5)=DC(I1:I1+4)

AL(BOUNDS) = ( QQ(X_1_LF_OF_X)+QQ(BOUNDS) )*0.5 &
           + ( DC(X_1_LF_OF_X)-DC(BOUNDS) )*(1./3.)
AL(I1-5:I1-1)=AL(IN-4:IN)
AL(IN+1:IN+5)=AL(I1:I1+4)

AR(BOUNDS) = AL( X_1_RT_OF_X )
AR(I1-5:I1-1)=AR(IN-4:IN)
AR(IN+1:IN+5)=AR(I1:I1+4)

A6(BOUNDS) = 3.0*( 2*QQ(BOUNDS) - AL(BOUNDS) - AR(BOUNDS) )
A6(I1-5:I1-1)=A6(IN-4:IN)
A6(IN+1:IN+5)=A6(I1:I1+4)


if ( .NOT. PRESENT(UNCONSTRAINED) ) then
  call LMTPPM( QQ , DC , AL , AR, A6, I1, IN )
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WHERE(  U >= 0.)
!!    FF(BOUNDS)=UU(BOUNDS)*QQ(Q_1_LF_OF_F)
!!ELSEWHERE
!!    FF(BOUNDS)=UU(BOUNDS)*QQ(Q_1_RT_OF_F)
!!END WHERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


WHERE(  U >= 0.)
   UPOS = 1.0
ELSEWHERE
   UPOS = 0.0
END WHERE



    FF(BOUNDS) =  UPOS(BOUNDS) *                          &
             (   AR(Q_1_LF_OF_F) + 0.5 * UU(BOUNDS) *     &
               ( AL(Q_1_LF_OF_F) - AR(Q_1_LF_OF_F) +      &
                 (1.-R23*UU(BOUNDS) )*A6(Q_1_LF_OF_F) ) ) &

             +  (1.0-UPOS(BOUNDS))*                       &
             (   AL(Q_1_RT_OF_F) - 0.5 * UU(BOUNDS) *     &
               ( AR(Q_1_RT_OF_F) - AL(Q_1_RT_OF_F) +      &
                 (1.+R23*UU(BOUNDS) )*A6(Q_1_RT_OF_F) ) )


!WHERE(  U >= 0.)
!    FF(BOUNDS) = AR(Q_1_LF_OF_F) + 0.5 * UU(BOUNDS) * &
!               ( AL(Q_1_LF_OF_F) - AR(Q_1_LF_OF_F) +  &
!                 (1.-R23*UU(BOUNDS) )*A6(Q_1_LF_OF_F) )
!ELSEWHERE
!    FF(BOUNDS) = AL(Q_1_RT_OF_F) - 0.5 * UU(BOUNDS) * &
!               ( AR(Q_1_RT_OF_F) - AL(Q_1_RT_OF_F) +  &
!                 (1.+R23*UU(BOUNDS) )*A6(Q_1_RT_OF_F) )
!END WHERE



!! PRINT *, FF

FF(I1-5:I1-1)=FF(IN-4:IN)
FF(IN+1:IN+5)=FF(I1:I1+4)

FF=FF*UU

Q = Q - ( FF(F_1_RT_OF_Q)- FF(F_1_LF_OF_Q) )

END SUBROUTINE fxppm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE upst1( U , Q , DX, DT, I1 , IN  )

IMPLICIT NONE

INTEGER, INTENT(IN)         :: I1, IN

REAL, INTENT(IN)            :: U(I1:IN+1)
REAL, INTENT(IN)            :: DX(I1:IN)
REAL, INTENT(IN)            :: DT
REAL, INTENT(INOUT)         :: Q(I1:IN)

REAL :: FF(I1:IN+1)

integer :: I

do i=i1+1,iN
   if ( u(i) <= 0. ) then
      ff(i) = u(i)*q(i)
   else
      ff(i) = u(i)*q(i-1)
   endif
end do 

ff(i1)=0.
ff(iN+1)=0.

do i=i1,iN
   q(I) = q(i) - ( ff(i+1)-ff(i) )*dt/dx(i)
end do



END SUBROUTINE upst1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE fxppm2( U , Q , DX, DT, I1 , IN  & 
                , UNCONSTRAINED )

IMPLICIT NONE

INTEGER, INTENT(IN)         :: I1, IN

#include "ppms.h"

REAL, INTENT(IN)            :: U(BOUNDS)
REAL, INTENT(IN)            :: DX(BOUNDS)
REAL, INTENT(IN)            :: DT
REAL, INTENT(INOUT)         :: Q(BOUNDS)
INTEGER, INTENT(IN),OPTIONAL   :: UNCONSTRAINED


REAL              :: UU(PADDED_BOUNDS)
REAL              :: FF(PADDED_BOUNDS)
REAL              :: QQ(PADDED_BOUNDS)
REAL              :: D4(PADDED_BOUNDS)
REAL              :: DC(PADDED_BOUNDS)
REAL              :: DXX(PADDED_BOUNDS)
REAL              :: PMAX(PADDED_BOUNDS)
REAL              :: PMIN(PADDED_BOUNDS)
REAL              :: TEST(PADDED_BOUNDS)

REAL              :: AR(PADDED_BOUNDS)
REAL              :: AL(PADDED_BOUNDS)
REAL              :: A6(PADDED_BOUNDS)
REAL              :: UPOS(BOUNDS)

INTEGER           :: I

REAL              :: R23=2.0/3.0

UU(I1:IN)=U
QQ(I1:IN)=Q
DXX(I1:IN)=DX

UU(I1-5:I1-1)=U(IN-4:IN)
UU(IN+1:IN+5)=U(I1:I1+4)

QQ(I1-5:I1-1)=0. !Q(IN-4:IN)
QQ(IN+1:IN+5)=0. !Q(I1:I1+4)

DXX(I1-5:I1-1)=DX(IN-4:IN)
DXX(IN+1:IN+5)=DX(I1:I1+4)

!! PRINT *, ' UU ========= '
!! PRINT *, UU 


!! 4th ORDER DERIVATIVE OF Q AT Q POINTS

D4(BOUNDS) = (  8.0 * (  QQ(X_1_RT_OF_X) - QQ(X_1_LF_OF_X) )    &
         -    (  QQ(X_2_RT_OF_X) - QQ(X_2_LF_OF_X) )  )/24.0


DO I=I1,IN
    PMAX(I) = AMAX1( QQ(I-1),QQ(I),QQ(I+1) ) - QQ(I)
    PMIN(I) = QQ(I) - AMIN1( QQ(I-1),QQ(I),QQ(I+1) ) 
    DC(I)   = SIGN( AMIN1( ABS(D4(I)),PMIN(I),PMAX(I) ), D4(I) )
END DO

DC(I1-5:I1-1)=DC(IN-4:IN)
DC(IN+1:IN+5)=DC(I1:I1+4)

AL(BOUNDS) = ( QQ(X_1_LF_OF_X)+QQ(BOUNDS) )*0.5 &
           + ( DC(X_1_LF_OF_X)-DC(BOUNDS) )*(1./3.)
AL(I1-5:I1-1)=AL(IN-4:IN)
AL(IN+1:IN+5)=AL(I1:I1+4)

AR(BOUNDS) = AL( X_1_RT_OF_X )
AR(I1-5:I1-1)=AR(IN-4:IN)
AR(IN+1:IN+5)=AR(I1:I1+4)

A6(BOUNDS) = 3.0*( 2*QQ(BOUNDS) - AL(BOUNDS) - AR(BOUNDS) )
A6(I1-5:I1-1)=A6(IN-4:IN)
A6(IN+1:IN+5)=A6(I1:I1+4)


if ( .NOT. PRESENT(UNCONSTRAINED) ) then
!  call LMTPPM( QQ , DC , AL , AR, A6, I1, IN )
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WHERE(  U >= 0.)
    FF(BOUNDS)=QQ(Q_1_LF_OF_F)
ELSEWHERE
    FF(BOUNDS)=QQ(Q_1_RT_OF_F)
END WHERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


WHERE(  U >= 0.)
   UPOS = 1.0
ELSEWHERE
   UPOS = 0.0
END WHERE



!    FF(BOUNDS) =  UPOS(BOUNDS) *                          &
!             (   AR(Q_1_LF_OF_F) + 0.5 * UU(BOUNDS) *     &
!               ( AL(Q_1_LF_OF_F) - AR(Q_1_LF_OF_F) +      &
!                 (1.-R23*UU(BOUNDS) )*A6(Q_1_LF_OF_F) ) ) &
!
!             +  (1.0-UPOS(BOUNDS))*                       &
!             (   AL(Q_1_RT_OF_F) - 0.5 * UU(BOUNDS) *     &
!               ( AR(Q_1_RT_OF_F) - AL(Q_1_RT_OF_F) +      &
!                 (1.+R23*UU(BOUNDS) )*A6(Q_1_RT_OF_F) ) )


!WHERE(  U >= 0.)
!    FF(BOUNDS) = AR(Q_1_LF_OF_F) + 0.5 * UU(BOUNDS) * &
!               ( AL(Q_1_LF_OF_F) - AR(Q_1_LF_OF_F) +  &
!                 (1.-R23*UU(BOUNDS) )*A6(Q_1_LF_OF_F) )
!ELSEWHERE
!    FF(BOUNDS) = AL(Q_1_RT_OF_F) - 0.5 * UU(BOUNDS) * &
!               ( AR(Q_1_RT_OF_F) - AL(Q_1_RT_OF_F) +  &
!                 (1.+R23*UU(BOUNDS) )*A6(Q_1_RT_OF_F) )
!END WHERE

FF(I1-5:I1-1)=FF(IN-4:IN)
FF(IN+1:IN+5)=FF(I1:I1+4)

FF=FF*UU

Q = Q - ( FF(F_1_RT_OF_Q)- FF(F_1_LF_OF_Q) )*DT/DX(BOUNDS)

END SUBROUTINE fxppm2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE lmtppm( QQ , DC , AL , AR, A6, I1, IN )

IMPLICIT NONE
INTEGER, INTENT(IN)         :: I1, IN

#include "ppms.h"

REAL, INTENT(IN)             :: QQ(PADDED_BOUNDS)
REAL, INTENT(IN)             :: DC(PADDED_BOUNDS)

REAL, INTENT(INOUT)             :: AR(PADDED_BOUNDS)
REAL, INTENT(INOUT)             :: AL(PADDED_BOUNDS)
REAL, INTENT(INOUT)             :: A6(PADDED_BOUNDS)

REAL              :: DA1(PADDED_BOUNDS)
REAL              :: DA2(PADDED_BOUNDS)

INTEGER           :: I

REAL              :: R23=2.0/3.0

DA1 = AR - AL
DA2 = DA1**2

WHERE ( DC == 0.0 )
   AR = QQ
   AL = QQ
   A6 = A6*0.0
ELSEWHERE

   WHERE(  A6*DA1+DA2 < 0.0)
      A6 = 3.*(AL-QQ)
      AR = AL-A6
   END WHERE
   WHERE(  A6*DA1-DA2 > 0.0)
      A6 = 3.*(AR-QQ)
      AL = AR-A6
   END WHERE

END WHERE



END SUBROUTINE lmtppm



SUBROUTINE upstad3( U , Q , X, DT, I1 , IN  )

IMPLICIT NONE

INTEGER, INTENT(IN)         :: I1, IN

REAL, INTENT(IN)            :: U(I1:IN)
REAL, INTENT(IN)            :: X(I1:IN)
REAL, INTENT(IN)            :: DT
REAL, INTENT(INOUT)         :: Q(I1:IN)

REAL :: FF(I1-2:IN+2), UU(I1-2:IN+2), QQ(I1-2:IN+2), XX(I1-2:IN+2)

REAL :: XDP, QDP 

integer :: I

qq(I1:IN) = q(I1:IN)
xx(I1:IN) = x(I1:IN)
!uu(I1:IN) = u(I1:IN)

qq(i1-2:i1-1) = q(i1)  ! q(iN-1:iN)
xx(i1-1) = xx(i1) -   (xx(i1+1)-xx(i1))
xx(i1-2) = xx(i1) - 2*(xx(i1+1)-xx(i1))



qq(iN+1:iN+2) = q(iN)  ! q(i1:i1+1)
xx(iN+1) = xx(iN) +   (xx(iN)-xx(iN-1))
xx(iN+2) = xx(iN) + 2*(xx(iN)-xx(iN-1))

do i=i1,iN

   xdp = xx(i) - u(i)*dt

   if ( u(i) <= 0. )  then 

       qdp  =                                                                              &
               qq(i) * ( xdp - xx(i-2) )*( xdp - xx(i-1) )*( xdp - xx(i+1) )               &
                  /  ( ( xx(i) - xx(i-2) )*( xx(i) - xx(i-1) )*( xx(i) - xx(i+1) ) )           &

            +  qq(i-1) * ( xdp - xx(i-2) )*( xdp - xx(i) )*( xdp - xx(i+1) )               &
                  /  ( ( xx(i-1) - xx(i-2) )*( xx(i-1) - xx(i) )*( xx(i-1) - xx(i+1) ) )       &

            +  qq(i+1) * ( xdp - xx(i-2) )*( xdp - xx(i) )*( xdp - xx(i-1) )               &
                  /  ( ( xx(i+1) - xx(i-2) )*( xx(i+1) - xx(i) )*( xx(i+1) - xx(i-1) ) )      &

            +  qq(i-2) * ( xdp - xx(i-1) )*( xdp - xx(i) )*( xdp - xx(i+1) )               &
                  /  ( ( xx(i-2) - xx(i-1) )*( xx(i-2) - xx(i) )*( xx(i-2) - xx(i+1) ) )     


     else 

       qdp   =                                                                              &
               qq(i) * ( xdp - xx(i+2) )*( xdp - xx(i-1) )*( xdp - xx(i+1) )               &
                   / ( ( xx(i) - xx(i+2) )*( xx(i) - xx(i-1) )*( xx(i) - xx(i+1) ) )           &

               +  qq(i-1) * ( xdp - xx(i+2) )*( xdp - xx(i) )*( xdp - xx(i+1) )               &
                  / ( ( xx(i-1) - xx(i+2) )*( xx(i-1) - xx(i) )*( xx(i-1) - xx(i+1) ) )       &

               +  qq(i+1) * ( xdp - xx(i+2) )*( xdp - xx(i) )*( xdp - xx(i-1) )     &
                  / ( ( xx(i+1) - xx(i+2) )*( xx(i+1) - xx(i) )*( xx(i+1) - xx(i-1) ) )      &

               +  qq(i+2) * ( xdp - xx(i-1) )*( xdp - xx(i) )*( xdp - xx(i+1) )     &
                  / ( ( xx(i+2) - xx(i-1) )*( xx(i+2) - xx(i) )*( xx(i+2) - xx(i+1) ) )       

      endif

      q(i) = qdp     

end do 



END SUBROUTINE upstad3




end module ppm
