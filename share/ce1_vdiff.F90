!  $Id: ce1_vdiff.F90 created 2013/10/31 NCAR $
!-----------------------------
MODULE CE1_VDIFF


use CE1_TYPES
use CE1_UTILS
use CE1_INFORM
use CE1_CONSTS
!------------------------------------
! !DESCRIPTION:

 IMPLICIT NONE
 PRIVATE

 PUBLIC VDIFF
 PUBLIC VERTICAL_DIFFUSION
 PUBLIC VERTICAL_DIFFUSION_TEST

   character*72, parameter :: I_am_module = "ce1_prsolv02" 



contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE VERTICAL_DIFFUSION (  PCECLST  , &
                                   DT ,       &
                                   LP    )

   type (T_ONE_PCE) ,    intent(inout) , dimension(:)  ::     PCECLST
   !!!type (T_BCKG)    ,    intent(inout)                 ::     BCKG
   REAL,                 INTENT(IN   )                 ::     DT
   integer,              INTENT(IN   )                 ::     LP
   type (T_ONE_PCE)    ::     PCE

   REAL, DIMENSION(LP)    :: RHL,ZL, AR,X,RHA,AM
   REAL, DIMENSION(LP+1)  :: ZE, KZZ
  
   REAL   :: wdamp

   INTEGER :: I,L,N, MYUP, IMYUP,NCL 
 
   NCL = SIZE( PCECLST%IPC )

   
   do I=1,NCL
      PCE=PCECLST(I)

      if (IS_ACTIVE_4_DYNMX(PCE) ) THEN       !

#if 1
         ZL  = C_ZGL( PCE )
         ZE  = X_ZGE( PCE )

         RHL = C_RHL( PCE )

#if 0
         KZZ =  C_KZZ( PCE )   !CONTROL%KZZ_BCKG ! C_KZZ( PCE )
         X   = PCE%A
         call VDIFF(  X, KZZ, ZE, ZL, RHL , LP, DT )
         PCE%A  = X
#endif

         RHA = MAX( C_RHL( PCE )*PCE%A, 0.0001 )
         !!RHA = C_RHL( PCE )

         KZZ = C_KZZ( PCE )
         X   = PCE%THETA + PCE%THBCK
         call VDIFF(  X, KZZ, ZE, ZL, RHA , LP, DT )
         PCE%THETA  = X -  PCE%THBCK

         X   = PCE%W
         call VDIFF(  X, KZZ, ZE, ZL, RHA , LP, DT )
         PCE%W  = X
#if 1
         X   = PCE%Q
         call VDIFF(  X, KZZ, ZE, ZL, RHA , LP, DT )
         PCE%Q  = X
         X   = PCE%QL
         call VDIFF(  X, KZZ, ZE, ZL, RHA , LP, DT )
         PCE%QL  = X
         X   = PCE%QI
         call VDIFF(  X, KZZ, ZE, ZL, RHA , LP, DT )
         PCE%QI  = X
#endif

         if (control%WDAMP_TAU > 0.0 ) then
         do L=1,LP
            wdamp = DT /  control % WDAMP_TAU
            wdamp = wdamp*(PCE%A0/PCE%A(L) )
            wdamp = MIN( wdamp , 0.1 )
            !!pce%w(L)     = (1.-wdamp)*pce%w(L)
            pce%theta(L) = (1.-wdamp)*pce%theta(L)
         end do
         endif
#endif


      endif
      
      PCECLST(I) = PCE

   end do

 end SUBROUTINE VERTICAL_DIFFUSION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE VERTICAL_DIFFUSION_TEST(  PCECLST  , &
                                   DT ,       &
                                   LP    )

   type (T_ONE_PCE) ,    intent(inout) , dimension(:)  ::     PCECLST
   !!!type (T_BCKG)    ,    intent(inout)                 ::     BCKG
   REAL,                 INTENT(IN   )                 ::     DT
   integer,              INTENT(IN   )                 ::     LP
   type (T_ONE_PCE)    ::     PCE

   REAL, DIMENSION(LP)    :: RHL,ZL, AR,X
   REAL, DIMENSION(LP+1)  :: ZE, KZZ
  

   INTEGER :: I,L,N, MYUP, IMYUP,NCL 
 

   NCL = SIZE( PCECLST%IPC )

   
   do I=1,NCL
      PCE=PCECLST(I)

      ZL  = C_ZGL( PCE )
      ZE  = X_ZGE( PCE )

      RHL = MAX( C_RHL( PCE )*PCE%A, 0.0001 )

      KZZ = C_KZZ( PCE )

      X   = PCE%Q

      call VDIFF(  X, KZZ, ZE, ZL, RHL , LP, DT )

   end do

 end SUBROUTINE VERTICAL_DIFFUSION_TEST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE VDIFF( X, KZZ, ZE, ZL, RHO , LP, DT )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sets up and solves tridiagonal matrix representation of
!
!         d_t( X ) = d_z( rho * K d_z(X) )   
!
! Problem is restated as 
!
!    a(i-1)*X(i-1) + b(i)*X(i) + c(i)*X(i+1) = X_t(i) for 1<i<LP
!
! where a,b, and c are sub-diagonal, diagonal, and super-diagonal
! of a tridiagonal matrix


   INTEGER,               INTENT(in   )   :: LP
   REAL,                  INTENT(in   )   :: DT
   REAL, DIMENSION(LP)  , INTENT(inout)   :: X
   REAL, DIMENSION(LP)  , INTENT(in   )   :: RHO, ZL
   REAL, DIMENSION(LP+1), INTENT(in   )   :: KZZ, ZE

   REAL              ::  RHE(1:LP+1)

   real :: S(0:LP+1)
   real :: A(0:LP), B(0:LP+1),C(0:LP), DD(LP+2)
   real :: DZL(1:LP),DZE(1:LP+1)
   
   integer :: L,I,J,INFO

   DO L=1,LP
      DZL(L)   = ZE(L)-ZE(L+1)
   END DO
   DO L=2,LP
      DZE(L)   = ZL(L-1)-ZL(L)
   END DO
   DZE(1)=DZE(2)
   DZE(LP+1)=DZE(LP)

   DO L=2,LP
      RHE(L) = 0.5*(RHO(L)+RHO(L-1))
   END DO
   RHE(1) = RHE(2)
   RHE(LP+1)=RHE(LP)


   DO L=1,LP
      A(L-1)  =   - KZZ(L)*RHE(L)*DT  / ( RHO(L)*DZL(L)*DZE(L) )
   END DO

   DO L=1,LP
      B(L)  =   1.0 +  KZZ(L)*RHE(L)*DT  / ( RHO(L)*DZL(L)*DZE(L) )  +  KZZ(L+1)*RHE(L+1)*DT  / ( RHO(L)*DZL(L)*DZE(L+1) )
   END DO

   DO L=1,LP
      C(L)  =     - KZZ(L+1)*RHE(L+1)*DT  / ( RHO(L)*DZL(L)*DZE(L+1) )
   END DO


! Set zero derivative BC's 
!--------------------------
   L=LP
   B(L+1) =  KZZ(L+1)*RHE(L+1)*DT  / ( RHO(L)*DZL(L)*DZE(L+1) ) 
   A(L)   = -B(L+1)

   L=1
   B(L-1) =  KZZ(L)*RHE(L)*DT  / ( RHO(L)*DZL(L)*DZE(L) )
   C(L-1) = -B(L-1)

   S(1:LP) = X(1:LP)
   S(0)    = 0.
   S(LP+1) = 0.


         call sgtsv(  LP+2 , 1, A , B , C , S , LP+2, INFO )

   X(1:LP) = S(1:LP)

 end SUBROUTINE VDIFF



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_KZZ( PCE ) RESULT(KZZ)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L+1)       :: KZZ
  REAL :: LAMBDA, AA,WW
  INTEGER :: LP,L
      LP=PCE%L
      LAMBDA = control % VDIFF_LAMBDA

      do L=2,LP

         AA     = 0.5*(PCE%A(L)+PCE%A(L-1) )
         WW     = SQRT(  0.5*( (PCE%W(L))**2 + (PCE%W(L-1))**2 ) )
         IF (AA > 0.01*PCE%A0 ) THEN
            KZZ(L)=MAX( WW*LAMBDA, control%kzz_bckg )
         ELSE
            KZZ(L)=control%kzz_bckg
         ENDIF

      end do

      KZZ(1)=KZZ(2)
      KZZ(LP+1)=KZZ(LP) 

    end function C_KZZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end MODULE CE1_VDIFF
