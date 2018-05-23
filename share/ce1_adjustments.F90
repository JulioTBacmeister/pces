!  $Id: ce1_adjustments.F90,v 1.25 2008/04/14 15:38:09 bacmj Exp $
!===================================================
! Mostly routines associated w/ interactions of pces
! with their environment 
!===================================================

!===================================================
! Log started 8/4/11 at NCAR
! 
!    08/04/11 - Changing BCKG structure.
!---------------------------------------------------

#define SKIPSUBSI

MODULE CE1_ADJUSTMENTS

use CE1_TYPES
use CE1_DIAGS
use CE1_UTILS
use CE1_INFORM
use CE1_UPHYS
use PPM, only : FXPPM,FXPPM2,UPST1,UPSTAD3
use CE1_CONSTS
use CE1_DYNAMICS, only : FLUX1D

 IMPLICIT NONE
 PRIVATE

 PUBLIC SUBSIDENCE_XYMOTION

 character*72, parameter :: I_am_module = "ce1_adjustments" 


CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  SUBROUTINE SUBSIDENCE_XYMOTION (  PCES    , &
                                    DT        )


   type (T_ONE_PCE),   intent(inout), dimension(:)    ::     PCES
   REAL,               INTENT(IN   )                  ::     DT

   !!type (T_ONE_PCE), dimension(:),allocatable :: PCEij
   type (T_ONE_PCE)                           :: PCE

   INTEGER :: LP,I,L,N,K,idx,IM,JM,J,NoInBox,NP,IJM,LM,NN,evin,nenvs

   LOGICAL :: DISABLE
   
   REAL, dimension( PCES(1)%L )   :: ZOPCE 
   REAL, dimension( PCES(1)%L )   :: AEN0,WSUB
   REAL, dimension( PCES(1)%L )   :: WE,RHO,DZ,X,DUM0,DUM1
   REAL, dimension( PCES(1)%L+1 ) :: WEE,FX

   character*72, parameter :: I_am = "adjust_bckg_in_"//trim(I_am_module)

   NP = SIZE( PCES%IPC )
   LP = PCES(1)%L

   !!evin = indx_my_enviro( PCES(1) , PCES ) 
   nenvs = 0
   do n=1,np
      if (IS_ENVIRONMENT(PCES(n))) then
         evin  = n
         nenvs = nenvs + 1
      endif
   end do
   
   
   if (nenvs /= 1 ) then
      write(*,*) " STOP : dont know what to do with more or less than 1 enviro "
      STOP
   end if


   ! Initialize WSUB diagnostic
   WSUB = 0.
   
   ! Now calculate effects of environmental subsidence induced by each
   ! PCE in the "cluster" PCES
   do n=1,np

      if ( (CAN_INDUCE_SUBSI(PCES(N)) )  &
           .and.(CONTROL% DO_ENVIRO_SUBSI) ) THEN
         
         RHO  = C_RHL(PCES(N))
         DZ   = C_LAYER_THICKNESS(PCES(N))
 
#if 1  
         ! Calculate enviro subsidence induced by PCES(n).
         ! Currently (12/22/15) assume full-compensation.
         !!WE       = - PCES(N)%A * PCES(N)%W(1:LP) / PCES(EVIN)%A
         WE       =    PCES(N)%W(1:LP)
         WEE(2:LP) =  (WE(2:LP)+WE(1:LP-1))*0.5
         WEE(LP+1) =  0. !WEE(LP)
         WEE(1)    =  0.

         do L=2,LP
            if(WEE(L)>=0.) then
               FX(L)=WEE(L)*PCES(n)%A_NM1(L)*RHO(L)
            else
               FX(L)=WEE(L)*PCES(n)%A_NM1(L-1)*RHO(L-1)
            end if
         end do
         FX(1)=0.
         FX(LP+1)=0.

         do L=2,LP
            if(FX(L)>=0.) then
               WEE(L) = -FX(L)/(RHO(L-1)*PCES(evin)%A(L-1))
            else
               WEE(L) = -FX(L)/(RHO(L)*PCES(evin)%A(L))
            end if
         end do
#else
         WE        = - PCES(N)%A * PCES(N)%W / PCES(EVIN)%A
         WEE(2:LP) =  (WE(2:LP)+WE(1:LP-1))*0.5
         WEE(LP+1) =  0. !WEE(LP)
         WEE(1)    =  0.
#endif

         
         ! store off previous time step enviro area 
         AEN0     =  PCES(EVIN)%A


         dum1(:)=1.0
         dum0(:)=1.0

         X=PCES(evin)%A
           call FLUX1D( X, dum0, dum1, RHO , WEE, DZ, DT, LP )
         PCES(evin)%A=X

         X=PCES(evin)%THBCK
           call FLUX1D( X, AEN0, PCES(evin)%A , RHO , WEE, DZ, DT, LP )
         PCES(evin)%THBCK=X

         X=PCES(evin)%U
           call FLUX1D( X, AEN0, PCES(evin)%A , RHO , WEE, DZ, DT, LP )
         PCES(evin)%U=X

         X=PCES(evin)%V
           call FLUX1D( X, AEN0, PCES(evin)%A , RHO , WEE, DZ, DT, LP )
         PCES(evin)%V=X

         X=PCES(evin)%Q
           call FLUX1D( X, AEN0, PCES(evin)%A , RHO , WEE, DZ, DT, LP )
         PCES(evin)%Q=X

         X=PCES(evin)%QL
           call FLUX1D( X, AEN0, PCES(evin)%A , RHO , WEE, DZ, DT, LP )
         PCES(evin)%QL=X

         X=PCES(evin)%QI
           call FLUX1D( X, AEN0, PCES(evin)%A , RHO , WEE, DZ, DT, LP )
         PCES(evin)%QI=X

         X=PCES(evin)%QR
           call FLUX1D( X, AEN0, PCES(evin)%A , RHO , WEE, DZ, DT, LP )
         PCES(evin)%QR=X

         X=PCES(evin)%QS
           call FLUX1D( X, AEN0, PCES(evin)%A , RHO , WEE, DZ, DT, LP )
         PCES(evin)%QS=X

         X=PCES(evin)%QH
           call FLUX1D( X, AEN0, PCES(evin)%A , RHO , WEE, DZ, DT, LP )
         PCES(evin)%QH=X

         WSUB = WSUB + 0.5*( WEE(1:LP) + WEE(2:LP+1) )
         
      endif ! not dormant
   END DO

   PCES(evin)%W = WSUB

   do n=1,np
        PCE=PCES(N)
        !!if ( (.NOT. IS_DORMANT(PCE) ).and.(IS_UPDRAFT(PCE) ) ) THEN
        if (.NOT. IS_DORMANT(PCE) ) THEN
            CALL XYMOM_CE1    (  PCE, PCES(evin) , DT ) 
        endif
        PCES(N)=PCE
   end do


      do n=1,np
         if(.not.IS_DORMANT(pces(n)) ) then
            PCES(n)%THBCK = PCES(evin)%THBCK
         end if
      end do

   
   !!CALL CONSERVATION( PCES , DT, NP, LP )


!!   deallocate( actg )

END SUBROUTINE SUBSIDENCE_XYMOTION



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE XYMOM_CE1(  PCE, PCENV, DT )
   type (T_ONE_PCE)  ,       intent(inout)  ::  PCE
   type (T_ONE_PCE)  ,       intent(inout)  ::  PCENV
   REAL,                     INTENT(in)     ::  DT

   REAL  :: RHL(PCE%L) , RAD(PCE%L), PDYN(PCE%L), VMAG(PCE%L)
   REAL  :: DU(PCE%L)  , DV(PCE%L), PX(PCE%L)  , PY(PCE%L)
   REAL  :: DZ(PCE%L)  , MASS(PCE%L)  
   REAL  :: TIME_IN_LAYER(PCE%L) 
   REAL  :: YMEAN, XMEAN, WGT, CDM
 
   INTEGER :: L,I,L2,LP

   CDM  = control%CD_HORZDRAG

   LP=PCE%L
   DZ   = C_LAYER_THICKNESS(PCE)
   MASS = C_LAYER_MASS(PCE)

   RHL =  C_RHL( PCE) ! 1.0 * C_PPL(PCE)/1000.
 
   RAD = C_RADIUS( PCE )

   DU  = PCE%U - PCEnv%U 
   DV  = PCE%V - PCEnv%V 

   VMAG = ( DU*DU + DV*DV ) 
   PDYN =  RHL* VMAG /2.0 
   VMAG = SQRT( VMAG )

   XMEAN = 0.0
   YMEAN = 0.0
   WGT   = 0.0

   DO L=1,LP
   !!if (PCE%A(L) > ModestFraction*PCE%A0 ) then
   if (PCE%A(L) > SmallFraction*PCE%A0 ) then
      WGT    = WGT   + MASS(L) * PCE%A(L)
      XMEAN  = XMEAN + MASS(L) * PCE%A(L) * PCE%X(L)
      YMEAN  = YMEAN + MASS(L) * PCE%A(L) * PCE%Y(L)
   end if
   END DO
   if ( WGT > 0 ) then
      XMEAN = XMEAN/WGT
      YMEAN = YMEAN/WGT
   end if

   WHERE( VMAG < 0.01 )
     PDYN = 0.0
     VMAG = 100000.
   ENDWHERE

   !-------------------------------
   ! Horizontal pressure force
   !-------------------------------
   !     dynamic press    unit
   !      * length        vector 
   !      * coef          in (U,V)
   PX = CDM*RAD*PDYN    * DU/VMAG
   PY = CDM*RAD*PDYN    * DV/VMAG

   !  Save U,V before drag
   DU    = PCE%U
   DV    = PCE%V

   WHERE( PCE%A > SmallFraction*PCE%A0 )
     PCE%U = PCE%U - PX*DT /(RHL*PCE%A )
     PCE%V = PCE%V - PY*DT /(RHL*PCE%A )
   END WHERE   

   !  Get velocity change by drag
   DU    = PCE%U - DU
   DV    = PCE%V - DV

   ! need to accelerate/decelerate environment
   !--------------------------------------------
   
   PCenv%U =  PCenv%U  -  PCE%A*DU / PCenv%A
   PCenv%V =  PCenv%V  -  PCE%A*DV / PCenv%A

   !!WHERE( PCE%A > ModestFraction*PCE%A0 )
   WHERE( PCE%A > SmallFraction*PCE%A0 )
     PCE%X = PCE%X + PCE%U*DT
     PCE%Y = PCE%Y + PCE%V*DT
   END WHERE   

   PCE%XMEAN = XMEAN
   PCE%YMEAN = YMEAN

   
   !!WHERE( PCE%A <= ModestFraction*PCE%A0 )
   WHERE( PCE%A <= SmallFraction*PCE%A0 )
     PCE%X = XMEAN
     PCE%Y = YMEAN
   end where
 
  END SUBROUTINE XYMOM_CE1 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  ! after elminating bck need to rethink this
  ! procedure (12/25/15)
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CONSERVATION (  PCECLST  , DT, NCL, LP ) 
!!-----------------------------------------
!! Repairs conservation errors arising from
!! inconsistent numerics
!!------------------------------------------
   type (T_ONE_PCE)   ,    intent(inout)    ::     PCECLST(NCL)
   real               ,    intent(in   )    ::     DT
   integer            ,    intent(in   )    ::     LP,NCL
   INTEGER :: I,J,L,N,NN,LM,LL,LI,LJ,NOOG
   REAL :: APCE(LP),ABK(LP),DAB(LP),AGRID,RHL(LP),DZ(LP)
   REAL :: QTBKEN, QTBKDE
   INTEGER :: OOG(NCL)

   ! Stuff to ensure conservation ....
   
end SUBROUTINE CONSERVATION
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE CE1_ADJUSTMENTS





  
