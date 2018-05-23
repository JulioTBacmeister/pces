#define GRIDPSOLV
#undef GRIDDIAG
#undef OLDPSOLV
#undef SLCPSOLV
!  $Id: ce1_dynamics.F90,v 1.39 2008/02/14 22:07:29 bacmj Exp $
!-----------------------------

MODULE CE1_DYNAMICS

use CE1_TYPES
use CE1_DIAGS
use CE1_UTILS
use CE1_INFORM
use PPM, only : FXPPM,FXPPM2,UPST1,UPSTAD3
use CE1_CONSTS

  use CE1_PRSOLVPH2,  only: PSOLVPH
  use CE1_PR3D

 IMPLICIT NONE
 PRIVATE

 PUBLIC DYNAMICS
 PUBLIC CFL_BUDDY
 PUBLIC FLUX1D

!Some arbitrary parameters that control numerical-type 
!aspects of plume/bubble simulation
   REAL, PARAMETER :: SMOOTHING_TIME     = 1.0e2 ! 1.0e2 ! (s)



   character*72, parameter :: I_am_module = "ce1_dynamics" 

   ! Energy Checker variables
   !------------------------------------------------------------------
   real, allocatable :: HTH(:), HTHZ(:)
   real, allocatable :: BT0(:),BT1(:),W0(:),W1(:)
   real, allocatable :: FUBDT(:),FUWDT(:)
   real :: PE0,PE1,PE11,KE0,KE1,KE2,PESRCb,PESRCx,KESRCb,PESRCt,KESRCt
   !------------------------------------------------------------------
   
   
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DYNAMICS (  PCECLST  , &
                         DT       )



   !!! INTEGER,                INTENT(IN   )                ::     NCL
   type (T_ONE_PCE)   ,    intent(inout), dimension(:)  ::     PCECLST
   REAL,                   INTENT(IN   )                ::     DT
   !!!type (T_BCKG)      ,    intent(INout)                ::     BCKG

   type (T_ONE_PCE)            ::     PCE

   INTEGER :: LP,I,L,N,NCL,NXC,NYC

   real, allocatable :: PPC(:,:),PPCZ(:,:),PPCX(:,:),PPCY(:,:),PPCR(:,:),BYCC(:,:),PHYC(:,:)

   NCL = SIZE( PCECLST%IPC )
   LP = PCECLST(1)%L

   allocate( PPC(LP+1,NCL) )
   allocate( PPCX(LP+1,NCL) )
   allocate( PPCY(LP+1,NCL) )
   allocate( PPCR(LP+1,NCL) )
   allocate( PPCZ(LP,NCL) )
   allocate( BYCC(LP,NCL) )
   allocate( PHYC(LP+1,NCL) )
   PPC=0.
   PPCX=0.
   PPCY=0.
   PPCR=0.
   PPCZ=0.
   PHYC=0.

      !! Moved here from ce1 (08/01/14)
      !!---------------------------------------------
   NXC=100 
   NYC=100

   if ((control%PSOLV == 3) .or. (control%PSOLV == 4)) then
      call GRIDOVERLAY( PCECLST, DT, PPC, PPCX, PPCY, PPCR, PPCZ, PHYC, BYCC, NCL , NXC , NYC , LP )  ! in ce1_pr3d
   endif


                 ! IS_ACTIVE_4_DYNMX is .TRUE. if PCE is (not latent) AND (status > 0)
                 ! (02/11/08)
   DO N=1,NCL

        PCE = PCECLST(N)

        if(allocated(PCEDIAGS%RHOL)) then 
          PCEDIAGS%RHOL(:,N) =  C_RHL( PCE )
        endif
        if(allocated(PCEDIAGS%MASS)) then 
          PCEDIAGS%MASS(:,N) =  C_LAYER_MASS( PCE )*100./grav
        endif
                                                       ! This function is in module
        if (IS_ACTIVE_4_DYNMX(PCECLST(N)) ) THEN       ! ce1_inform

         CALL UPDATE_ONE_CE1   (    PCE              , &
                                    PCECLST          , &
                                    PPC , PPCZ       , &
                                    PPCX, PPCY, PPCR , & 
                                    PHYC, BYCC       , &
                                    DT , NCL, LP        )
        endif

        PCECLST(N) = PCE

   end do

   deallocate( PPC )
   deallocate( PPCX )
   deallocate( PPCY )
   deallocate( PPCR )
   deallocate( PPCZ )
   deallocate( BYCC )
   deallocate( PHYC )

   end SUBROUTINE DYNAMICS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE UPDATE_ONE_CE1 (   PCE      , &
                                PCECLST  , &
                                PPC,PPCZ , &
                                PPCX,PPCY,PPCR, &
                                PHYC, BYCC    , &
                                DT,NCL,LP         )

   type (T_ONE_PCE)   ,    intent(inout)                   ::     PCE
   type (T_ONE_PCE)   ,    intent(inout)  ,  dimension(:)  ::     PCECLST
   REAL,                   INTENT(IN   )                   ::     DT
   REAL,                   INTENT(IN   )                   ::     PPC(LP+1,NCL),PPCZ(LP,NCL),BYCC(LP,NCL)
   REAL,                   INTENT(IN   )                   ::     PPCX(LP+1,NCL),PPCY(LP+1,NCL),PPCR(LP+1,NCL),PHYC(LP+1,NCL)
   integer,                INTENT(IN   )                   ::     LP,NCL


   REAL              ::  BYCZ( PCE%L+1 ),BYC( PCE%L ), PP(PCE%L+1), PZ(PCE%L+1), MASS( PCE%L ), PHYD(PCE%L+1),PZL(PCE%L),PPL(PCE%L)
   REAL              ::  RCDOT(PCE%L),subclad(PCE%L), FX(PCE%L+1),AA(PCE%L),FXTERM(PCE%L),DFTERM(PCE%L),WL(PCE%L),DZ(PCE%L)
   REAL              ::  DENS(PCE%L)

   INTEGER :: I,IPC,IPCS(NCL),IdPC,LUN_ERG
   logical :: energy_checking

        !-----------------
         energy_checking = .true.
         lun_erg=511   
         MASS       =  C_LAYER_MASS( PCE )*100./grav
         DZ         =  C_LAYER_THICKNESS(PCE)
         DENS       =  C_RHL(PCE)
        !------------

         PCE%A_NM1  = PCE%A
         IPCS       = PCECLST%IPC
         IdPC       = PCE%IPC
         IPC        = FIND_IDX_OF(IPCS, IdPC )


   PHYD = 0. ! initialize possible diag

 
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Advect most things in the vertical. Precip
           ! not here, currently in ce1_precipitation.
           ! (07/29/11)
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !++++++++++++++++++++++++++++

      CALL W_AREA_CE1 ( PCE  , &
                        DT , IPC    )

      call ENERGY_CHECKER( PCE, PCECLST, IPC, DT, LP, 1)

      CALL W_QS_CE1    ( PCE  , &
                         DT     )

      CALL W_UR_CE1    ( PCE  , &
                         DT   )

      CALL W_UV_CE1    ( PCE  , &
                         DT   )
 
      CALL W_THETA_CE1 ( PCE  ,   PCECLST, &
                         DT, IPC   )

      call ENERGY_CHECKER( PCE, PCECLST, IPC, DT, LP, 2 )

     ! Kinetic Energy budget pre-calculations
     !------------------------------------------
      FX = C_EDGE_WS( PCE ) 
   
      WL = PCE%W
      if (allocated(PCEDIAGS%KEI))    & 
              PCEDIAGS%KEI(ipc)     = SUM( MASS * PCE%A_NM1 * PCE%W * PCE%W , 1 )
      if (allocated(PCEDIAGS%DKEI1))  &
              PCEDIAGS%DKEI1(ipc)   = SUM( MASS * PCE%A_NM1 * PCE%W * PCE%W , 1 )

      CALL W_W_CE1    ( PCE  ,   &
                        DT   )

      call ENERGY_CHECKER( PCE, PCECLST, IPC, DT, LP, 3 )

      if ((IPC == 2).and.energy_checking) then
      call FLUX1D_KE_CHECK( WL, PCE%W, PCE%A_Nm1, PCE%A, DENS , FX, DZ, DT, LP, lun_erg )
      end if
      
      ! Kinetic Energy budget post-calculations
      !------------------------------
      if (allocated(PCEDIAGS%DKEI1))  PCEDIAGS%DKEI1(ipc)   =  ( PCEDIAGS%DKEI1(ipc)  & 
                             - SUM( MASS * PCE%A * PCE%W * PCE%W , 1 ) )/DT


      if(control%do_subcloud_flux) then
         CALL SUBCLOUD_FLUX ( PCE , DT )
         write(*,*) " ...... *** USED subcloud fluxes ***** "
      else
         write(*,*) " ...... skipped subcloud fluxes "
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  06/26/16:
      !  Two options for pressure soln. (1) is 1-D elliptic soln.
      !  (4) uses 3-D reconstruction of bubble to force fully 3D
      !  elliptic pressure solution. Both use Davies-Jones HYD/NOHYD
      !  decomp (I think). Good anelastic form since ...
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      if (control%PSOLV == 1) then
      !--elliptic 1D solution style
         write(*,*) " option  1  1  1  1 psolvPH (Belljar Style) "
         CALL PSOLVPH       ( PCE , PCECLST, PP, BYC, PHYD, DT, IPC ) 
         call ENERGY_CHECKER( PCE, PCECLST, IPC, DT, LP, 4 )
         call ENERGY_CHECKER( PCE, PCECLST, IPC, DT, LP, 5 )
      endif

      if (control%PSOLV == 4) then
      !--gridded 3d solution style
         PP=PPC(:,IPC)
         PZ(1:LP)=PPCZ(:,IPC)
         BYC=BYCC(:,IPC)
         PHYD=PHYC(:,IPC)
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Calculate vertical acceleration
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CALL ZMOMG3_CE1    (  PCE, PP , PZ(1:LP), BYC, DT )
      endif

           if (allocated(PCEDIAGS%PP))    PCEDIAGS%PP(:,ipc)   = PP
           if (allocated(PCEDIAGS%PHYD))  PCEDIAGS%PHYD(:,ipc) = PHYD
           !!if (allocated(PCEDIAGS%RCDOT)) PCEDIAGS%RCDOT(:,ipc) = RCDOT

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Diagnose radial motion from PCE area change
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL RMOM_CE1    (  PCE, PP , PHYD, DT )


         !!call IMPOSE_W_AT_SFC( PCE, DT )

      if (control%use_prescribed_plume ) then
         call PRESC_PLUME_CE1( PCE , DT )
      endif

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! A little more ad-hoc smoothing
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (control%adhoc_smoothing_after_dyn) then 
          CALL SMOOTH_CE1  ( PCE , DT )
      else 
          write(*,*) "                                    NO AD HOC spatial smoothing !!!!! "
      endif
         
      ! CLEAN UP
      if (control%cleanup_after_dyn) then 
          CALL CLEANUP_CE1  ( PCE , DT )
      else
          write(*,*) "                                    NO KLUDGY cleanup !!!!! "
      endif
 
      ! Time Smoothing
      if (control%time_smoothing) then 
          CALL TIME_SMOOTH  ( PCE , DT )
      else
          write(*,*) "                                    NO TIME-smoothing !!!!! "
      endif
 
      ! CFL_BUDDY
      if (control%cfl_limiter) then 
          CALL CFL_BUDDY  ( PCE , DT )
      else
          write(*,*) "                                    NO CFL limiter !!!!! "
      endif
 


  END SUBROUTINE UPDATE_ONE_CE1 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CFL_BUDDY(    PCE , DT   )
!--------------------------------------------------
! added on April 11,2012
!--------------------------------------------------
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   real, intent(in) :: DT

   real :: FACTS,WCFL(PCE%L) ,DZ(PCE%L)

   INTEGER :: L,I

   FACTS = 0.01

   L=PCE%L   
   DZ   = C_LAYER_THICKNESS(PCE)

   WCFL = DZ/DT

   where(  PCE%W >= 0.75*WCFL )
        PCE%W = 0.75*WCFL
   end where
    where( -PCE%W >= 0.75*WCFL )
        PCE%W = -0.75*WCFL
   end where
  

   PCE%W(1:1)=0.
   

end SUBROUTINE CFL_BUDDY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE TIME_SMOOTH (    PCE , DT   )
!--------------------------------------------------
! added on April 11,2012
!--------------------------------------------------
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   real, intent(in) :: DT

   real :: FACTS

   INTEGER :: L,I

   FACTS = 0.01

   L=PCE%L

   PCE%A  = PCE%A*(1.-FACTS) + FACTS*PCE%A_NM1
   

   PCE%W(1:1)=0.
   

 end SUBROUTINE TIME_SMOOTH


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CLEANUP_CE1 (    PCE , DT   )
!--------------------------------------------------
! added on April 11,2012
!--------------------------------------------------
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   real, intent(in) :: DT

   real :: GDAMP(PCE%L)

   INTEGER :: L,I

   L=PCE%L

   !TDAMP = 100000.*DT
   GDAMP = 1.e-5  ! 1.0/100000.
   where( PCE%A <= control%cleanup_maxarea_ratio *PCE%A0 )
      GDAMP = 1.e-5 + ( control%cleanup_maxarea_ratio *PCE%A0 - PCE%A ) * &
                      ( 0.2-1.e-5 ) / &
              ( control%cleanup_maxarea_ratio *PCE%A0 -  control%cleanup_minarea_ratio *PCE%A0 ) 
   endwhere

   where( PCE%A <= control%cleanup_minarea_ratio *PCE%A0 )
      GDAMP = 0.2
   endwhere

   GDAMP(1:5) = 0.2
  
   GDAMP = 1.-GDAMP

      PCE%W=PCE%W * GDAMP !  (1.0-DT/TDAMP)
      PCE%THETA=PCE%THETA* GDAMP  ! (1.0-DT/TDAMP)
      !!->  PCE%Q=  (PCE%Q-Q_enviro)* GDAMP  + Q_enviro

   

   PCE%W(1:1)=0.
   

 end SUBROUTINE CLEANUP_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IMPOSE_W_AT_SFC (    PCE , &
                                  DT   )
!--------------------------------------------------
! added on April 11,2012
!--------------------------------------------------
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT

   INTEGER :: L,I

   L=PCE%L

      

IF (IS_UPDRAFT(PCE)) THEN   

   PCE%A(L) = PCE%A0
   PCE%W(L) = 5.0

ENDIF




END SUBROUTINE IMPOSE_W_AT_SFC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PRESC_PLUME_CE1 (    PCE , &
                                  DT   )
!--------------------------------------------------
! added on April 11,2012
!--------------------------------------------------
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT

   INTEGER :: L,I

   L=PCE%L

      

IF (IS_UPDRAFT(PCE)) THEN   
   WHERE( PCE%ZGE(1:L) < control%prescribed_plume_height )
        pce%w = control%prescribed_plume_w
   elsewhere
        pce%w=0.
   endwhere

   if ( control%prescribed_plume_baser > 0. ) then 
      ! write(*,*) " ADJUSTING BASE RADIUS "
      !WHERE( PCE%ZGE(1:L) < 1000. )
      !  pce%a = PI * control%prescribed_plume_baser**2
      !endwhere
      !pce%a(L-1:L) = PI * control%prescribed_plume_baser**2
   endif
    write(*,*) " USING A *Prescribed* PLUME !!!!!"
ELSE
    write(*,*) " Prescribed dynamics encountered non-updraft "
      pce%w = 0.
    write(*,*) "      ... set W=0 "

ENDIF




  END SUBROUTINE PRESC_PLUME_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE W_UR_CE1 (    PCE , &
                           DT  )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT

   real :: DENS(PCE%L) , DZ(PCE%L), X(PCE%L),RM1(PCE%L) ,WF(PCE%L+1)
   real :: WFX(PCE%L+1), Y(PCE%L),WL(PCE%L),ZL(PCE%L),DUMMY,alsm
   INTEGER :: L,I,LP

   LP=PCE%L

   DENS = C_RHL(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)
   DUMMY=100.
   WL = C_W(PCE)
   ZL = C_ZGL(PCE)

   WFX = C_EDGE_WS( PCE ) 

   WF(1:LP) = -1.* WFX(1:LP) * ( DT / DZ )
   WF(LP+1) =  WFX(LP+1)   

     ! Advect UR vertically
     !---------------------
     X=PCE%UR
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%UR=X



  END SUBROUTINE W_UR_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE W_W_CE1 (    PCE , &
                           DT   )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT
   !REAL,                  INTENT(IN)       ::     SUBCLAD(PCE%L)

   real :: DENS(PCE%L) , DZ(PCE%L), X(PCE%L),RM1(PCE%L), WL(PCE%L), ZL(PCE%L)
   real :: WFX(PCE%L+1), Y(PCE%L),DUMMY,alsm
   INTEGER :: L,I,LP,IUNC

   LP=PCE%L

   IUNC =1
   DENS = C_RHL(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)
   DUMMY=0. ! 100.

   WL = C_W(PCE)
   ZL = C_ZGL(PCE)

   WFX = C_EDGE_WS( PCE ) 

     ! Advect vertically
     X=PCE%W
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%W=X
     !PCE%W(LP)=0.
          write(*,*) 'option - WADV Simple Upwind'



if ( control%w_nu_del4 >0.0 ) then 
   alsm = control%w_nu_del4
   Y = PCE%W ! * PCE%A
   X = Y
   do L = 3, LP-2 
      X(L)   =  (1.0-6*alsm)*Y(L)   & 
               +  2*alsm*( Y(L-1) +  Y(L+1) ) &
               +  1*alsm*( Y(L-2) +  Y(L+2) ) 
   end do
      L=LP-1
      X(L)   =  (1.0-6*alsm)*Y(L)   & 
               +  2*alsm*( Y(L-1) +  Y(L+1) ) &
               +  1*alsm*( Y(L-2) +  Y(L+1) ) 
      L=LP
      X(L)   =  (1.0-6*alsm)*Y(L)   & 
               +  2*alsm*( Y(L-1) +  Y(L) ) &
               +  1*alsm*( Y(L-2) +  Y(L) ) 

   PCE%W = X ! /PCE%A
else
   write(*,*)" Skip del4 smoothing of W W  W W W !!! !! ! "
endif


  END SUBROUTINE W_W_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE W_QS_CE1 (    PCE , &
                           DT  )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT

   real :: DENS(PCE%L) , DZ(PCE%L), X(PCE%L),RM1(PCE%L) ,WF(PCE%L+1)
   real :: WFX(PCE%L+1), Y(PCE%L),DUMMY,DENSQ0,DENSQ1
   INTEGER :: L,I,LP

   L=PCE%L
   LP=PCE%L

   DENS = C_RHL(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)
   DUMMY=100.

   WFX = C_EDGE_WS( PCE ) 

   WF(1:L) = -1.* WFX(1:L) * ( DT / DZ )
   WF(L+1) =  WFX(L+1)   

     X=PCE%Q
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%Q=X

     X=PCE%QL
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%QL=X

     X=PCE%QI
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%QI=X

     X=PCE%QR
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%QR=X

     X=PCE%QS
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%QS=X

     X=PCE%QH
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%QH=X

    if (control%use_morrison_uphys) then 
      ! Move number concentrations if using
      ! Morrison microphysics. These are already
      ! denisties so DENS factor is not needed
         DENS(1:LP) = 1.0  !***
     X=PCE%NL
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%NL=X

     X=PCE%NI
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%NI=X

     X=PCE%NR
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%NR=X

     X=PCE%NS
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%NS=X

     X=PCE%NH
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%NH=X
   endif

  END SUBROUTINE W_QS_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE W_AREA_CE1 (  PCE , &
                           DT  , IPC )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT
   integer,               intent(in)       ::     IPC

   real :: DENS(PCE%L) , DZ(PCE%L), X(PCE%L),RM1(PCE%L)
   real :: WFX(PCE%L+1), Y(PCE%L),alsm,XX(PCE%L),dum0(PCE%L),dum1(PCE%L)
   INTEGER :: L,I,LP   !,IPC


   LP=PCE%L
   !IPC=PCE%IPC

   DENS = C_RHL(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)

   !if(allocated(PCEDIAGS%RHOL)) PCEDIAGS%RHOL(:,IPC)=DENS

   WFX = C_EDGE_WS( PCE ) 

   dum1(:)=1.0
   dum0(:)=1.0

     X=PCE%A
       call FLUX1D( X, dum0, dum1, DENS, WFX, DZ, DT, LP )
     PCE%A=X

  END SUBROUTINE W_AREA_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE W_THETA_CE1 (  PCE   , PCECLST, &
                            DT    , IPC  )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   type (T_ONE_PCE)   ,    intent(inout)  ,  dimension(:)  ::     PCECLST
   REAL,                  INTENT(IN)       ::     DT
   integer,               INTENT(IN)       ::     IPC
   !!REAL,                  INTENT(IN)       ::     SUBCLAD(PCE%L)

   real :: DENS(PCE%L) , DZ(PCE%L), X(PCE%L),RM1(PCE%L)
   real :: WFX(PCE%L+1), THD(PCE%L), ZGL(PCE%L), WL(PCE%L),Y(PCE%L),DTx,THDs(PCE%L)
   INTEGER :: LP,I,Nsubdiv,l,IUNC

   real :: DUMMY,alsm

   DUMMY=0. !100.
   Nsubdiv=5
   IUNC=1

   LP   = PCE%L

   DENS = C_RHL(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)

   WFX = C_EDGE_WS( PCE ) 
   ZGL = C_ZGL( PCE )
   WL  = C_W( PCE )

select case(control%theta_advection)

case(-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1st order upwind flux form applied to perturbation THETA field.
! Then, background THETA is advected using centered 2nd order scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------   
! Advection of deviation THETA using
! 1st order flux scheme
!-------------------------------------   
     X=PCE%THETA / pce%thbck
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     THD=X * pce%thbck
          write(*,*) 'option -1 THADV Simple UPWIND '
     PCE%THETA = THD
          
         call ENERGY_CHECKER( PCE, PCECLST, IPC, DT, LP, 11 )

          
!-----------------------------------------
! Add advection of background THETA by W
! 2nd order centered diffs          
!-----------------------------------------          
#if 1
   X(2:LP-1) = WL(2:LP-1)*(PCE%THBCK(1:LP-2) - PCE%THBCK(3:LP) ) &
                       /(        ZGL(1:LP-2) -       ZGL(3:LP) )
   X(1)  = WL(1) * (PCE%THBCK(1) - PCE%THBCK(2) ) &
                 /(  ZGL(1) -       ZGL(2) )
   X(LP) = WL(LP) * (PCE%THBCK(LP-1) - PCE%THBCK(LP) ) &
                 /(  ZGL(LP-1) -       ZGL(LP) )
   !X(1)  = X(2)
   !X(LP) = X(LP-1) 
#else
! Idealized term 
! Works only when thbck =  TH00 * EXP(z/H_th)
! and H_th=100000.
!---------------------------   
   X = WL*PCE%THBCK/100000.
#endif
   
!----------------------------
!  Add tendency from BCKG adv
!----------------------------
   THD   = THD - X*DT   
   PCE%THETA = THD 

case(0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1st order upwind flux form applied to total THETA field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     X=PCE%THETA + PCE%THBCK 
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%THETA = X - PCE%THBCK 
          write(*,*) 'option 0 THADV Simple UPWIND '
          
end select


if (control%theta_nu_del4>0.0) then
   alsm = control%theta_nu_del4
   Y  = PCE%THETA +  PCE%THBCK ! * PCE%A
   X  = Y
   do l = 3, lp-2 
      X(l)   =  (1.0-6*alsm)*Y(L)   & 
               +  2*alsm*( Y(L-1) +  Y(L+1) ) &
               +  1*alsm*( Y(L-2) + Y(L+2) ) 
   end do

      L=LP-1
      X(L)   =  (1.0-6*alsm)*Y(L)   & 
               +  2*alsm*( Y(L-1) +  Y(L+1) ) &
               +  1*alsm*( Y(L-2) +  Y(L+1) ) 
      L=LP
      X(L)   =  (1.0-6*alsm)*Y(L)   & 
               +  2*alsm*( Y(L-1) +  Y(L) ) &
               +  1*alsm*( Y(L-2) +  Y(L) ) 

   PCE%THETA = X -  PCE%THBCK  ! / PCE%A
else
   write(*,*) " ******  *** * SKIPPED del4 Theta smooth "
endif

END SUBROUTINE W_THETA_CE1 





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE W_UV_CE1 (    PCE , &
                           DT  )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT

   real :: DENS(PCE%L) , DZ(PCE%L), X(PCE%L),X0(PCE%L),RM1(PCE%L) ,WF(PCE%L+1),WL(PCE%L),ZGL(PCE%L)
   real :: WFX(PCE%L+1), Y(PCE%L),DUMMY,DENSQ0,DENSQ1,ALSM,DELX(PCE%L),DELXS(PCE%L),ALSZ(PCE%L)
   INTEGER :: LP,L,I

   LP=PCE%L

   DENS = C_RHL(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)

   WFX = C_EDGE_WS( PCE ) 


     X=PCE%U
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%U=X
     X=PCE%V
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%V=X

     !---------------------------------------------------------
     ! Horizontal position of clouds is a conserved property 
     ! under vertical advection. Think about stair step cloud
     ! shifted vertically (07/26/16)
     !---------------------------------------------------------
     X=PCE%X
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%X=X
     X=PCE%Y
       call FLUX1D( X, PCE%A_NM1, PCE%A , DENS, WFX, DZ, DT, LP )
     PCE%Y=X


  END SUBROUTINE W_UV_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SMOOTH_CE1 (  PCE, DT )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT

   REAL    :: AXT,A1,A0,B1,B0,ALPHA_SCALING
   INTEGER :: LP,I
   real :: MASS(PCE%L) , DZ(PCE%L), X(PCE%L),RM1(PCE%L) ,WF(PCE%L+1)
   real :: WFX(PCE%L+1), THD(PCE%L), ZGL(PCE%L), WL(PCE%L),Y(PCE%L)
   real :: AR(PCE%L)

#ifdef NOSMOOTH
    write(*,*) " ------------------   **  ** no no No NO nO SmOootHHH HthhinGGG  "
    RETURN
#endif



   LP   = PCE%L
   MASS = C_LAYER_MASS(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)

   ALPHA_SCALING = (150./MAXVAL(DZ))**2

   A0   = 0.5
   A1   = (1.-A0)*0.5
   B1   = ALPHA_SCALING*DT/SMOOTHING_TIME
   B0   = 1.0-B1


   AR   = PCE%A*MASS
   X    = PCE%THETA
   Y    = X
   Y(2:LP-1) = (A1*AR(1:LP-2)*X(1:LP-2) + A0*X(2:LP-1)*AR(2:LP-1) + A1*AR(3:LP)*X(3:LP)) &
               /( A1*AR(1:LP-2) + A0*AR(2:LP-1) + A1*AR(3:LP) )
   X    = B0 * X + B1*Y
   PCE%THETA = X

   X    = PCE%W
   Y    = X
   Y(2:LP-1) = (A1*AR(1:LP-2)*X(1:LP-2) + A0*X(2:LP-1)*AR(2:LP-1) + A1*AR(3:LP)*X(3:LP)) &
               /( A1*AR(1:LP-2) + A0*AR(2:LP-1) + A1*AR(3:LP) )
   X    = B0 * X + B1*Y
   PCE%W = X

   X    = PCE%UR
   Y    = X
   Y(2:LP-1) = (A1*AR(1:LP-2)*X(1:LP-2) + A0*X(2:LP-1)*AR(2:LP-1) + A1*AR(3:LP)*X(3:LP)) &
               /( A1*AR(1:LP-2) + A0*AR(2:LP-1) + A1*AR(3:LP) )
   X    = B0 * X + B1*Y
   PCE%UR = X


END SUBROUTINE SMOOTH_CE1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE RMOM_CE1(  PCE, PP,PHYD, DT )
   type (T_ONE_PCE)  ,       intent(inout)  ::  PCE
   REAL, DIMENSION(PCE%L+1), INTENT(in)     ::  PP,PHYD
   REAL,                     INTENT(in)     ::  DT

   REAL  :: RAD(PCE%L), RADM(PCE%L)
   REAL  :: RHL(PCE%L),ZGL(PCE%L)
   REAL  :: TDAMP_UR,RADMAX
   INTEGER :: L,I,L2,LP

      !-------------------------------------
      !-------------------------------------


   LP=PCE%L

   RHL   = C_RHL( PCE )
   ZGL   = C_ZGL(PCE)


   RAD   = SQRT( PCE%A      / PI )
   RADM  = SQRT( PCE%A_NM1  / PI )

   RADMAX = maxval( RAD )

   PCE%UR = ( RAD - RADM ) / DT

#if 1
   do L=2,LP-1
   if (PCE%W(L) >= 0) then
       PCE%UR(L) = PCE%UR(L) + PCE%W(L)*(RADM(L)-RADM(L+1))/(ZGL(L)-ZGL(L+1))
   else
       PCE%UR(L) = PCE%UR(L) + PCE%W(L)*(RADM(L-1)-RADM(L))/(ZGL(L-1)-ZGL(L))
   end if
   end do 
#endif


#if 0
   do L=1,LP
      PCE%UR(L) = PCE%UR(L)+ 1*DT*( PP(L)+PHYD(L)+PP(L+1)+PHYD(L+1) )/( RHL(L)*RADMAX )
   end do 
#endif

  ! where( PCE%A <= 0.01*PCE%A0 )
  !      PCE%UR = 0.
  ! end where   
  END SUBROUTINE RMOM_CE1 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE RCDOT_CE1(  PCE, DT, RCDOT )  

   type (T_ONE_PCE)  ,       intent(in)     ::  PCE
   REAL, DIMENSION(PCE%L),   INTENT(inout)  ::  RCDOT
   REAL,                     INTENT(in)     ::  DT

   REAL  :: DZ(PCE%L),W(PCE%L)
   REAL  :: RM(PCE%L), R0(PCE%L)

   INTEGER :: L,I,L2,LP

   LP = PCE%L

   DZ  = C_LAYER_THICKNESS(PCE)
   R0  = SQRT( PCE%A      / PI )
   RM  = SQRT( PCE%A_NM1  / PI )

   do L=2,LP-1

     RCDOT(L) = (R0(L)-RM(L))/DT + PCE%W(L)*( R0(L-1)-R0(L+1) )/ (2.*DZ(L) )

   end do
   L=1  
   RCDOT(L) = (R0(L)-RM(L))/DT + PCE%W(L)*( R0(L)-R0(L+1) )/ (1.*DZ(L) )
   L=LP     
   RCDOT(L) = (R0(L)-RM(L))/DT + PCE%W(L)*( R0(L-1)-R0(L) )/ (1.*DZ(L) )
     


end SUBROUTINE RCDOT_CE1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE ENERGY_CHECKER(  PCE, PCECLST, IPC,DT,LP, STAGE )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  This subroutine is supposed to track energy 
   !!  (kinetic and potential) through dynamics loop
   !!  Needs some documentation and notes (4/25/17)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   type (T_ONE_PCE)  ,    intent(in)       ::     PCE
   type (T_ONE_PCE)   ,   intent(inout), dimension(:)  ::     PCECLST
   REAL,                  INTENT(IN)       ::     DT
   integer,               INTENT(IN)       ::     IPC,STAGE,LP

   real :: ZL(LP),RHL(LP),DZ(LP),BETA(LP),KEi,KEf,PEi,PEf,KESRCv,PESRCv,TESRCb
   real :: QD(LP),QCD(LP)
   integer :: L,myEV,veb,clb
   !------------------------------------------------
   ! Many variables are declared in preamble of this module
   ! e.g., HTH ...
   
   

   ZL   = C_ZGL(PCE)
   RHL  = C_RHL(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)
   myEV = indx_my_enviro(PCE,PCECLST)
   veb=CONTROL % BUOY_VIRTUAL
   clb=CONTROL % BUOY_CONDLOAD

   if (IS_UPDRAFT(PCE)) then
      ! Top of dynamics
      ! before any calls to w, q, or theta evolution routines
      !--------------------------------------------------
      if (STAGE==1) then
         allocate(HTH(LP),HTHZ(LP),BT0(LP),BT1(LP),W0(LP),W1(LP))         
         allocate(FUWDT(LP),FUBDT(LP))         
         QD   =  veb*(PCE%Q- PCECLST(myev)%Q)        
         QCD  =  clb*( ( PCE%QL - PCECLST(myev)%QL )  & 
                   +   ( PCE%QI - PCECLST(myev)%QI )  &      
                   +   ( PCE%QR - PCECLST(myev)%QR )  &      
                   +   ( PCE%QS - PCECLST(myev)%QS )  &      
                   +   ( PCE%QH - PCECLST(myev)%QH )  )      
         BT0  =  PCE%THETA/PCE%THBCK   + 0.61*QD - QCD
         BETA =  PCE%THETA/PCE%THBCK   + 0.61*QD - QCD
         W0   =  PCE%W
         do L=2,LP-1
            HTH(L) = ( PCE%THBCK(L-1)-PCE%THBCK(L+1) )/( ZL(L-1)-ZL(L+1) )
         end do
         L=1
            HTH(L) = ( PCE%THBCK(L)-PCE%THBCK(L+1) )/( ZL(L)-ZL(L+1) )
         L=LP
            HTH(L) = ( PCE%THBCK(L-1)-PCE%THBCK(L) )/( ZL(L-1)-ZL(L) )
         HTH = PCE%THBCK / HTH
         do L=2,LP-1
            HTHZ(L) = ( HTH(L-1) - HTH(L+1) )/( ZL(L-1)-ZL(L+1) )
         end do
         L=1
            HTHZ(L) = ( HTH(L) - HTH(L+1) )/( ZL(L)-ZL(L+1) )
         L=LP
            HTHZ(L) = ( HTH(L-1) - HTH(L) )/( ZL(L-1)-ZL(L) )

         KE0 = 0.5 * SUM( PCE%W * PCE%W * PCE%A_NM1 *RHL * DZ , 1)
         PE0 = 0.5* grav* SUM( BETA  * BETA  * HTH * PCE%A_NM1 * RHL * DZ, 1 )
      end if

      ! After nonlinear advection of q and qc 
      ! and of theta/beta
      ! d_z(r*a*w*b)
      ! ---------------------------------------
      if (STAGE==11) then
         QD   =  veb*(PCE%Q- PCECLST(myev)%Q)        
         QCD  =  clb*( ( PCE%QL - PCECLST(myev)%QL )  & 
                   +   ( PCE%QI - PCECLST(myev)%QI )  &      
                   +   ( PCE%QR - PCECLST(myev)%QR )  &      
                   +   ( PCE%QS - PCECLST(myev)%QS )  &      
                   +   ( PCE%QH - PCECLST(myev)%QH )  )      
         BETA = PCE%THETA/PCE%THBCK  + 0.61*QD - QCD  ! "beta star"
         PE11   = 0.5* grav* SUM( BETA  * BETA  * HTH * PCE%A * RHL * DZ, 1 )
         FUBDT  = -( RHL*PCE%A*BETA - RHL*PCE%A_NM1*BT0 )
      end if

      ! After linear advection of background
      ! (1/Th)*w*d_z(Th) ~ w*H(z) 
      !----------------------------------------
      if (STAGE==2) then
         QD   =  veb*(PCE%Q- PCECLST(myev)%Q)        
         QCD  =  clb*( ( PCE%QL - PCECLST(myev)%QL )  & 
                   +   ( PCE%QI - PCECLST(myev)%QI )  &      
                   +   ( PCE%QR - PCECLST(myev)%QR )  &      
                   +   ( PCE%QS - PCECLST(myev)%QS )  &      
                   +   ( PCE%QH - PCECLST(myev)%QH )  )      
         BETA = PCE%THETA/PCE%THBCK + 0.61*QD - QCD
         BT1  = PCE%THETA/PCE%THBCK + 0.61*QD - QCD
         PE1   = 0.5 * grav* SUM( BETA  * BETA  * HTH * PCE%A * RHL * DZ, 1 )
         PESRCb = 0.5 * grav* SUM( RHL * PCE%A * W0 * (BT1+BT0)  * DZ, 1 )*DT
         PESRCx = 0.5 * grav* SUM(  PCE%W * BETA  * BETA  * HTHz * PCE%A * RHL * DZ, 1 )*DT
         PESRCt = 0.5 * grav* SUM(  BT0  * BT1  * HTH * (PCE%A-PCE%A_NM1) * RHL * DZ, 1 )
      end if

      ! After advection (nonlinear) of w
      ! d_z(r*a*w*w)
      ! ---------------------------------------
      if (STAGE==3) then
         KE1 = SUM( 0.5* PCE%W * PCE%W * PCE%A * RHL * DZ , 1)
         !PCE%W is "w star" at this stage
         FUWDT = -( RHL*PCE%A*PCE%W - RHL*PCE%A_NM1*W0 )
      end if

      ! After buoyancy force added to w
      ! g*b. Note: pce%a is a(t+1)
      !----------------------------------------
      if (STAGE==4) then
         W1  = PCE%W
         KESRCb = -0.5 * grav * SUM( RHL * PCE%A * BT1 * (W1+W0)  * DZ, 1 )*DT
         KESRCt = 0.5* SUM(  W0  * W1  * (PCE%A-PCE%A_NM1) * RHL * DZ, 1 )
         KE2 = SUM( 0.5* PCE%W * PCE%W * PCE%A * RHL * DZ , 1)
      end if
      
      if (STAGE==5) then
         KEi = 0.5 * SUM( W0 * W0 * PCE%A_NM1 *RHL * DZ , 1)
         PEi = 0.5* grav* SUM( BT0  * BT0  * HTH * PCE%A_NM1 * RHL * DZ, 1 )
         KEf = 0.5 * SUM( W1 * W1 * PCE%A *RHL * DZ , 1)
         PEf = 0.5* grav* SUM( BT1  * BT1  * HTH * PCE%A * RHL * DZ, 1 )

         TESRCb = 0.5 * grav * SUM( RHL * PCE%A * (BT1*W1-BT0*W0)  * DZ, 1 )*DT

         PESRCv = 0.5 * grav * SUM( FUBDT * HTH *( BT1+BT0) * DZ, 1)
         KESRCv = 0.5 * SUM( FUWDT *( W1+W0) * DZ, 1)
         
         write(114) pce%time,ke0,ke1,ke2,pe0,pe1,kesrcb,pesrcb,pesrcx,pe11,kesrct,pesrct,kei,pei,kef,pef,tesrcb,kesrcv,pesrcv
         deallocate(HTH,HTHZ)
         deallocate(BT0,BT1,W0,W1)         
         deallocate(FUWDT,FUBDT)         
      end if
         
end if
   
 end SUBROUTINE ENERGY_CHECKER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SUBCLOUD_FLUX(  PCE, DT )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT

   REAL    :: AXT,A1,A0,B1,B0,ALPHA_SCALING
   INTEGER :: LP,I,L
   real :: RHL(PCE%L) , DZ(PCE%L), X(PCE%L),RM1(PCE%L) ,FXT(PCE%L+1),FXW(PCE%L+1)
   real :: WFX(PCE%L+1), THD(PCE%L), ZGL(PCE%L), WL(PCE%L),Y(PCE%L)
   real :: AR(PCE%L),SC1(PCE%L),sbw(PCE%L),sbt(PCE%L),SC0


   real, parameter :: theta_bump=1.0


   LP   = PCE%L
   RHL  = C_RHL(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)
   ZGL  = C_ZGL(PCE)


   
   WFX = C_EDGE_WS( PCE ) 

   do L=2,LP-1
      SC1(L) = (PCE%A(L-1)-PCE%A(L+1) ) /(ZGL(L-1) - ZGL(L+1) )
   end do
   L=LP
    SC1(L) = (PCE%A(L-1)-PCE%A(L) ) /(ZGL(L-1) - ZGL(L) ) 
   L=1
    SC1(L) = (PCE%A(L)-PCE%A(L+1) ) /(ZGL(L) - ZGL(L+1) )

   SC0 = PCE%a0 / 2500.
   
   SC1 = control%subcloud_adv*abs(SC1) / SC0
   where( SC1 > control%subcloud_adv )
      SC1 = control%subcloud_adv 
   endwhere

      sbt = SC1
      sbw =  control%subcloud_adv*PCE%A/PCE%A0 
   where( SBW > control%subcloud_adv )
      SBW = control%subcloud_adv 
   endwhere


      FXT(1)=0.
      FXT(LP+1)=0.
      FXW(1)=0.
      FXW(LP+1)=0.
      do L=2,LP
         if ( wfx(L) < 0 ) then
            FXT(L)=  wfx(L)*sbt(L-1)*pce%A(L-1)*RHL(L-1)*pce%THETA(L-1)  !++dbg
            FXW(L)=  wfx(L)*sbw(L-1)*pce%A(L-1)*RHL(L-1)*pce%W(L-1)
         else
            FXT(L) = wfx(L)*sbt(L)*pce%A(L)*RHL(L)*pce%THETA(L)
            FXW(L) = wfx(L)*sbw(L)*pce%A(L)*RHL(L)*pce%W(L)
         endif     
      end do

      X = RHL*PCE%A*PCE%THETA
      do L=1,LP
         X(L) = X(L) + DT*(FXT(L+1)-FXT(L))/DZ(L)
      end do
      PCE%THETA = X/(RHL*PCE%A)
 
#if 1
      X = RHL*PCE%A*PCE%W
      do L=1,LP
         X(L) = X(L) + DT*(FXW(L+1)-FXW(L))/DZ(L)
      end do
      PCE%W = X/(RHL*PCE%A)
#endif


! Subcloud fluxes of water quantities
!-------------------------------------
      FXT(1)=0.
      FXT(LP+1)=0.
      do L=2,LP
         if ( wfx(L) < 0 ) then
            FXT(L)=  wfx(L)*sbt(L-1)*pce%A(L-1)*RHL(L-1)*pce%Q(L-1)  !++dbg
         else
            FXT(L) = wfx(L)*sbt(L)*pce%A(L)*RHL(L)*pce%Q(L)
         endif     
      end do

      X = RHL*PCE%A*PCE%Q
      do L=1,LP
         X(L) = X(L) + DT*(FXT(L+1)-FXT(L))/DZ(L)
      end do
      PCE%Q = X/(RHL*PCE%A)
! 
      FXT(1)=0.
      FXT(LP+1)=0.
      do L=2,LP
         if ( wfx(L) < 0 ) then
            FXT(L)=  wfx(L)*sbt(L-1)*pce%A(L-1)*RHL(L-1)*pce%QL(L-1)  !++dbg
         else
            FXT(L) = wfx(L)*sbt(L)*pce%A(L)*RHL(L)*pce%Q(L)
         endif     
      end do

      X = RHL*PCE%A*PCE%QL
      do L=1,LP
         X(L) = X(L) + DT*(FXT(L+1)-FXT(L))/DZ(L)
      end do
      PCE%QL = X/(RHL*PCE%A)
! 
      FXT(1)=0.
      FXT(LP+1)=0.
      do L=2,LP
         if ( wfx(L) < 0 ) then
            FXT(L)=  wfx(L)*sbt(L-1)*pce%A(L-1)*RHL(L-1)*pce%QI(L-1)  !++dbg
         else
            FXT(L) = wfx(L)*sbt(L)*pce%A(L)*RHL(L)*pce%Q(L)
         endif     
      end do

      X = RHL*PCE%A*PCE%QI
      do L=1,LP
         X(L) = X(L) + DT*(FXT(L+1)-FXT(L))/DZ(L)
      end do
      PCE%QI = X/(RHL*PCE%A)

     write(*,*) " SUBCLOUD advection ctl parameter is ",control%subcloud_adv
     write(*,*) " SUBCLOUD advection vector MAX ",maxval( SC1)
     !write(*,*) " SUBCLOUD advection theta bump ",theta_bump

 end SUBROUTINE SUBCLOUD_FLUX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FLUX1D( X, A0, A1, RHL, WFX, DZ, DT, LP )

   INTEGER,               INTENT(IN)       ::     LP
   REAL,                  INTENT(IN)       ::     DT, RHL(LP), A0(LP) , A1(LP) , DZ(LP) , WFX(LP+1)
   REAL,                  INTENT(INOUT)    ::     X(LP)

   REAL    :: FX(LP+1)
   INTEGER :: I,L


      FX(1)=0.
      FX(LP+1)=0.

      do L=2,LP
         if ( wfx(L) < 0 ) then
            FX(L)=  wfx(L)*A0(L-1)*RHL(L-1)*X(L-1)
         else
            FX(L) = wfx(L)*A0(L)*RHL(L)*X(L)
        endif     
      end do

      X = RHL*A0*X
      do L=1,LP
         X(L) = X(L) - DT*(FX(L)-FX(L+1))/DZ(L)
      end do
      X = X/(RHL*A1)
 

   end SUBROUTINE FLUX1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FLUX1D_KE_CHECK( WT0, WT1, A0, A1, RHL, FX, DZ, DT, LP, lun_erg )

   INTEGER,               INTENT(IN)       ::     LP , lun_erg
   REAL,                  INTENT(IN)       ::     RHL(LP), A1(LP) , DT , DZ(LP) , FX(LP+1), A0(LP)
   REAL,                  INTENT(IN)       ::     WT0(LP), WT1(LP) 
   REAL   :: TTW(LP), RHS1(LP),RHS2(LP),FAX(LP+1),A0X(LP),RHS3(LP),KRHS(LP),dke(lp)
   REAL   :: klhs1(lp), klhs2(lp), klhs3(lp),flx(lp+1),plx(lp+1)

   INTEGER :: I,K


   RHS1=0.
   RHS2=0.
   RHS3=0.
   KRHS=0.
   DKE =0.
   A0X = RHL * A0
   do K=1,LP
      TTW(k) = RHL(k) * A1(k) * (WT1(k)-WT0(k))/DT
   end do

   ! 1st order upwind fluxes of (rho a) and (rho a w)
   !-------------------------------------------------
   plx(1)=0.
   plx(lp+1)=0.
   do K=2,LP
      plx(k) =       0.5*     FX(k)  * A0x(k-1)    &
                   - 0.5*abs(FX(k))  * A0x(k-1)    &
                   + 0.5*     FX(k)  * A0x(k)      &
                   + 0.5*abs(FX(k))  * A0x(k)
   end do
   flx(1)=0.
   flx(lp+1)=0.
   do K=2,LP
      flx(k) =       0.5*     FX(k)  * A0x(k-1)*WT0(k-1)    &
                   - 0.5*abs(FX(k))  * A0x(k-1)*WT0(k-1)    &
                   + 0.5*     FX(k)  * A0x(k)*WT0(k)    &
                   + 0.5*abs(FX(k))  * A0x(k)*WT0(k) 
   end do
   !--------------------------------------------------

   
   do K=2,LP-1
      RHS1(k) =   ( &
                 (   0.5*     FX(k)  * A0x(k-1)*( WT0(k-1)-WT0(k) )    &
                   - 0.5*abs(FX(k))  * A0x(k-1)*( WT0(k-1)-WT0(k) ) ) &
              -  (   0.5*   FX(k+1)  * A0x(k+1)*( WT0(k+1)-WT0(k) )    &
                   - 0.5*abs(FX(k+1))  * A0x(k+1)*( WT0(k)-WT0(k+1) ) ) )/DZ(k)
   end do
   do K=2,LP-1
      RHS2(k) =   ( &
                    0.5*     FX(k)  * A0x(k-1)*( WT0(k-1)-WT0(k) )    &
                   - 0.5*abs(FX(k))  * A0x(k-1)*( WT0(k-1)-WT0(k) )  &
                   +   0.5*   FX(k+1)  * A0x(k+1)*( WT0(k)-WT0(k+1) )    &
                   + 0.5*abs(FX(k+1))  * A0x(k+1)*( WT0(k)-WT0(k+1) )  )/DZ(k)
   end do
   do K=2,LP-1
      RHS3(k) =   ( &
                    0.5*wt0(k-1)*( fx(k)*a0x(k-1) - abs(fx(k))*a0x(k-1) ) &
                    +  0.5*wt0(k)  *( abs(fx(k))*a0x(k-1) - fx(k)*a0x(k-1) &
                    + fx(k+1)*a0x(k+1) + abs(fx(k+1))*a0x(k+1) ) &
                    - 0.5*wt0(k+1)*( fx(k+1)*a0x(k+1) + abs(fx(k+1))*a0x(k+1) ) &
                    )/DZ(k)
   end do

      k=1
      RHS3(k) =   ( &
                     0.5*wt0(k)  *(  &
                    fx(k+1)*a0x(k+1) + abs(fx(k+1))*a0x(k+1) ) &
                    - 0.5*wt0(k+1)*( fx(k+1)*a0x(k+1) + abs(fx(k+1))*a0x(k+1) ) &
                    )/DZ(k)

      k=LP
      RHS3(k) =   ( &
                    0.5*wt0(k-1)*( fx(k)*a0x(k-1) - abs(fx(k))*a0x(k-1) ) &
                    +  0.5*wt0(k)  *( abs(fx(k))*a0x(k-1) - fx(k)*a0x(k-1) ) &
                    )/DZ(k)


   do K=1,LP
      KRHS(k) =  wt0(k) * rhs3(k)
   end do
   do K=1,LP
      KLHS1(k) =  wt0(k) * ttw(K)
   end do
   do K=1,LP
      KLHS2(k) = ( 0.5*rhl(k)*a1(k)*(wt1(k)-wt0(k))**2 )/DT
   end do
   do K=1,LP
      KLHS3(k) = ( 0.5*rhl(k)*(wt0(k)**2)*(a1(k)-a0(k)) )/DT
   end do
   do K=1,LP
      DKE(k) =  0.5*( rhl(k)*a1(k)*wt1(k)**2 - rhl(k)*a0(k)*wt0(k)**2 )/DT
   end do
   
   write(lun_erg) TTW,RHS3,KRHS,dke,klhs1,klhs2,klhs3,flx,plx,rhl,wt0,wt1,a0,a1

 end SUBROUTINE FLUX1D_KE_CHECK


end module ce1_dynamics
