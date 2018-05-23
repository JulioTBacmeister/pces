!----------------------------------------------------------------
! These are mostly highly-kludged ways of representing interactions
! of elements with each other, e.g. sucking each other in .. .
!-----------------------------------------------------------------
! Extensively modified between 09-2015 and 12-2015.
! Mutual entrainment finally implemented.
!-----------------------------------------------------------------
MODULE CE1_INTERACTIONS

use CE1_CONSTS
use CE1_TYPES
use CE1_INFORM
use CE1_UTILS
use CE1_uPHYS
use NUMERICAL_UTILITIES, only : QSAT, DQSAT
use CE1_CONSTS, only: PI
use PPM, only : FXPPM,FXPPM2,UPST1,UPSTAD3
use CE1_GCM_CPLR
use CE1_DIAGS, only : PCEDIAGS

 IMPLICIT NONE
 PRIVATE


 PUBLIC INTERACTIONS
 PUBLIC REFRESH_AENV

 character*72, parameter :: I_am_module = "ce1_interactions" 


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !  Some control/tuning parameters that probably should
 !  go into some form of init subroutine
 !-----------------------------------------------------
! SEQ_ENTRAIN
 real :: UER_FIXED ! = -0.01    ! m s-1 
 real :: W_2_UER   ! = 0.01    ! UER = -W_2_UER * W
 real :: UR_2_UER_ENT  ! = 0.01    ! UER = -W_2_UER * W
 real :: UR_2_UER_DET  ! = 0.01    ! UER = -W_2_UER * W


 logical :: USE_FIXED_UE !  = .true.
 logical :: USE_W_2_UE   ! = .false.
 logical :: USE_UR_2_UE  ! = .false.


 CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE INTERACTIONS (  PCECLST  , &
                             DT       )

   type (T_ONE_PCE) ,    intent(inout) , dimension(:)  ::     PCECLST
   REAL,                 INTENT(IN   )                 ::     DT

   type (T_ONE_PCE)                         ::     PCE 

   INTEGER :: LP,I,L,N,MYUP,NGO,NCL,NUSE

   logical :: usevector(16)

   usevector = .false.

   USE_FIXED_UE = CONTROL%USE_FIXED_UE
   USE_W_2_UE   = CONTROL%USE_W_2_UE
   W_2_UER      = CONTROL%W_2_UER
   USE_UR_2_UE  = CONTROL%USE_UR_2_UE
   UR_2_UER_ENT     = CONTROL%UR_2_UER_ENT
   UR_2_UER_DET     = CONTROL%UR_2_UER_DET
   UER_FIXED    = CONTROL%UER_FIXED  


   usevector(1)  = CONTROL%USE_FIXED_UE
   usevector(2)  = CONTROL%USE_W_2_UE
   usevector(3)  = CONTROL%USE_UR_2_UE

   nuse=count( usevector )
   if (nuse>1) THEN
     write(*,*) " too many entrainment options specified "
     STOP
   endif

    !!!pause

   LP = PCECLST(1)%L 

   NCL = SIZE( PCECLST%IPC )


                       ! Entraimnents
                       ! 
   !! Temporary Kluge - really need to rethink PCE%THBCK thing
   DO N=1,NCL
      PCECLST(N)%THETA = PCECLST(N)%THETA  +  PCECLST(N)%THBCK
   end do
   DO N=1,NCL
        PCE = PCECLST(N)
        !!if (.NOT. IS_DORMANT(PCE) ) THEN
        if (IS_ACTIVE_4_ENTR(PCE) ) THEN
          CALL CE_CE_ENTRAIN   (   PCE     , &
                                   NCL , N , & 
                                   PCECLST , &
                                   DT     ) 
        endif
        PCECLST(N) = PCE
   end do
   !! Temporary Kluge - really need to rethink PCE%THBCK thing
   !! For now PCE%THBCK remains, but other *BCKs are
   !! being absorbed by environment. - 12/14/15
   DO N=1,NCL
      if(is_environment(PCECLST(n))) PCECLST(n)%THBCK= PCECLST(n)%THETA
      PCECLST(N)%THETA = PCECLST(N)%THETA  -  PCECLST(N)%THBCK
   end do


   end SUBROUTINE INTERACTIONS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CE_CE_ENTRAIN (       PCE     , & 
                                   NCL,NPC , &
                                   PCECLST , &
                                   DT     ) 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  HISTORY.
!    09/26/13: Copied from grid_ce_entrain
!------------------------------------------------------------------

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   integer           ,    intent(in   )    ::     NCL,NPC
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCECLST(:)
   REAL , INTENT(IN)                       ::     DT
   
   REAL    :: AXT,E(PCE%L),E2(PCE%L),FIL(PCE%L),UR_MUNCH, Efx, DA, DPCEA(PCE%L)
   real    :: PCEAchk(PCE%L)
   INTEGER :: LP,I, TypeOfCE, NP, NumITouch,N,L, KP,Npo,INEST,IENP,LL,LM,LC,IE,JE,EVIN

   REAL    :: THBCK(PCE%L),AEN0(PCE%L),UR2(PCE%L),MxFx(PCE%L)

   REAL, parameter    :: ZERO_TKE_BCK = 0.0

   real, dimension( NCL ) :: XX, YY, AA 
   real, dimension( 0:NCL,PCE%L ) :: FracE,MutE
   real    :: mra
  
   LP   = PCE%L
   KP   = 1
   !!NCL  = SIZE( PCECLST%IPC )
   
   evin = INDX_MY_ENVIRO( PCE, PCECLST )



! Temporary until we figure out what to do 
! about PCES mutually entraining
!------------------------------------------
   !FracE(0  , : )  = 1.0
   !FracE(1: , : )  = 0.0
   FracE = Get_FracE ( PCE, PCECLST, NCL , LP ) 


   ! calculating a kinematic radial velocity
   ! i.e., dr_dt + w*dr_dz
   ! This is for new (12/2014) entrainment rules
   !--------------------------------------------
   PCEPHYS%UR(:,NPC) = DIAG_UR( PCE , DT  )
   PCEDIAGS%RCDOT(:,NPC) = PCEPHYS%UR(:,NPC)
 


   TypeOfCE  = PCE%TYPE_OF_CE

   SELECT CASE (TypeOfCE)

     CASE (UpdraftCode) 

       if (USE_W_2_UE) then 
          E  = C_ENTR_FLUX_W_2_UE( PCE , DT ) ! 
       endif 

       if (USE_UR_2_UE) then 
          !!E  = C_ENTR_FLUX_UR_2_UE( PCE ,  PCEPHYS%UR(:,NPC), DT ) ! 
          E  = C_ENTR_FLUX_SHAPE_PRES( PCE , DT ) ! 
       endif 

       if (USE_FIXED_UE)   E = C_ENTR_FLUX_FIXED_UE( PCE , DT ) ! 

       if (control%use_prescribed_plume ) then 
          E2 = C_ENTR_FLUX_SFC_RAD( PCE , E, DT )
          E=E2
       endif

       if (control%limit_ce_area ) then 
          E2 = C_DETR_AT_LARGE_RAD( PCE , E, DT )
          E=E2
       endif

     CASE (PrecipShaftCode)
       if (control%share_precip_v > 0 ) then
          if (control%top_down_shaft) then
             E = C_ENTR_FLUX_FIXED_UE( PCE , DT )
          else
             E = C_ENTR_FLUX_SHAFT( PCE , DT )
          endif
       else
          E = C_ENTR_SIMPLE_SHAFT( PCE , DT )
       endif
       
     CASE default
       E = 0.0

   END SELECT

       if (control%use_sponge_lyr ) then 
          E2 = C_ENTR_FLUX_SPONGE( PCE , E, DT )
          E=E2
       endif


! The following is an attempt to repsresent the detrainment that 
! occurs when PCEs "splat" as they lose their buoyancy and spread
! rapidly. This prevents PCEs from becoming ridiculously 
! large with supersonic radial winds ... (12/02/2014)
!----------------------------------------------------------------
   if(control%do_gust_munching) then
     !!WHERE( ABS(PCE%W) < 0.5 .and. PCE%UR > control%gust_munch_ur )
     !!WHERE( ABS(PCE%W) <  control%gust_munch_ur  .and. PCEPHYS%UR(:,NPC) > control%gust_munch_ur )
     WHERE( PCEPHYS%UR(:,NPC) > control%gust_munch_ur )
      E        =  E - 2.0 * PI * C_RADIUS(PCE)   & 
        * (PCEPHYS%UR(:,NPC)-control%gust_munch_ur)
     ENDWHERE
   endif

! Munch around the edges of the cloud
!---------------------------------------
   UR_MUNCH =  CONTROL%MUNCH_VEL  ! 0.1 ! 1.0  ! (m/s) 
   E        =  E - 2.0 * PI * C_RADIUS(PCE) * UR_MUNCH



  where( PCE%A <= control%shutoff_ent_area_ratio*PCE%A0 )
       E=0.
  endwhere

   PCEAchk    = PCE%A  + E*DT
  
          ! prevent detraiment from creating negative areas
          ! in detraining cloud and/or fill the cell 
          ! if it gets too thin
   where( PCEAchk <= control%minallowedarea_ratio*PCE%A0 )
        E = -( PCE%A- control%minallowedarea_ratio*PCE%A0 )/DT
   endwhere

! WRITE ENTRAINMENT DIAGNOSTIC 
   if( ALLOCATED(PCEDIAGS%E)) PCEDIAGS%E(:,NPC) = E  !PCE%UR 
  MutE = Mutual_Ent(  E, PCE, PCECLST, NCL , LP , DT ) 
  
   !! --- Id myself
          print *," used new ce ce ",pce%ipc 
    

  MUTUAL_PCEN: do n=1,NCL
     MUTUAL_PCEL : do LL=1,LP
         !Efx = E(LL)*FracE(n,LL)
         Efx = MutE(n,LL)
         IF ( Efx >= 0. ) then
            DA               =  PCE%A(LL) + Efx*DT
            PCE%THETA(LL)    = (PCE%A(LL)*PCE%THETA(LL)  +  Efx*PCECLST(n)%THETA(LL) *DT ) / DA
            PCE%U(LL)        = (PCE%A(LL)*PCE%U(LL)      +  Efx*PCECLST(n)%U(LL) *DT     ) / DA
            PCE%V(LL)        = (PCE%A(LL)*PCE%V(LL)      +  Efx*PCECLST(n)%V(LL) *DT     ) / DA
            !++jtb testing 09/17
            !!PCE%W(LL)        = (PCE%A(LL)*PCE%W(LL)      +  Efx*PCECLST(n)%W(LL) *DT     ) / DA
            PCE%Q(LL)        = (PCE%A(LL)*PCE%Q(LL)      +  Efx*PCECLST(n)%Q(LL) *DT     ) / DA
            PCE%QL(LL)       = (PCE%A(LL)*PCE%QL(LL)     +  Efx*PCECLST(n)%QL(LL) *DT    ) / DA
            PCE%QI(LL)       = (PCE%A(LL)*PCE%QI(LL)     +  Efx*PCECLST(n)%QI(LL) *DT    ) / DA
            PCE%QR(LL)       = (PCE%A(LL)*PCE%QR(LL)     +  Efx*PCECLST(n)%QR(LL) *DT    ) / DA
            PCE%QS(LL)       = (PCE%A(LL)*PCE%QS(LL)     +  Efx*PCECLST(n)%QS(LL) *DT    ) / DA
            PCE%QH(LL)       = (PCE%A(LL)*PCE%QH(LL)     +  Efx*PCECLST(n)%QH(LL) *DT    ) / DA
            PCE%A(LL)        =  PCE%A(LL) + Efx*DT
            PCECLST(n)%A(LL) =  PCECLST(n)%A(LL) - Efx*DT
         end if
         IF ( Efx < 0. ) then
            DA                   = PCECLST(n)%A(LL) - Efx*DT
            PCECLST(n)%THETA(LL) = (PCECLST(n)%THETA(LL)*PCECLST(n)%A(LL) - Efx*PCE%THETA(LL) *DT  ) / DA
            PCECLST(n)%U(LL)     = (PCECLST(n)%U(LL)*PCECLST(n)%A(LL)     - Efx*PCE%U(LL)*DT       ) / DA
            PCECLST(n)%V(LL)     = (PCECLST(n)%V(LL)*PCECLST(n)%A(LL)     - Efx*PCE%V(LL)*DT       ) / DA
            !++jtb testing 09/17
            !!PCECLST(n)%W(LL)     = (PCECLST(n)%W(LL)*PCECLST(n)%A(LL)     - Efx*PCE%W(LL)*DT       ) / DA
            PCECLST(n)%Q(LL)     = (PCECLST(n)%Q(LL)*PCECLST(n)%A(LL)     - Efx*PCE%Q(LL)*DT       ) / DA
            PCECLST(n)%QL(LL)    = (PCECLST(n)%QL(LL)*PCECLST(n)%A(LL)    - Efx*PCE%QL(LL)*DT      ) / DA
            PCECLST(n)%QI(LL)    = (PCECLST(n)%QI(LL)*PCECLST(n)%A(LL)    - Efx*PCE%QI(LL)*DT      ) / DA
            PCECLST(n)%QR(LL)    = (PCECLST(n)%QR(LL)*PCECLST(n)%A(LL)    - Efx*PCE%QR(LL)*DT      ) / DA
            PCECLST(n)%QS(LL)    = (PCECLST(n)%QS(LL)*PCECLST(n)%A(LL)    - Efx*PCE%QS(LL)*DT      ) / DA
            PCECLST(n)%QH(LL)    = (PCECLST(n)%QH(LL)*PCECLST(n)%A(LL)    - Efx*PCE%QH(LL)*DT      ) / DA
            PCECLST(n)%A(LL)     = PCECLST(n)%A(LL) - Efx*DT
            PCE%A(LL)            =  PCE%A(LL) + Efx*DT
         end if
     end do MUTUAL_PCEL
  end do MUTUAL_PCEN

       
#if 1
  MxFx = (PCE%A0/PCE%A) * (ABS(PCE%W)/5.)   * DT* 0.01*ABS(PCE%W) * 2.0 * PI * C_RADIUS(PCE) / PCE%A
  WHERE( MxFx > 0.999 )
    MxFx = 0.999
  endWHERE
  MIXITUP: do LL=LP,1,-1 
            PCE%W(LL)      = PCE%W(LL)     * (1.0-  MxFx(LL)) + MxFx(LL)*PCECLST(evin)%W(LL)   
            PCE%THETA(LL)  = PCE%THETA(LL) * (1.0-  MxFx(LL)) + MxFx(LL)*PCECLST(evin)%THBCK(LL)    
  end do MIXITUP
#endif





END SUBROUTINE CE_CE_ENTRAIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !MutE = Mutual_Ent(  E, PCE, PCECLST, NCL , LP , DT ) 
  FUNCTION Mutual_Ent ( E, PCE, PCECLST, NCL , LP , DT) RESULT(MutE)
   type (T_ONE_PCE)  ,    intent(in   )    ::     PCE
   integer           ,    intent(in   )    ::     NCL,LP
   type (T_ONE_PCE)  ,    intent(in   )    ::     PCECLST(NCL)
   real              ,    intent(in   )    ::     E(LP), DT
   real :: MutE(0:NCL, LP)

   integer n,myup,mysh,L,myev
   !---------------------------------
   ! Entrainment decision tree.
   ! 1) Is PCECLST(n) the same as PCE?
   !    Yes -> dont entrain
   !    No  ->
   !      2) Is PCE an updraft?
   !         No  -> go to (3)
   !         Yes ->
   !           2a) Is PCECLST(n) the precip shaft of PCE?
   !             Yes -> entrain some
   !             No ->
   !           2b) Is PCECLST(n) the enviro of PCE?
   !             Yes -> entrain some
   !             No  ->
   !           2c) Other reason to entrain some of PCECLST(n)?
   !             Yes -> entrain some
   !             No  ->
   !      3) Is PCE a precip shaft
   !         No  -> go to (4)
   !         Yes ->
   !           3a) Is PCECLST(n) the enviro of PCE?
   !             Yes -> entrain some
   !             No  ->
   !      4) n++ go to (1)
   !-------------------------------------

   
      
     ! Initialize mutual entrainment matrix.
     ! To start only entrain from background
     ! i.e., n=0
      MutE(0    , 1:LP )=E(1:LP)
      MutE(1:NCL, 1:LP )=0.

#if 1
     ! Me = PCE , This = PCECLST(n) 
     do n=1,NCL
       if (PCE%IPC==PCECLST(N)%IPC) then ! This is Me
          MutE(n,1:LP) = 0.              ! dont entrain myself
       else
             ! This is NOT me. So decide whether to entrain 
             ! some of this PCECLST(n).
          if(IS_UPDRAFT(PCE)) then   ! I am an UPDRAFT

            mysh = PCE%MY_PRSHFT     ! ID of my precip shaft
            myev = PCE%MY_ENV        ! ID of my environment

            ! 2a) Deal with precip shaft
            if (  (PCECLST(n)%IPC==mysh) .and. (IS_PASSIVE_4_ENTR( PCECLST(n) ) ) ) then 
                            ! This is my precip shaft, so I entrain a bit of it
                do L=1,LP
                   MutE(n,L) = 0.0*E(L) 
                end do
                            ! Ensure I dont drain the precip shaft
                do L=1,LP
                   MutE(n,L) = MIN( MutE(n,L), 0.99*PCECLST(n)%A(L)/DT )
                end do
            endif

            ! 2b) Deal with environment
            if (  (PCECLST(n)%IPC==myev) .and. (IS_PASSIVE_4_ENTR( PCECLST(n) ) ) ) then 
                            ! This is my precip shaft, so I entrain a bit of it
                do L=1,LP
                   MutE(n,L) = 1.00*E(L) 
                end do
                            ! Ensure I dont drain the environment
                do L=1,LP
                   MutE(n,L) = MIN( MutE(n,L), 0.99*PCECLST(n)%A(L)/DT )
                end do
            endif

          endif

          if(IS_SHAFT(PCE)) then  ! I am a precip shaft

            myev = PCE%MY_ENV     ! ID of my environment

            ! 3a) Deal with environment
            if (  (PCECLST(n)%IPC==myev) .and. (IS_PASSIVE_4_ENTR( PCECLST(n) ) ) ) then 
                            ! This is my environment, so I entrain a bit of it
                do L=1,LP
                   MutE(n,L) = 1.00*E(L) 
                end do
                            ! Ensure I dont drain the environment
                do L=1,LP
                   MutE(n,L) = MIN( MutE(n,L), 0.99*PCECLST(n)%A(L)/DT )
                end do
            endif

             
          endif
       endif 
     end do ! n=1,NCL loop

             ! Now adjust MutE so that 
             ! sum across n is unchanged
     do n=1,NCL
        MutE( 0 , 1:LP) =  MutE( 0 , 1:LP)- MutE( n , 1:LP)
     end do

#endif


  end FUNCTION Mutual_Ent


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Get_FracE ( PCE, PCECLST, NCL , LP ) RESULT(FracE)
   type (T_ONE_PCE)  ,    intent(in   )    ::     PCE
   integer           ,    intent(in   )    ::     NCL,LP
   type (T_ONE_PCE)  ,    intent(in   )    ::     PCECLST(NCL)
   real :: FracE(0:NCL, LP)

    
      FracE(0    , 1:LP )=1.00
      FracE(1:NCL, 1:LP )=0.

  end FUNCTION Get_FracE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 FUNCTION DIAG_UR(  PCE , DT, UR2  ) RESULT(DGUR)

   type (T_ONE_PCE)  ,       intent(inout)  ::  PCE

   real, intent(in)  :: DT   
   real  :: DGUR(PCE%L)   
   real, intent(out), optional :: UR2(PCE%L)   

   REAL  :: RAD(PCE%L), RADM(PCE%L)
   REAL  :: ZGL(PCE%L)
   INTEGER :: L,I,L2,LP

      !-------------------------------------
      !-------------------------------------


   LP=PCE%L

   ZGL   = C_ZGL(PCE)


   RAD   = SQRT( PCE%A      / PI )
   RADM  = SQRT( PCE%A_NM1  / PI )


   DGUR = ( RAD - RADM ) / DT

   if(PRESENT(UR2)) UR2 = DGUR

#if 0
   do L=2,LP-1
   if (PCE%W(L) >= 0) then
       DGUR(L) = DGUR(L) + PCE%W(L)*(RADM(L)-RADM(L+1))/(ZGL(L)-ZGL(L+1))
   else
       DGUR(L) = DGUR(L) + PCE%W(L)*(RADM(L-1)-RADM(L))/(ZGL(L-1)-ZGL(L))
   end if
   end do 
#endif


END FUNCTION DIAG_UR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION CE_SURROUNDED_AT_L ( X , Y, A, XC, YC, AC, NCL, BY_CE ) RESULT(IS_SURROUNDED)
   integer, intent(in) ::  NCL
   real,    intent(in) ::  X , Y, A, XC(NCL), YC(NCL), AC(NCL)
   integer, intent(out) :: BY_CE
   logical             ::  IS_SURROUNDED

      IS_SURROUNDED = .false.
      BY_CE         = -99

  end FUNCTION CE_SURROUNDED_AT_L 




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE REFRESH_AENV (  PCECLST  , NCL, LP ) 
   type (T_ONE_PCE)   ,    intent(inout)    ::     PCECLST(NCL)
   integer            ,    intent(in   )    ::     LP,NCL
   INTEGER :: I,J,L,N,LM,LL,LI,LJ,evin
   REAL :: APCE(LP),A0,A1,DLA

! this may not work with shared enviro's (12/21/15)

!! Loop over clouds 
APCE=0.
do N=1,NCL
 do L=1,LP
    APCE(L)  = APCE(L) +PCECLST(N)%A(L)
 end do
end do


do N=1,NCL
  !!--bck PCECLST(evin)%A = PCECLST(N)%AGRID - APCE
#if 0
   do L=1,LP
      A0 =  PCECLST(evin)%A(L)
      A1 =  PCECLST(N)%AGRID - APCE(L)
      PCECLST(evin)%A(L) = PCECLST(evin)%AGRID - APCE
   end do
#endif
end do 

end SUBROUTINE REFRESH_AENV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_ENTR_FLUX_FIXED_UE( PCE , DT ) RESULT(F)
!==================================================
!
! initial             08/15/11
!==================================================


  type (T_ONE_PCE), intent(INout) :: PCE
  real,             intent(in)    :: DT
  REAL, DIMENSION(PCE%L)       :: F

  REAL  :: ZGL(PCE%L) , TAU(PCE%L), UER(PCE%L), RAD(PCE%L), RWK(PCE%L), MINSHEAR , CON1, CON2

  INTEGER :: LP

      LP=PCE%L

         !--------------------------------------------------------------------
      UER     = UER_FIXED ! 0.1e-5 ! m s-1 

      F = -2.0 * PI * C_RADIUS(PCE) * UER

      IF ( PCE%STATUS < 0. ) F = -2.0 * PI * C_RADIUS(PCE) * 1.0

  end function C_ENTR_FLUX_FIXED_UE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_ENTR_FLUX_W_2_UE( PCE , DT ) RESULT(F)
!==================================================

!==================================================


  type (T_ONE_PCE), intent(INout) :: PCE
  real,             intent(in)    :: DT
  REAL, DIMENSION(PCE%L)       :: F

  REAL  :: UER(PCE%L)

      UER     = - W_2_UER * ABS(PCE%W)

      F = -2.0 * PI * C_RADIUS(PCE) * UER

      IF ( PCE%STATUS < 0. ) F = -2.0 * PI * C_RADIUS(PCE) * 1.0

  end function C_ENTR_FLUX_W_2_UE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_ENTR_FLUX_UR_2_UE( PCE , UR, DT ) RESULT(F)
!==================================================

!==================================================


  type (T_ONE_PCE), intent(INout) :: PCE
  real,             intent(in)    :: DT
  real, DIMENSION(PCE%L), intent(in)    :: UR
  REAL, DIMENSION(PCE%L)       :: F

  REAL  :: UER(PCE%L),DGUR(PCE%L)
  INTEGER :: LP,L

             ! UER<0 means entrainment
             ! UER>0 means detrainment
      where( UR<=0 )
      UER     = UR_2_UER_ENT * UR
      elsewhere
      UER     = UR_2_UER_DET * UR
      end where

      F = -2.0 * PI * C_RADIUS(PCE) * UER


      IF ( PCE%STATUS < 0. ) F = -2.0 * PI * C_RADIUS(PCE) * 1.0

   end function C_ENTR_FLUX_UR_2_UE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_ENTR_FLUX_SHAPE_PRES( PCE , DT ) RESULT(F)
!==================================================

!==================================================


  type (T_ONE_PCE), intent(INout) :: PCE
  real,             intent(in)    :: DT
  !!real, DIMENSION(PCE%L), intent(in)    :: UR
  REAL, DIMENSION(PCE%L)       :: F

  REAL  :: FX(PCE%L),DENS(PCE%L),ZL(PCE%L)
  INTEGER :: LP,L

      DENS = C_RHL(PCE)
      !DZ   = C_LAYER_THICKNESS(PCE)
      ZL   = C_ZGL(PCE)

      LP=PCE%L

      do L=2,LP-1
         FX(L) = PCE%A(L) * (DENS(L-1)*PCE%W(L-1)-DENS(L+1)*PCE%W(L+1) ) / ( ZL(L-1)-ZL(L+1) )
      end do
      L=1
         FX(L) = PCE%A(L) * (DENS(L)*PCE%W(L)-DENS(L+1)*PCE%W(L+1) ) / ( ZL(L)-ZL(L+1) )
      L=LP
         FX(L) = PCE%A(L) * (DENS(L-1)*PCE%W(L-1)-DENS(L)*PCE%W(L) ) / ( ZL(L-1)-ZL(L) )
     
      F = FX/DENS


      where( F>0.)
         !F=0.1*F
         F=UR_2_UER_ENT*F
      elsewhere
         F=UR_2_UER_DET*F
      end where



      !IF ( PCE%STATUS < 0. ) F = -2.0 * PI * C_RADIUS(PCE) * 1.0

   end function C_ENTR_FLUX_SHAPE_PRES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_ENTR_FLUX_UR( PCE ) RESULT(F)
  type (T_ONE_PCE), intent(INout) :: PCE
  REAL, DIMENSION(PCE%L)       :: F

  REAL  :: ZGL(PCE%L) , UER(PCE%L), A1

  INTEGER :: LP

      LP=PCE%L
      ZGL  = C_ZGL( PCE )

       
      F = -2.0 * PI * C_RADIUS(PCE) * PCE%UR


      IF ( PCE%STATUS < 0. ) F = -2.0 * PI * C_RADIUS(PCE) * 1.0

  end function C_ENTR_FLUX_UR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_ENTR_FLUX_SHAFT( PCE, DT ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, INTENT(IN)             :: DT
  REAL, DIMENSION(PCE%L)       :: F
  REAL, DIMENSION(PCE%L)       :: G


  INTEGER :: LP
      
       select case (control%shaft_entrainment)
       case (1)  
          F =  ( PCE%ASH -  PCE%A ) /DT  ! maintain shaft area at PCE%ASH

       case (2)  
          write(*,*)" no longer avail."
          STOP                           ! Used top-down overlaps to define shaft

       end select       

       ! F =  MAX( F , 0.0 )  ! Allow spreading but not contraction

       where( F < 0.0 )       ! Slow down contraction relative to spreading
          !!F = 1.0 * F 
          F = 0.1 * F 
       endwhere


  end function C_ENTR_FLUX_SHAFT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_ENTR_FLUX_SFC_RAD( PCE, E, DT ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, INTENT(IN)             :: DT
  REAL, INTENT(INOUT),DIMENSION(PCE%L)       :: E
  REAL, DIMENSION(PCE%L)       :: F
  REAL                         :: ASFC


  INTEGER :: LP

       LP = PCE%L      
       F =  E
       if ( control%prescribed_plume_baser > 0. ) then 
           write(*,*) " ADJUSTING BASE RADIUS thru ENTRAINMENT"
          ASFC = PI * control%prescribed_plume_baser**2
          F(LP-1:LP) = ( ASFC -  PCE%A(LP-1:LP) ) /DT  ! maintain shaft area near sfc
       endif
   
  end function C_ENTR_FLUX_SFC_RAD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_ENTR_GEOMETRIC( PCE, DT ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, INTENT(IN)             :: DT
  REAL, DIMENSION(PCE%L)       :: F
  REAL                         :: ASFC


  INTEGER :: LP

       LP = PCE%L      
       F =  0.
       where( PCE%A < PCE%A0  )
          F    = ( PCE%A0 -  PCE%A ) /DT  ! 
       endwhere
       where( PCE%A >= PCE%A0  )
          F    = ( PCE%A0 -  PCE%A ) /DT  ! 
       endwhere
   
     end function C_ENTR_GEOMETRIC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_ENTR_FLUX_SPONGE( PCE, E, DT ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, INTENT(IN)             :: DT
  REAL, INTENT(INOUT),DIMENSION(PCE%L)       :: E
  REAL, DIMENSION(PCE%L)       :: F,ZGL
  REAL                         :: ASFC


  INTEGER :: LP

       LP = PCE%L      
       F =  E
       ZGL = C_ZGL(PCE)
       ASFC = 0.001*PCE%A0
       where( ZGL>control%sponge_lyr_base )
          F    = ( ASFC -  PCE%A ) /DT  ! 
       endwhere
   
    end function C_ENTR_FLUX_SPONGE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_DETR_AT_LARGE_RAD( PCE, E, DT ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, INTENT(IN)             :: DT
  REAL, INTENT(INOUT),DIMENSION(PCE%L)       :: E
  REAL, DIMENSION(PCE%L)       :: F
  REAL                         :: AMAX,maxarf


  INTEGER :: LP

       LP = PCE%L      
       F =  E

       maxarf=control%max_ce_areal_fraction

       AMAX = maxarf * PCE%AGRID

       if ( maxval( PCE%A ) > maxarf * PCE%AGRID ) then

             print *," about to clip "

       endif
       
       where ( PCE%A > maxarf * PCE%AGRID )
          F    = ( AMAX -  PCE%A ) /DT  ! maintain shaft area near sfc
       endwhere


  end function C_DETR_AT_LARGE_RAD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_ENTR_SIMPLE_SHAFT( PCE , DT ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, INTENT(IN)             :: DT
  REAL, DIMENSION(PCE%L)       :: F


      F = -( PCE%A - control%SIMPLE_SHARE_AREA_FCTR * PCE%A0 )/DT

      !------------------------------------------------
      ! 07/21/16
      ! What if we allow shaft to grow but not shrink?
      where (F<0.)
        F = 0.
      endwhere
      !------------------------------------------------

  end FUNCTION C_ENTR_SIMPLE_SHAFT

!======================================
END MODULE CE1_INTERACTIONS
