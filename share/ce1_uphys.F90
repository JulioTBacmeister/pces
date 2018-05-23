!  $Id: ce1_uphys.F90,v 1.16 2008/02/11 23:06:40 bacmj Exp $
!
!-----------------------------
#define SKIPREVAP
!!!!#define SKIPAUTOC

MODULE CE1_UPHYS

use CE1_TYPES
use CE1_UTILS
use CE1_INFORM
use NUMERICAL_UTILITIES, only : QSAT, DQSAT
use CE1_CONSTS
use PPM, only : FXPPM,FXPPM2,UPST1,UPSTAD3
use MODULE_MP_MORR_TWO_MOMENT
use CE1_DIAGS, only : PCEDIAGS

 IMPLICIT NONE
 PRIVATE

 PUBLIC UPHYS
 PUBLIC RAIN_TVEL
 PUBLIC RAIN_TVEL2
 !PUBLIC CONDENSE_CE1

 character*72, parameter :: I_am_module = "ce1_uphys" 

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE UPHYS (     PCECLST  , &
                         DT        )

   type (T_ONE_PCE), dimension(:),  intent(inout)    ::     PCECLST
   REAL,                            INTENT(IN   )    ::     DT

   type (T_ONE_PCE)            ::     PCE

   INTEGER :: LP,I,L,N,NCL

   NCL = SIZE( PCECLST%IPC )

                 ! IS_ACTIVE_4_UPHYS is .TRUE. if PCE is (not latent)
                 ! (02/11/08)
   DO N=1,NCL

     PCE = PCECLST(N)

     if ( IS_ACTIVE_4_UPHYS(PCE) ) THEN

     !CALL FILLQS ( PCE , N )

      if(control%use_morrison_uphys) then 

       CALL MORRISON_WRAPPER( PCE , DT, N )

      else
       CALL CONDENSE_CE1 ( PCE    , &
                           DT    , & 
                           CP,GRAV,RKAP,ALHL,ALHS     ) 

#ifdef SKIPAUTOC
             write(*,*) "SKIPPING SKIPPING SKIPPING AUTOCON "
#else
       CALL AUTOCONV_CE1 ( PCE  , &
                           DT     )
#endif

! Should put vertical motion of PRECIP here
       CALL W_QPREC_CE1 (  PCE ,  & 
                           DT )

#ifdef SKIPREVAP
             write(*,*) "SKIPPING SKIPPING SKIPPING REVAP "
#else
       CALL REVAP_CE1    ( PCE ,  &
                           DT       )
#endif

    
       CALL DRAIN_CE1    ( PCE  )

       endif !whicu uphys


     endif
   
     PCECLST(N) = PCE

   END DO


  !!CALL PRSHAFTS_NEEDED( PCECLST , NCL )

  END SUBROUTINE UPHYS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MORRISON_WRAPPER ( PCE ,  DT , IPC   )

   type (T_ONE_PCE)  , intent(inout)    ::   PCE
   !type (T_ONE_PCE)  , intent(inout)    ::   PCENV
   REAL,               INTENT(IN)       ::   DT
   integer, intent(IN)                  ::   IPC

   integer ::   IDS,IDE, JDS,JDE, KDS,KDE               & ! domain dims
               ,IMS,IME, JMS,JME, KMS,KME               & ! memory dims
               ,ITS,ITE, JTS,JTE, KTS,KTE                 ! tile   dims            )
   integer ::   ITIMESTEP, iinum, LP,L,LI

   REAL :: PRECRT,SNOWRT, QCOL_tot_0, QCOL_tot_1, rQCOL_tot_0, rQCOL_tot_1, riQCOL_tot_0

   real, dimension(PCE%L) :: MASS, PRES,DZ,QRCUTEN,QSCUTEN,QICUTEN,WVAR
   real, dimension(PCE%L) :: QVTEN,TETEN,TE,TE0,RHOA
   real, dimension(PCE%L) :: EFFL,EFFI,EFFS,EFFR,EFFH
   real, dimension(PCE%L) :: QRTEN,QSTEN,QITEN,QHTEN,QLTEN 
   real, dimension(PCE%L) :: NRTEN,NSTEN,NITEN,NHTEN,NLTEN 
   real, dimension(PCE%L) :: QRFTEN,QSFTEN,QIFTEN,QHFTEN,QLFTEN ! "fall"==>Sedimentation tends 
   REAL, DIMENSION(PCE%L) :: C2PREC,CSED,ISED,SSED,GSED,RSED    

   REAL, DIMENSION(PCE%L) :: PCQ, PCTH,PCQL,PCQI,PCQR,PCQS,PCQH
   REAL, DIMENSION(PCE%L) :: PCPR,PCTE,PCNL,PCNI,PCNR,PCNS,PCNH,PCDZ,PCW,PCA,PCRHO,RHOL


   LP   = PCE%L
         itimestep = INT(pce%time / dt )
   DZ   = C_LAYER_THICKNESS(PCE)
   TE   = C_TEMPERATURE( PCE, RKAP)
   MASS = C_LAYER_MASS( PCE )*100./grav ! layer mass is in hPa
   TE0  = TE

   WVAR = 0.5
   QRCUTEN=0.
   QSCUTEN=0. 
   QICUTEN=0.
   PRES =  PCE%PPL * 100.         
   RHOA =  PRES/(RDRY*TE)
   RHOL =  C_RHL( PCE)

          ! if(allocated(PCEDIAGS%RHOL)) PCEDIAGS%RHOL(:,IPC)=RHOL       
          if(allocated(PCEDIAGS%RHOA)) PCEDIAGS%RHOA(:,IPC)=RHOA       
          ! if(allocated(PCEDIAGS%MASS)) PCEDIAGS%MASS(:,IPC)=MASS       
 
   !!RHOA =  PCE%A * PRES/(RDRY*TE)
 
! not sure what this is doing
   iinum=1
   PCE%NL=0.


   QCOL_tot_0 = sum( MASS *PCE%A * ( PCE%Q+ PCE%QL+ PCE%QI+  PCE%QR+  PCE%QS+  PCE%QH ) , 1 )
   !rQCOL_tot_0 = sum( RHOA * ( PCE%Q+ PCE%QL+ PCE%QI+  PCE%QR+  PCE%QS+  PCE%QH )*DZ , 1 )

   write(*,*)  "  "
   write(*,*)  "CALLING mORR MUPHYS, PCE idx " ,IPC
   write(*,*)  "Start Mass(dm)  in Morr Wrapper ",QCOL_tot_0
   !write(*,*)  "Start Mass(rho) in Morr Wrapper ",rQCOL_tot_0
 



QLTEN=0.
QITEN=0.
QRTEN=0.
QSTEN=0.     
QHTEN=0.     
       
NLTEN=0.
NITEN=0.
NRTEN=0.
NSTEN=0.
NHTEN=0.   

TETEN=0.
QVTEN=0.

! invert height index
do L=1,LP
   LI      = LP+1-L
   PCRHO(L)= RHOL(Li)
   PCTE(L) = TE(Li)
   PCPR(L) = PRES(Li)
   PCDZ(L) = DZ(Li)
   PCA(L)  = PCE%A(Li)
   PCW(L)  = PCE%W(Li)
   PCQ(L)  = PCE%Q(Li)
   PCQL(L) = PCE%QL(Li)
   PCQI(L) = PCE%QI(Li)
   PCQR(L) = PCE%QR(Li)
   PCQS(L) = PCE%QS(Li)
   PCQH(L) = PCE%QH(Li)
   PCNL(L) = PCE%NL(Li)
   PCNI(L) = PCE%NI(Li)
   PCNR(L) = PCE%NR(Li)
   PCNS(L) = PCE%NS(Li)
   PCNH(L) = PCE%NH(Li)
end do

   !riQCOL_tot_0 = sum( PCRHO*PCA * ( PCQ + PCQL+ PCQI+  PCQR+  PCQS+  PCQH )*PCDZ , 1 )
   !write(*,*)  "Inv Start Mass(rho) in Morr Wrapper ",riQCOL_tot_0

   
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call MORR_TWO_MOMENT_MICRO(QLTEN,QITEN,QSTEN,QRTEN,         &
       NITEN,NSTEN,NRTEN, &
       PCQL,PCQI,PCQS,PCQR, & 
       PCNI,PCNS,PCNR,              &
       TETEN,QVTEN,PCTE,PCQ,PCPR,PCDZ,PCW,WVAR,PRECRT,SNOWRT,            &
       EFFL,EFFI,EFFS,EFFR,DT,                                                   &
!                                            IMS,IME, JMS,JME, KMS,KME,           &
!                                            ITS,ITE, JTS,JTE, KTS,KTE,           & ! 
                                              1,  1 ,  1 , 1,  1,PCE%L,           &
                                              1,  1 ,  1 , 1,  1,PCE%L,           &
                        QHTEN,NHTEN,PCQH,PCNH,EFFH,   & ! graupel vars
                        qrcuten,qscuten, qicuten,    & ! convective tendencies; ZERO here
                        QHFTEN,QRFTEN,QIFTEN,QSFTEN,QLFTEN, &
                        PCNL,NLTEN,iinum, & ! wrf-chem
				c2prec,CSED,ISED,SSED,GSED,RSED  &  ! hm added, wrf-chem
                        ,IPC, PCA , PCRHO )! IPC,PCA added ++jtb

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
do L=1,LP
   Li        = LP+1-L
   TE(Li)    = PCTE(L)
   PCE%Q(Li) = PCQ(L) 
   PCE%QL(Li)= PCQL(L)  
   PCE%QI(Li)= PCQI(L)
   PCE%QR(Li)= PCQR(L)
   PCE%QS(Li)= PCQS(L)
   PCE%QH(Li)= PCQH(L)
   PCE%NL(Li)= PCNL(L)
   PCE%NI(Li)= PCNI(L)
   PCE%NR(Li)= PCNR(L)
   PCE%NS(Li)= PCNS(L)
   PCE%NH(Li)= PCNH(L)
end do


   !!RHOA =  PCE%A * PRES/(RDRY*TE)
   !!rQCOL_tot_1 = sum( RHOA * ( PCE%Q+ PCE%QL+ PCE%QI+  PCE%QR+  PCE%QS+  PCE%QH )*DZ , 1 )


   QCOL_tot_1 = sum( MASS *PCE%A * ( PCE%Q+ PCE%QL+ PCE%QI+  PCE%QR+  PCE%QS+  PCE%QH ) , 1 )


   write(*,*)  "End Mass(dm)  in Morr Wrapper   ",QCOL_tot_1
   !write(*,*)  "End Mass(rho) in Morr Wrapper   ",rQCOL_tot_1
 
   write(*,*)  "Mass change(dm)  in Morr Wrapper ",QCOL_tot_1-QCOL_tot_0

   write(*,*)  "Precip flux         :",PRECRT
   write(*,*)  "Imbalance           :",PRECRT+SNOWRT + QCOL_tot_1-QCOL_tot_0

   write(*,*) " "

   PCE%THETA =  PCE%THETA + (TE-TE0)*( (1000./C_PPL(PCE))**RKAP )

   PCE%RAIN_GAUGE = PCE%RAIN_GAUGE + PRECRT
   PCE%SNOW_GAUGE = PCE%SNOW_GAUGE + SNOWRT


       end  SUBROUTINE MORRISON_WRAPPER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CONDENSE_CE1     ( PCE ,         & !  PCENV , &
                                DT          , & 
                                CP,GRAV,RKAP,ALHL,ALHS     )

   type (T_ONE_PCE)  , intent(inout)    ::   PCE
   !type (T_ONE_PCE)  , intent(inout)    ::   PCENV
   REAL,               INTENT(IN)       ::   DT
   REAL,               INTENT(IN)       ::   CP, GRAV, RKAP,ALHL,ALHS

   real :: TE(PCE%L) , DZ(PCE%L), X(PCE%L), PPL(PCE%L)
   real :: QST(PCE%L), DQ(PCE%L),DTE(PCE%L),DQST(PCE%L)
   INTEGER :: LP,I

   LP=PCE%L

   DZ   = C_LAYER_THICKNESS(PCE)

   QST =   QSAT  (  C_TEMPERATURE( PCE, RKAP)   & 
                  , C_PPL(PCE)                  )
   DQST=   DQSAT (  C_TEMPERATURE( PCE, RKAP)   & 
                  , C_PPL(PCE)                  )


 
      ! Save diagnostics
      !------------------------------------

   DQ  =   MAX( PCE%Q - QST , 0. )

   DQ  =   DQ / ( 1. + (ALHL/CP)*DQST )

   PCE%Q  = PCE%Q  - DQ
   PCE%QL = PCE%QL + DQ

   DTE    = ALHL*DQ / CP
   
   PCE%THETA =  PCE%THETA + DTE*( (1000./C_PPL(PCE))**RKAP )

  END SUBROUTINE CONDENSE_CE1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE AUTOCONV_CE1     ( PCE ,  &
                                DT       )
    

   type (T_ONE_PCE)  , intent(inout)    ::   PCE
   REAL,               INTENT(IN)       ::   DT

   real    :: DQ(PCE%L)
   INTEGER :: LP,I

   real, parameter    ::  QL_CRIT=1.0E-3
   real, parameter    ::  C_0    =3.0E-3


   LP=PCE%L

   DQ  =   MAX(  PCE%QL -  QL_CRIT  , 0. )

   DQ  =   C_0*DQ*DT

   PCE%QL  = PCE%QL - DQ
   PCE%QR  = PCE%QR + DQ


  END SUBROUTINE AUTOCONV_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE REVAP_CE1        ( PCE ,  &
                                DT       )
    

   type (T_ONE_PCE)  , intent(inout)    ::   PCE
   REAL,               INTENT(IN)       ::   DT

   real    :: DQ(PCE%L),QST(PCE%L),RHc(PCE%L),Tc(PCE%L),PL(PCE%L), AA, BB
   real    :: RnTvel(PCE%L),DQST(PCE%L),DZ(PCE%L),RnDiam(PCE%L)
   INTEGER :: LP,I,L

   real    :: EVAP, FallSpLiq, DTE, efactor, efactor2, TinLayer, DropRad, Ventil

   real, parameter    ::  QL_CRIT=1.0E-3
   real, parameter    ::  C_0    =1.0E-3


   LP     =  PCE%L
   RnTvel =  RAIN_TVEL( PCE )
   RnDiam =  RAINDROPSIZE(PCE)
   Tc     =  C_TEMPERATURE( PCE, RKAP)
   PL     =  C_PPL( PCE )
   QST    =  QSAT  (Tc, PL )
   DQST   =  DQSAT (Tc, PL )
   DZ     =  C_LAYER_THICKNESS(PCE)
   RHc    =  PCE%Q / QST
 
      ! Save diagnostics
      !------------------------------------

   Ventil =  3.     ! <<ventilation factor>> - ....  ha, ha, ha 
  
    do L=1,LP


         call MICRO_AA_BB(Tc(L),PL(L),QST(L),AA,BB)

         Efactor = ( 1.00 - RHc(L) ) / ( CE1_RHO_W * ( AA  + BB  ))  ! / (1.00 - RHc(L) )
         Efactor = MAX( Efactor, 0.0 )

             ! Get properties of typical rain drop
             !---------------------------------------

         FallSpLiq = RnTvel(L)
         DropRad   = RnDiam(L)/2.0
      


         
              ! Evaporate some rain into ddf
              !------------------------------------
                   ! First find right time scale
         TinLayer =  MIN( DT , DZ(L)/FallSpLiq )

         Efactor2 =   Efactor * Ventil / ( 1.0 + DQST(L)*GEOS_ALHL/GEOS_CP )/ (DropRad**2)
         EVAP     =   PCE%QR(L)*(1.0 - EXP( -TinLAYER*Efactor2  ) )

              ! Add effects of phase changes into in-cloud ddf quantities
         PCE%QR(L)    = PCE%QR(L)  - EVAP
         PCE%Q(L)     = PCE%Q(L)   + EVAP
         
         DTE          =  - EVAP * GEOS_ALHL / GEOS_CP

         PCE%THETA(L) = PCE%THETA(L) + DTE*( (1000./PL(L) )**RKAP )

         ! Save diagnostics
         !------------------------------------

    end do


  END SUBROUTINE REVAP_CE1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DRAIN_CE1     ( PCE )
    type (T_ONE_PCE)  , intent(inout)    ::   PCE
    INTEGER :: LP,L1
    real :: DENS(PCE%L) , WFX(PCE%L+1)
 
    !!MASS = C_LAYER_MASS(PCE)
    DENS = C_RHL(PCE)

    LP=PCE%L
    L1=PCE%L-1

    WFX = C_EDGE_WS( PCE ) !! - RAIN_TVEL( PCE )
    WFX(2:LP+1) = WFX(2:LP+1) - RAIN_TVEL( PCE )

    !!PCE%RAIN_GAUGE  = PCE%RAIN_GAUGE + MASS(LP)*PCE%A(LP)*PCE%QR(LP)      
    !!PCE%SNOW_GAUGE  = PCE%SNOW_GAUGE + MASS(LP)*PCE%A(LP)*PCE%QS(LP)      
    !!PCE%HAIL_GAUGE  = PCE%HAIL_GAUGE + MASS(LP)*PCE%A(LP)*PCE%QH(LP)      

    if ( WFX(L1) < 0. ) then 
       PCE%RAIN_GAUGE  = -DENS(L1)*PCE%A(L1)*WFX(L1)*PCE%QR(L1)      
       PCE%SNOW_GAUGE  = -DENS(L1)*PCE%A(L1)*WFX(L1)*PCE%QS(L1)      
       PCE%HAIL_GAUGE  = -DENS(L1)*PCE%A(L1)*WFX(L1)*PCE%QH(L1)      
    else
       PCE%RAIN_GAUGE  = 0.     
       PCE%SNOW_GAUGE  = 0.      
       PCE%HAIL_GAUGE  = 0.
    endif      
 


    PCE%QR(LP) = 0.00
    PCE%QS(LP) = 0.00
    PCE%QH(LP) = 0.00

  end SUBROUTINE DRAIN_CE1 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE W_QPREC_CE1 ( PCE , &
                           DT  )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCE
   REAL,                  INTENT(IN)       ::     DT

   real :: DENS(PCE%L) , DZ(PCE%L), X(PCE%L),WF(PCE%L+1),WT(PCE%L)
   real :: WFX(PCE%L+1)
   INTEGER :: L,I

   L=PCE%L

   DENS = C_RHL(PCE)
   DZ   = C_LAYER_THICKNESS(PCE)

       ! Precipitation falling. RAIN_TVEL is
       ! in ce1_uphys.F90. 
       ! Currently (9/5/13) RAIN_TVEL = 10. ms-1
       ! always.
       
#if 0
   WFX = C_EDGE_WS( PCE )
   WFX(2:L+1) = WFX(2:L+1) - RAIN_TVEL( PCE )
   WF(1:L) = -1.* WFX(1:L) * ( DT / DZ )
   WF(L+1) =  0. ! WFX(L+1)   
#endif


   WT  = RAIN_TVEL(PCE)
   WF(1:L) =  WT(1:L) * ( DT / DZ )
   WF(L+1) =  WF(L)   


   X=DENS*PCE%A_NM1*PCE%QR
   CALL FXPPM ( WF , X , 1 , L )
   PCE%QR=X/(PCE%A*DENS)

   

  END SUBROUTINE W_QPREC_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION RAIN_TVEL( PCE ) RESULT(WF)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)     :: WF
  INTEGER :: LP,L



     LP=PCE%L

     WF       = 10.0


  end FUNCTION RAIN_TVEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION RAIN_TVEL2( IPC , LP ) RESULT(WF)
  integer, intent(IN) :: IPC,LP
  REAL, DIMENSION(LP)     :: WF,WFS,WMX2
  INTEGER :: L,IPASS

     do L = 1,LP
       WF(L) = MAXVAL( (/ 0., PCEPHYS%RF(L,IPC), PCEPHYS%RF(L,IPC), PCEPHYS%RF(L,IPC) /)  )
     end do


     Wmx2=WF

     DO L=2,LP
        if( WF(L) > WMX2(L-1) ) then 
            WMX2(L)=WF(L)
        else 
            WMX2(L)=WMX2(L-1)
        end if
     END DO
     

  end FUNCTION RAIN_TVEL2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION RAINDROPSIZE( PCE ) RESULT(Diam)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)       :: Diam ! size in meters 
  INTEGER :: LP,L

     Diam = 0.001    
     !Diam = 0.00033    
          
   
  end FUNCTION RAINDROPSIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE FILLQS ( PCE , N  )
  type (T_ONE_PCE), intent(INout) :: PCE
  integer, intent(IN) :: N
  integer :: L,LP

  LP=PCE%L
  do L=1,LP
     IF( PCE%Q(L) < 0.0 ) PCE%Q(L)=0.
  end do  
   


  end SUBROUTINE FILLQS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MICRO_AA_BB(TEMP,PR,Q_SAT,AA,BB)

real, intent(in ) :: TEMP,PR,Q_SAT
real, intent(out) :: AA,BB

real  :: E_SAT

real, parameter  :: K_COND  =  2.4e-2        ! J m**-1 s**-1 K**-1
real, parameter  :: DIFFU   =  2.2e-5        ! m**2 s**-1
real, parameter  :: epsilon =  GEOS_H2OMW/GEOS_AIRMW

E_SAT = 100.* PR * Q_SAT /( EPSILON + (1.0-EPSILON)*Q_SAT )  ! (100 converts 
                                                             ! from mbar to Pa)

 AA  = ( GEOS_ALHL**2 ) / ( K_COND*GEOS_RVAP*(TEMP**2) )

 BB  = GEOS_RVAP*TEMP / ( DIFFU*(1000./PR)*E_SAT )

end subroutine MICRO_AA_BB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE CE1_UPHYS
