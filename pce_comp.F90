!  11/27/2015

module pce_comp_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use CE1_TYPES


IMPLICIT NONE
PRIVATE

PUBLIC PCES
PUBLIC pce_init_nml
PUBLIC pce_init_pces_0
PUBLIC pce_run_sngl

type (T_ONE_PCE), allocatable, save, dimension(:) ::  PCES

contains
  
!=====================================================
subroutine pce_init_nml
use ce1_diags, only : DIAGCTL
  
real                            :: DIAG_DT=50.
real                            :: CreatePrecShaftValue = 1.0e-4
logical                         :: USE_SIMPLE_SHARE_PREC = .false.
real                            :: SIMPLE_SHARE_PREC_FRAC = 1.0
real                            :: SIMPLE_SHARE_AREA_FCTR = 1.0
real                            :: UER_FIXED,W_2_UER,PRESCRIBED_PLUME_HEIGHT,PRESCRIBED_PLUME_W,PRESCRIBED_PLUME_BASER &
                                  ,MUNCH_VEL,UR_2_UER_ENT,UR_2_UER_DET
logical                         :: USE_FIXED_UE,USE_W_2_UE,USE_PRESCRIBED_PLUME,USE_UR_2_UE
logical                         :: NEVER_KILL_UPDS = .false.
real                            :: TINY_CONDENSATE_MR = 1.e-6
integer                         :: SHARE_PRECIP_V = 2
integer                         :: SHARE_PRECIP_COMPLEXITY = 1
integer                         :: INEST_FOR_ENTR_ENV = 1

logical                         :: DO_GUST_MUNCHING = .false.
real                            :: GUST_MUNCH_UR = 50.
!logical                         :: LIMIT_CE_AREA = .false.
!real                            :: MAX_CE_AREAL_FRACTION = 1.00
integer                         :: SHAFT_ENTRAINMENT = 1
logical                         :: SHAFT_PASS4DYNMX = .false. 
logical                         :: SHARE_THROUGH_COLUMN = .false. 
logical                         :: TOP_DOWN_SHAFT = .false.
logical                         :: UPDFT_PRECIP_COLLAR = .true.

logical                         :: DO_ENVIRO_SUBSI = .true.

logical                         :: USE_MORRISON_UPHYS = .false.

logical                         :: USE_SPONGE_LYR  = .false.
real                            :: SPONGE_LYR_BASE = 50000. ! meters

real                            :: KZZ_BCKG      = 10.   ! m+2 s-1
real                            :: VDIFF_LAMBDA  = 10.   ! m 
real                            :: WDAMP_TAU     = 1000. ! s 

integer                         :: BUOY_CONDLOAD = 1
integer                         :: BUOY_VIRTUAL  = 1
!!integer                         :: W_ADVECTION   = 1 ! 0=PPM,1=3rd UPST ADV
integer                         :: THETA_ADVECTION   = 1 ! -1=Split- BCKG(2nd Order), pert(PPM)
                                                         ! 0=PPM,1=3rd UPST ADV
real                            :: W_NU_DEL4 = 0.003
real                            :: THETA_NU_DEL4 = 0.003

real                            :: INITBBL_RADMAX = 10000.
real                            :: INITBBL_ZCEN   = 300.
real                            :: INITBBL_THP0   = 1.0
real                            :: INITBBL_W0     = 0.5
real                            :: INITBBL_DPTH   = 500.


real                            :: CD_HORZDRAG = 1.


integer                         :: PSOLV=1

integer :: UNIT

    logical      ::  do_dynamics      = .true.
    logical      ::     do_subcloud_flux      = .true.
    logical      ::  do_vdiff         = .true.
    logical      ::  do_uphys         = .true.
    logical      ::  do_precipitation = .true.
    logical      ::  do_interactions  = .true.
    logical      ::  do_subsidence    = .true.


namelist /diags/ DIAG_DT

namelist /entparams/ UER_FIXED,W_2_UER,UR_2_UER_ENT,UR_2_UER_DET,USE_FIXED_UE,USE_W_2_UE,USE_UR_2_UE,MUNCH_VEL,INEST_FOR_ENTR_ENV &
                    , DO_GUST_MUNCHING, GUST_MUNCH_UR   &
                    , SHAFT_ENTRAINMENT

namelist /pplumeparams/ USE_PRESCRIBED_PLUME,PRESCRIBED_PLUME_HEIGHT & 
                      , PRESCRIBED_PLUME_W,PRESCRIBED_PLUME_BASER

namelist /shaftparams/ CreatePrecShaftValue,USE_SIMPLE_SHARE_PREC,SHAFT_PASS4DYNMX   & 
                     , SHARE_PRECIP_V,SIMPLE_SHARE_PREC_FRAC,SIMPLE_SHARE_AREA_FCTR  & 
                     , TOP_DOWN_SHAFT,SHARE_THROUGH_COLUMN,SHARE_PRECIP_COMPLEXITY &
                     , UPDFT_PRECIP_COLLAR 

namelist /mgmtparams/ NEVER_KILL_UPDS,TINY_CONDENSATE_MR

namelist /uphysparams/ USE_MORRISON_UPHYS

namelist /spongeparams/ USE_SPONGE_LYR,SPONGE_LYR_BASE

namelist /adjparams/ DO_ENVIRO_SUBSI

namelist /dynparams/ BUOY_CONDLOAD,BUOY_VIRTUAL,THETA_ADVECTION,THETA_NU_DEL4,W_NU_DEL4

namelist /vdiffparams/ KZZ_BCKG, VDIFF_LAMBDA, WDAMP_TAU

namelist /initbblparams/ INITBBL_RADMAX,INITBBL_ZCEN,INITBBL_THP0,INITBBL_W0,INITBBL_DPTH

namelist /psolvparams/ PSOLV

namelist /processcntrl/ do_dynamics,do_vdiff,do_uphys,do_precipitation,do_interactions,do_subsidence,do_subcloud_flux

OPEN( UNIT=UNIT, FILE="control.nml" ) !, NML =  cntrls )
READ( UNIT=UNIT, NML=diags)
READ( UNIT=UNIT, NML=adjparams)
READ( UNIT=UNIT, NML=entparams)
READ( UNIT=UNIT, NML=dynparams)
READ( UNIT=UNIT, NML=vdiffparams)
READ( UNIT=UNIT, NML=spongeparams)
READ( UNIT=UNIT, NML=pplumeparams)
READ( UNIT=UNIT, NML=shaftparams)
READ( UNIT=UNIT, NML=uphysparams)
READ( UNIT=UNIT, NML=mgmtparams)
READ( UNIT=UNIT, NML=initbblparams)
READ( UNIT=UNIT, NML=psolvparams)
READ( UNIT=UNIT, NML=processcntrl)
CLOSE(UNIT=UNIT)

CONTROL% W_2_UER      =  W_2_UER
CONTROL% UR_2_UER_ENT     =  UR_2_UER_ENT
CONTROL% UR_2_UER_DET     =  UR_2_UER_DET
CONTROL% UER_FIXED    =  UER_FIXED
CONTROL% USE_FIXED_UE =  USE_FIXED_UE
CONTROL% USE_W_2_UE   =  USE_W_2_UE
CONTROL% USE_UR_2_UE  =  USE_UR_2_UE
CONTROL% MUNCH_VEL    =  MUNCH_VEL
CONTROL% INEST_FOR_ENTR_ENV      =  INEST_FOR_ENTR_ENV
CONTROL% DO_GUST_MUNCHING        =  DO_GUST_MUNCHING
CONTROL% GUST_MUNCH_UR           =  GUST_MUNCH_UR
CONTROL% SHAFT_ENTRAINMENT       =  SHAFT_ENTRAINMENT 

CONTROL% BUOY_VIRTUAL            =  BUOY_VIRTUAL
CONTROL% BUOY_CONDLOAD           =  BUOY_CONDLOAD
CONTROL% THETA_ADVECTION         =  THETA_ADVECTION
CONTROL% W_NU_DEL4               =  W_NU_DEL4
CONTROL% THETA_NU_DEL4           =  THETA_NU_DEL4

CONTROL% USE_PRESCRIBED_PLUME    =  USE_PRESCRIBED_PLUME
CONTROL% PRESCRIBED_PLUME_HEIGHT =  PRESCRIBED_PLUME_HEIGHT
CONTROL% PRESCRIBED_PLUME_W      =  PRESCRIBED_PLUME_W
CONTROL% PRESCRIBED_PLUME_BASER  =  PRESCRIBED_PLUME_BASER

CONTROL% CreatePrecShaftValue    =  CreatePrecShaftValue
CONTROL% USE_SIMPLE_SHARE_PREC   =  USE_SIMPLE_SHARE_PREC
CONTROL% SHARE_PRECIP_V          =  SHARE_PRECIP_V
CONTROL% SHARE_PRECIP_COMPLEXITY =  SHARE_PRECIP_COMPLEXITY
CONTROL% SIMPLE_SHARE_PREC_FRAC  =  SIMPLE_SHARE_PREC_FRAC
CONTROL% SIMPLE_SHARE_AREA_FCTR  =  SIMPLE_SHARE_AREA_FCTR
CONTROL% SHAFT_PASS4DYNMX        =  SHAFT_PASS4DYNMX
CONTROL% SHARE_THROUGH_COLUMN    =  SHARE_THROUGH_COLUMN
CONTROL% TOP_DOWN_SHAFT          =  TOP_DOWN_SHAFT  
CONTROL% UPDFT_PRECIP_COLLAR     =  UPDFT_PRECIP_COLLAR
!!! overide choice of top down if version=4
if (CONTROL% SHARE_PRECIP_V==4) CONTROL% TOP_DOWN_SHAFT=.true.  

CONTROL% NEVER_KILL_UPDS         =  NEVER_KILL_UPDS
CONTROL% TINY_CONDENSATE_MR      =  TINY_CONDENSATE_MR

CONTROL% DO_ENVIRO_SUBSI         =  DO_ENVIRO_SUBSI
CONTROL% USE_MORRISON_UPHYS      =  USE_MORRISON_UPHYS

CONTROL% USE_SPONGE_LYR          =  USE_SPONGE_LYR
CONTROL% SPONGE_LYR_BASE         =  SPONGE_LYR_BASE

CONTROL% KZZ_BCKG                =  KZZ_BCKG
CONTROL% VDIFF_LAMBDA            =  VDIFF_LAMBDA
CONTROL% WDAMP_TAU               =  WDAMP_TAU

CONTROL% CD_HORZDRAG             =  CD_HORZDRAG

CONTROL% PSOLV                   =  PSOLV

CONTROL% DO_DYNAMICS             = DO_DYNAMICS
  CONTROL% DO_SUBCLOUD_FLUX        =   DO_SUBCLOUD_FLUX
CONTROL% DO_UPHYS                = DO_UPHYS
CONTROL% DO_VDIFF                = DO_VDIFF
CONTROL% DO_PRECIPITATION        = DO_PRECIPITATION
CONTROL% DO_SUBSIDENCE           = DO_SUBSIDENCE
CONTROL% DO_INTERACTIONS         = DO_INTERACTIONS

CONTROL %INITBBL_RADMAX        = INITBBL_RADMAX    
CONTROL %INITBBL_DPTH          = INITBBL_DPTH
CONTROL %INITBBL_ZCEN          = INITBBL_ZCEN
CONTROL %INITBBL_THP0          = INITBBL_THP0
CONTROL %INITBBL_W0            = INITBBL_W0

DIAGCTL%WRITE_INTERVAL         = DIAG_DT


end subroutine pce_init_nml

!============================================================
subroutine pce_init_pces_0(LP,NP)
use ce1_utils, only :  create_ce1
  
 integer, intent(in) :: LP,NP

 integer :: KP,N
 
 allocate( PCES(NP) )

 write(*, *) " INITIALIZING ",NP," elements !! "

   ! Initial Allocation of memory-space for PCES. 
   ! All fields allocd and set to zero or prespec flags
   !-----------------------------------------------------
   KP=1
   do N=1,NP
      call CREATE_CE1 ( LP , KP , PCES(N) )  ! in ce1_utils
   end do
 write(*, *) " DONE INITIALIZING ",NP," elements !! "


end subroutine pce_init_pces_0
!============================================================

subroutine pce_run_sngl(dt,nsub)
use ce1, only : run_ce1
  
  REAL, intent(in)    :: DT
  INTEGER, intent(in) :: nsub
  
      REAL   , PARAMETER    :: GRAV   =  9.81
      real   , parameter    :: CP     =  1003.5
      REAL   , PARAMETER    :: RKAP   =  0.286
      REAL   , PARAMETER    :: ALHL   =  2.4548E6
      REAL   , PARAMETER    :: ALHS   =  2.8368E6

      integer :: i
      
   do i=1,nsub
      CALL RUN_CE1 (   PCES     , &      ! in module ce1 
                       DT       , &
                       CP       , & 
                       GRAV     , & 
                       RKAP     , & 
                       ALHL     , & 
                       ALHS       )

       write(*,*) " did step " , I, PCES(1)%MODELDATE
   end do

      

      
end subroutine pce_run_sngl
!=============================================================
end module pce_comp_mod
