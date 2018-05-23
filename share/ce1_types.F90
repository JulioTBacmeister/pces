!  $Id: ce1_types.F90,v 1.30 2008/07/03 20:37:11 bacmj Exp $
!-----------------------------
MODULE CE1_TYPES
!-----------------------------------------------------------
! !DESCRIPTION:
!
! This module contains derived types and some global integer
! codes and some non-physical constatnts used by pces
!-----------------------------------------------------------


 IMPLICIT NONE   
 PRIVATE

 PUBLIC T_ONE_PCE
 PUBLIC T_PCE_CONTROLS
 PUBLIC T_BCKG
 PUBLIC T_DIAG_CTL
 PUBLIC T_AGCM_GRID
 PUBLIC T_PCE_DIAGS 
 PUBLIC T_PCE_PHYS

 PUBLIC PCEPHYS
 PUBLIC CONTROL

 PUBLIC EnvironmentCode
 PUBLIC UpdraftCode
 PUBLIC MergedUpdraftCode
 PUBLIC DowndraftCode
 PUBLIC PrecipShaftCode
 PUBLIC AnvilCode
 PUBLIC LatentCode
 PUBLIC SenescentCode
 PUBLIC TrueCode
 PUBLIC FalseCode
 PUBLIC MissingCode
 PUBLIC MissingValue
 PUBLIC BadValue

 PUBLIC TinyFraction
 PUBLIC SmallFraction
 PUBLIC ModestFraction


 PUBLIC DormantValue
 PUBLIC DormantStatusCode
 !!PUBLIC CreatePrecShaftValue

 character*72, parameter :: I_am_module = "ce1_types" 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                     C O D E S   A N D   C O N S T A N T S 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real    , parameter :: TinyFraction    =  1.e-9
 real    , parameter :: SmallFraction   =  1.e-2
 real    , parameter :: ModestFraction  =  1.e-1
 real    , parameter :: TinyArea        =  100.  ! m+2


 integer , parameter :: EnvironmentCode =  0
 integer , parameter :: UpdraftCode     =  1
 integer , parameter :: DowndraftCode   =  2
 integer , parameter :: PrecipShaftCode =  3
 integer , parameter :: AnvilCode       =  10
 integer , parameter :: MergedUpdraftCode =  11
 integer , parameter :: LatentCode      = -1
 integer , parameter :: SenescentCode   = -100
 integer , parameter :: FalseCode       =  0
 integer , parameter :: TrueCode        =  1
 integer , parameter :: MissingCode     = -999
 real    , parameter :: MissingValue    = -999e9
 real    , parameter :: BadValue        = -999e9
 
 real    , parameter :: DormantValue         = -9999.9
 real    , parameter :: DormantStatusCode    = -999.9

 real    , parameter :: CreatePrecShaftValue = 1.0e-4 ! kg kg-1: Value of maxval (PCE%condensate) at which 
                                                      !   precip shaft is turned on.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                      D E R I V E D   D A T A   T Y P E S
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Basic "state" for prognostic convective elements.
!! Cleaned up to contain profiles of only (almost)  
!! prognostic quantities,as well as,  some other key 
!! non-profile quantities . Could argue with PPL,PPE
!! and ZGE, which are vertical grid descriptions and
!! probably redundant, or with A_NM1 which holds 
!! area profile of PCEs from last time-step, but is
!! used throughout ce1_dynamics (09/30/13)
!!------------------------------------------------------ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTERNAL CONVECTIVE ELEMENT PROPERTIES (ONE TIME LEVEL except for A)
  TYPE T_ONE_PCE
!-------------------------------
! ALLOCATABLE PROFILE QUANTITIES
!-------------------------------
! vertical grid
    REAL, DIMENSION(:), ALLOCATABLE                 :: ZGE   ! geometric height above local surface (m)
    REAL, DIMENSION(:), ALLOCATABLE                 :: PPL   ! pressure at mid-levels (hPa NOT Pa !!)
    REAL, DIMENSION(:), ALLOCATABLE                 :: PPE   ! pressure at edges (hPa NOT Pa !!)


! "dynamics" state
    REAL, DIMENSION(:), ALLOCATABLE                 :: U     ! Zonal wind (m/s)
    REAL, DIMENSION(:), ALLOCATABLE                 :: V     ! Meridional wind (m/s)
    REAL, DIMENSION(:), ALLOCATABLE                 :: W     ! Vertical wind (m/s)
    REAL, DIMENSION(:), ALLOCATABLE                 :: UR    ! Radial wind (m/s)
    REAL, DIMENSION(:), ALLOCATABLE                 :: THETA ! Potential Temp. (deviation from BCKG in K)
  !---------------------------------------------------------------------
  ! Needed for cleaner Anelastic solution strategy
  !---------------------------------------------------------------------
    REAL, DIMENSION(:), ALLOCATABLE                 :: THBCK ! background pot. temp. (K)

! water substance
    REAL, DIMENSION(:), ALLOCATABLE                 :: Q     ! water vapor mixing ratio (kg/kg) 
    REAL, DIMENSION(:), ALLOCATABLE                 :: QL    ! cloud liquid water mixing ratio (kg/kg) 
    REAL, DIMENSION(:), ALLOCATABLE                 :: QI    ! cloud ice water mixing ratio (kg/kg) 
    REAL, DIMENSION(:), ALLOCATABLE                 :: QR    ! rain mixing ratio (kg/kg) 
    REAL, DIMENSION(:), ALLOCATABLE                 :: QS    ! snow mixing ratio (kg/kg) 
    REAL, DIMENSION(:), ALLOCATABLE                 :: QH    ! hail/graupel mixing ratio (kg/kg) 

! second moments for uphys
    REAL, DIMENSION(:), ALLOCATABLE                 :: NL    ! cloud liquid particle number density (#/m^3) (I think)
    REAL, DIMENSION(:), ALLOCATABLE                 :: NI    ! cloud ice particle number density (#/m^3) (I think)
    REAL, DIMENSION(:), ALLOCATABLE                 :: NR    ! rain drop number density (#/m^3) (I think)
    REAL, DIMENSION(:), ALLOCATABLE                 :: NS    ! snowflake number density (#/m^3) (I think)
    REAL, DIMENSION(:), ALLOCATABLE                 :: NH    ! hail/graupel number density (#/m^3) (I think)

! areas and location
    REAL, DIMENSION(:), ALLOCATABLE                 :: A     ! cross-sctional area of PCE cloud obj (m^2)
    REAL, DIMENSION(:), ALLOCATABLE                 :: A_NM1 ! cross-sctional area of PCE cloud obj (m^2) at PCE previous timestep (here for numerical reasons)
    REAL, DIMENSION(:), ALLOCATABLE                 :: X     ! x-location of cloud obj (relative to something in m)
    REAL, DIMENSION(:), ALLOCATABLE                 :: Y     ! y-location of cloud obj (relative to something in m)



!-------------------------------
! TIME STUFF
!-------------------------------
    REAL, DIMENSION(6) :: MODELDATE ! year,mon,day,hour,min,sec
    REAL    :: TIME
    REAL    :: BIRTHTIME
    REAL    :: DEATHTIME
    REAL    :: AGE

!-------------------------------
! OTHER 0-Dim STUFF, PARAMETERS, CODES
!-------------------------------
    REAL    :: RAIN_GAUGE
    REAL    :: SNOW_GAUGE
    REAL    :: HAIL_GAUGE

    REAL    :: Z0
    REAL    :: PS
    REAL    :: A0
    REAL    :: ASH
    REAL    :: SHTOP
    REAL    :: AGRID
    REAL    :: XMEAN
    REAL    :: YMEAN


             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !! Below are used in CHECK_STATUS to track PCEs.
             !! MAXWTF is a "reddened" trace of a MAXVAL of 
             !! W (for updrafts) or Prec. Cond. (for shafts)
             !! This procedure needs to be cleaned up
             !! in check_status -01/08/15
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL    :: STATUS
    REAL    :: MAXW00
    REAL    :: MAXWTF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!   TYPE_OF_CE =  EnvironmentCode : Environment                       !!  
    !!              =  UpdraftCode     : Updraft                           !!
    !!              =  DowndraftCode   : Downdraft                         !!
    !!              =  LatentCode      : "Latent" (waiting for activation) !!
    !!              =  AnvilCode       : Anvil scud passively moving along !!
    !!              =  SenescentCode   : Dead (waiting for destruction)    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER :: TYPE_OF_CE


            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!  PCE Relationship Codes 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!   
    INTEGER :: MY_ENV
    INTEGER :: MY_DNDRFT
    INTEGER :: MY_PRSHFT
    INTEGER :: MY_UPDRFT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!  PCE Management Codes 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!   
    INTEGER :: NEED_PRSHFT
    INTEGER :: NEED_DNDRFT


            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!  PCE Identifiers 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!   
    INTEGER :: IPC
    INTEGER :: L 
    INTEGER :: K
    INTEGER :: I
    INTEGER :: J
    INTEGER :: NPOP
    CHARACTER*32 :: CONFIG
  END TYPE T_ONE_PCE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





! BACKGROUND/GRIDBOX PROPERTIES (ONE TIME LEVEL)
  TYPE T_AGCM_GRID
               ! Grid type 0=unstruc., 2=latlon or xy
    INTEGER                                         :: GRID_TYPE
               ! Logically-rectangular 2D horz
    REAL, DIMENSION(:,:),   ALLOCATABLE                 :: AREA2
    REAL, DIMENSION(:,:),   ALLOCATABLE                 :: XX2  ! corner locs
    REAL, DIMENSION(:,:),   ALLOCATABLE                 :: YY2  !     "
    REAL, DIMENSION(:,:),   ALLOCATABLE                 :: DX 
    REAL, DIMENSION(:,:),   ALLOCATABLE                 :: DY
               ! Unstructured
    REAL, DIMENSION(:),     ALLOCATABLE                 :: AREA1
    REAL, DIMENSION(:),     ALLOCATABLE                 :: XX1
    REAL, DIMENSION(:),     ALLOCATABLE                 :: YY1

   end TYPE T_AGCM_GRID


! BACKGROUND/GRIDBOX PROPERTIES (ONE TIME LEVEL)
! HERE A LOGICALLY-REC. XY GRID IS ASSUMED
  TYPE T_BCKG
    REAL, DIMENSION(:,:),   ALLOCATABLE                 :: AREA
    REAL, DIMENSION(:,:),   ALLOCATABLE                 :: XX
    REAL, DIMENSION(:,:),   ALLOCATABLE                 :: YY
    REAL, DIMENSION(:,:),   ALLOCATABLE                 :: ZSF    ! Sfc topo ht (m)
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: U
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: V
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: TH
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: W
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: Q
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: QI
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: QL
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: QR
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: QS
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: QH
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: PPL
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: PPE
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: ZGE
    REAL, DIMENSION(:,:,:), ALLOCATABLE                 :: AA
 
    INTEGER                                         :: IMS
    INTEGER                                         :: JMS
    INTEGER                                         :: LMS
   
  END TYPE T_BCKG

! TIME TENDENCIES FROM MODEL
  TYPE T_BCKG_TNDS
    REAL, DIMENSION(:), ALLOCATABLE                 :: TH
    REAL, DIMENSION(:), ALLOCATABLE                 :: U
    REAL, DIMENSION(:), ALLOCATABLE                 :: V
    REAL, DIMENSION(:), ALLOCATABLE                 :: W
    REAL, DIMENSION(:), ALLOCATABLE                 :: Q
    REAL, DIMENSION(:), ALLOCATABLE                 :: QL
    REAL, DIMENSION(:), ALLOCATABLE                 :: QI
    REAL, DIMENSION(:), ALLOCATABLE                 :: QR
    REAL, DIMENSION(:), ALLOCATABLE                 :: QS
    REAL, DIMENSION(:), ALLOCATABLE                 :: QH
  END TYPE T_BCKG_TNDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 TYPE T_PCE_CONTROLS

                                ! Basic process controls for ce1.F90
                                !-----------------------------------
    logical      ::  do_dynamics
    logical      ::     do_subcloud_flux  
    logical      ::  do_vdiff
    logical      ::  do_uphys
    logical      ::  do_precipitation
    logical      ::  do_interactions
    logical      ::  do_subsidence

                                ! Numerical/stability parameters. Added 11/25/2014.
                                ! These are primarily numerical band-aids and
                                ! there doesn't seem to be a good reason to tune
                                ! them extensively
                                !---------------------------------------------
                      ! *FUNDAMENTAL PARAMETER* 
                      ! minimum allowed area of PCES (invoked in ent and precip)
                      ! but has global impact on PCE-shape
    real         ::  minallowedarea_ratio=1.0e-2 ! 1.0e-2
                      !
                      ! turns on/off cleanup routine in dyn
    logical      ::  cleanup_after_dyn=.false.
                       ! area at which cleanup routine kicks in
    real         ::  cleanup_maxarea_ratio=1.e-1
                       ! area at which cleanup routine maxes
    real         ::  cleanup_minarea_ratio=1.e-2
                       ! turns on/off adhoc smoothing routine (dyn, not used since 2013 or so)
    logical      ::  adhoc_smoothing_after_dyn=.false.
                       ! turns on/off time smoothing routine
    logical      ::  time_smoothing =.false.
                       ! minimum radius at which 3D reconstruction is done (3D pressure)
    real         ::  pr3d_recon_radius=-1. ! (m)
                       ! minimum area for FV type z-momentum eq's (3D pressure)
    real         ::  zmom_minarea_ratio=1.e-6 ! 1.e-4
                       ! area* at which w damping kicks in (3D pressure)
    real         ::  zmom_wdamp_areax = 1.e-6 ! 1.e-2
                       ! area* at which w damping time (3D pressure)
    real         ::  zmom_wdamp_tau = 1000000. ! 10.
                       ! not sure why this is here
    real         ::  shutoff_ent_area_ratio=5.0e-6 ! 5.0e-2
                       ! turns on/off limit on PCE area
    LOGICAL      ::  LIMIT_CE_AREA=.false.          
                       ! max allowed frac w/resp to grid box area
    REAL         ::  MAX_CE_AREAL_FRACTION=0.25
                       ! turns on/off CFL limiter on PCE W
    LOGICAL      ::  CFL_LIMITER =.true.          

                           
 
                                ! Namelist parameters for various processes
                                !---------------------------------------------
    CHARACTER*32 ::  CONFIG
    REAL         ::  INITIAL_RADIUS
    REAL         ::  INITIAL_W
    REAL         ::  THETA_PERT
    REAL         ::  CENTER
    REAL         ::  SPREAD
                                ! esponge layer parameters  11/01/13
    LOGICAL      ::  USE_SPONGE_LYR
    REAL         ::  SPONGE_LYR_BASE
                                ! entrainment parameters (noted) 4/11/12
    LOGICAL      ::  USE_UR_2_UE
    REAL         ::  UR_2_UER_ENT
    REAL         ::  UR_2_UER_DET
    LOGICAL      ::  USE_W_2_UE
    REAL         ::  W_2_UER
    LOGICAL      ::  USE_FIXED_UE
    REAL         ::  UER_FIXED
    REAL         ::  MUNCH_VEL
    INTEGER      ::  INEST_FOR_ENTR_ENV
    LOGICAL      ::  DO_GUST_MUNCHING
    REAL         ::  GUST_MUNCH_UR
    INTEGER      ::  SHAFT_ENTRAINMENT
                                ! dynamical parameters
    INTEGER      ::  BUOY_VIRTUAL
    INTEGER      ::  BUOY_CONDLOAD
    INTEGER      ::  THETA_ADVECTION
    REAL         ::  W_NU_DEL4
    REAL         ::  THETA_NU_DEL4
    REAL         ::  SUBCLOUD_ADV =  0.666 ! 0.333 ! 0.1666 ! 0.333 ! 0.666  ! 0.333
                                ! precribed plume parameters 4/11/12
    LOGICAL      ::  USE_PRESCRIBED_PLUME
    REAL         ::  PRESCRIBED_PLUME_HEIGHT
    REAL         ::  PRESCRIBED_PLUME_W
    REAL         ::  PRESCRIBED_PLUME_BASER
                                ! precipitation shaft parameters 4/14/12
    REAL         ::  CreatePrecShaftValue
    LOGICAL      ::  USE_SIMPLE_SHARE_PREC
    LOGICAL      ::  SHARE_THROUGH_COLUMN
    LOGICAL      ::  SHAFT_PASS4DYNMX
    INTEGER      ::  SHARE_PRECIP_V
    INTEGER      ::  SHARE_PRECIP_COMPLEXITY
    REAL         ::  SIMPLE_SHARE_PREC_FRAC !added 8/26/13
    REAL         ::  SIMPLE_SHARE_AREA_FCTR !added 8/30/13
    LOGICAL      ::  TOP_DOWN_SHAFT ! 9/07/13
    LOGICAL      ::  UPDFT_PRECIP_COLLAR ! 10/16/13
                               ! pce management parameters 4/14/12
    LOGICAL      ::  NEVER_KILL_UPDS
    REAL         ::  TINY_CONDENSATE_MR
                               ! uphys parameters
    LOGICAL      ::  USE_MORRISON_UPHYS
                                ! enviro. adjustment parameters 9/28/12
    LOGICAL      ::  DO_ENVIRO_SUBSI
                                ! diffusion/turbulence parameters 11/08/2013
    REAL         ::  KZZ_BCKG
    REAL         ::  VDIFF_LAMBDA
    REAL         ::  WDAMP_TAU
                                ! some more updraft parameters 07/10/2014
    REAL         ::  CD_HORZDRAG

                                ! pressure solution code
                                !   1 => old elliptic style
    INTEGER      ::  PSOLV

                                ! Initial bubble params
    REAL         :: INITBBL_RADMAX
    REAL         :: INITBBL_DPTH
    REAL         :: INITBBL_ZCEN
    REAL         :: INITBBL_THP0
    REAL         :: INITBBL_W0

 END TYPE T_PCE_CONTROLS

 TYPE T_DIAG_CTL
    CHARACTER*32 ::  CONFIG
    REAL         ::  LAST_WRITE = -99999.9
    REAL         ::  WRITE_INTERVAL = -999.99
    REAL         ::  LAST_WRITE_BCKG = -99999.9
    REAL         ::  WRITE_INTERVAL_BCKG = -999.99
 END TYPE T_DIAG_CTL


!================================================
! CONVECTIVE ELEMENT DIAGNOSTIC QUANTITIES.
!================================================ 
  TYPE T_PCE_DIAGS
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: DHN
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: MUTT
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: PRACS
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: TMUPH
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: QMUPH
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: QVSL
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: QVSI

    REAL, DIMENSION(:,:), ALLOCATABLE                 :: E
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: RHOA
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: RHOL
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: MASS
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: DZRW
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: PFZ

    REAL, DIMENSION(:,:), ALLOCATABLE                 :: PP
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: PHYD
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: RCDOT

    REAL, DIMENSION(:  ), ALLOCATABLE                 :: QTBKEN
    REAL, DIMENSION(:  ), ALLOCATABLE                 :: QTBKDE
    REAL, DIMENSION(:  ), ALLOCATABLE                 :: ABKEN
    REAL, DIMENSION(:  ), ALLOCATABLE                 :: ABKDE
    REAL, DIMENSION(:  ), ALLOCATABLE                 :: SHEXTOP
    REAL, DIMENSION(:  ), ALLOCATABLE                 :: SHEXBOT
    REAL, DIMENSION(:  ), ALLOCATABLE                 :: KEI
    REAL, DIMENSION(:  ), ALLOCATABLE                 :: DKEI1
    
 END TYPE T_PCE_DIAGS



!=================================
! DECLARE a CONTROL vector 
!=================================
 type (T_PCE_CONTROLS) , save  ::  CONTROL



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE T_PCE_PHYS
!-------------------------------
! ALLOCATABLE PROFILE QUANTITIES
!-------------------------------
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: HF
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: SF
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: RF
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: UR
    REAL, DIMENSION(:,:), ALLOCATABLE                 :: WDZA
 end type T_PCE_PHYS

!=================================
! DECLARE a PHYS coupler
!=================================
 type (T_PCE_PHYS) , save  ::  PCEPHYS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
END MODULE CE1_TYPES
