MODULE CE1_CONSTS

 IMPLICIT NONE
 PRIVATE


 PUBLIC PI

 PUBLIC GRAV 
 PUBLIC CP  
 PUBLIC RKAP 
 PUBLIC ALHL 
 PUBLIC ALHS 
 PUBLIC RGAS
 PUBLIC RDRY
 PUBLIC RVIR
 PUBLIC EP_2

 PUBLIC FILLING_CRITERION 

 PUBLIC CE1_GRAV 
 PUBLIC CE1_CP  
 PUBLIC CE1_RKAP 
 PUBLIC CE1_ALHL 
 PUBLIC CE1_ALHS 
 PUBLIC CE1_PI   
 PUBLIC CE1_FILLING_CRITERION 
 PUBLIC CE1_MIN_AREA_FAC


 character*72, parameter :: I_am_module = "ce1_consts" 

      real,  parameter, public :: CE1_RHO_W  =  1.0e3  ! Density of liquid water in kg/m^3
      real,  parameter, public :: Cnv2Pascal =  100.   ! Conversion of grid box pressure if needed

      logical, parameter, public :: pressure_by_slice = .FALSE.  ! do slice-wise pressure soln
                                                                 ! (see pp_ce1_in_ce1_dynamics)

      REAL,    PARAMETER :: PI=3.1415927


      REAL   , PARAMETER    :: CE1_GRAV   =  9.81
      real   , parameter    :: CE1_CP     =  1003.5
      REAL   , PARAMETER    :: CE1_RKAP   =  0.286
      REAL   , PARAMETER    :: CE1_ALHL   =  2.4548E6
      REAL   , PARAMETER    :: CE1_ALHS   =  2.8368E6
      REAL   , PARAMETER    :: CE1_PI     =  3.1415927

      REAL   , PARAMETER    :: CE1_MIN_AREA_FAC = 1.0e-3     ! Min allowed ratio of area to A0

      REAL   , PARAMETER    :: CE1_FILLING_CRITERION  = 10.0 ! Min allowed area

      REAL   , PARAMETER    :: FILLING_CRITERION  = 10.0 ! Min allowed area


      REAL   , PARAMETER    :: GRAV   =  9.81
      real   , parameter    :: CP     =  1003.5
      REAL   , PARAMETER    :: RKAP   =  0.286
      REAL   , PARAMETER    :: ALHL   =  2.4548E6
      REAL   , PARAMETER    :: ALHS   =  2.8368E6
      real   , PARAMETER    :: RGAS   =  CP*RKAP
      real   , PARAMETER    :: RDRY   =  CP*RKAP
      real   , PARAMETER    :: RVIR   =  461.6
      real   , PARAMETER    :: EP_2   =  Rdry/Rvir
   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constants from GEOS-5 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



real(kind=8), parameter, public :: GEOS_PI_R8     = 3.14159265358979323846
real, parameter, public :: GEOS_PI     = GEOS_PI_R8
real, parameter, public :: GEOS_GRAV   = 9.80                   ! m^2/s
real, parameter, public :: GEOS_RADIUS = 6376.0E3               ! m
real, parameter, public :: GEOS_OMEGA  = 2.0*GEOS_PI/86164.0    ! 1/s
real, parameter, public :: GEOS_ALHL   = 2.4548E6               ! J/kg
real, parameter, public :: GEOS_ALHS   = 2.8368E6               ! J/kg
real, parameter, public :: GEOS_ALHF   = GEOS_ALHS-GEOS_ALHL    ! J/kg
real, parameter, public :: GEOS_STFBOL = 5.6734E-8              ! W/(m^2 K^4)
real, parameter, public :: GEOS_AIRMW  = 28.97                  ! kg/Kmole
real, parameter, public :: GEOS_H2OMW  = 18.01                  ! kg/Kmole
real, parameter, public :: GEOS_O3MW   = 47.9982                ! kg/Kmole
real, parameter, public :: GEOS_RUNIV  = 8314.3                 ! J/(Kmole K)
real, parameter, public :: GEOS_KAPPA  = 2.0/7.0                ! --
real, parameter, public :: GEOS_RVAP   = GEOS_RUNIV/GEOS_H2OMW  ! J/(kg K)
real, parameter, public :: GEOS_RGAS   = GEOS_RUNIV/GEOS_AIRMW  ! J/(kg K)
real, parameter, public :: GEOS_CP     = GEOS_RGAS/GEOS_KAPPA   ! J/(kg K)
real, parameter, public :: GEOS_P00    = 100000.0               ! Pa
real, parameter, public :: GEOS_CAPWTR = 4218.                  ! J/(K kg)
real, parameter, public :: GEOS_RHOWTR = 1000.                  ! kg/m^3
real, parameter, public :: GEOS_NUAIR  = 1.533E-5               ! m^2/S (@ 18C)
real, parameter, public :: GEOS_TICE   = 273.16                 ! K
real, parameter, public :: GEOS_UNDEF  = 1.0e15                 ! --
real, parameter, public :: GEOS_SRFPRS = 98470                  ! Pa
real, parameter, public :: GEOS_KARMAN = 0.40                   ! --
real, parameter, public :: GEOS_USMIN  = 1.00                   ! m/s
real, parameter, public :: GEOS_VIREPS = GEOS_AIRMW/GEOS_H2OMW-1.0   ! --

integer,parameter, public :: GEOS_R8 = selected_real_kind(12) ! 8 byte real
integer,parameter, public :: GEOS_R4 = selected_real_kind( 6) ! 4 byte real
integer,parameter, public :: GEOS_RN = kind(1.0)              ! native real
integer,parameter, public :: GEOS_I8 = selected_int_kind (13) ! 8 byte integer
integer,parameter, public :: GEOS_I4 = selected_int_kind ( 6) ! 4 byte integer
integer,parameter, public :: GEOS_IN = kind(1)                ! native integer



end MODULE CE1_CONSTS
