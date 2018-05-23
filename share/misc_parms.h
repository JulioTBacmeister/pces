      real, parameter          :: ZERO         = 0.0
      real, parameter          :: ONE          = 1.0
      real, parameter          :: TWO          = 2.0
      real, parameter          :: THREE        = 3.0
      real, parameter          :: FOUR         = 4.0
      real, parameter          :: FIVE         = 5.0
      real, parameter          :: SIX          = 6.0
      real, parameter          :: HALF         = ONE/TWO
      real, parameter          :: THIRD        = ONE/THREE
      real, parameter          :: FOURTH       = ONE/FOUR
      real, parameter          :: TIMFILT      = 0.05
      integer, parameter       :: MAX_STR      = 132
      real, parameter          :: PI           = 3.1415926535898
      real, parameter          :: DTOR         = PI/180.0
      integer, parameter       :: SECONDS_PER_DAY = 86400
      REAL   , PARAMETER       :: GRAV   = 9.81
      real, parameter          :: CP = 1003.5
      real, parameter          :: EARTH_RADIUS = 6.371E6
      REAL   , PARAMETER      :: ALHL   = 2.4548E6
      REAL   , PARAMETER      :: ALHS   = 2.8368E6
      REAL   , PARAMETER      :: STFBOL = 5.67E-8
      REAL   , PARAMETER      :: AIRMW  = 28.97
      REAL   , PARAMETER      :: H2OMW  = 18.01
      REAL   , PARAMETER      :: RUNIV  = 8314.3
      REAL   , PARAMETER      :: RVAP   = RUNIV/H2OMW
      REAL   , PARAMETER      :: RKAP   = 0.286
      REAL(kind=8)   , PARAMETER      :: RKAP_8   = 0.286
      REAL   , PARAMETER      :: RGAS   = CP*RKAP
      REAL   , PARAMETER      :: P00    = 1000.0
      REAL   , PARAMETER      :: HEATW  = 597.2
      REAL   , PARAMETER      :: HEATI  = 680.0
      REAL   , PARAMETER      :: TICE   = 273.16
      REAL   , PARAMETER      :: UNDEF  = -999.0
      REAL   , PARAMETER      :: SRFPRS = 984.7