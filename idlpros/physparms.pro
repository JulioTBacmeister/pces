pro physparms,radius_e=radius_e  $
             ,grav=grav  $
             ,cpd=cp $
             ,rkap=rkap $
             ,ALHL=alhl $
             ,alhs=alhs $
             ,alhf=alhf $
             ,rgas=rgas

;      REAL   , PARAMETER    :: GRAV   =  9.81
;      real   , parameter    :: CP     =  1003.5
;      REAL   , PARAMETER    :: RKAP   =  0.286
;      REAL   , PARAMETER    :: ALHL   =  2.4548E6
;      REAL   , PARAMETER    :: ALHS   =  2.8368E6
;      real   , PARAMETER    :: RGAS   =  CP*RKAP
;      real   , PARAMETER    :: RDRY   =  CP*RKAP
;      real   , PARAMETER    :: RVIR   =  461.6
;      real   , PARAMETER    :: EP_2   =  Rdry/Rvir

RADIUS_E = 6376.0E3  
GRAV   =  9.81
CP     =  1003.5
RKAP   =  0.286
ALHL   =  2.4548E6
ALHS   =  2.8368E6
ALHF   =  ALHS - ALHL
RGAS   =  CP*RKAP
RDRY   =  CP*RKAP
RVIR   =  461.6
EP_2   =  Rdry/Rvir


return
end
