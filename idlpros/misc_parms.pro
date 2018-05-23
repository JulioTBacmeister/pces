pro misc_parms
common misc_parms, alhl,alhs,cp,rgas,airmw,h2omw,grav,rkap,EARTH_RADIUs $
                 ,tice

            PI           = 3.1415926535898
            DTOR         = PI/180.0
           SECONDS_PER_DAY = 86400
            GRAV   = 9.81
            CP = 1003.5
            EARTH_RADIUS = 6.371E6
 

           ALHL   = 2.4548E6
           ALHS   = 2.8368E6
           STFBOL = 5.67E-8
           AIRMW  = 28.97
           H2OMW  = 18.01
           RUNIV  = 8314.3
           RVAP   = RUNIV/H2OMW
           RKAP   = 0.286
           RKAP_8   = 0.286
           RGAS   = CP*RKAP
           P00    = 1000.0
           HEATW  = 597.2
           HEATI  = 680.0
           TICE   = 273.16
           UNDEF  = -999.0
           SRFPRS = 984.7


return
end
