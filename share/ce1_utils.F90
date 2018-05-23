!  $Id: ce1_utils.F90,v 1.53 2008/10/08 20:32:56 bacmj Exp $
!
!-----------------------------
MODULE CE1_UTILS

use CE1_TYPES
use CE1_INFORM
use NUMERICAL_UTILITIES, only : QSAT, DQSAT
use CE1_CONSTS

 IMPLICIT NONE
 PRIVATE

 PUBLIC ADVANCE_TIME

 PUBLIC LOCATE_PCENV
 PUBLIC CREATE_CE1
 PUBLIC DESTROY_CE1
 PUBLIC COPY_PCE
 PUBLIC MERG_PCE
 PUBLIC INIT_ONE_UPD
 PUBLIC INIT_ONE_ENVIRO

 PUBLIC STARTBUBBLE

 PUBLIC SEPXY
 PUBLIC OVRCIRCLE

 PUBLIC C_PPL
 PUBLIC X_ZGE
 PUBLIC C_ZGL
 PUBLIC X_PPE

 PUBLIC C_THETA_D
 PUBLIC SET_PCE_ZGRID_1
 PUBLIC C_RADIUS
 PUBLIC C_TEMPERATURE
 PUBLIC C_LAYER_MASS
 PUBLIC C_LAYER_THICKNESS
 PUBLIC C_RHL
 PUBLIC C_W

 PUBLIC C_EDGE_WS
 PUBLIC C_MASS_XSCT2
 PUBLIC C_TYPICAL_RADIUS

 !PUBLIC C_RADIAL_AVG_4_BUOY

 PUBLIC C_TYPICAL_AREA

 PUBLIC ZINTRP
 
 character*72, parameter :: I_am_module = "ce1_utils" 

 CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION ZINTRP( AIN, ZIN, ZTG ) RESULT(ATG)
  REAL, DIMENSION(:)         , intent(IN)  ::  AIN, ZIN, ZTG
  REAL, DIMENSION(size(ZTG))               ::  ATG

  real :: ZMAX, ZMIN
  integer :: L ,LIN, LTG, L2

  LTG = size(ZTG)
  LIN = size(ZIN)
 
  do L=1,LTG

     if ( ZTG(L) >= ZIN(1)   ) THEN
       ATG(L) = AIN(1)
       !exit
     endif
     if ( ZTG(L) <= ZIN(LIN) ) THEN
       ATG(L) = AIN(LIN)
       !exit
     endif

     do L2=2,LIN-1

        if ( ( ZIN(L2) >= ZTG(L) ) .AND. ( ZTG(L) > ZIN(L2+1) ) ) then        
           ATG(L)  = AIN(L2+1) + (ZTG(L)-ZIN(L2+1))*(AIN(L2)-AIN(L2+1))/(ZIN(L2)-ZIN(L2+1))
        endif

     end do

  end do

  end function ZINTRP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION LAT_FRM_Y( Y,A_EARTH ) RESULT(LAT)
  REAL, intent(IN) :: Y, A_EARTH
  REAL             :: LAT, RLAT


       !! DX =  PI *A_EARTH * COS(RLAT0) * DX / 180.
  RLAT= PI*LAT/180.
  LAT = 180.* Y / ( PI *A_EARTH )
       
   

  end function LAT_FRM_Y


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION SEPXY( R1,R2 ) RESULT(DR)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INPUTS (from PCEs) : 
  !   R[1,2](1,2)        = Position coordinates in m
  ! OUTPUT:
  !   DR                 = separation in m
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL, DIMENSION(2), intent(IN) ::  R1,R2
  REAL                           ::  DR
  REAL                           ::  DX,DY


        DX = ( R1(1) - R2(1) )

        DY = ( R1(2) - R2(2) )

        DR = SQRT( DY**2 + DX**2 )

  end function SEPXY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION OVRCIRCLE( R1,R2,D ) RESULT(AR)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INPUTS (from PCEs) : 
  !   R[1,2](1,2)        = Position coordinates in m
  ! OUTPUT:
  !   DR                 = separation in m
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL,               intent(IN) ::  R1,R2,D
  REAL                           ::  AR
  REAL                           ::  RR,r,cos1,cos2,ang1,ang2,ar1,ar2

   RR = r1
   r  = r2
   if (r1.lt.r2) then
      RR = r2
      r  = r1
   endif

   ! no separation of centers
   if (d/RR .lt. 1e-4) then
      AR = pi*r**2
      return
   endif

   ! no overlap
   if (d .gt. (r+RR)) then
      AR = 0.
      return
   endif

   !smaller circle totally contained
   if ((d+r) .lt. RR) then
      ar=pi*r**2
      return
   endif


   cos1 = (d**2+r1**2-r2**2) / (2*d*r1 ) 
   cos2 = (d**2+r2**2-r1**2) / (2*d*r2 )

   if (cos1.lt.-0.99999) cos1=-0.99999
   if (cos2.lt.-0.99999) cos2=-0.99999
   if (cos1.gt. 0.99999) cos1= 0.99999
   if (cos2.gt. 0.99999) cos2= 0.99999


   ang1 = 2*acos( cos1 )
   ang2 = 2*acos( cos2 )


   ar1 = 0.5 * r1**2 *( ang1 - sin(ang1) )
   ar2 = 0.5 * r2**2 *( ang2 - sin(ang2) )

   ar  = ar1 + ar2

  end function OVRCIRCLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION X_ZGE( PCE ) RESULT(FQ)
  type (T_ONE_PCE), intent(IN) ::  PCE
  REAL, DIMENSION(PCE%L+1)     :: FQ
      FQ = PCE%ZGE
  end function X_ZGE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION X_W( PCE ) RESULT(W)
  type (T_ONE_PCE), intent(IN) ::  PCE
  REAL, DIMENSION(PCE%L)       :: W
      W = PCE%W
  end function X_W

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION X_PPE( PCE ) RESULT(FQ)
  type (T_ONE_PCE), intent(IN) ::  PCE
  REAL, DIMENSION(PCE%L+1)     :: FQ
      FQ = PCE%PPE
  end function X_PPE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_ZGL( PCE ) RESULT(FQ)
  type (T_ONE_PCE), intent(IN) ::  PCE
  REAL, DIMENSION(PCE%L)     :: FQ
  INTEGER :: LP
      LP=PCE%L
      FQ = ( PCE%ZGE(1:LP)+PCE%ZGE(2:LP+1) )/2.0
  end function C_ZGL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_PPL( PCE ) RESULT(FQ)
  type (T_ONE_PCE), intent(IN) ::  PCE
  REAL, DIMENSION(PCE%L)     :: FQ
  INTEGER :: LP
      LP=PCE%L
      FQ = ( PCE%PPE(1:LP)+PCE%PPE(2:LP+1) )/2.0
  end function C_PPL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  FUNCTION C_RADIUS( PCE ) RESULT(FQ)
  type (T_ONE_PCE), intent(IN) ::  PCE
  REAL, DIMENSION(PCE%L)       :: FQ
      where( PCE%A > 0.) 
        FQ = SQRT ( PCE%A / PI )
      elsewhere
        FQ = 1.
      endwhere
  end function C_RADIUS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_W( PCE ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)       :: F
  REAL, DIMENSION(PCE%L)       :: A
         F = PCE%W
  end function C_W

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_TEMPERATURE( PCE, RKAP ) RESULT(FQ)
  type (T_ONE_PCE), intent(IN) :: PCE
  !!! type (T_ONE_PCE), intent(IN) :: PCENV
  REAL,             intent(IN) :: RKAP
  REAL, DIMENSION(PCE%L)       :: FQ
  INTEGER :: LP
      LP = PCE%L
      IF ( PCE%TYPE_OF_CE .NE. EnvironmentCode ) THEN
          FQ = PCE%THETA+PCE%THBCK
      ELSE 
          FQ = PCE%THETA+PCE%THBCK
      ENDIF     
      FQ = FQ*( ( C_PPL(PCE) /1000.)**RKAP )
  end function C_TEMPERATURE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_BACK_TO_THETA( PCE, TE , RKAP ) RESULT(FQ)
  type (T_ONE_PCE), intent(INout) :: PCE
  !!! type (T_ONE_PCE), intent(IN) :: PCENV
  REAL,             intent(IN) :: RKAP,TE(PCE%L)
  REAL, DIMENSION(PCE%L)       :: FQ
  INTEGER :: LP
      LP = PCE%L
  
      FQ  = TE*( ( 1000. / C_PPL(PCE))**RKAP )
      IF ( PCE%TYPE_OF_CE .NE. EnvironmentCode ) THEN
          PCE%THETA = FQ - PCE%THBCK
      ELSE 
          PCE%THETA = FQ - PCE%THBCK
      ENDIF     
 
    end function C_BACK_TO_THETA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_Q( PCE ) RESULT(FQ)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)       :: FQ
      FQ = PCE%Q
  end function C_Q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_QI( PCE ) RESULT(FQ)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)       :: FQ
      FQ = PCE%QI
  end function C_QI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_QL( PCE ) RESULT(FQ)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)       :: FQ
      FQ = PCE%QL
  end function C_QL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_LAYER_MASS( PCE ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)       :: F
  INTEGER :: LP
      LP=PCE%L
      F = PCE%PPE(2:LP+1)-PCE%PPE(1:LP)
  end function C_LAYER_MASS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_RHL( PCE ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)       :: F
  INTEGER :: LP
      LP=PCE%L
      F = Cnv2Pascal*( PCE%PPE(2:LP+1)-PCE%PPE(1:LP)   )   &
        / ( PCE%ZGE(1:LP)  -PCE%ZGE(2:LP+1) ) / CE1_GRAV
  end function C_RHL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_LAYER_THICKNESS( PCE ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)       :: F
  INTEGER :: LP
      LP=PCE%L
      F = PCE%ZGE(1:LP)  -PCE%ZGE(2:LP+1)
  end function C_LAYER_THICKNESS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_TYPICAL_RADIUS( PCE ) RESULT(R)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL                         :: R,A

     A = C_MASS_XSCT2( PCE )
     R = SQRT( MAX(A,PI) /PI )

  end function C_TYPICAL_RADIUS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_TYPICAL_AREA( PCE ) RESULT(R)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL                         :: R,A

  real :: MASS(PCE%L),ARMN,MASSx

     ARMN = control%minallowedarea_ratio*PCE%A0 
     MASS = C_LAYER_MASS(PCE) 
     
     where( PCE%A <= 1.0001*ARMN )
          MASS = 0.00
     endwhere

     MASSx = SUM( MASS )
     if ( MASSx > 0. ) then 
         R = SUM( MASS*PCE%A )/ MASSx
     else
         R = ARMN
     endif   


  end function C_TYPICAL_AREA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION C_MASS_XSCT2( PCE ) RESULT(R)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL                         :: R

  real :: MASS(PCE%L) , MASK(PCE%L) 

    MASS = C_LAYER_MASS(PCE) 
    MASK = 0.
    where( PCE%A > (1.0e-5)*PCE%A0 )
        MASK=1.0
    endwhere

    IF (MAXVAL(MASK) > 0.) THEN
      R = SUM( MASS*MASK*PCE%A**2 ) / SUM(MASS*MASK )
    ELSE
      R = 1.E-6
    ENDIF

    R = SQRT( MAX(R,0.) )

  end function C_MASS_XSCT2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  FUNCTION C_EDGE_WS( PCE ) RESULT(WF)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L+1)     :: WF
  REAL, DIMENSION(PCE%L+1)     :: AF,DA
  INTEGER :: LP,L
  REAL    :: AREA_MAX,WMAX

     LP=PCE%L
     
     WF(2:LP) = ( PCE%W(2:LP)+PCE%W(1:LP-1) )*0.5

     WF(LP+1) = 0.
     WF(1)    = 0.
     WMAX     = 500.

      WHERE( WF < -WMAX )
           WF=-WMAX
      ENDWHERE
      WHERE( WF > WMAX )
           WF= WMAX
      ENDWHERE
       

  end function C_EDGE_WS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION C_THETA_D( PCE ) RESULT(F)
  type (T_ONE_PCE), intent(IN) :: PCE
  REAL, DIMENSION(PCE%L)       :: F
  INTEGER :: LP,L
  REAL    :: AREA_MAX,TMAX
     LP=PCE%L
      TMAX=1000.
      F = PCE%THETA
      WHERE( F < -TMAX )
           F=-TMAX
      ENDWHERE
      WHERE( F > TMAX )
           F= TMAX
      ENDWHERE

     
     
  end function C_THETA_D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! S O M E   S U B R O U T I N E S !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ADVANCE_TIME ( PCES, DT )

   type (T_ONE_PCE)  ,  intent(inout)     :: PCES(:)
   real              ,  intent(in)        :: DT
   real, dimension(6)  :: PCEDATE

   INTEGER   ::  IP , NP , DaysThisMonth

   DaysThisMonth=30

   NP = SIZE( PCES%IPC )

   DO IP=1,NP
      PCES(IP)%TIME = PCES(IP)%TIME + DT
      if (.NOT. IS_DORMANT(PCES(IP)) ) THEN
         PCES(IP)%AGE  = PCES(IP)%AGE  + DT
      endif

      PCEDATE = PCES(IP)%MODELDATE

      PCEDATE(6) = PCEDATE(6)+DT

      IF (PCEDATE(6) > 60.) THEN
         PCEDATE(5) = PCEDATE(5)+1.0
         PCEDATE(6) = PCEDATE(6)-60.
      endif
      IF (PCEDATE(5) > 60.) THEN
         PCEDATE(4) = PCEDATE(4)+1.0
         PCEDATE(5) = PCEDATE(5)-60.
      endif
      IF (PCEDATE(4) > 24.) THEN
         PCEDATE(3) = PCEDATE(3)+1.0
         PCEDATE(4) = PCEDATE(4)-24.
      endif
      IF (PCEDATE(3) > DaysThisMonth) THEN
         PCEDATE(2) = PCEDATE(2)+1.0
         PCEDATE(3) = 1.0
      endif

      PCES(IP)%MODELDATE = PCEDATE
      

   END DO


   if ( PCES(1)%TIME >= 5299. ) THEN
      write(*,*) "STOP POINT "
   end if


  END SUBROUTINE ADVANCE_TIME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DESTROY_CE1 ( PCE )   
                       
 
  
   type (T_ONE_PCE)  , intent(inout)    ::     PCE

   logical :: NULLFIRST = .FALSE.
   logical :: NULLAFTER = .FALSE.
   

           ! Deallocate space 
    DEALLOCATE(PCE%PPL      )
    DEALLOCATE(PCE%PPE      )
    DEALLOCATE(PCE%ZGE      )
    DEALLOCATE(PCE%UR       )
    DEALLOCATE(PCE%U        )
    DEALLOCATE(PCE%V        )
    DEALLOCATE(PCE%W        )
    DEALLOCATE(PCE%THETA    )
    DEALLOCATE(PCE%THBCK    )
    DEALLOCATE(PCE%Q        )
    DEALLOCATE(PCE%QL       )
    DEALLOCATE(PCE%QI       )
    DEALLOCATE(PCE%QS       )
    DEALLOCATE(PCE%QR       )
    DEALLOCATE(PCE%QH       )
    DEALLOCATE(PCE%NL       )
    DEALLOCATE(PCE%NI       )
    DEALLOCATE(PCE%NS       )
    DEALLOCATE(PCE%NR       )
    DEALLOCATE(PCE%NH       )
    DEALLOCATE(PCE%A        )
    DEALLOCATE(PCE%A_NM1    )
    DEALLOCATE(PCE%X        )
    DEALLOCATE(PCE%Y        )

  end SUBROUTINE DESTROY_CE1 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FREE_PCE ( PCE )   
                       
 
   type (T_ONE_PCE)  , intent(inout)    ::     PCE
   

  

! RESET CONFIGURATION
    PCE%CONFIG=""

! RESET SOME OF THE SCALAR COMPONENTS OF PCE

    PCE%IPC  = -999
    PCE%NPOP = 0

    PCE%Z0   = 0.00 
    PCE%PS   = 0.00 ! PS
    PCE%AGE  = 0.00
    PCE%A0   = 0.00
    PCE%AGRID= 0.00

    PCE%XMEAN  = MissingValue
    PCE%YMEAN  = MissingValue
 
    PCE%BIRTHTIME = -999.0
    PCE%DEATHTIME = -999.0

    PCE%MAXW00 = 0.00
    PCE%MAXWTF = 0.00
    PCE%STATUS = DormantValue
    PCE%TYPE_OF_CE  = LatentCode    

    PCE%RAIN_GAUGE= 0.0 
    PCE%HAIL_GAUGE= 0.0 
    PCE%SNOW_GAUGE= 0.0


    PCE%NEED_DNDRFT = FalseCode
    PCE%NEED_PRSHFT = FalseCode
    PCE%MY_DNDRFT   = MissingCode
    PCE%MY_PRSHFT   = MissingCode
    PCE%MY_ENV      = MissingCode

    
! RESET PROFILE COMPONENTS OF PCE

   ! location profile
    PCE%X    = 0. 
    PCE%Y    = 0.

   ! pressure grid (this will 
   ! be updated during run)
    PCE%PPL  = 0. ! 
    PCE%PPE  = 0. ! 

   ! zero height grid initially,
   ! replace after theta_env is 
   ! known 
    PCE%ZGE  = 0.

   ! motionless/empty ICs
    PCE%UR     = 0.
    PCE%U      = 0.
    PCE%V      = 0.
    PCE%W      = 0.
    PCE%A      = 0.
    PCE%THETA  = 0.
    PCE%THBCK  = 0.
    PCE%Q      = 0.
    PCE%QL     = 0.
    PCE%QI     = 0.
    PCE%QR     = 0.
    PCE%QS     = 0.
    PCE%QH     = 0.
    PCE%NL     = 0.
    PCE%NI     = 0.
    PCE%NR     = 0.
    PCE%NS     = 0.
    PCE%NH     = 0.

  END SUBROUTINE FREE_PCE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CREATE_CE1 ( L , K , PCE )   
                       
 
   INTEGER, INTENT(IN) ::  L , K
   type (T_ONE_PCE)  , intent(out)    ::     PCE
   

   logical :: write_pce_precshare = .true.
   logical :: write_pce_precfall  = .true.
   logical :: write_pce_shaftovr  = .true.
   logical :: write_pce_phyd_tot  = .true.
   logical :: write_pce_phyd_th   = .true.
   logical :: write_pce_phyd_qv   = .true.
   logical :: write_pce_phyd_qprc = .true.
   logical :: write_pce_phyd_qcld = .true.
   logical :: write_pce_xxbyc     = .true.
   logical :: write_pce_entr_th   = .true.
   logical :: write_pce_entr_qv   = .true.

    ALLOCATE(PCE%PPL   (L)     )
    ALLOCATE(PCE%PPE   (L+1)   )
    ALLOCATE(PCE%ZGE   (L+1)   )
    ALLOCATE(PCE%UR    (L)     )
    ALLOCATE(PCE%U     (L)     )
    ALLOCATE(PCE%V     (L)     )
    ALLOCATE(PCE%W     (L)     )
    ALLOCATE(PCE%THETA (L)     )
    ALLOCATE(PCE%THBCK (L)     )
    ALLOCATE(PCE%Q     (L)     )
    ALLOCATE(PCE%QL    (L)     )
    ALLOCATE(PCE%QI    (L)     )
    ALLOCATE(PCE%QS    (L)     )
    ALLOCATE(PCE%QR    (L)     )
    ALLOCATE(PCE%QH    (L)     )
    ALLOCATE(PCE%NL    (L)     )
    ALLOCATE(PCE%NI    (L)     )
    ALLOCATE(PCE%NS    (L)     )
    ALLOCATE(PCE%NR    (L)     )
    ALLOCATE(PCE%NH    (L)     )
    ALLOCATE(PCE%A     (L)     )
    ALLOCATE(PCE%A_NM1 (L)     )
    ALLOCATE(PCE%X     (L)     )
    ALLOCATE(PCE%Y     (L)     )
    
! INITIALIZE CONFIGURATION
    PCE%CONFIG=""

! INITIALIZE SOME OF THE SCALAR COMPONENTS OF PCE

    PCE%I    = -1
    PCE%J    = -1
    PCE%L    = L
    PCE%K    = K     
    PCE%IPC  = 0 
    PCE%NPOP = 0

    PCE%Z0   = 0.00 
    PCE%PS   = 0.00 ! PS
    PCE%AGE  = 0.00
    PCE%TIME = 0.00
    PCE%A0   = 0.00
    PCE%AGRID= 0.00

    PCE%XMEAN  = MissingValue
    PCE%YMEAN  = MissingValue
 
    PCE%BIRTHTIME = -999.0
    PCE%DEATHTIME = -999.0
    PCE%MODELDATE = (/ 0. , 0., 0. , 0. , 0., 0. /)-999.

    PCE%MAXW00 = 0.00
    PCE%MAXWTF = 0.00
    PCE%STATUS = DormantValue
    PCE%TYPE_OF_CE  = LatentCode    

    PCE%RAIN_GAUGE= 0.0 
    PCE%HAIL_GAUGE= 0.0 
    PCE%SNOW_GAUGE= 0.0
    PCE%ASH  = 0.00
    PCE%SHTOP= 0.00


    PCE%NEED_DNDRFT = FalseCode
    PCE%NEED_PRSHFT = FalseCode
    PCE%MY_DNDRFT   = MissingCode
    PCE%MY_PRSHFT   = MissingCode
    PCE%MY_ENV      = MissingCode
    
! INITIALIZE PROFILE COMPONENTS OF PCE

   ! location profile
    PCE%X    = 0. 
    PCE%Y    = 0.

   ! pressure grid (this will 
   ! be updated during run)
    PCE%PPL  = 0. ! 
    PCE%PPE  = 0. ! 

   ! zero height grid initially,
   ! replace after theta_env is 
   ! known 
    PCE%ZGE  = 0.

   ! motionless/empty ICs
    PCE%UR     = 0.
    PCE%U      = 0.
    PCE%V      = 0.
    PCE%W      = 0.
    PCE%A      = 0.
    PCE%THETA  = 0.
    PCE%THBCK  = 0.
    PCE%Q      = 0.
    PCE%QL     = 0.
    PCE%QI     = 0.
    PCE%QR     = 0.
    PCE%QS     = 0.
    PCE%QH     = 0.
    PCE%NL     = 0.
    PCE%NI     = 0.
    PCE%NR     = 0.
    PCE%NS     = 0.
    PCE%NH     = 0.

    
  END SUBROUTINE CREATE_CE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COPY_PCE ( PCE, PCETRG )   
                       
   type (T_ONE_PCE)  , intent(in   )    ::     PCE
   type (T_ONE_PCE)  , intent(inout)    ::     PCETRG
   
    PCETRG%PPL   = PCE%PPL
    PCETRG%PPE   = PCE%PPE
    PCETRG%ZGE   = PCE%ZGE
    PCETRG%UR    = PCE%UR 
    PCETRG%U     = PCE%U 
    PCETRG%V     = PCE%V 
    PCETRG%W     = PCE%W 
    PCETRG%THETA = PCE%THETA 
    PCETRG%THBCK = PCE%THBCK
    PCETRG%Q     = PCE%Q 
    PCETRG%QL    = PCE%QL 
    PCETRG%QI    = PCE%QI 
    PCETRG%QS    = PCE%QS 
    PCETRG%QR    = PCE%QR 
    PCETRG%QH    = PCE%QH 
    PCETRG%NL    = PCE%NL 
    PCETRG%NI    = PCE%NI 
    PCETRG%NS    = PCE%NS 
    PCETRG%NR    = PCE%NR 
    PCETRG%NH    = PCE%NH 
    PCETRG%A     = PCE%A
    PCETRG%A_NM1 = PCE%A_NM1 
    PCETRG%X     = PCE%X 
    PCETRG%Y     = PCE%Y 

    
!  CONFIGURATION
    PCETRG%CONFIG=PCE%CONFIG

! INITIALIZE SOME OF THE SCALAR COMPONENTS OF PCE

    PCETRG%I     = PCE%I    
    PCETRG%J     = PCE%J    
    PCETRG%L     = PCE%L   
    PCETRG%K     = PCE%K       
    PCETRG%IPC   = PCE%IPC  
    PCETRG%NPOP  = PCE%NPOP 

    PCETRG%Z0     = PCE%Z0    
    PCETRG%PS     = PCE%PS   ! PS
    PCETRG%AGE    = PCE%AGE  
    PCETRG%TIME   = PCE%TIME 
    PCETRG%A0     = PCE%A0   
    PCETRG%AGRID  = PCE%AGRID   

    PCETRG%MODELDATE = PCE%MODELDATE
    PCETRG%BIRTHTIME = PCE%BIRTHTIME 

    PCETRG%MAXW00    = PCE%MAXW00 
    PCETRG%MAXWTF    = PCE%MAXWTF  
    PCETRG%XMEAN     = PCE%XMEAN 
    PCETRG%YMEAN     = PCE%YMEAN  
    PCETRG%ASH       = PCE%ASH  
    PCETRG%SHTOP     = PCE%SHTOP  
        
    PCETRG%STATUS        = PCE%STATUS 
    PCETRG%RAIN_GAUGE    = PCE%RAIN_GAUGE 
    PCETRG%HAIL_GAUGE    = PCE%HAIL_GAUGE 
    PCETRG%SNOW_GAUGE    = PCE%SNOW_GAUGE 


    PCETRG%NEED_DNDRFT   = PCE%NEED_DNDRFT 
    PCETRG%NEED_PRSHFT   = PCE%NEED_PRSHFT 
    PCETRG%MY_DNDRFT     = PCE%MY_DNDRFT 
    PCETRG%MY_PRSHFT     = PCE%MY_PRSHFT
    PCETRG%MY_UPDRFT     = PCE%MY_UPDRFT
    PCETRG%MY_ENV        = PCE%MY_ENV 
    PCETRG%TYPE_OF_CE    = PCE%TYPE_OF_CE
 
  END SUBROUTINE COPY_PCE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MERG_PCE ( PCE, PCETRG )   
                       
   type (T_ONE_PCE)  , intent(in   )    ::     PCE
   type (T_ONE_PCE)  , intent(inout)    ::     PCETRG
   
    PCETRG%PPL   = PCE%PPL
    PCETRG%PPE   = PCE%PPE
    PCETRG%ZGE   = PCE%ZGE

    where( (PCE%A  + PCETRG%A) > 0. )
      PCETRG%UR    = ( PCE%A*PCE%UR     +  PCETRG%A*PCETRG%UR    )/( PCE%A  + PCETRG%A )
      PCETRG%U     = ( PCE%A*PCE%U      +  PCETRG%A*PCETRG%U     )/( PCE%A  + PCETRG%A )
      PCETRG%V     = ( PCE%A*PCE%V      +  PCETRG%A*PCETRG%V     )/( PCE%A  + PCETRG%A )
      PCETRG%W     = ( PCE%A*PCE%W      +  PCETRG%A*PCETRG%W     )/( PCE%A  + PCETRG%A )
      PCETRG%THETA = ( PCE%A*PCE%THETA  +  PCETRG%A*PCETRG%THETA )/( PCE%A  + PCETRG%A )
      PCETRG%THBCK = ( PCE%A*PCE%THBCK  +  PCETRG%A*PCETRG%THBCK )/( PCE%A  + PCETRG%A )
      PCETRG%Q     = ( PCE%A*PCE%Q      +  PCETRG%A*PCETRG%Q     )/( PCE%A  + PCETRG%A )
      PCETRG%QL    = ( PCE%A*PCE%QL     +  PCETRG%A*PCETRG%QL    )/( PCE%A  + PCETRG%A )
      PCETRG%QI    = ( PCE%A*PCE%QI     +  PCETRG%A*PCETRG%QI    )/( PCE%A  + PCETRG%A )
      PCETRG%QS    = ( PCE%A*PCE%QS     +  PCETRG%A*PCETRG%QS    )/( PCE%A  + PCETRG%A ) 
      PCETRG%QR    = ( PCE%A*PCE%QR     +  PCETRG%A*PCETRG%QR    )/( PCE%A  + PCETRG%A )
      PCETRG%QH    = ( PCE%A*PCE%QH     +  PCETRG%A*PCETRG%QH    )/( PCE%A  + PCETRG%A )
      PCETRG%X     = ( PCE%A*PCE%X      +  PCETRG%A*PCETRG%X     )/( PCE%A  + PCETRG%A )
      PCETRG%Y     = ( PCE%A*PCE%Y      +  PCETRG%A*PCETRG%Y     )/( PCE%A  + PCETRG%A )
   elsewhere
      PCETRG%UR    = 0.
      PCETRG%U     = 0. 
      PCETRG%V     = 0. 
      PCETRG%W     = 0. 
      PCETRG%THETA = 0. 
      PCETRG%THBCK = 0.
      PCETRG%Q     = 0. 
      PCETRG%QL    = 0. 
      PCETRG%QI    = 0. 
      PCETRG%QS    = 0. 
      PCETRG%QR    = 0. 
      PCETRG%QH    = 0. 
      PCETRG%X     = 0.5*( PCE%X  + PCETRG%X )  
      PCETRG%Y     = 0.5*( PCE%Y  + PCETRG%Y )
   endwhere


    PCETRG%A_NM1 = PCE%A_NM1 
   
    IF( IS_UPDRAFT(PCE) ) PCETRG%A     = PCE%A  + PCETRG%A  
    IF( IS_SHAFT(PCE) )   then 
       where( (PCE%QR + PCE%QS + PCE%QH) > 1.0e-5 )
         PCETRG%A     = PCE%A  + PCETRG%A  
       endwhere
    endif


!  CONFIGURATION
    PCETRG%CONFIG="merged"
    PCETRG%TYPE_OF_CE    = MergedUpdraftCode
    PCETRG%BIRTHTIME     = PCE%TIME 
    PCETRG%STATUS        = 1.00
    PCETRG%MAXW00        = 0.00
    PCETRG%MAXWTF        = 0.00  
    PCETRG%MODELDATE     = PCE%MODELDATE

! INITIALIZE SOME OF THE SCALAR COMPONENTS OF PCE


    PCETRG%XMEAN     = PCE%XMEAN 
    PCETRG%YMEAN     = PCE%YMEAN  
    
    PCETRG%RAIN_GAUGE    = 0.00
    PCETRG%HAIL_GAUGE    = 0.00 
    PCETRG%SNOW_GAUGE    = 0.00 


    PCETRG%NEED_DNDRFT   = FalseCode 
    PCETRG%NEED_PRSHFT   = FalseCode  
    PCETRG%MY_DNDRFT     = MissingCode
    PCETRG%MY_PRSHFT     = MissingCode 
    PCETRG%MY_UPDRFT     = MissingCode 
    PCETRG%MY_ENV        = MissingCode 
 

 
  END SUBROUTINE MERG_PCE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE INIT_ONE_UPD( PCE,PCENV )   
  !------------------------------------------------                     
   
   type (T_ONE_PCE)  , intent(inout)     ::     PCE
   type (T_ONE_PCE)  , intent(inout)     ::     PCENV
   real  :: RMODE 

    integer :: N,L,LP,KP,I,J,KK,NP

    real    :: AREA0, BDPTH0, W0, BCEN0, THP0, RAD0

    RMODE = control%INITBBL_RADMAX
  
        PCE%NPOP = 1
        PCE%I    = 1 ! these should depend sensibly on some BCKG grid
        PCE%J    = 1 ! these should depend sensibly on some BCKG grid
    
        PCE%TYPE_OF_CE = UpdraftCode
        PCE%STATUS     = -9999.9

        !Inherent from enviro

        PCE%THBCK  = PCENV%THBCK
        PCE%Q      = PCENV%Q

        PCE%QL     = PCENV%QL
        PCE%QI     = PCENV%QI
        PCE%QR     = PCENV%QR
        PCE%QS     = PCENV%QS
        PCE%QH     = PCENV%QH
        
        PCE%NL     = PCENV%NL
        PCE%NI     = PCENV%NI
        PCE%NR     = PCENV%NR
        PCE%NS     = PCENV%NS
        PCE%NH     = PCENV%NH

        PCE%W      = PCENV%W
        PCE%V      = PCENV%V
        PCE%U      = PCENV%U

        PCE%PPL    = PCENV%PPL
        PCE%PPE    = PCENV%PPE

        
       RAD0    =  RMODE
       AREA0   =  PI   * RAD0**2
       BDPTH0  =  control%INITBBL_DPTH !  00. * TANH( RAD0 / RMODE ) * 2
       BCEN0   =  control%INITBBL_ZCEN !! 0.4  * BDPTH0                   * 2
       THP0    =  control%INITBBL_THP0 !!  1.0  * TANH( RAD0 / RMODE ) * 4    !* 2
       W0      =  control%INITBBL_W0   !! 0.5  * TANH( RAD0 / RMODE ) * 4    !* 0

     
       call STARTBUBBLE    ( PCE , AREA0, THP0, W0, BCEN0 , BDPTH0  )
  

         PCE%X          = 0. 
         PCE%Y          = 0.
         PCE%XMEAN      = 0.
         PCE%YMEAN      = 0.
         PCE%BIRTHTIME  = 0.

         PCENV%A  = PCENV%A - PCE%A
         PCENV%A0 = PCE%A0 ! not really relevant but needed
                           ! for now by analysis routines
         

  END SUBROUTINE INIT_ONE_UPD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE INIT_ONE_ENVIRO( PCE )   
  !------------------------------------------------                     
   
   type (T_ONE_PCE)  , intent(inout)     ::     PCE
   real  :: RMODE 

    integer :: N,L,LP,KP,I,J,KK,NP

    real    :: AREA0, BDPTH0, W0, BCEN0, THP0, RAD0

  
        PCE%NPOP = 1
        PCE%I    = 1 ! these should depend sensibly on some BCKG grid
        PCE%J    = 1 ! these should depend sensibly on some BCKG grid

            ! Each environment is its own environment
        PCE%MY_ENV = PCE%IPC


        
        PCE%TYPE_OF_CE = EnvironmentCode
        PCE%STATUS     = -9999.9

       

         PCE%X          = 0. 
         PCE%Y          = 0.
         PCE%XMEAN      = 0.
         PCE%YMEAN      = 0.
         PCE%BIRTHTIME  = 0.

       PCE%THETA     = 0.0

       PCE%U     = 0.0 
       PCE%V     = 0.0
       PCE%W     = 0.0

       PCE%Q     = 0.0  
       PCE%QL    = 0.0   
       PCE%QI    = 0.0 
       PCE%QR    = 0.0
       PCE%QS    = 0.0
       PCE%QH    = 0.0

       PCE%NL    = 0.0
       PCE%NI    = 0.0
       PCE%NR    = 0.0
       PCE%NS    = 0.0
       PCE%NH    = 0.0


      END SUBROUTINE INIT_ONE_ENVIRO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SET_PCE_ZGRID_1    ( PCES, LP, ZEPCE )

    type (T_ONE_PCE)  ,  dimension(:) , intent(inout)     :: PCES
    integer, intent(in)                                   :: LP
    real, intent(in) , dimension(LP+1)                    :: ZEPCE




   INTEGER  ::   NP
   
   integer :: N,L,KP,I,J,IM,JM


   NP = SIZE( PCES%IPC )


   do N=1,NP
      PCES(N)%ZGE(:)  =  ZEPCE(:)
   end do


 end SUBROUTINE SET_PCE_ZGRID_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE LOCATE_PCENV(  PCE      , &
                            PCENV    , &
                            LATS, LONS, DXYS ,ICS, JCS, WGTS, NPS )

   type (T_ONE_PCE)   ,    intent(in)    ::     PCE
   type (T_ONE_PCE)   ,    intent(IN)    ::     PCENV

   REAL, DIMENSION (:,:),   INTENT(IN)    ::  LONS
   REAL, DIMENSION (:,:),   INTENT(IN)    ::  LATS
   REAL, DIMENSION (:,:),   INTENT(IN)    ::  DXYS

   INTEGER, DIMENSION(9), INTENT(OUT)     ::  ICS,JCS
   REAL   , DIMENSION(9), INTENT(OUT)     ::  WGTS

   INTEGER, INTENT(OUT) :: NPS

   NPS=1
   ICS=1
   JCS=1

   WGTS = PCENV%A0 / DXYS(1,1)

 
  END SUBROUTINE LOCATE_PCENV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE STARTBUBBLE    ( PCE , AR0, TD0, W0, Z0, ZSIG )

   type (T_ONE_PCE)  ,     intent(inout) ::  PCE
 
   REAL,                   INTENT(IN)    ::  AR0, TD0, Z0, ZSIG, W0
   REAL :: THD(PCE%L),ZGL(PCE%L),EE(PCE%L),EE2(PCE%L)
   INTEGER :: L,I

   IF ( PCE%AGE > 1.e-5 ) RETURN

   L = PCE%L
   ZGL = C_ZGL(PCE)

   EE  = EXP (  -( (ZGL-Z0)/ZSIG)**2 )
   EE2 = EE

   THD        = TD0*EE

   !!PCE%A      = AR0 !! <<bubble-in-Sleeve>>
   PCE%A      = MAX( AR0*EE2 , AR0*0.001 ) !! Bubble
   PCE%A0     = AR0

   PCE%W(:)   = W0*EE !0.0
   PCE%THETA  = THD
   !!PCE%Q      = PCENV%Q

  END SUBROUTINE STARTBUBBLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE CE1_UTILS
