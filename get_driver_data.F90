!  $Id: get_driver_data.F90,v 1.3 2008/07/03 20:37:11 bacmj Exp $
!
!-----------------------------
MODULE GET_DRIVER_DATA

use CE1_CONSTS


 IMPLICIT NONE
 PRIVATE

PUBLIC GET_THE_DATA



 character*72, parameter :: I_am_module = "get_driver_data" 

CONTAINS

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE GET_THE_DATA(  DATA$ , ModelTime  ,  Prefs , Q, T, TH, U, V, PLE, PLO, ZLE, ZLO )      

!=======================================================================
    character(len=*)   , intent(in)                :: DATA$
    real, dimension(6),  intent(in)                :: ModelTime
    real, dimension(0:), intent(in)                :: Prefs 

    !!real, dimension(:,:,:), intent(out)            :: Q,T,TH,U,V,PLE,PLO,ZLE,ZLO
    real, dimension(:), intent(out)            :: Q,T,TH,U,V,PLE,PLO,ZLE,ZLO

!!!    real,                intent(in)                :: DT 
 
  ! temporary garbage dump for profile data

  ! Prefs come in in Pa. 


      real, ALLOCATABLE, SAVE, DIMENSION(:       ) :: time        !Calenday day
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::  yy        !Year
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::  mo        !Month
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::  dd        !Day
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::  hh        !Hour
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::  mm        !Minutes
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: LHF
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: SHF 
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: PSFC
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: PCP_OBS 
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: TS_AIR
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: TG_SOIL
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: TSKIN
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: QSKIN
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: QSFCAIR
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: tt
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: qq
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: uu
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: vv
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: T_H_adv 
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: T_V_adv
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q_H_adv
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q_V_adv
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q1 
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q2 
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: PLE_DATA
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: OMEGA


      INTEGER :: NT, NLEVEL,I,J,VERTADV,II,L

      LOGICAL :: DATA_DRIVER, AT_START
      LOGICAL, SAVE :: ALREADY_HAVE_DATA
      LOGICAL, allocatable, dimension(:), SAVE :: INITIALIZED_PCES
      integer, save :: I_time_step

      integer :: N_PCES,LM,STATUS
      real    :: DT_PCES,FAC0, FAC1

      real    :: PKE( size(PLE,1) )
      !!real    :: PKE( size(PLE,1) , size(PLE,2) , size(PLE,3) )



      LM = size( PREFS )-1



     select CASE(trim(DATA$))
       case("b500rpsw" )
         NT = 124
	 NLEVEL = 72
	 DATA_DRIVER=.true.
       case("CPT1")
         NT = 1440
	 NLEVEL = 26
	 DATA_DRIVER=.true.
       case("CPT2")
         NT = 1440
	 NLEVEL = 26
	 DATA_DRIVER=.true.
       case("CPT3")
         NT = 1440
	 NLEVEL = 26
	 DATA_DRIVER=.true.
       case("BOMEX")
         NT = 120
	 NLEVEL = 75
	 DATA_DRIVER=.true.
       case("SCSMEX_SESA")
         NT = 244
	 NLEVEL = 40
	 DATA_DRIVER=.true.
       case("SCSMEX_NESA")
         NT = 244 
	 NLEVEL = 40
	 DATA_DRIVER=.true.
       case("TOGA_COARE")
         NT = 85
	 NLEVEL = 38
	 DATA_DRIVER=.true.
       case("arm_95jul")
         NT = 141
	 NLEVEL = 18
	 DATA_DRIVER=.true.
       case("arm_97jul")
         NT = 233
	 NLEVEL = 18
	 DATA_DRIVER=.true.
       case("arm_97sep")
         NT = 167
	 NLEVEL = 18
	 DATA_DRIVER=.true.
       case("arm_98may")
         NT = 166
	 NLEVEL = 18
	 DATA_DRIVER=.true.
       case("arm_99jan")
         NT = 166
	 NLEVEL = 18
	 DATA_DRIVER=.true.
       case("arm_99jan_m5K")
         NT = 166
	 NLEVEL = 18
	 DATA_DRIVER=.true.
       case("arm_99mar")
         NT = 166
	 NLEVEL = 18
	 DATA_DRIVER=.true.
       case("arm_99jul")
         NT = 78
	 NLEVEL = 18
	 DATA_DRIVER=.true.
       case("arm_00mar")
         NT = 166
	 NLEVEL = 35
	 DATA_DRIVER=.true.
       case("arm_kwjx")
         NT = 210
	 NLEVEL = 37
	 DATA_DRIVER=.true.
       case("arm_scmx")
         NT = 182
	 NLEVEL = 37
	 DATA_DRIVER=.true.
       case("ARM95")
         NT=141
         DATA_DRIVER=.true.
       case("GATE")
         NT=160
         DATA_DRIVER=.true.
       case("FUNNY")
         NT=9
         DATA_DRIVER=.true.
       case("RADCONVEQ")
         NT=1
         DATA_DRIVER=.false.
       case("STRATUS")
         NT=1
         DATA_DRIVER=.false.
       case("NEURAL")
         NT=1
         DATA_DRIVER=.false.
       case("FORECAST")
         NT=1
         DATA_DRIVER=.false.
       case default
         NT=1
         DATA_DRIVER=.false.
     end select




! Allocate arrays for driver data
!----------------------------------------
 ! Here key off of whether TIME is allocated to determine
 ! whether data has been alloc'd and read in or not. Data arrays
 ! have "SAVE" attribute, so are allocated once and not deallocated.
 ! Kludgey but seems to work - JTB 7/21/04
if (.not.(ALLOCATED(TIME)) ) THEN 
      allocate(time (NT) )       
      allocate(yy (NT) )        !Year
      allocate(mo (NT) )        !Month
      allocate(dd (NT) )        !Day
      allocate(hh (NT) )        !Hour
      allocate(mm (NT) )        !Minutes
      ALLOCATE(PSFC (NT) )
      ALLOCATE(PCP_OBS (NT) ) 
      ALLOCATE(TS_AIR (NT) )
      ALLOCATE(TG_SOIL (NT) )
      ALLOCATE(TSKIN (NT) )
      ALLOCATE(QSKIN (NT) )
      ALLOCATE(QSFCAIR(NT) )
      ALLOCATE(LHF(NT) )
      ALLOCATE(SHF(NT) )
      ALLOCATE(tt (NT , LM ) )
      ALLOCATE(qq (NT , LM ) )
      ALLOCATE(uu (NT , LM ) )
      ALLOCATE(vv (NT , LM ) )
      ALLOCATE(T_H_adv (NT , LM ) ) 
      ALLOCATE(T_V_adv (NT , LM ) )
      ALLOCATE(Q_H_adv (NT , LM ) )
      ALLOCATE(Q_V_adv (NT , LM ) )
      ALLOCATE(Q1 (NT , LM ) ) 
      ALLOCATE(Q2 (NT , LM ) ) 
      ALLOCATE(OMEGA (NT , 0:LM ) )
      ALLOCATE(PLE_DATA (NT , 0:LM ) )

     
      I_time_step=0
      ALREADY_HAVE_DATA=.FALSE.
      
endif

select CASE(trim(DATA$))

     case("arm_00mar","arm_95jul","arm_97jul","arm_97sep","arm_98may",  &
          "arm_99jan","arm_99jan_m5K","arm_99jul","arm_99mar","arm_kwjx","arm_scmx",    & 
          "TOGA_COARE","SCSMEX_NESA","SCSMEX_SESA" )

          CALL ARM2(trim(DATA$)//'.dat', NT, NLEVEL,     &
                     LM ,                                &
                     PREFs,                            &
                     time,                               &
                     yy,                                 &
                     mo,                                 &
                     dd,                                 &
                     hh,                                 &
                     mm,                                 &
                     PCP_OBS,                            &
                     TS_AIR,                             &
                     TG_SOIL,                            &
                     TSKIN,                              &
                     QSFCAIR,                            &
                     QSKIN,                              &
                     PSFC,                               &
		     LHF,                                &
		     SHF,                                &
                     tt,                                 &
                     qq,                                 &
                     uu,                                 &
                     vv,                                 &
                     T_H_adv,                            &
                     T_V_adv,                            &
                     Q_H_adv,                            &
                     Q_V_adv,                            &
                     Q1,                                 &
                     Q2,                                 & 
                     PLE_DATA                             ) 

                 OMEGA = 0.


        case DEFAULT
         write(*,*) " No good data "
         STATUS=1

       
      end SELECT
      ALREADY_HAVE_DATA=.TRUE.





            ii=WhereInTime ( yy,mo,dd,hh,mm, nt, MODELtime, Fac0, Fac1 , RC=STATUS)



      do l=1,lm+1
         PLE(L) = Prefs(L-1) / 100.0
      end do
      do l=1,lm
         PLO(L) = ( Prefs(L-1)+Prefs(L) ) / 2.0 / 100.0
      end do

      do l=1,lm
        Q(l)   =( Fac0*QQ(ii,l) + Fac1*QQ(ii+1,l) )/1000.
        T(l)   =( Fac0*TT(ii,l) + Fac1*TT(ii+1,l) )
        U(l)   =( Fac0*UU(ii,l) + Fac1*UU(ii+1,l) )
        V(l)   =( Fac0*VV(ii,l) + Fac1*VV(ii+1,l) )
      end do

      TH     = T    * ( ( 1000. / PLO )**CE1_RKAP )
      PKE    = ( ( PLE / 1000. )**CE1_RKAP )

      ZLE(LM+1) = 0.0

      do l=LM+1,2,-1
         ZLE(L-1) = ZLE(L) + CE1_CP * TH(L-1)*(PKE(L)-PKE(L-1))/CE1_GRAV
      end do


end SUBROUTINE GET_THE_DATA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function WhereInTime ( yy,mo,dd,hh,mm, nt,SCMtime,Fac0,Fac1, RC ) result(ii)
integer,                intent(in) :: nt
real,    dimension(nt), intent(in) :: yy,mo,dd,hh,mm
real,    dimension(6),  intent(in) :: SCMtime
!!!real,                   intent(in) :: DT
real,                   intent(out):: Fac0,Fac1
integer, optional,      intent(out):: RC
integer                            :: ii

integer :: i

real    ::  SCM_SecOfDay,Da1_SecOfDay,Da0_SecOfDay,hh_i,hh_i1
real    ::  SCM_SecN,Da1_SecN,Da0_SecN


ii=-1
Fac0=0.
Fac1=0.
RC=0

SCM_SecN = SecOfYear( SCMtime(1), SCMtime(2), SCMtime(3), SCMtime(4), SCMtime(5), SCMtime(6) )

if ( nt == 1) then
   i=1
   Da0_SecN = SecOfYear( yy(i),   mo(i),   dd(i),   hh(i)  , mm(i)  , 0. ) 
   if (ABS(Da0_SecN-SCM_SecN) < 60.) then
      ii=1
   else
      ii=-99
   endif
   RETURN
end if

do i=1,nt-1

   if ( yy(i) .ne. SCMtime(1) ) cycle

   Da0_SecN = SecOfYear( yy(i),   mo(i),   dd(i),   hh(i)  , mm(i)  , 0. ) 
   Da1_SecN = SecOfYear( yy(i+1), mo(i+1), dd(i+1), hh(i+1), mm(i+1), 0. ) 
       

   if ( ( Da0_SecN .gt. SCM_SecN ) .or. ( Da1_SecN .le. SCM_SecN ) ) cycle

         !! now we should be in the correct Data interval
         !! Lets stop here
       ii=i

  
       Fac0 = ( Da1_SecN - SCM_SecN ) / ( Da1_SecN - Da0_SecN ) 
       Fac1 = 1.0-Fac0
  
end do

IF (II.lt.1) THEN
  RC=-1
  print *, ' BAD Time in SCM ',SCMtime
         print *," This data starts on "
         print *, yy(1),mo(1),dd(1),hh(1),mm(1)
         print *," This data ends on "
         print *, yy(nt),mo(nt),dd(nt),hh(nt),mm(nt)
         print *, Da0_SecN, SCM_SecN, Da1_SecN

         !do i=1,nt
         !print *, yy(i),mo(i),dd(i),hh(i),mm(i)
         !enddo 

ENDIF

end function WhereInTime 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function SecOfYear(yr,mo,dy,hh,mm,ss) result(SecN)
real, intent(in) :: yr,mo,dy,hh,mm,ss
real             :: SecN 
integer, dimension(12)      :: monlen
integer, dimension(0:12)    :: DaysSoFar
integer :: i, DayN, iyr, imo

imo=INT(mo)
iyr=INT(yr)

if (MOD(iyr,4).ne.0) then
   !!        J   F   M   A   M   J   J   A   S   O   N   D
   monlen=(/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
ELSE
   monlen=(/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
endif

DaysSoFar(0)=0
do i=1,12
   DaysSoFar(i)=DaysSoFar(i-1)+monlen(i)
end do

DayN=DaysSoFar(imo-1)+dy

SecN = DayN*24.*3600. + hh*3600. + mm*60. + ss

end function SecOfYear


end MODULE GET_DRIVER_DATA
