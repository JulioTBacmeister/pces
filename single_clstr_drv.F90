!  Created 11/27/2015 to replace fake_gcm and PceGridComp
!  with simpler more transparent driver for a simple test 
!  configuration of PCEs.
!  Initial cut was to make run binarily identical w/ 
!  original "pces.x" driver.
!-----------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program MNPRG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use pce_comp_mod, only : pces, pce_init_nml , pce_init_pces_0, pce_run_sngl
  use NUMERICAL_UTILITIES, only : QSAT, DQSAT
  use GET_DRIVER_DATA, only : get_the_data


  use ce1_utils, only :  set_pce_zgrid_1, init_one_upd, init_one_enviro
  use ce1_gcm_cplr, only : set_pce_enviro
  use MODULE_MP_MORR_TWO_MOMENT, only : morr_two_moment_init
  use ce1_diags, only : init_diags


IMPLICIT NONE

integer, parameter :: IM=10
integer, parameter :: JM=10
integer, parameter :: LM=72
integer, parameter :: NWAT=1
integer, parameter :: indxqv=1

real, parameter    :: ZTOP=30000.

logical  :: INITIALIZE
real :: DT
integer :: NSTEPS,I,L,J,UNIT

 
real, allocatable, dimension(:) :: pGCM_U
real, allocatable, dimension(:) :: pGCM_V
real, allocatable, dimension(:) :: pGCM_T
real, allocatable, dimension(:) :: pGCM_TH
real, allocatable, dimension(:) :: pGCM_Q
real, allocatable, dimension(:) :: pGCM_ZE
real, allocatable, dimension(:) :: pGCM_PE
real, allocatable, dimension(:) :: pGCM_PO
real, allocatable, dimension(:) :: pGCM_ZO
 
character(len=12)               :: DATAFILE = "arm_jjjj"    
real, dimension(6)              :: ModelTime = (/ 1998. , 5.0 , 25. , 0.0 , 0.0 , 0.0 /)
real, dimension(0:LM)           :: Prefs
real                            :: PCE_DT=10.
real                            :: DIAG_DT=50.
integer                         :: LenSeg =180
integer                         :: NSEGS=4
integer                         :: NUPDS=2
real                            :: DOMX = 1.e6
real                            :: DOMY = 1.e6
real                            :: SHRDPTH = -9999.
real                            :: USHR = 0.0
real                            :: VSHR = 0.0
real                            :: U0   = 0.0
real                            :: V0   = 0.0
logical                         :: uniform_strat=.false.

real :: strat0,bckgarea
integer :: LP,NP,N

real, allocatable, dimension(:) :: ZEPCE

namelist /cntrls/ ModelTime,DATAFILE,PCE_DT,LenSeg,NSEGS,NUPDS,DOMX,DOMY,SHRDPTH,USHR,VSHR,U0,V0,uniform_strat

allocate( pGCM_PE(1:LM+1) )
allocate( pGCM_ZE(1:LM+1) )
allocate( pGCM_ZO(1:LM) )
allocate( pGCM_Q(1:LM) )
allocate( pGCM_V(1:LM) )
allocate( pGCM_U(1:LM) )
allocate( pGCM_T(1:LM) )
allocate( pGCM_TH(1:LM) )
allocate( pGCM_PO(1:LM) )

UNIT=222
OPEN( UNIT=UNIT, FILE="control.nml" ) !, NML =  cntrls )
READ( UNIT=UNIT, NML=cntrls)

CLOSE(UNIT=UNIT)


LP=100
NP=3
        ! 73-level GEOS-5 edge pressures (MERRA vintage)
        !-------------------------------------------------
prefs = & 
(/       1.00293  ,    2.00293,      3.27293,      4.76143,      6.60293,      8.93743,      11.9732,      15.9524,      21.1378,  &
        27.8555   ,   36.5070 ,     47.5835 ,     61.6808 ,     79.5163 ,     101.947 ,     130.054 ,     165.082 ,     208.500 ,  & 
       262.024    ,  327.646  ,    407.660  ,    504.683  ,    621.683  ,    761.987  ,    929.297  ,    1127.69  ,    1364.34  ,  &
      1645.71     , 1979.16   ,   2373.04   ,   2836.78   ,   3381.00   ,   4017.54   ,   4764.39   ,   5638.79   ,   6660.34   ,  &
      7851.23     , 9236.58   ,   10866.3   ,   12783.7   ,   15039.3   ,   17693.0   ,   20825.3   ,   24528.2   ,   28900.9   ,  &
      33915.1     , 37727.1   ,   41541.2   ,   45358.7   ,   49177.9   ,   52999.1   ,   56821.1   ,   60645.0   ,   64468.9   ,  &
      68294.1     , 70844.6   ,   73395.1   ,   75945.8   ,   78497.0   ,   81048.3   ,   83089.3   ,   84620.2   ,   86151.2   ,  &
      87682.2     , 89213.2   ,   90744.3   ,   92275.4   ,   93806.6   ,   95337.9   ,   96869.0   ,   98400.3   ,   99931.7   ,  &
      101454. /)


         ! Prefs into GET_THE_DATA are supplied in Pa.  Within GET_THE_DATA 
         ! Prefs is indexed as (0:LM).  Pressures coming out of GET_THE_DATA 
         ! are in hPa as desired by PCES and are indexed  PLE(:,:,1:LM+1) and
         ! PLO(:,:,1:LM).  Subroutine GET_THE_DATA lives in source file
         ! "get_driver_data.F90" 
         !-------------------------------------------------------------------

call GET_THE_DATA(  DATA$=TRIM(DATAFILE) , ModelTime=ModelTime  ,  Prefs=Prefs &
              , U   = pGCM_U                          &
              , V   = pGCM_V                          &
              , T   = pGCM_T                          &
              , TH  = pGCM_TH                         &
              , Q   = pGCM_Q                          &
              , PLE = pGCM_PE                         &
              , PLO = pGCM_PO                         &
              , ZLE = pGCM_ZE                         &
              , ZLO = pGCM_ZO                         &
                      ) 


where( pGCM_ZE(1:LM) < SHRDPTH )
pGCM_U=U0 + USHR * pGCM_ZE(1:LM)/SHRDPTH
pGCM_V=V0 + VSHR * pGCM_ZE(1:LM)/SHRDPTH
elsewhere
pGCM_U=U0 + USHR
pGCM_V=V0 + VSHR
endwhere


if (uniform_strat) then
#if 0
 strat0 = 9.825e-5
 pGCM_ZO(1:LM) = 0.5 *( pGCM_ZE(2:LM+1)+pGCM_ZE(1:LM) )
 pGCM_TH(LM)=300. - ( strat0 * 300.  /9.81 ) * (150. - pGCM_ZO(LM))
      write(*,*) LM,pGCM_TH(LM),pGCM_ZO(LM)
 do L=LM,2,-1
    pGCM_TH(L-1) = ( strat0 * pGCM_TH(L)  /9.81 ) * (pGCM_ZO(L-1)- pGCM_ZO(L)) + pGCM_TH(L)

      write(*,*) L-1,pGCM_TH(L-1),pGCM_ZO(L-1)

 end do
#endif

 pGCM_ZO(1:LM) = 0.5 *( pGCM_ZE(2:LM+1)+pGCM_ZE(1:LM) )
 do L=LM,1,-1
    !pGCM_TH(L) = 300.* EXP(  pGCM_ZO(L) /( (7./2.)*7000.) )
    pGCM_TH(L) = 300.* EXP(  pGCM_ZO(L) / 100000. )
 end do


endif



#if 0
pGCM_U=5.+15.*(pGCM_ZE(1:LM)/15000.)
where ( pGCM_U > 20. )
      pGCM_U=20.
endwhere
#endif


!Impose a dry layer between 2000 and 6000 m altitude. Added long ago 
!Noted - 4/3/12
where ( (pGCM_ZE(1:LM) > 2000.).and.(pGCM_ZE(1:LM) < 10000.))
      pGCM_q = 0.
      !pGCM_q = 0.333*pGCM_q
endwhere


BCKGAREA =  ( DOMY/(1.*JM) ) * ( DOMX/(1.*IM) )


!************************************************************
! Data for run has been read in/created. Now set-up PCES
!*************************************************************

! Read in remaining namelist variables
!-------------------------------------
call pce_init_nml


! Create and allocate space for PCES structure
!---------------------------------------------
call pce_init_pces_0( LP, NP )

! Create and allocate PCEDIAGS structure
!---------------------------------------------
call INIT_DIAGS(LP,NP)


! Give each pce a unique id number
!---------------------------------
   do N=1,NP
     PCES(N)%IPC = N + 199
   end do


! Initialize model date fields
!-----------------------------
   do n=1,np 
      pces( n )%MODELDATE = (/ 1.00, 1.00, 1.00, 0.0,0.0,0.0 /)
   end do


! Create and give pce's their vertical grids
!-------------------------------------------
   allocate( ZEPCE(1:LP+1) )

   ZEPCE(LP+1) = 0.0 ! surface height
   DO L=LP,1,-1 
      ZEPCE(L) = ZEPCE(L+1) + ZTOP/LP
   END DO
   call set_pce_zgrid_1( PCES, LP, ZEPCE )   


! Set backgrounds fields for pce's 
!---------------------------------
!call SET_PCE_BCKG ( LP, LM, NP, PCES &
!     , pGCM_U, pGCM_V, pGCM_TH , pGCM_Q &
!     , pGCM_PE, pGCM_ZE, BCKGAREA )

!
! 
!---------------------------------------------
write(*,*) "at init one enviro"
call init_one_enviro( PCES(1) )
pces(1)%status=1.0

write(*,*) "at set pce enviro"
call SET_PCE_ENVIRO ( LP, LM, PCES(1) &
     , pGCM_U, pGCM_V, pGCM_TH , pGCM_Q &
     , pGCM_PE, pGCM_ZE, BCKGAREA )

! Create a bubble with parameters in namelist.
!---------------------------------------------
call init_one_upd( PCES(2) , PCES(1) )
pces(2)%birthtime=0.
pces(2)%my_env = pces(1)%ipc

  ! kluge needed for diagnostics package
pces(1)%a0 = pces(2)%a0

! Initialize reams of physical constants etc. for Morrison microphysics 
!----------------------------------------------------------------------
call MORR_TWO_MOMENT_INIT



! Run pces
!---------
do i=1,NSEGS
call pce_run_sngl(PCE_DT, LenSeg )
end do

end program MNPRG
