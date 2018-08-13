!  $Id: ce1.F90,v 1.35 2008/02/21 22:14:55 bacmj Exp $
!-----------------------------
!!!!#define SKIPINT
!!!!#define SKIPUPH
!!!!#define SKIPADJENV
!!!!#define SKIPPREC

!======================================================
! Log started 8/4/11 at NCAR
! 
!    08/04/11 - Eliminating BCKG structure from all 
!               processes except ADJUST_BCKG. Changing
!               BCKG/PCE structure to reflect likelihood
!               of sparse PCES in a model.
!
!    06/07/18 - Github version
!-------------------------------------------------------

MODULE CE1

use CE1_CONSTS, only: grav
use CE1_TYPES
use CE1_DIAGS
use CE1_INFORM
use CE1_UTILS
use CE1_UPHYS
use CE1_DYNAMICS
use CE1_VDIFF
use CE1_ADJUSTMENTS
use CE1_PRECIPITATION
use CE1_INTERACTIONS,   only: INTERACTIONS
use CE1_CREATE_DESTROY
!use CE1_GCM_CPLR,       only: REFRESH_PCEBCKG

 IMPLICIT NONE
 PRIVATE

 PUBLIC RUN_CE1
 PUBLIC FINALIZE_CE1

 character*72, parameter :: I_am_module = "ce1" 

CONTAINS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE RUN_CE1    (    PCES       , &
                             DT         , & 
                             CP         , & 
                             GRAV       , & 
                             RKAP       , & 
                             ALHL       , & 
                             ALHS          )




   !!! type (T_ONE_PCE)  ,    pointer    ::  PCES(:)
   type (T_ONE_PCE)  ,    intent(inout)    ::  PCES(:)

   REAL, INTENT( IN    )  :: DT,CP,GRAV,RKAP,ALHL,ALHS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !! type (T_ONE_PCE) ,    pointer      ::  PCE
       !! type (T_BCKG) ,       pointer      ::  BCKG
   type (T_ONE_PCE)          ::  PCE
       ! Structure to record incoming state of PCES 05/03/12
   type (T_ONE_PCE), allocatable, dimension(:)            ::  PCES0


   REAL    :: ADUM, QTOTAL1, QTOTAL2, QTOTAL3,DUM1 , wlp,masskids
   REAL    :: ourmass,tmass,tmass0,tmass_ini,tmass_fin
   integer :: NSUB,LP,LUNW,IP,L,IC,i,KP

   integer :: NP,N,NXC,NYC,evin

   integer :: SEVERITY=1

   logical :: report=.true.

   np = size( PCES%IPC)
   LP = PCES(1)%L
   KP = 1
       ! Structure to record incoming state of PCES 05/03/12

   LUNW=111

      call PHYSINI(LP,NP)

      do n=1,np
         PCES(n)%A_NM1 = PCES(n)%A
         !call refresh_pcebckg(PCES(n),BCKG)  ! in ce1_gcm_coupler
      end do

       call CHECK_MASS(  PCES , NP ,  LP  ,  " top o' time step" , tmass0 )

      tmass_ini = tmass0


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Adjust PCE environment/gridbox, e.g.,
      ! through subsidence. Also calculate 
      ! movement of PCEs through grid-scale
      ! Done first to ensure strict constancy
      ! of TOTAL( PCE%A ) -01/01/16
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(control%do_subsidence) then
      CALL SUBSIDENCE_XYMOTION(  PCES     , &           ! in ce1_adjustments
                                 DT        )
             call CHECK_CEs (       PCES , NP ,  SEVERITY  ,  " after subsidence " )
             call CHECK_WATER(  PCES , NP ,  LP  ,  " after subsidence" )
      else
      if(REPORT) write(*,*) "     skipped subsidence ... "
         do n=1,np
         if(.not.IS_DORMANT(pces(n)).and..not.IS_ENVIRONMENT(pces(n)) ) then
            evin = indx_my_enviro( PCES(n) , PCES )
            PCES(n)%THBCK = PCES(evin)%THBCK
         end if
         end do
      endif



       call CHECK_MASS(  PCES , NP ,  LP  ,  " before dynamics" , tmass0 )
     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Do the 1-D "vertical" dynamics on the 
      ! clouds (and shafts). This includes
      ! advection of area, and non-precip
      ! species. Also, calcuate vertical 
      ! pressure forces. Only operate on
      ! elements with IS_ACTIVE_4_DYNMX()=.true..
      ! This function is in ce1_inform.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(control%do_dynamics) then
      CALL DYNAMICS   (          PCES        , &          ! in ce1_dynamics
                                 DT    )
             call CHECK_CEs (       PCES , NP ,  SEVERITY  ,  " after dynamics " )
             call CHECK_WATER(  PCES , NP ,  LP  ,  " after dynamics" )
      else
      if(REPORT) write(*,*) "skipped dynamics ... "
      endif

      call CHECK_MASS(  PCES , NP ,  LP  ,  " after dynamics" , tmass )

      write(*,*) " mass change across dynamics ",tmass - tmass0


      if(control%do_vdiff) then
      CALL VERTICAL_DIFFUSION (  PCES        , &          ! in ce1_dynamics
                                 DT , LP    )
             call CHECK_CEs (       PCES , NP ,  SEVERITY  ,  " after vertical diffusion " )
             call CHECK_WATER(  PCES , NP ,  LP  ,  " after vertical diffusion" )
      else
      if(REPORT) write(*,*) " skipped v diffusion ... "
      endif

      call CHECK_MASS(  PCES , NP ,  LP  ,  " after vdiff " , tmass )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Do actual microphysics - condensation
      ! re-evaporation, autoconversion etc..
      ! Operates on all PCEs where 
      ! IS_ACTIVE_4_UPHYS(PCE)=.true.. Look 
      ! in ce1_inform for codes.  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(control%do_uphys) then
      CALL UPHYS   (             PCES     , &           ! in ce1_uphys
                                 DT        )
             call CHECK_CEs (       PCES , NP ,  SEVERITY  ,  " after uphys " )
             call CHECK_WATER(  PCES , NP ,  LP  ,  " after uphys " )
      else
      if(REPORT) write(*,*) "  skipped uphys ... "
      endif
      call CHECK_MASS(  PCES , NP ,  LP  ,  " after uphys " , tmass0 )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate interactions within the group of PCEs.
      ! Includes entrainment and "mutual attraction". This
      ! code is a bear. You'd better have a good look at it.
      ! (07/29/11) 
      ! Works on all non dormant PCEs
      ! IS_DORMANT(PCE)=.false. 
      !----------------
      ! Rethinking entrainment (10/09/12). Nested entrainment
      ! was not really working anyway.  Probably more impt to 
      ! entraiment from grid scale better. Restored BCKG to
      ! args to allow this to be done
      ! ---------------
      ! Rethinking entrainment yet again (9/26/13) (annual cycle?).
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(control%do_interactions) then
      CALL INTERACTIONS   (      PCES     , &            ! in ce1_interactions
                                 DT        )
             call CHECK_CEs (       PCES , NP ,  SEVERITY  ,  " after interactions " )
             call CHECK_WATER(  PCES , NP ,  LP  ,  " after interactions " )
      else
      if(REPORT) write(*,*) "   skipped interactions ... "
      endif
      call CHECK_MASS(  PCES , NP ,  LP  ,  " after interactions " , tmass )
      write(*,*) " mass change across interactions ",tmass - tmass0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This code calculates precipitation 
      ! MACROphysics, i.e., shaft dynamics,
      ! vertical movement of precipitating
      ! condensate q_{r,s,h}, transfer of
      ! q_{r,s,h} from updrafts to shafts.  
      ! Loops over shafts not updrafts.    
      ! (08/04/11)            
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(control%do_precipitation) then
      CALL PRECIPITATION   (     PCES     , &            ! in ce1_precipitation
                                 DT        )
             call CHECK_CEs (       PCES , NP ,  SEVERITY  ,  " after precipitation " )
             call CHECK_WATER(  PCES , NP ,  LP  ,  " after precipitation " )
      else
      if(REPORT) write(*,*) "    skipped precipittation ... "
      endif
      call CHECK_MASS(  PCES , NP ,  LP  ,  " after precipitation " , tmass )


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Turn PCEs on and off etc.. Not very well
      ! developed as of 10/01/13
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      call MANAGE_PCES( PCES , DT )                   ! in ce1_create_destroy
      write(*,*) " after manage pces"

      ! Manually sync up THBCK with enviro -12/14/15
      ! This needs to move somewhere else.
      !---------------------------------------------
  
      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Advance everyones time.  Here DT should be the model << interval
  ! of consciousness >> , i.e., the shortest meaningful timestep
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL ADVANCE_TIME( PCES , DT )                       ! in ce1_utils
  !!CALL SIMPLE_WRITE_BCKG( 111 , PCES,  BCKG , DT , NP )          ! in ce1_diags


  CALL SIMPLE_WRITE( 112 , PCES, DT , NP )          ! in ce1_diags

  call CHECK_MASS(  PCES , NP ,  LP  ,  " end of step " , tmass_fin )

       if( abs(tmass_fin - tmass_ini)/tmass_fin > 1.e-12 ) then
           write(*,*) tmass_fin,tmass_ini , tmass_fin - tmass_ini
           write(*,*) abs(tmass_fin - tmass_ini)/tmass_fin
           write(*,*) " Mass violation "
           STOP
       endif

  call PHYSDESTROY

  END SUBROUTINE RUN_CE1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CHECK_MASS (     PCES           , &
                               NP             , &
                               LP    , text  , TMASS )

   INTEGER         ,    INTENT(IN)       :: NP,LP

   type (T_ONE_PCE),    intent(in)        ::  PCES(NP)
   character(LEN=*),    intent(in )       :: text
   real , intent(out) :: tmass

   real :: MASS(LP),QVT,QVB
   integer :: L,N,I

   TMASS = 0.
   do N=1,NP
      MASS   = C_LAYER_MASS(PCES(N))*100./grav
      TMASS  = SUM( MASS * PCES(n)%A , 1) + TMASS
   end do
#if 0
   write(161,*) " Mass at:"//trim(text)
   write(161,*) "TMASS =", TMASS
#endif

end SUBROUTINE CHECK_MASS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CHECK_WATER (     PCES           , &
                               NP             , &
                               LP    , text   )

   INTEGER         ,    INTENT(IN)       :: NP,LP

   type (T_ONE_PCE),    intent(in)        ::  PCES(NP)
   character(LEN=*),    intent(in )       :: text

   real :: MASS(LP),QVT,QVB
   integer :: L,N,I

   QVT = 0.
   write(*,*) " Water at:"//trim(text)
   do N=1,NP
      MASS = C_LAYER_MASS(PCES(N))*100./grav
      QVT  = SUM( MASS * PCES(n)%A * PCES(n)%Q , 1) + QVT
      !!write(*,*) N, QVT
   end do
   QVB  = SUM( MASS * PCES(1)%A * PCES(1)%Q , 1)

   write(*,*) "TPW =", QVB+QVT


 end SUBROUTINE CHECK_WATER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CHECK_CEs (       PCES           , &
                               NP             , &
                               SEVERITY       , &
                               text             )

   INTEGER         ,    INTENT(IN)       :: NP,SEVERITY
   character(LEN=*),    intent(in )      :: text

   type (T_ONE_PCE),    intent(in)        ::  PCES(NP)

   type (T_ONE_PCE)      ::  PCE
   integer :: N

     if ( severity < 0 ) RETURN

     do n = 1,NP 
         if ( MINVAL( PCES(N)%A) < 0 ) then 
            write(*,*) trim(text)//":: negative Area in PCE" ,N
            write(*,*) PCES(N)%A
            if ( SEVERITY >= 1 )  STOP
         end if
     end do

     end subroutine CHECK_CEs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE FINALIZE_CE1 (    PCES           , &
                               RESTART_FILE   , &
                               NP                 )

   INTEGER         ,    INTENT(IN)       :: NP
   character(LEN=*),    intent(in )      :: RESTART_FILE

   type (T_ONE_PCE),    intent(in)        ::  PCES(NP)

   type (T_ONE_PCE)      ::  PCE
   integer :: lun , ip

   lun=100

!#include  "finalize.code"
      
   END SUBROUTINE FINALIZE_CE1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PHYSINI ( LP, NP )
!=================================
! DECLARE a PHYS coupler
!=================================
! type (T_PCE_PHYS) , save  ::  PCEPHYS
integer, intent(IN) :: LP,NP

 allocate( PCEPHYS%RF(LP,NP) )
 allocate( PCEPHYS%SF(LP,NP) )
 allocate( PCEPHYS%HF(LP,NP) )
 allocate( PCEPHYS%UR(LP,NP) )
 !allocate( PCEPHYS%WDZA(LP,NP) )

 if (allocated(PCEPHYS%RF))  PCEPHYS%RF   = -999.
 if (allocated(PCEPHYS%SF))  PCEPHYS%SF   = -999.
 if (allocated(PCEPHYS%HF))  PCEPHYS%HF   = -999.
 if (allocated(PCEPHYS%UR))  PCEPHYS%UR   = -999999.
 if (allocated(PCEPHYS%WDZA))  PCEPHYS%WDZA = 0.0
 

end SUBROUTINE PHYSINI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PHYSDESTROY
!=================================
! DECLARE a PHYS coupler
!=================================
! type (T_PCE_PHYS) , save  ::  PCEPHYS

 if (allocated(PCEPHYS%RF))   deallocate( PCEPHYS%RF )
 if (allocated(PCEPHYS%SF))   deallocate( PCEPHYS%SF )
 if (allocated(PCEPHYS%HF))   deallocate( PCEPHYS%HF )
 if (allocated(PCEPHYS%UR))   deallocate( PCEPHYS%UR )
 if (allocated(PCEPHYS%WDZA)) deallocate( PCEPHYS%WDZA )

end SUBROUTINE PHYSDESTROY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




END MODULE CE1
