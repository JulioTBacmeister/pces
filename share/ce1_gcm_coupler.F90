!  $Id: ce1_gcm_coupler.F90,v 1.14 2008/07/01 21:08:03 bacmj Exp $
!
!-----------------------------
MODULE CE1_GCM_CPLR

use CE1_TYPES
use CE1_UTILS
use CE1_INFORM
use NUMERICAL_UTILITIES, only : QSAT


 IMPLICIT NONE
 PRIVATE

 !public set_pce_bckg
 public set_pce_enviro
 !PUBLIC REFRESH_PCEBCKG


      REAL   , PARAMETER       :: GRAV   = 9.81
      real, parameter          :: CP = 1003.5
      REAL   , PARAMETER       :: RKAP   = 0.286
      REAL   , PARAMETER      :: ALHL   =  2.4548E6
      REAL   , PARAMETER      :: ALHS   =  2.8368E6

 character*72, parameter :: I_am_module = "ce1_gcm_coupler" 

CONTAINS

 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vertical mapping of AGCM profiles to PCE z-grid.
!
! History: 
!    12/15/2015 - Branched from set_pce_bckg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SET_PCE_ENVIRO ( LP, LM, PCE , U, V, TH , Q,  PGE, ZGE, BCKGAREA )

 integer, intent(in) :: LP,LM
  type (T_ONE_PCE)  ,  intent(inout)     :: PCE
 real, intent(in), dimension(LM)   :: U, V, TH , Q  ! QL, QI,
 real, intent(in), dimension(LM+1) :: PGE, ZGE 
 real, intent(in) :: BCKGAREA

 ! local
 real, dimension( LM )  :: ZOGCM , RHOGCM
 real, dimension( LP ) :: ZOPCE , RHOPCE,TE,QST,QST2,TE2
 integer :: L,I,J,IW,JW,N


    PCE%A         = BCKGAREA
    PCE%AGRID     = BCKGAREA
    ZOPCE         = 0.5 * ( PCE%ZGE(1:LP)  + PCE%ZGE(2:LP+1)  )
    ZOGCM         = 0.5 * ( ZGE(1:LM) + ZGE(2:LM+1) )
    PCE%U         = ZINTRP( U(:) , ZOGCM , ZOPCE )
    PCE%V         = ZINTRP( V(:) , ZOGCM , ZOPCE )
    PCE%THBCK     = ZINTRP( TH(:) , ZOGCM , ZOPCE )
    PCE%THETA     = 0. 
    PCE%Q         = ZINTRP( Q(:) , ZOGCM , ZOPCE )
    PCE%QL        = 0. 
    PCE%QI        = 0. 
    PCE%QR        = 0. 
    PCE%QS        = 0. 
    PCE%QH        = 0. 
    PCE%W         = 0. 

    RHOGCM          = ( PGE(2:LM+1) - PGE(1:LM) )/( ZGE(2:LM+1) - ZGE(1:LM) )
    RHOPCE          = ZINTRP( RHOGCM , ZOGCM , ZOPCE )

    PCE%PPE(LP+1)  =  PGE(LM+1) !PE(I,J,LM+1)
    DO L=LP,1,-1 
       PCE%PPE(L)  =  PCE%PPE(L+1) - RHOPCE(L)*(PCE%ZGE(L+1)-PCE%ZGE(L))
    end do
    PCE%PPL(1:LP)  =  (   PCE%PPE(1:LP)  +  PCE%PPE(2:LP+1) )/2.
 
end SUBROUTINE SET_PCE_ENVIRO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE CE1_GCM_CPLR
