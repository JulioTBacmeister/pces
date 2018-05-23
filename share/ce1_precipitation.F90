!  $Id: ce1_precipitation.F90,v 1.2 2008/02/12 23:22:42 bacmj Exp $
!------------------------------------------------------------------
! Routines in charge of disposing of precipitation that
! has been created elsewhere.  First stage (as of 02/11/08)
! is to hand some precip from updraft to associated "shaft".
!------------------------------------------------------------------


!===================================================
! CHANGE Log started 8/4/11 at NCAR
! 
!    08/04/11 - Eliminating BCKG structure.
!      ...
!    11/15/15 - Added super simple precipitation
!---------------------------------------------------

MODULE CE1_PRECIPITATION

use CE1_CONSTS
use CE1_TYPES
use CE1_INFORM
use CE1_UTILS
use CE1_DIAGS
use CE1_uPHYS
use NUMERICAL_UTILITIES, only : QSAT, DQSAT
use CE1_CONSTS, only: PI
use PPM, only : FXPPM,FXPPM2,UPST1,UPSTAD3

 IMPLICIT NONE
 PRIVATE


 PUBLIC PRECIPITATION

 character*72, parameter :: I_am_module = "ce1_precipitation" 


 CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PRECIPITATION (  PCECLST  , &
                              DT       )

   type (T_ONE_PCE) ,    intent(inout) , dimension(:)  ::     PCECLST
   !!!type (T_BCKG)    ,    intent(inout)                 ::     BCKG
   REAL,                 INTENT(IN   )                 ::     DT

   INTEGER, DIMENSION( SIZE( PCECLST%IPC ) ) :: IDXS

   INTEGER :: LP,I,L,N, MYUP, IMYUP,NCL 
 
   LP = PCECLST(1)%L 

   NCL = SIZE( PCECLST%IPC )

   
   do I=1,NCL
      IF(  IS_SHAFT( PCECLST(I) ) .and. ( .NOT. IS_DORMANT(PCECLST(I) ) )   ) THEN
           myup             = pceclst(i)%my_updrft
           idxs             = pceclst%ipc
           IMYUP            = -999
           DO N=1,NCL
              IF( IDXS(N) .EQ. MYUP ) IMYUP=N 
           END DO

           SELECT CASE (control%share_precip_v)
  
             CASE(0)
               call SUPER_SIMPLE_PRECIP( pceclst(imyup) , pceclst(i), DT,I )
             CASE(5)
                 write(*,*) " need to retink with enviros "
              ! call SHARE_PRECIP_5( pceclst(imyup) , pceclst(i), DT,I )

           END SELECT

      END IF


   end do


   end SUBROUTINE PRECIPITATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SUPER_SIMPLE_PRECIP (  PCUPDFT ,   &
                                    PCSHAFT ,   &
                                    DT ,II         ) 

   REAL,                  INTENT(IN   )    ::     DT
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCUPDFT
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCSHAFT
   integer           ,    intent(in   )    ::     II     

   real, dimension(pcupdft%l) :: X,Y,RAD,QR,A,WT,WU,DZ,MASS,RAD2,MASS2,HN,DHN,PSRC,QPCU,QPCS,QPCU0,DQPCU,ZO,DLA
   real :: DD,DV,AV,afctr,OFR,FALL,WREL,NEWASH,QVS,QVB
   real :: RP_U_S, TauShaftArea, cval,zmean_sh,delz_sh
   integer :: L,LP,LSHAFTBOT

   LP = PCUPDFT%L 

   
! Shaft is locked to X,Y of its Updraft
!-------------------------------------------------------------------------
#if 0
           PCSHAFT%X = PCUPDFT%X
           PCSHAFT%Y = PCUPDFT%Y
           PCSHAFT%XMEAN = PCUPDFT%XMEAN
           PCSHAFT%YMEAN = PCUPDFT%YMEAN
#endif

           PCSHAFT%ASH   = MAXVAL( PCSHAFT%A)

           write(*,*) "INside super-simple prec "

           MASS = C_LAYER_MASS(PCUPDFT)
           
     do L=1,LP-1               
              FALL = CONTROL%SIMPLE_SHARE_PREC_FRAC ! originally  0.1                  
              ! Mass cons addition of precip to shaft at L+1
              !----------------------------------------------
              ! Rain
              PCSHAFT%QR(L+1)    = PCSHAFT%QR(L+1) + FALL * PCUPDFT%QR(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L+1)*PCSHAFT%A(L+1))
              PCUPDFT%QR(L)      = (1.0-FALL) * PCUPDFT%QR(L)  ! Removal from upd at L
              PCSHAFT%NR(L+1)    = PCUPDFT%NR(L) 
              ! Snow
              PCSHAFT%QS(L+1)    = PCSHAFT%QS(L+1) + FALL * PCUPDFT%QS(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L+1)*PCSHAFT%A(L+1))
              PCUPDFT%QS(L)      = (1.0-FALL) * PCUPDFT%QS(L)  ! Removal from upd at L
              PCSHAFT%NS(L+1)    = PCUPDFT%NS(L) 
              ! Hail/Graupel
              PCSHAFT%QH(L+1)    = PCSHAFT%QH(L+1) + FALL * PCUPDFT%QH(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L+1)*PCSHAFT%A(L+1))
              PCUPDFT%QH(L)      = (1.0-FALL) * PCUPDFT%QH(L)  ! Removal from upd at L
              PCSHAFT%NH(L+1)    = PCUPDFT%NH(L) 
              ! X and Y
              !PCSHAFT%X(L+1)    = PCSHAFT%X(L+1) + FALL * PCUPDFT%X(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L+1)*PCSHAFT%A(L+1))
              !PCSHAFT%Y(L+1)    = PCSHAFT%Y(L+1) + FALL * PCUPDFT%Y(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L+1)*PCSHAFT%A(L+1))
              FALL              = FALL * MIN( PCUPDFT%A(L)/PCSHAFT%A(L+1) , 1.0 )
              PCSHAFT%X(L+1)    = (1.0-FALL)*PCSHAFT%X(L+1) + FALL * PCUPDFT%X(L)
              PCSHAFT%Y(L+1)    = (1.0-FALL)*PCSHAFT%Y(L+1) + FALL * PCUPDFT%Y(L) 
     end do

           

   end SUBROUTINE SUPER_SIMPLE_PRECIP

#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SHARE_PRECIP_5 (  PCUPDFT ,   &
                               PCSHAFT ,   &
                               DT ,II         ) 

   REAL,                  INTENT(IN   )    ::     DT
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCUPDFT
   type (T_ONE_PCE)  ,    intent(inout)    ::     PCSHAFT
   integer           ,    intent(in   )    ::     II     

   real, dimension(pcupdft%l) :: X,Y,RAD,QR,A,WT,WU,DZ,MASS,RAD2,MASS2,HN,DHN,PSRC,QPCU,QPCS,QPCU0,DQPCU,ZO,DLA
   real :: DD,DV,AV,afctr,OFR,FALL,WREL,NEWASH,QVS,QVB
   real :: RP_U_S, TauShaftArea, cval,zmean_sh,delz_sh
   integer :: L,LP,LSHAFTBOT


           afctr=PI
           LP  = PCUPDFT%L
           
           X   = PCUPDFT%X
           Y   = PCUPDFT%Y
           WU  = PCUPDFT%W
           QR  = PCUPDFT%QR
           QPCU= PCUPDFT%QR+PCUPDFT%QS+PCUPDFT%QH
           A   = PCUPDFT%A
           RAD = C_RADIUS(PCUPDFT)
           RAD2= C_RADIUS(PCSHAFT)
           MASS = C_LAYER_MASS(PCUPDFT)*100./grav

! Diagnose overhanging areas in updraft
!----------------------------------
           do l=1,lp-1
              DD    = SEPXY( (/ X(L+1),Y(L+1) /) , &     ! separation of centers of adj upd lyrs in the vertical 
                          (/ X(L),  Y(L)   /)  )             
              DHN(L) = OVRCIRCLE( RAD(L),RAD(L+1),DD )   ! area of overlap 
              DHN(L) = A(L)-DHN(L)                       ! subtract overlap from upper lyr of updraft (cant be <0)
           end do

           where( DHN < 0.)
              DHN = 0.
           endwhere

           DHN(lp) =  0.
           HN(1)   =  0.
          
           if(allocated(PCEDIAGS%DHN))  PCEDIAGS%DHN(:,II) =DHN       


               ! To make this correct maximal shaft as diagnosed by
               ! projecting from above, DHN at L should be minimum overhang
               ! w/resp to *every* layer below, i.e., L+1->LM
               ! (07/03/2014)
         call MAXIMAL_SHAFT( PCSHAFT ,PCUPDFT, DHN,DT,LP, LSHAFTBOT )
 

! Shaft is locked to X,Y of its Updraft
!-------------------------------------------------------------------------
           PCSHAFT%X = PCUPDFT%X
           PCSHAFT%Y = PCUPDFT%Y
           PCSHAFT%XMEAN = PCUPDFT%XMEAN
           PCSHAFT%YMEAN = PCUPDFT%YMEAN

           PCSHAFT%ASH   = MAXVAL( PCSHAFT%A)


             
                  ! Find top of shaft
              PCSHAFT%SHTOP=0.  
           
              do l=1,lp-1
                 if ( QPCU(L)>control%TINY_CONDENSATE_MR ) then
                      PCSHAFT%SHTOP =  PCSHAFT%ZGE(L)
                      exit
                 end if
              end do


end SUBROUTINE SHARE_PRECIP_5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MAXIMAL_SHAFT ( PCSHAFT , PCUPDFT, PCENV, DHN, DT, LP , LSHAFT_BOT )

   type (T_ONE_PCE)  ,    intent(inout)    ::     PCSHAFT,PCUPDFT,PCENV
   integer,               INTENT(IN)       ::     LP
   integer,               INTENT(  out)    ::     LSHAFT_BOT
   REAL,                  INTENT(IN)       ::     DHN(LP),DT
   integer :: i,j,l,n,m,l0,l1,lmx(1),ldz
   real :: zmean,delz,dzmx,zmx,tau
   real :: AA,FALL,A0,A1,DLA,SCALE_EXCH_LYR
   real :: mass(LP),qpcu(LP),ZO(LP),AX(LP),APRV(LP),DA(LP),AEN0(LP),GEOM1
   !! real, save :: QRX(LP)


              ! This is a brute-force way of determining
              ! how much of the "potential" shaft cross-
              ! section from overhanging updraft is actually
              ! used (0<x<1). Nothing prevents >1 but
              ! what would that mean? Used with "complexity"=2
              ! - 10/15/13 (Taipei)              
      GEOM1 = CONTROL% SIMPLE_SHARE_AREA_FCTR ! 1.0 ! 0.333
      SCALE_EXCH_LYR = 1.0

      write(*,*) "INside Maximal Shaft "


     APRV = PCSHAFT%A

           QPCU  =  PCUPDFT%QS + PCUPDFT%QR + PCUPDFT%QH
           MASS = C_LAYER_MASS(PCUPDFT)
           ZO   = 0.5 * ( PCUPDFT%ZGE(1:LP)  + PCUPDFT%ZGE(2:LP+1)  )

     TAU = DT

     PCSHAFT%A = MAX( PCSHAFT%A , control%minallowedarea_ratio*PCSHAFT%A0  )
     PCENV%A   = PCENV%A - (PCSHAFT%A-APRV)
     APRV = PCSHAFT%A
     
     lmx = MAXLOC (PCUPDFT%A ) !( dhn )
     l1=lmx(1)



     if (ALLOCATED(PCEDIAGS%SHEXTOP)) PCEDIAGS%SHEXTOP(:) = -999. 
     if (ALLOCATED(PCEDIAGS%SHEXBOT)) PCEDIAGS%SHEXBOT(:) = -999.

         ! New *target* Maximal shaft (not really the same 
         ! "maximal" determined by looking from above)
         ! 06/20/2014
     AX =0.
     do l=1,LP-1
        AX(L+1) = AX(L)+GEOM1*DHN(L)
     end do

     DA   = 0.

        ! Shaft expansion/contraction is an en/de-trainment process. Need to address conservation here
     do l=1,LP
        DA(L)=AX(L)-PCSHAFT%A(L)
     end do

     if (PCSHAFT%TIME < PCSHAFT%BIRTHTIME+Tau) then
        PCSHAFT%A = PCSHAFT%A + DA      
     else
        ! Growth (DA>0) is fast, but shrinking (DA<0) is slow
        where( DA > 0 )
           PCSHAFT%A  = PCSHAFT%A + DA
        elsewhere
            PCSHAFT%A  = PCSHAFT%A + (DT/5000.)*DA
        endwhere
     endif



     PCSHAFT%A = MAX( PCSHAFT%A , control%minallowedarea_ratio*PCSHAFT%A0  )

     AEN0=PCENV%A




!!===============================================================
!! NOTE!! Shaft can be shrinking as well as growing
!! So real en/de/trainment has to be done.
!!===============================================================
     PCSHAFT%THETA=PCSHAFT%THETA+PCSHAFT%THBCK
     do L=1,LP
        A0  = APRV(L)
        A1  = PCSHAFT%A(L)
        DLA = A1 - A0 
        PCENV%A(L)  = PCENV%A(L)  - DLA
        if (DLA >=0.) then
              ! EN-trainment calc
          A0  = A0/A1
          DLA = DLA/A1
          PCSHAFT%THETA(L) = PCSHAFT%THETA(L)*A0 + DLA*PCENV%THBCK(L)  
          PCSHAFT%U(L)  = PCSHAFT%U(L)*A0  + DLA*PCENV%U(L)  
          PCSHAFT%V(L)  = PCSHAFT%V(L)*A0  + DLA*PCENV%V(L)  
          PCSHAFT%Q(L)  = PCSHAFT%Q(L)*A0  + DLA*PCENV%Q(L)  
          PCSHAFT%QL(L) = PCSHAFT%QL(L)*A0 + DLA*PCENV%QL(L)  
          PCSHAFT%QI(L) = PCSHAFT%QI(L)*A0 + DLA*PCENV%QI(L)  
          PCSHAFT%QR(L) = PCSHAFT%QR(L)*A0 + DLA*PCENV%QR(L)  
          PCSHAFT%QS(L) = PCSHAFT%QS(L)*A0 + DLA*PCENV%QS(L)  
          PCSHAFT%QH(L) = PCSHAFT%QH(L)*A0 + DLA*PCENV%QH(L)  
              ! end of EN-trainment calc
        else  
              ! Now detrainment into BCKG
          A1  = PCENV%A(L)-DLA
          A0  = PCENV%A(L)/A1
          DLA = DLA/A1
          PCENV%THBCK(L) = PCENV%THBCK(L)*A0 - DLA*PCSHAFT%THETA(L)  
          PCENV%U(L)  = PCENV%U(L)*A0  - DLA*PCSHAFT%U(L)  
          PCENV%V(L)  = PCENV%V(L)*A0  - DLA*PCSHAFT%V(L)  
          PCENV%Q(L)  = PCENV%Q(L)*A0  - DLA*PCSHAFT%Q(L)  
          PCENV%QL(L) = PCENV%QL(L)*A0 - DLA*PCSHAFT%QL(L)  
          PCENV%QI(L) = PCENV%QI(L)*A0 - DLA*PCSHAFT%QI(L)  
          PCENV%QR(L) = PCENV%QR(L)*A0 - DLA*PCSHAFT%QR(L)  
          PCENV%QS(L) = PCENV%QS(L)*A0 - DLA*PCSHAFT%QS(L)  
          PCENV%QH(L) = PCENV%QH(L)*A0 - DLA*PCSHAFT%QH(L)  
        end if
     end do
     PCSHAFT%THETA=PCSHAFT%THETA-PCSHAFT%THBCK



     do L=1,LP-1
        
              if ( PCUPDFT%A(L)>=ModestFraction * PCUPDFT%A0 ) then
                 select CASE ( control%share_precip_complexity )
                 case(0) 
                     FALL = CONTROL%SIMPLE_SHARE_PREC_FRAC
                 case(1) 
                     FALL = CONTROL%SIMPLE_SHARE_PREC_FRAC * MIN( PCSHAFT%A(L+1)/ PCUPDFT%A(L), 1.00 )
                 case(2)  
                     FALL =  CONTROL%SIMPLE_SHARE_PREC_FRAC * GEOM1 *DHN(L)/ PCUPDFT%A(L)
                 end select
              else
                 FALL = 0.
              endif
                                                               ! Mass cons addition of precip to shaft at L+1
! Rain
              PCSHAFT%QR(L+1)    = PCSHAFT%QR(L+1) + FALL * PCUPDFT%QR(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L+1)*PCSHAFT%A(L+1))
              PCUPDFT%QR(L)      = (1.0-FALL) * PCUPDFT%QR(L)  ! Removal from upd at L
              PCSHAFT%NR(L+1)    = PCUPDFT%NR(L) 
! Snow
              PCSHAFT%QS(L+1)    = PCSHAFT%QS(L+1) + FALL * PCUPDFT%QS(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L+1)*PCSHAFT%A(L+1))
              PCUPDFT%QS(L)      = (1.0-FALL) * PCUPDFT%QS(L)  ! Removal from upd at L
              PCSHAFT%NS(L+1)    = PCUPDFT%NS(L) 
! Hail/Graupel
              PCSHAFT%QH(L+1)    = PCSHAFT%QH(L+1) + FALL * PCUPDFT%QH(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L+1)*PCSHAFT%A(L+1))
              PCUPDFT%QH(L)      = (1.0-FALL) * PCUPDFT%QH(L)  ! Removal from upd at L
              PCSHAFT%NH(L+1)    = PCUPDFT%NH(L) 
     end do

     
#if 1     
     !! clean up precip
     do L=1,LP 
        if ( PCUPDFT%A(L) < ModestFraction * PCUPDFT%A0 ) then
           if ( PCSHAFT%A(L) <  ModestFraction * PCUPDFT%A0 ) then
              PCENV%QR(L)= PCENV%QR(L) + PCUPDFT%A(L)*PCUPDFT%QR(L)/(PCENV%A(L))
              PCENV%QS(L)= PCENV%QS(L) + PCUPDFT%A(L)*PCUPDFT%QS(L)/(PCENV%A(L))
              PCENV%QH(L)= PCENV%QH(L) + PCUPDFT%A(L)*PCUPDFT%QH(L)/(PCENV%A(L))
           else
              PCSHAFT%QR(L)    = PCSHAFT%QR(L) + PCUPDFT%QR(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L)*PCSHAFT%A(L))
              PCSHAFT%QS(L)    = PCSHAFT%QS(L) + PCUPDFT%QS(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L)*PCSHAFT%A(L))
              PCSHAFT%QH(L)    = PCSHAFT%QH(L) + PCUPDFT%QH(L)*MASS(L)*PCUPDFT%A(L)/(MASS(L)*PCSHAFT%A(L))
              PCSHAFT%NR(L)    = PCSHAFT%NR(L) 
              PCSHAFT%NS(L)    = PCSHAFT%NS(L)
              PCSHAFT%NH(L)    = PCSHAFT%NH(L)
           endif
           PCUPDFT%QR(L) = 0.
           PCUPDFT%QS(L) = 0.
           PCUPDFT%QH(L) = 0.
           PCUPDFT%NR(L) = 0.
           PCUPDFT%NS(L) = 0.
           PCUPDFT%NH(L) = 0.
        endif
     end do
#endif


end SUBROUTINE MAXIMAL_SHAFT
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end MODULE CE1_PRECIPITATION
