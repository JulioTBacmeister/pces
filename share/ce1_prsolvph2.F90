#undef FVPRESSURE
! 
!-----------------------------
MODULE CE1_PRSOLVPH2


use CE1_DIAGS, only : PCEDIAGS
use CE1_TYPES
use CE1_UTILS
use CE1_INFORM
use CE1_CONSTS
!------------------------------------
! DESCRIPTION:
! This will solve for mean pressure in a "bell jar" with radius
! constant in height.  
!

 IMPLICIT NONE
 PRIVATE

 PUBLIC PSOLVPH

   character*72, parameter :: I_am_module = "ce1_prsolvph2" 

    ! These parameters are supposed to control imagined 
    ! geometry of bubble relevant to pressure sol'n 
   REAL :: MU_P = 1.1 ! 1.
   REAL :: MU_H = 1.1 ! 1. !1.
   REAL :: ABELL_FAC = 1.5 ! 1.0 ! 0.5 !0.5
   REAL :: SUBCLD_FAC = 1.0 


! Jun-Ichi says you should wind up with
!   dzz (P) - (2/Rc)**2 P  
!--------------------------------

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PSOLVPH( PCE, PCECLST, PP, BYC, PHYD, DT, IPC )

   type (T_ONE_PCE)   ,    intent(inout)                   ::     PCE
   !!type (T_BCKG)      ,    intent(IN   )                   ::     BCKG
   type (T_ONE_PCE)   ,    intent(inout)  ,  dimension(:)  ::     PCECLST


   REAL, DIMENSION(PCE%L)  , INTENT(out)   ::  BYC
   REAL, DIMENSION(PCE%L+1), INTENT(out)   ::  PP
   REAL, DIMENSION(PCE%L+1), INTENT(out)   ::  PHYD
   REAL,                     INTENT(IN)    ::  DT

   integer, INTENT(IN)                     ::  IPC

   REAL              ::  RHS(PCE%L+1), KE(PCE%L+1) , ABELL, rcount

   INTEGER           ::  LP, L

      LP = PCE%L
      ABELL  = 0.
      rcount = 0.

      ABELL = ABELL_FAC * MAXVAL( PCE%A )  

      CALL PRHS      ( PCE , PCECLST, ABELL, PHYD, RHS, BYC, KE )

      CALL ELLPSOLV  ( PCE , RHS , ABELL, PP , KE )


            PP=SUBCLD_FAC*PP ! "downscale" PP

      CALL ZMOM_CE1    (  PCE, BYC, PP , PHYD, ABELL, DT, IPC )


end  SUBROUTINE PSOLVPH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE PRHS    ( PCE , PCECLST, ABELL, PHYD, RHS, BYC, KE)
      !------------------------------------------------------------
      ! Calculates hydrostatic pressure $p_h$ and $\nabla^2_H p_h$ 
      ! for RHS of elliptic p-eq
      !------------------------------------------------------------ 

   type (T_ONE_PCE)  ,       intent(inout) ::  PCE
   type (T_ONE_PCE)   ,    intent(inout)  ,  dimension(:)  ::     PCECLST ! added 09/01/11
   !!! type (T_BCKG)     ,       intent(in)    ::  PCENV ! Note 6/7/06 PCENV=>BCKG ! Eliminate 8/4/11
   REAL, DIMENSION(PCE%L+1), INTENT(out)   ::  PHYD
   REAL, DIMENSION(PCE%L+1), INTENT(out)   ::  RHS,KE
   REAL, DIMENSION(PCE%L),   INTENT(out)   ::  BYC
   REAL,                     INTENT(in)    ::  ABELL


   INTEGER :: L,LP,NCL,CLB, VEB, ILOC(1),N,NSMTH,MYEV

   real :: RHO(PCE%L), THD(PCE%L),DZ(PCE%L),ZL(PCE%L),W2(PCE%L)
   real :: PPL(PCE%L),PJCT(PCE%L)
   real :: TH00,ATYPCL,BYCX(PCE%L),RHD(PCE%L),AFAC  !,KRADIUS2
   real :: QD(PCE%L), QLD(PCE%L),QID(PCE%L),RHSX(PCE%L+1)
   real :: QRD(PCE%L),QHD(PCE%L),QSD(PCE%L),DKE(PCE%L),DDKE(PCE%L+1)

   type (T_ONE_PCE)                             :: PCEX

   myEV=indx_my_enviro(PCE,PCECLST)

   
   LP  = PCE%L
   RHO = C_RHL(PCE)   !*0.0+1.
   PPL = C_PPL(PCE)   !*0.0+1.
   DZ  = C_LAYER_THICKNESS( PCE )   
   NCL = SIZE( PCECLST%IPC )
   ZL  = C_ZGL(PCE)

   NSMTH=10

   TH00=PCECLST(myev)%THBCK(LP)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! All of these routines to find 
        ! 'perturbations' for buoyancy
        ! are found ce1_utils.F90  
        ! 
        ! Of these only C_QL_D and C_QI_D
        ! explicitly use 'PCENV' 
        ! quantities.
        ! 
        ! PCE%THETA is already a 'perturbation'
        ! quantity. See advection in 
        ! ce1_dynamics_w_theta_ce1.                                     
        !               07/15/2011
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

veb=CONTROL % BUOY_VIRTUAL
clb=CONTROL % BUOY_CONDLOAD

   THD =  PCE%THETA !- PCEX %THETA
   QD  =  veb*(PCE%Q- PCECLST(myev)%Q)     !- PCEX %Q
   QLD =  clb*(PCE%QL-PCECLST(myev)%QL)    !- PCEX %QL
   QID =  clb*(PCE%QI-PCECLST(myev)%QI)    !- PCEX %QI
   QSD =  clb*(PCE%QS-PCECLST(myev)%QS)    !- PCEX %QS
   QRD =  clb*(PCE%QR-PCECLST(myev)%QR)    !- PCEX %QR  
   QHD =  clb*(PCE%QH-PCECLST(myev)%QH)    !- PCEX %QH

   !---------------------------------------------
   ! Relate peak theta pert to mean in 
   ! cloud disk of radius R:
   !     TH^ = THD/( 1.0 - 2/[n+2] )
   !         = THD*( [n+2]/n )
   ! where TH^ is peak perturbation, THD is cynlinder-mean,
   ! and n>0 is exponent in
   !     TH(r) = 1.0-(r/R)^n
   !--------------------------------------------
   !thd = thd*2.0


   RHD   = - THD / PCECLST(myev)%THBCK 

   DO L=1,LP                     !
      BYC(L)    = -1.0*GRAV *                 &
         ( RHD (L) - 0.61*QD(L) + ( QLD(L) + QID(L) + QRD(L) + QSD(L) + QHD(L) ) ) 
   END DO


   PHYD(1) = 0.0
   DO L=1,LP                     !
      PHYD(L+1)  = PHYD(L) - DZ(L) * BYC(L)  
   END DO
   !-------------------------------------------
   ! FYI: At this point PHYD at z=0 is actually
   ! negative if BYC is positive everywhere
   !-------------------------------------------

   DO L=2,LP                     !
      RHS(L) =     (2*PI*MU_H/ABELL) * PHYD(L)
   END DO


#if 1
   DO L=1,LP                     !
      AFAC   =  MIN ( PCE%A(L) / ABELL , 1.00  )
      KE(L)  = 0.5 * AFAC * ( PCE%W(L)**2 )
   END DO
   KE(LP+1) = KE(LP)
   KE(1:LP) = 0.5 * ( KE(1:LP) + KE(2:LP+1) )

   DDKE = 0.0

   DO L=1,LP                     !
      DKE(L) = ( KE(L)-KE(L+1) )/DZ(L)
   END DO
   DO L=2,LP                     !
      DDKE(L) = ( RHO(L-1)*DKE(L-1) - RHO(L)*DKE(L) ) / ( ZL(L-1) - ZL(L) )
   END DO

   DO N=1,nsmth
   RHSX = DDKE
      DO L=2,LP                     !
         DDKE(L) = 0.5 * RHSX(L) +  0.25 *RHSX(L-1) + 0.25 * RHSX(L+1)
      END DO
   END DO
#endif



#if 0
   DO L=1,LP                     !
      AFAC   =  MIN ( PCE%A(L) / ABELL , 1.00  )
      W2(L)  =  AFAC * PCE%W(L)
   END DO

   DO L=2,LP                     !
      DDKE(L) = ( W2(L-1)-W2(L) ) / ( ZL(L-1) - ZL(L) )
   END DO
   DO L=2,LP                     !
      DDKE(L) = 0.5 * ( RHO(L-1)+RHO(L) )*  ( DDKE(L)**2 )
   END DO
   DDKE(1)    = DDKE(2)
   DDKE(LP+1) = DDKE(LP)
#endif



   DO L=2,LP                     !
      !RHS(L) = RHS(L)  ! - DELTA_K * DDKE(L)  + (2*PI*MU_K/ABELL) * 0.5* (RHO(L)+RHO(L-1))*KE(L)
      RHS(L) = RHS(L)   ! - DDKE(L)  + (2*PI*MU_H/ABELL) * 0.5* (RHO(L)+RHO(L-1))*KE(L)
   END DO

   ! dP/dz = rho*beta at z=0,H set here
   ! rescaled at 0 and H in ELLPSOLV
   RHS(LP+1) = 0.  ! RHO(LP)* BYC(LP)
   RHS(1)    = 0.  ! RHO(1) * BYC(1)
   !RHS(LP+1) =  RHO(LP)* BYC(LP)
   !RHS(1)    =  RHO(1) * BYC(1)

    

  END SUBROUTINE PRHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ELLPSOLV    ( PCE , B , ABELL, P, KE )
!------------------------------------------------------------
! Sets up and solves 
!              dd/dzz[p'] -k^2 p = +k^2 PH
! Subject to
!              d/dz[p'] = rho*B at z=0,H
!------------------------------------------------------------
!
! Inputs:
!  PCE                  A derived type containing convective 
!                       element data (see ce1_types.F90)
!
!  B         real[L+1]  RHS of p-eq 
!                         - B(2:L)= (d_x,d_y)*( [d_x,d_y]* PH ) ~ k^2 PH
! Output:
!  P         real[L+1]  Nonhydrostatic pressure (Pa)
!                       
!----------------------------------------------------------------------------------
!
! Problem setup:
!
!           D(1)  DU(1)                            /P(1)  \        /B(1) = 0\
!           DL(1) D(2) DU(2)                      |        |      | B(2)     |
!                 .    .    .                     | .      |      |   .      |
!                      .   .   .               *  | .      |   =  |   .      |
!                         .   .   .               | .      |      |   .      |
!                       DL(L-1) D(L) DU(L)        |        |      | B(L)     |
!                              DL(L) D(L+1)        \P(L+1)/        \B(L+1)=0/ 
!-----------------------------------------------------------------------------------



   type (T_ONE_PCE)  ,       intent(in)      ::  PCE
   REAL, DIMENSION(PCE%L+1), INTENT(IN)      ::  B,KE
   REAL, DIMENSION(PCE%L+1), INTENT(OUT)     ::  P
   REAL, INTENT(IN)      ::  ABELL

   
   real :: BB(PCE%L+1)
   real :: DU(PCE%L), D(PCE%L+1),DL(PCE%L), DD(PCE%L+1)
   real :: ZGL(PCE%L), ZGE(PCE%L+1),RHO(PCE%L),RHZN(PCE%L+1)
   INTEGER :: L,I,INFO,L2,LP,LPP,K
   !!REAL :: KRADIUS2, ERADIUS , NORM
   REAL :: NORM
   
  

#ifdef SKIPPFRC
!   P=0.0 
    write(*,*) " PFORCE=0 PFORCE=0  PFORCE=0  PFORCE=0  PFORCE=0  PFORCE=0 "
    RETURN
#endif

   LP  = PCE%L
   LPP = PCE%L+1
   
   ZGL   = C_ZGL(PCE)
   ZGE   = X_ZGE(PCE)
   RHO = C_RHL(PCE)   !*0.0+1.
     

                  ! Anelastic needs (1/rho)*(drho/dz)*dP/dz term
    do k=2,LP 
       rhzn(k) = (2./ ( rho(k-1)+rho(k)) ) * ( rho(k-1)-rho(k) ) / (zgl(k-1)-zgl(k)  )
    end do
    k=0
    !!rhzn(k)=0.
    k=1
    rhzn(k)=(1./rho(k)) * ( rho(k)-rho(k+1) ) / (zgl(k)-zgl(k+1) )
    k=LP+1
    rhzn(k)=rhzn(k-1) ! (1./rho(k)) * ( rho(k-1)-rho(k) ) / (zgl(k-1)-zgl(k) )



   BB =  B
   ! dP/dz = rho*beta at z=0,H set in PRHS
   ! rescaled at 0 and H below
  

   DO L=2,LP
      DL(L-1)  =   1.0  / ( ( ZGE(L-1) - ZGE(L) )*( ZGL(L-1) - ZGL(L) ) )   &
                      + rhzn(L) /( ZGE(L+1)- ZGE(L-1) )
   END DO
   DO L=2,LP
      DU(L)   =    1.0  / ( ( ZGE(L) - ZGE(L+1) )*( ZGL(L-1) - ZGL(L) ) )  &
                      - rhzn(L) /( ZGE(L+1)- ZGE(L-1) )
   END DO

   DO L=2,LP
      D(L)   =           1.0  /  ( ZGE(L-1) - ZGE(L)   )

      D(L)   =   D(L) +  1.0  /  ( ZGE(L)   - ZGE(L+1) )

      D(L)   =   D(L) / ( ZGL(L-1) - ZGL(L) ) 

      D(L)   =   -1.0 * D(L) 
   END DO

   D(LP+1)  =  -1./( ZGE(LP)   - ZGE(LP+1) )
   DL(LP)   =  -D(LP+1)
   
   D(1)    =  1./( ZGE(1)   - ZGE(2) )
   DU(1)   = -D(1)  
   !DU(1)   =  0.   


   ! rescale top and bottom rows of matrix calc
     BB(LP+1) =  BB(LP+1)/ ( ZGL(LP-1) - ZGL(LP) ) 
     D(LP+1)  =  D(LP+1) / ( ZGL(LP-1) - ZGL(LP) ) 
     DL(LP)   =  DL(LP)  / ( ZGL(LP-1) - ZGL(LP) ) 

     BB(1)    =  BB(1) / ( ZGL(1) - ZGL(2) ) 
     D(1)     =  D(1)  / ( ZGL(1) - ZGL(2) ) 
     DU(1)   =   DU(1) / ( ZGL(1) - ZGL(2) ) 
   
              ! p=0 at H
     !BB(1)    =  0.  
     !D(1)     =  D(1)  / ( ZGL(1) - ZGL(2) ) 
     !DU(1)   =   0.
   


          ! Do a single pressure solution for entire column. Requires
          ! deciding on a single effective cloud radius
          ! KRADIUS2 now a module parameter calculated in PRHS
          !----------------------------------------------------------

      D(2:LP)  =  D(2:LP) - 2*PI*MU_P/ABELL

      call sgtsv(  LP+1 , 1, DL, D, DU, BB, LP+1, INFO )

      P = BB(1:LP+1)

      P = P ! - KE

  END SUBROUTINE ELLPSOLV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE ZMOM_CE1(  PCE, BYC, PP , PHYD, ABELL, DT , IPC )  

   type (T_ONE_PCE)  ,       intent(inout)  ::  PCE
   REAL, DIMENSION(PCE%L)  , INTENT(in)     ::  BYC
   REAL, DIMENSION(PCE%L+1), INTENT(in)     ::  PP,PHYD
   REAL,                     INTENT(in)     ::  DT, ABELL
   integer, INTENT(IN)                     ::  IPC

   
   REAL  :: P(PCE%L+1), BY(PCE%L), AX(PCE%L+1), DAZ(PCE%L), DAT(PCE%L)
   REAL  :: PPE(PCE%L+1) , ZGE(PCE%L+1),  PFZ(PCE%L)
   REAL  :: PPL(PCE%L)   , ZGL(PCE%L), TDAMP(PCE%L), RHL(PCE%L), RAW(PCE%L)

   REAL  :: KYY   

   INTEGER :: L,I,L2,LP
   LP=PCE%L

   KYY = 1.E3
   KYY = 1.E-3

   ZGL = C_ZGL(PCE)
   PPL = C_PPL(PCE)  
   ZGE = X_ZGE(PCE)
   PPE = X_PPE(PCE)   
   RHL  = C_RHL(PCE)
  
   ! Mapp belljar avg's back to cloud

   P  = 0.
   BY = 0.
   AX(2:LP) = 0.5 * (PCE%A(1:LP-1)+PCE%A(2:LP))
   AX(LP+1) = AX(LP)
   AX(1)    = AX(2)  

   P = PP
   where( AX > 0 )
      P  = PP ! * ( 1.0 + SQRT(ABELL/AX)) / 2.0
   end where

   where( PCE%A > 0. )
      BY = BYC * ABELL/PCE%A
   end where
    
   DAT = ( PCE%A - PCE%A_NM1) / DT

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Vertical gradient of pressure
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
   PFZ(1:LP) = ( P(2:LP+1)-P(1:LP) ) / ( ZGE(2:LP+1)-ZGE(1:LP) ) ! / RHL
       write(*,*) " ZMOM: nonhydrostatic pressure gradient force ======<><>V^V^V^<>><><> "
#else
   PFZ(1:LP) = -BYC(1:LP) 
       write(*,*) " ZMOM: Idealized Testing with buoyancy ONLY !!!!!B !! B !! B !! B !!!!!! ..... "          
#endif

   DAZ(1:LP) = 0.5* ( P(2:LP+1)+P(1:LP) ) * ( AX(2:LP+1)-AX(1:LP) ) & 
                                          / ( ZGE(2:LP+1)-ZGE(1:LP) )   !/ RHL
  where( PCE%A > 0. )
      DAZ = DAZ/PCE%A
  elsewhere
      DAZ = 0.
  end where

  PCEDIAGS%PFZ(:,IPC) = PFZ
  PCEDIAGS%PP(:,IPC)  = P



#ifndef FVPRESSURE
   PCE%W = PCE%W - PFZ*DT   !+ .05*DAZ*DT
               write(*,*) " DOING straight DW/Dt in prsolvph2 !!!!!!!! ********* "
#else
   where( PCE%A > control%zmom_minarea_ratio * PCE%A0 )
      RAW   = PCE%A * PCE%W
      RAW   = RAW   - PCE%A_NM1 * PFZ*DT  ! + (GAMMA_P-1.) * DAZ  !!- PCE%W*DT/TDAMP   !! + PCE%A_NM1 * BY*DT
      PCE%W = RAW / PCE%A
   elsewhere
      PCE%W = 0.
   end where   
               write(*,*) " DOING FV-zmom in prsolvph2  !!!!!!!! ********* "
#endif
 
   !if ( PCE%W(LP) < 0. ) PCE%W(LP)  = (1.- DT / 100. )*PCE%W(LP) ! If W is downward in lowest layer add extra damping 
   
         ! Save diagnostics
         !------------------------------------

 end subroutine ZMOM_CE1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end MODULE CE1_PRSOLVPH2
