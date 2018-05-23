!  $Id: ce1_prsolv01.F90,v 1.4 2007/02/02 19:21:47 bacmj Exp $
!-----------------------------
MODULE CE1_PRSOLV01


use CE1_TYPES
use CE1_UTILS
use CE1_CONSTS
!------------------------------------
! !DESCRIPTION:
!
!  Approximate form of diagnostic Poisson equation
!  for pressure is solved.  This code uses traditional
!  anelastic form e.g Bannon, Rottunno in which a 
!  pressure $p'$ perturbation from a horizontally-independent
!  mean is found. This $p'$ may be hydrostatic or not. Buoyancy 
!  related RHS term is $\partial_z \beta$
!  

 IMPLICIT NONE
 PRIVATE

 PUBLIC PSOLV

   character*72, parameter :: I_am_module = "ce1_prsolv01" 


CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PSOLV( PCE, BCKG, PP, BYC, PHYD )

   type (T_ONE_PCE)   ,    intent(in   )                   ::     PCE
   type (T_BCKG)      ,    intent(IN   )                   ::     BCKG

   REAL, DIMENSION(PCE%L)  , INTENT(out)   ::  BYC
   REAL, DIMENSION(PCE%L+1), INTENT(out)   ::  PP
   REAL, DIMENSION(PCE%L+1), INTENT(out)   ::  PHYD

   REAL              ::  BYCZ( PCE%L+1 )

        write(*,*) " new psolv, old decomposition "
      CALL PRHS      ( PCE , BCKG , BYC, BYCZ )

      CALL ELLPSOLV  ( PCE , BYCZ , PP )


      PHYD = 0.

  end  SUBROUTINE PSOLV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE PRHS    ( PCE , PCENV , BYC, BYCZ  )
      !-----------------------------------------------------
      ! Calculates buoyancy and d/dz [rho*buoyancy] RHS
      ! for elliptic p-eq
      !----------------------------------------------------- 

   type (T_ONE_PCE)  ,       intent(in)    ::  PCE
   type (T_BCKG)     ,       intent(in)    ::  PCENV ! Note 6/7/06 PCENV=>BCKG
   REAL, DIMENSION(PCE%L)  , INTENT(out)   ::  BYC
   REAL, DIMENSION(PCE%L+1), INTENT(out)   ::  BYCZ


   INTEGER :: L,LP

   real :: RHO(PCE%L), THD(PCE%L), ZGE(PCE%L+1),DZ(PCE%L)
   real :: ZGL(PCE%L),DZL(PCE%L-1),PPL(PCE%L),PJCT(PCE%L)
   real :: TH00,ATYPCL,BYCX(PCE%L) 
   real :: QD(PCE%L), QLD(PCE%L),QID(PCE%L)
   real :: QRD(PCE%L),QHD(PCE%L),QSD(PCE%L)

   LP  = PCE%L
   ZGL = C_ZGL(PCE)
   ZGE = X_ZGE(PCE)
   RHO = C_RHL(PCE)   !*0.0+1.
   PPL = C_PPL(PCE)   !*0.0+1.

   TH00=PCE%THBCK(LP)

   THD =  C_THETA_D(PCE,PCENV)  
   QD  =  C_Q_D(PCE,PCENV)  
   QLD =  C_QL_D(PCE,PCENV)  
   QID =  C_QI_D(PCE,PCENV)  
   QSD =  C_QS_D(PCE,PCENV)  
   QRD =  C_QR_D(PCE,PCENV)  
   QHD =  C_QH_D(PCE,PCENV)  

   BYC  = GRAV * ( THD / TH00 + 0.61*QD - QLD - QID - QRD - QSD - QHD ) 

   ATYPCL = PCE%A0

   PJCT     = 1.00   

   BYCX    =  PJCT * BYC

   BYCZ(2:LP) = ( RHO( 2:LP)*BYCX( 2:LP ) - RHO( 1:LP-1)*BYCX( 1:LP-1) ) & 
               /(            ZGL( 2:LP ) -              ZGL( 1:LP-1) )
  


         ! Set up BCs for p-eq
         !     d/dz[p'] = rho*B at z=0,H
         !----------------------------------------
   BYCZ(1)    = RHO(1)* BYC(1)   
   BYCZ(LP+1) = RHO(LP)*BYC(LP) 

  END SUBROUTINE PRHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ELLPSOLV    ( PCE , B , P )
!------------------------------------------------------------
! Sets up and solves 
!              dd/dzz[p'] -k^2 p = d/dz[rho*B]
! Subject to
!              d/dz[p'] = rho*B at z=0,H
!------------------------------------------------------------
! External controls:
!  PRESSURE_BY_SLICE   logical parameter set in ce1_consts.F90
!
! Inputs:
!  PCE                  A derived type containing convective 
!                       element data (see ce1_types.F90)
!
!  B         real[L+1]  RHS of p-eq 
!                         - B(2:L)=d/dz(rho*buoy)
!                         - B([1,L+1]) = rho*buoy|_[1,L]
!                         - rho is kg m-3. buoy is m s-2 
! Output:
!  P         real[L+1]  Perturbation pressure (Pa)
!                       due to buoyancy force.
!------------------------------------------------------------
!  Two solutions possible depending on 
!------------------------------------------------------------

   type (T_ONE_PCE)  ,       intent(in)      ::  PCE
   REAL, DIMENSION(PCE%L+1), INTENT(IN)      ::  B
   REAL, DIMENSION(PCE%L+1), INTENT(OUT)     ::  P

   
   real :: BB(PCE%L+1), BB0(PCE%L+1)
   real :: DU(PCE%L), D(PCE%L+1),DL(PCE%L), DD(PCE%L+1)
   real :: ZGL(PCE%L), ZGE(PCE%L+1),DDU(PCE%L),DDL(PCE%L)
   INTEGER :: L,I,INFO,L2,LP
   REAL :: KRADIUS2
  

#ifdef SKIPPFRC
   P=0.0
    write(*,*) " PFORCE=0 PFORCE=0  PFORCE=0  PFORCE=0  PFORCE=0  PFORCE=0 "
    RETURN
#endif

   LP=PCE%L
   
   ZGL = C_ZGL(PCE)
   ZGE = X_ZGE(PCE)
     
   DO L=2,LP
      DL(L-1)  =   1.0  / ( ( ZGE(L-1) - ZGE(L) )*( ZGL(L-1) - ZGL(L) ) )
   END DO
   DO L=2,LP
      DU(L)   =    1.0  / ( ( ZGE(L) - ZGE(L+1) )*( ZGL(L-1) - ZGL(L) ) )
   END DO

   DO L=2,LP
      D(L)   =           1.0  /  ( ZGE(L-1) - ZGE(L)   )

      D(L)   =   D(L) +  1.0  /  ( ZGE(L)   - ZGE(L+1) )

      D(L)   =   D(L) / ( ZGL(L-1) - ZGL(L) )  ! + KRADIUS2

      D(L)   =   -1.0 * D(L) 
   END DO



   D(LP+1) = -1.0 / ( ZGE(LP) - ZGE(LP+1) ) ! *( ZGL(LP-1) - ZGL(LP) ))
   DL(LP)  =  1.0 / ( ZGE(LP) - ZGE(LP+1) ) ! *( ZGL(LP-1) - ZGL(LP) ))
   


   D(1)    =  1.0 / ( ZGE(1) - ZGE(2) )    ! *( ZGL(1) - ZGL(2) ))
   DU(1)   = -1.0 / ( ZGE(1) - ZGE(2) )    ! *( ZGL(1) - ZGL(2) ))
  
   if ( .not. pressure_by_slice  ) THEN
          ! Do a single pressure solution for entire column. Requires
          ! deciding on a single effective cloud radius
          !----------------------------------------------------------
      KRADIUS2 = C_TYPICAL_RADIUS(PCE) ! Effective radius estimate.
      KRADIUS2 = ( 2.4 / KRADIUS2 )**2 ! Rescale based on root of J_0 Bessel Fun.
                                       ! (Holton, 1973, MWR 101, 201)
      D(2:LP)  =  D(2:LP) - KRADIUS2
      BB       =  B
      BB(1)    =  B(1)
      BB(LP+1) =  B(LP+1)

      call sgtsv(  LP+1 , 1, DL, D, DU, BB, LP+1, INFO )

      P = BB

   else
                                write(*,*) " SSSLIIIICCEEEE  PRESSURE "
           ! Do a slice-by-slice pressure soln, i.e., adding up 
           ! buoyancy contributions from each layer to total 
           ! pressure perturbation. Will allow hght-dep cloud 
           ! radii. Note sgtsv has side on DU, D, DL so these 
           ! need to be refreshed before each call
           !----------------------------------------------------- 
      P        = 0.0   !Initialize pressure pert to zero
      do l = 1,lp+1
         BB       =  0.    
         !!KRADIUS2 = C_TYPICAL_RADIUS(PCE) ! Column-wide effective radius estimate.
         KRADIUS2 = MAX( SQRT( PCE%A(l) / PI ) , 10.0 ) ! Local effective radius estimate.
         KRADIUS2 = ( 2.4 / KRADIUS2 )**2 ! Based on root of J_0 Bessel Fun.
                                          ! (Holton, 1973, MWR 101, 201)
         DD(2:LP) =  D(2:LP) - KRADIUS2
         DD(1)    =  D(1)
         DD(LP+1) =  D(LP+1)
         BB(L)    =  B(L)
         DDL      =  DL
         DDU      =  DU

         call sgtsv(  LP+1 , 1, DDL, DD, DDU, BB, LP+1, INFO )
         P = P + BB
      enddo
   endif

  END SUBROUTINE ELLPSOLV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end MODULE CE1_PRSOLV01
