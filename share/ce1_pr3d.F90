!  $Id: ce1_pr3d.F90 created 2014/06/11 NCAR $
!-------------------------------------------------------
!      2015/07/30: Needs to be "Anelasticized" according 
!                  to Lipps Hemler, EuLag etc.
!--------------------------------------------------------
MODULE CE1_PR3D


use CE1_TYPES
use CE1_UTILS
use CE1_INFORM
use CE1_CONSTS
!------------------------------------
! !DESCRIPTION:
!
!  Ultimate goal: Pressure solution based on 
!  3D reconstruction of source terms for elliptic
!  pressure equation. 

 IMPLICIT NONE
 PRIVATE

 PUBLIC GRIDOVERLAY
 PUBLIC ZMOMG3_CE1

   character*72, parameter :: I_am_module = "ce1_pr3d" 

  real :: pscl=1.0

  integer :: debug_wrt_intrv = 10

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define DEBUG
!#define CUBICWGHT
#define QUADWGHT
!!#define BELLJAR
!!#define BELLRECON
!!#define NODYNPRES
#define FVPRESSURE
!!#define SMOOTHDYNSRC
!!#define SMOOTHBYCSRC
#define SMOOTHTOTSRC


  SUBROUTINE GRIDOVERLAY( PCES, DT, PP, PPX, PPY, PPR, PPZ, PHYC, BYC, NCL, NXC, NYC, LP )

   type (T_ONE_PCE)  , intent(in   )    ::  PCES(NCL)
   integer, intent(in  ) :: NCL, LP,NXC,NYC
   real, intent(out)     :: PP(LP+1,NCL),PPZ(LP,NCL),BYC(LP,NCL),PHYC(LP+1,NCL)
   real, intent(out)     :: PPX(LP+1,NCL),PPY(LP+1,NCL),PPR(LP+1,NCL)
   real, intent(in)      :: DT
   



   real, allocatable  :: celli(:,:,:),xc(:),yc(:),dc(:),cell0(:,:,:)
   real, allocatable  :: dphyd(:,:,:),cellu(:,:,:),celld(:,:,:),phyd(:,:,:),pnhz(:,:,:) 
   real, allocatable  :: w3(:,:,:),ur3(:,:,:),d2ph(:,:,:),pnh2(:,:,:)
   real, allocatable  :: pnh(:,:,:)
   !real, allocatable, save  :: pnh(:,:,:)
   real, allocatable  :: p2(:,:,:), p2x(:,:,:), p2y(:,:,:), p2r(:,:,:),u3(:,:,:),v3(:,:,:),th3(:,:,:)
   real, allocatable  :: ztx3(:,:,:),zty3(:,:,:),ztz3(:,:,:)
   real :: xx, ys, yn, dx, dy, rad, xxl2, xw,xe,swght,xci,ycj, ARMAX
   
   real :: RHOL(LP),DZ(LP),THD(LP),RHD(LP)
   real :: QD(LP), QLD(LP),QID(LP)
   real :: QRD(LP),QHD(LP),QSD(LP)
   real :: KE(LP+1),WX(LP+1)
   real :: DUCZ_diag(nxc,lp),DKE_diag(nxc,lp),UR_diag(nxc,lp),ZAX_diag(nxc,lp),DIV_diag(nxc,lp),Exz_diag(nxc,lp)

   integer :: i,j,n,l,nshft,ncells,si

   logical :: do_shafts=.false. ! .true.

#ifdef DEBUG
   real :: THD_S(LP),RHD_S(LP)
   real :: QD_S(LP), QLD_S(LP),QID_S(LP)
   real :: QRD_S(LP),QHD_S(LP),QSD_S(LP)
   real :: THD_U(LP),RHD_U(LP)
   real :: QD_U(LP), QLD_U(LP),QID_U(LP)
   real :: QRD_U(LP),QHD_U(LP),QSD_U(LP)
   RHD_S = 0.  
   RHD_u = 0.
   THD_S = 0.  
   THD_u = 0.
   QD_S  = 0. 
   QD_u  = 0.
   QLD_S = 0.  
   QLD_u = 0.
   QID_S  = 0. 
   QID_u  = 0.
   QRD_S = 0.  
   QRD_u = 0.
   QSD_S  = 0. 
   QSD_u  = 0.
   QHD_S  = 0. 
   QHD_u  = 0.
#endif


   dy=1000.0 ! meters
   dx=1000.0

    allocate ( pnh2(nxc , nyc, 0:lp) )
    allocate ( pnh(nxc , nyc, 0:lp) )
 
  
    allocate ( celli(nxc , nyc, lp) )
    allocate ( cellu(nxc , nyc, lp) )
    allocate ( celld(nxc , nyc, lp) )
    allocate ( cell0(nxc , nyc, lp) )
    allocate ( th3(nxc , nyc, lp) )
    allocate ( dphyd(nxc , nyc, lp) )
    allocate ( pnhz(nxc , nyc, lp) )
    allocate ( phyd(nxc , nyc, 0:lp) )
    allocate ( d2ph(nxc , nyc, 0:lp) )
    allocate ( u3(nxc , nyc, lp) )
    allocate ( v3(nxc , nyc, lp) )
    allocate ( w3(nxc , nyc, lp) )
    allocate ( ztx3(nxc , nyc, lp) )
    allocate ( zty3(nxc , nyc, lp) )
    allocate ( ztz3(nxc , nyc, lp) )
    allocate ( ur3(nxc , nyc, lp) )
    allocate ( p2(nxc , nyc, 0:lp) )
    allocate ( p2x(nxc , nyc, 0:lp) )
    allocate ( p2y(nxc , nyc, 0:lp) )
    allocate ( p2r(nxc , nyc, 0:lp) )
    allocate ( xc(nxc) )
    allocate ( yc(nyc) )
    allocate ( dc(nyc) )


    pnh=0.0
    pnh2=0.0
    celli=0.0
    cellu=0.0
    celld=0.0
    cell0=0.0
    dphyd=0.0
    th3=0.0
    w3=0.0
    u3=0.0
    v3=0.0
    ztx3=0.0
    zty3=0.0
    ztz3=0.0

    BYC = 0.



    do n=1,ncl 

    if (IS_UPDRAFT( pces(n) ) ) THEN
       RHOL =  C_RHL( PCES(n))
       DZ   =  C_LAYER_THICKNESS( PCES(n) )   
       THD  =  PCES(n)%THETA !- PCEX %THETA
             ! need to rethink in light of 
             ! shared enviro -12/18/15
       QD   =  0. 
       QLD  =  0.
       QID  =  0. 
       QSD  =  0.
       QRD  =  0.
       QHD  =  0.

       !   - anelastic "p" is actually p/rho
       RHD  = -THD / PCES(n)%THBCK ! - 0.61*QD + QLD + QID + QRD + QSD + QHD
       !RHD  = -RHOL * THD / PCES(n)%THBCK ! - 0.61*QD + QLD + QID + QRD + QSD + QHD
      
       BYC(:,n) = -GRAV * RHD

#ifdef DEBUG
         RHD_u = RHD
         THD_u = THD
         QD_u  = QD
         QLD_u = QLD
         QID_u = QID
         QRD_u = QRD
         QSD_u = QSD
         QHD_u = QHD
#endif


       do i=1,nxc
          xc(i) = (i-nxc/2.)*dx +  pces(n)%xmean
       end do
       do i=1,nyc
          yc(i) = (i-nyc/2.)*dy +  pces(n)%ymean
       end do

       ARMAX = MAXVAL(PCES(n)%A)

       do L=1,lp
          rad=sqrt( pces(n)%a0 / PI )
          ys = pces(n)%y(l)-rad
          yn = pces(n)%y(l)+rad
          do j=1,nyc
             if ( (yc(j)<=yn) .and. (yc(j)>=ys) ) then
                xxl2 = sqrt( rad**2-( yc(j)-pces(n)%y(l) )**2 )
                xw   = pces(n)%x(l) - xxl2
                xe   = pces(n)%x(l) + xxl2
                dc   = sqrt( (xc - pces(n)%x(l))**2 + (yc(j) - pces(n)%y(l))**2 )
                where( (xc >=xw) .and. (xc<=xe) )
                     cell0(:,j,L)=1.0
                endwhere
             endif
           end do
        end do



       do L=1,lp
          rad=sqrt( pces(n)%a(l) / PI )
#ifdef BELLJAR
          rad=sqrt( ARMAX / PI )
#endif
          IF ( rad > control%pr3d_recon_radius ) then 
          ys = pces(n)%y(l)-rad
          yn = pces(n)%y(l)+rad
          do j=1,nyc
             if ( (yc(j)<=yn) .and. (yc(j)>=ys) ) then
                xxl2 = sqrt( rad**2-( yc(j)-pces(n)%y(l) )**2 )
                xw   = pces(n)%x(l) - xxl2
                xe   = pces(n)%x(l) + xxl2
                dc   = sqrt( (xc - pces(n)%x(l))**2 + (yc(j) - pces(n)%y(l))**2 )
                where( (xc >=xw) .and. (xc<=xe) )
                     celli(:,j,L)=1.0
                     cellu(:,j,L)= 1.0
#ifdef CUBICWGHT
                     cellu(:,j,L)= 1.0   - (dc/rad)**3
#endif
#ifdef QUADWGHT
                     cellu(:,j,L)= 1.0  - 1.*(dc/rad)**2 !++dbg
#endif
#ifdef BELLJAR
                     cellu(:,j,L)= pces(n)%a(L)/ARMAX
#endif 
                endwhere
             endif
          end do
          endif
       end do

       do L=1,lp
          ncells = count( celli(:,:,L) > 0. )
          swght  = sum( cellu(:,:,L ) )+1.e-6
          cellu(:,:,L)= (ncells/swght)*cellu(:,:,L)
       end do

#ifdef QUADWGHT
       where( cellu > 2.1 )
            cellu=2.1
       end where
#endif




       do L=1,lp
          where( celli(:,:,L) > 0.)
               dphyd(:,:,L) =  dphyd(:,:,L) + GRAV * cellu(:,:,L) * RHD(L)                              
               !dphyd(:,:,L) =  dphyd(:,:,L) + GRAV * celli(:,:,L) * RHD(L)                              
               th3(:,:,L)  =  th3(:,:,L) + cellu(:,:,L) * THD(L)                              
               w3(:,:,L)   =  w3(:,:,L) + cellu(:,:,L) * PCES(n)%W(L)                              
          endwhere
       end do
                      !   RHD   = -RHO * THD / PCE%THBCK

       nshft =   INDX_MY_PRSHAFT (  PCES(n) , PCES )

#ifdef DEBUG
   write(*,*) "Shaft index of this updraft is: ", nshft
#endif
       ! active prshaft exists for updraft 'n' and shafts 
       ! are used in perssure calc.
       !--------------------------------------------------
       if ( ( nshft > 0) .and. do_shafts )  then 

        RHOL =  C_RHL( PCES(nshft))
        DZ   =  C_LAYER_THICKNESS( PCES(nshft) )   
        THD  =  PCES(nshft)%THETA !- PCEX %THETA
             ! need to rethink in light of 
             ! shared enviro -12/18/15
        QD   =  0. 
        QLD  =  0. 
        QID  =  0.
        QSD  =  0.
        QRD  =  0. 
        QHD  =  0.

        RHD  = -THD / PCES(nshft)%THBCK - 0.61*QD + QLD + QID + QRD + QSD + QHD
        !RHD  = -RHOL * THD / PCES(nshft)%THBCK - 0.61*QD + QLD + QID + QRD + QSD + QHD

        BYC(:,nshft) = -GRAV * RHD

#ifdef DEBUG
         RHD_s = RHD
         THD_s = THD
         QD_s  = QD
         QLD_s = QLD
         QID_s = QID
         QRD_s = QRD
         QSD_s = QSD
         QHD_s = QHD
#endif

      ! Maximal shaft
       do L=2,lp
          where( (celli(:,:,L-1) >= 1.0) .and. celli(:,:,L) == 0.0 )
               celli(:,:,L) = 2.0
               celld(:,:,L) = 1.0
          endwhere 
       end do

       do L=1,lp
          ncells = count( celld(:,:,L) > 0. )
          swght  = sum( celld(:,:,L ) )+1.e-6
          celld(:,:,L)= (ncells/swght)*celld(:,:,L)
       end do

       do L=1,lp
          where( celld(:,:,L) > 0.)
               dphyd(:,:,L) =  dphyd(:,:,L) + GRAV * celld(:,:,L) * RHD(L)                              
          endwhere
       end do

      end if  ! prshaft exists and shafts are used prssure calc.

    end if ! is an updraft
    end do ! loop over clouds (only updrafts actually cause an action)
    
       phyd(:,:,0)=0.
       do L=1,lp
          phyd(:,:,L) = phyd(:,:,L-1) + dphyd(:,:,L)*dz(L)
       end do
       
       !if (.not. allocated(pnh)) then
       !   allocate ( pnh(nxc , nyc, 0:lp) )
       !   pnh=-phyd
       !end if
       !if (.not. allocated(pnh2)) then
       !   allocate ( pnh2(nxc , nyc, 0:lp) )
       !   pnh2=-phyd
       !end if


       call radial_wind_cart3d ( W3, RHOL, UR3, DX,DY,DZ(100),pces(1)%time, DT, NXc,NYc,LP )
          
            ! calculate pressure without dynamical forcing for diagnostic
            ! purposes. ZTX3 is =0.  
       call poiss3d_fft( Pnh2, Phyd , ztx3, d2PH, RHOL , DX,DY,DZ(100),NXc,NYc,LP ) 


       call poiss3d_fft( Pnh, Phyd , UR3, d2PH, RHOL , DX,DY,DZ(100),NXc,NYc,LP )


       
       do L=1,LP
          pnhz(:,:,L) =   ( pnh(:,:,L-1) - pnh(:,:,L) )/DZ(L)
       end do

    p2  = phyd+pnh
    p2x = 0.
    p2y = 0.
    p2r = 0.

          p2x(2:nxc-1,:,:) =   ( p2(3:nxc,:,:) - p2(1:nxc-2,:,:) )/DX
          p2y(:,2:nyc-1,:) =   ( p2(:,3:nyc,:) - p2(:,1:nyc-2,:) )/DY



    PP=0.
    PPZ=0.       



    do n=1,ncl 
    if (IS_UPDRAFT( pces(n) ) ) THEN

#ifdef FVPRESSURE
       write(*,*)" Area averaging P and Pforce ***FV-Pressure*** "
       do L=1,LP
          PPZ(L,n)  = sum( PNHZ(:,:,L)*celli(:,:,L) )/ ( sum( celli(:,:,L) )+1.e-4 )
          !PP(L,n)   = sum( PNH(:,:,L-1)*celli(:,:,L) )/ ( sum( celli(:,:,L) )+1.e-4 )
          PP(L,n)   = sum( PNH(:,:,L-1)*cell0(:,:,L) )/ ( sum( cell0(:,:,L) )+1.e-4 )
          PHYC(L,n) = sum( PHYD(:,:,L-1)*cell0(:,:,L) )/ ( sum( cell0(:,:,L) )+1.e-4 )
       end do
#else
       write(*,*)" Scaled ***Central profiles*** of P and Pforce: SCALE =  ",pscl
       do L=1,LP
          PPZ(L,n)  =  pscl*PNHZ(50,50,L)
          !PP(L,n)   =  pscl*PNH(50,50,L-1) 
          !PHYC(L,n) =  pscl*PHYD(50,50,L-1)
          PP(L,n)   = sum( PNH(:,:,L-1)*cell0(:,:,L) )/ ( sum( cell0(:,:,L) )+1.e-4 )
          PHYC(L,n) = sum( PHYD(:,:,L-1)*cell0(:,:,L) )/ ( sum( cell0(:,:,L) )+1.e-4 )
       end do
#endif

         ! Calculate radial pressure gradients for updraft.
         ! Not sure how to proceed for shafts
       do L=1,lp
        do j=1,nyc
         ycj = yc(j)-pces(1)%y(L)
         do i=1,nxc
            xci = xc(i)-pces(1)%x(L)
            rad = sqrt( xci**2 + ycj**2 )+1.e-5
            p2r(i,j,L) = (1./rad)*( (xci*p2x(i,j,L)) + (ycj*p2y(i,j,L)) )
         end do
        end do
       end do
       p2r(:,:,0)=p2r(:,:,1)

      if ( ( nshft > 0) .and. do_shafts )  then 
          do L=1,LP
             PPZ(L,nshft) = sum( PNHZ(:,:,L)*celld(:,:,L) )/ ( sum( celld(:,:,L) )+1.e-4 )
             PP(L,nshft) = sum( PNH(:,:,L-1)*celld(:,:,L) )/ ( sum( celld(:,:,L) )+1.e-4 )
         end do
       endif

    endif
    end do



    PP(LP+1,:)=PP(LP,:)
    PHYC(LP+1,:)=PHYC(LP,:)

   
#ifdef DEBUG
    if ( MOD( pces(1)%time, debug_wrt_intrv*DT ) .eq. 0. )  then
    !!if ( ( MOD( pces(1)%time, 100.) .eq. 0.) .or. (nshft>0) ) then
           write(*,*) " Writing from PR3d "
       write(311) nxc,nyc,lp,pces(1)%time
       write(311) RHOL,PCES(1)%THBCK,PCES(1)%THETA,PCES(1)%W,PCES(1)%A
       write(311) rhd_u,thd_u,qd_u,qld_u,qid_u,qrd_u,qsd_u,qhd_u
       write(311) rhd_s,thd_s,qd_s,qld_s,qid_s,qrd_s,qsd_s,qhd_s
       write(311) PP,PPZ
       write(311) celli
       write(311) cellu
       write(311) celld
       write(311) dphyd , th3
       write(311) phyd
       write(311) pnh
        write(311) pnh2
      !write(311) w3
       write(311) d2ph
       write(311) p2
       write(311) p2x
       write(311) p2y
       write(311) p2r
   endif
#endif


   
    deallocate ( cellu )
    deallocate ( celld )
    deallocate ( xc  )
    deallocate ( yc  )
    deallocate ( celli )
    deallocate ( cell0 )
    deallocate ( dphyd )
    deallocate ( pnhz )
    deallocate ( pnh2 )
    deallocate ( phyd )
    deallocate ( d2ph )
    deallocate ( u3 )
    deallocate ( v3 )
    deallocate ( w3 )
    deallocate ( th3 )
    deallocate ( ztx3 )
    deallocate ( zty3 )
    deallocate ( ztz3 )
    deallocate ( ur3 )
    deallocate ( p2 )
    deallocate ( p2x )
    deallocate ( p2y )
    deallocate ( p2r )
    deallocate ( pnh )
    deallocate ( dc )



   end SUBROUTINE GRIDOVERLAY




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine poiss3d_fft ( P, Ph, UR3, d2ph, RHO, DX,DY,DZ,NX,NY,NZ )

   integer, intent(in  ) :: NX,NY,NZ
   real, intent(inout)   :: P(nx,ny,0:nz)
   !real, intent(  out)   :: P2(nx,ny,0:nz)
   real, intent(in   )   :: Ph(nx,ny,0:nz)
   real, intent(in   )   :: dx,dy,dz
   real, intent(in)      :: RHO(nz)
   real, intent(in)      :: UR3(nx,ny,nz)
   real, intent(out)     :: d2ph(nx,ny,0:nz)

   real, dimension(NX,NY,0:NZ) :: rhsp
   real, dimension(NX,NY,0:NZ) :: ke,tmp

   integer :: i,j,k,n
 
   real :: A(0:NZ),B(0:NZ),C(0:NZ),C1,C2,RHZN(0:NZ)
   !real :: UR3(nx,ny,nz)

   integer :: xperod,yperod,zperod,ierror

! One source term for elliptic p_nh eq is HORIZONTAL del2 of p_h
!--------------------------------------------------------------- 

 !!!call radial_wind ( W3, RHO, UR3, DX,DY,DZ,NX,NY,NZ )

  rhsp=0.

  do k=0,nz
  do j=2,ny-1
  do i=2,nx-1
       

    d2ph(i,j,k) =     ( ph(i+1,j  ,k  ) + ph(i-1,j  ,k  ) - 2*ph(i,j,k)  ) / (dx**2)   &
                    + ( ph(i,  j+1,k  ) + ph(i,  j-1,k  ) - 2*ph(i,j,k)  ) / (dy**2)  

  end do
  end do
  end do
  

  d2ph(nx,:,: )     = 0. !d2ph(2,   :,:)
  d2ph(1, :,: )     = 0. !d2ph(nx-1,:,:)
  d2ph(:,ny,: )     = 0. !d2ph(:,2,   :)
  d2ph(:,1, : )     = 0. !d2ph(:,ny-1,:)

#ifdef SMOOTHBYCSRC
! Smooth 
   do n=1,5
     tmp = d2ph
     do k=1,nz-1
     do j=2,ny-1 
     do i=2,nx-1 
        d2ph(i,j,k) =  (1./7.)* (tmp(i,j,k) +    tmp(i+1,j,k)  +  tmp(i-1,j,k) &
                               + tmp(i,j+1,k)  + tmp(i,j-1,k)   + tmp(i,j,k+1) &
                               + tmp(i,j,k-1) )
     end do
     end do
     end do
   end do
#endif



#ifdef NODYNPRES
  rhsp = d2ph
    write(*,*) " excluding Dyn Press in PR3D !!!! "
#else
  ! ur3 same sign as e2ij which should be same sign as d2ph
  !---------------------------------------------------------
  rhsp(:,:,0   ) = d2ph(:,:, 0)
  rhsp(:,:,1:nz) = d2ph(:,:,1:nz) + ur3(:,:,1:nz)
    write(*,*) " 'Dyn Press' source **RETAINED** in PR3D !!!! "
#endif


#ifdef SMOOTHTOTSRC
! Smooth 
   do n=1,5
     tmp = rhsp
     do k=1,nz-1
     do j=2,ny-1 
     do i=2,nx-1 
        rhsp(i,j,k) =  (1./7.)* (tmp(i,j,k) + tmp(i+1,j,k)  + tmp(i-1,j,k)   & 
                               + tmp(i,j+1,k)  + tmp(i,j-1,k)   + tmp(i,j,k+1) &
                               + tmp(i,j,k-1) )
     end do
     end do
     end do
   end do
#endif



    do k=2,nz-1 
       rhzn(k) = (1./rho(k)) * ( rho(k-1)-rho(k+1) ) / (2*dz)
    end do
    k=0
    rhzn(k)=0.
    k=1
    rhzn(k)=(1./rho(k)) * ( rho(k)-rho(k+1) ) / (1*dz)
    k=nz
    rhzn(k)=(1./rho(k)) * ( rho(k-1)-rho(k) ) / (1*dz)




! Fishpack soln
!
!                        C1*(X(I-1,J,K) -2.*X(I,J,K) +X(I+1,J,K)) +
!                        C2*(X(I,J-1,K) -2.*X(I,J,K) +X(I,J+1,K)) +
!                        A(K)*X(I,J,K-1) +B(K)*X(I,J,K)+ C(K)*X(I,J,K+1)
!                        = F(I,J,K)
!

  C1 = 1./dx**2
  C2 = 1./dy**2

  A  = 1./dz**2 + rhzn/(2.*dz)  
  B  =-2./dz**2
  C  = 1./dz**2 - rhzn/(2.*dz)  


  xperod=1 ! zero BCs at sides
  yperod=1 !

#if 0
!--------------------------------------
!Zero BCs on P and both top and bottom
!--------------------------------------
  A(0)=0.
  C(NZ)=0.
  zperod=1
!--------------------------------------
#endif

#if 0
!--------------------------------------
! First order deriv w/ ghost point at NZ+1
! is zero
!--------------------------------------
  A(0)=0.
  C(NZ)=0.
  B(NZ)=-1./(DZ**2)
  zperod=1
!---------------------------------------
#endif

#if 1
!--------------------------------------
! Second order deriv w/ ghost point at NZ+1
! AND j=0 is zero
!--------------------------------------
  A(0)=0.
  C(0) =2./(DZ**2)
  C(NZ)=0.
  A(NZ)=2./(DZ**2)
  zperod=1
!---------------------------------------
#endif

#if 0
!--------------------------------------
! Second order deriv w/ ghost point at NZ+1
! is zero
!--------------------------------------
  A(0)=0.
  C(NZ)=0.
  A(NZ)=2./(DZ**2)
  zperod=1
!---------------------------------------
#endif

  ! RHS of P-eq is like  -d2ph
  p=-rhsp

write(*,*) " going in t pois3d from fishpack "

  call POIS3D(XPEROD, NX, C1, YPEROD, NY, C2, ZPEROD, NZ+1, A, B, C &
      , NX, NY, p , IERROR)



write(*,*) " got out pois3d from fishpack, exit condition: ", IERROR




end subroutine poiss3d_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE ZMOMG3_CE1(  PCE, P , PZ , BYC, DT )  

   type (T_ONE_PCE)  ,       intent(inout)  ::  PCE
   REAL, DIMENSION(PCE%L+1), INTENT(in)     ::  P
   REAL, DIMENSION(PCE%L)  , INTENT(in)     ::  PZ,BYC
   REAL,                     INTENT(in)     ::  DT

   REAL  :: BY(PCE%L), PFZ(PCE%L),RHL(PCE%L),RAW(PCE%L)
   REAL  :: ZGE(PCE%L+1)
   REAL  :: DZ(PCE%L), ZGL(PCE%L), TDAMP(PCE%L)

   REAL  :: KYY   

   INTEGER :: L,I,L2,LP
   LP=PCE%L

   KYY = 1.E3
   RHL  = C_RHL(PCE)
  

   ZGL = C_ZGL(PCE)
   ZGE = X_ZGE(PCE)

   DZ  = ZGE(1:LP)-ZGE(2:LP+1)


if (pce%ipc==1) then 
write(114) pce%time ,pce%a0
write(114) PCE%W 
write(114) PCE%A
write(114) RHL
write(114) P
write(114) PZ
write(114) BYC
endif


       !   - anelastic "p" is actually p/rho
#ifdef FVPRESSUREXXX
   where( ( PCE%A > control%zmom_minarea_ratio * PCE%A0 ).and.( PCE%A_NM1 > control%zmom_minarea_ratio * PCE%A0 ) )
      RAW = PCE%W * PCE%A !* RHL 
      RAW   = RAW - PCE%A * PZ * DT   ! + BYC * PCE%A *RHL *DT
      PCE%W = RAW  / ( PCE%A )  ! * RHL ) !  (PCE%A*RHL)
      !PCE%W = PCE%W - PCE%W*(PCE%A-PCE%A_NM1)*DT/PCE%A
   elsewhere
        PCE%W = 0.
   end where   
#else

      !PFZ = PZ
      PFZ = ( P(1:LP)-P(2:LP+1) ) / ( ZGE(1:LP)-ZGE(2:LP+1) ) 

      PCE%W  = PCE%W - PFZ * DT ! / RHL   ! + BYC * PCE%A *RHL *DT
              write(*,*) " DOING straight DW/Dt in PR3D !!!!!!!! ********* "
#endif
 


if (pce%ipc==1) then 
write(114) PCE%W 
endif

   if ( PCE%W(LP) < 0. ) PCE%W(LP)  = (1.- DT / 100. )*PCE%W(LP) ! If W is downward in lowest layer add extra damping 

   where( PCE%A < control%zmom_wdamp_areax * PCE%A0 )
       PCE%W  = PCE%W*( 1. - DT /    control%zmom_wdamp_tau )
   end where


 end subroutine ZMOMG3_CE1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine radial_wind_cart3d ( W3, RHO, UR3, DX,DY,DZ, Time, Dt, NX,NY,NZ )
    

   integer, intent(in  ) :: NX,NY,NZ
   real, intent(in   )   :: dx,dy,dz,time,dt
   real, intent(inout)   :: W3(nx,ny,nz),RHO(nz)
   real, intent(out)     :: UR3(nx,ny,nz)

   real   :: U3(nx,ny,nz),V3(nx,ny,nz),ztx3(nx,ny,nz),zty3(nx,ny,nz),ztz3(nx,ny,nz)
   real   :: DUCZ_diag(nx,nz),DKE_diag(nx,nz),UR_diag(nx,nz),ZAX_diag(nx,nz),DIV_diag(nx,nz),Exz_diag(nx,nz)


   real :: e2ij(nx,ny,nz),ztmg(nx,ny,nz)
   real :: ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz)
   real :: vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
   real :: wx(nx,ny,nz),wy(nx,ny,nz),wz(nx,ny,nz)

   integer :: i,j,k,n,l,i2
   real :: R2(nx,ny),XX(nx),YY(ny),ZZ(nz),UR2(nx,nz),UR2s(nx,nz),r0,ZAX(nx,nz),UCZR(nx,nz),UCZZ(nx,nz),Exz(nx,nz)
   real :: DRKE(nx,nz),DUCZ(nx,nz),KE(nx,nz),KEZ(nx,nz),KER(nx,nz),DKE(nx,nz)
   real :: tmp(nx,ny,nz),mfu(nz)
   real :: ARsubs,RADsubs

   RADsubs = NX*DX/8.
   UR3 =0.
   U3 =0.
   V3 =0.
   Ztx3 =0.
   Zty3 =0.
   Ztz3 =0.
   e2ij =0.

   do i=1,nx
      xx(i) = (i-nx/2)*DX
   end do
   do i=1,ny
      yy(i) = (i-ny/2)*DY
   end do
   do i=1,nz
      zz(i) = (NZ-i+0.5)*DZ 
   end do
   do j=1,ny
   do i=1,nx
      r2(i,j) = SQRT( xx(i)**2 + yy(j)**2 )
   end do
   end do

   mfu=0.
   do k=1,nz
   do j=1,ny
   do i=1,nx
      mfu(k)=mfu(k)+w3(i,j,k)*DX*DY ! *rho(k)
   end do
   end do
   end do

   ARsubs = pi *    (RADsubs**2)

   do k=1,nz
   do j=1,ny
   do i=1,nx
      if ( r2(i,j) < RADsubs ) w3(i,j,k)= w3(i,j,k) - mfu(k)/ARsubs
   end do
   end do
   end do

   


   



! Diagnose radial wind using continuity
   UR2 = 0.
   do j=2,nz-1
   do i=nx/2+1,nx-1
      UR2(i+1,j) = UR2(i-1,j) - (DX/DZ)*(1./rho(j))*( xx(i)*rho(j-1)*w3(i,ny/2,j-1) -   xx(i)*rho(j+1)*w3(i,ny/2,j+1)   )
   end do
   end do


   do j=1,nz
   do i=1,nx/2-1
      UR2(nx/2-i,j) = UR2(nx/2+i,j)
   end do
   end do

   do j=2,nz-1
   do i=nx/2+1,nx
      UR2(i,j) = UR2(i,j)/abs(xx(i))
   end do
   end do
   do j=2,nz-1
   do i=nx/2-1,1,-1
      UR2(i,j) = UR2(i,j)/abs(xx(i))
   end do
   end do

   ! at this point radial wind as function of x along y=0 (or as f[y] along x=0) is known



   do k=1,nz
   do j=ny/2,ny-1
      where(  (r2>=xx(j)).and.(r2<xx(j+1)) )
        !ur3(:,:,k) = ur2(j,k)
        ur3(:,:,k) = ur2(j,k)*(r2-xx(j+1))/(xx(j)-xx(j+1)) +  ur2(j+1,k)*(xx(j)-r2)/(xx(j)-xx(j+1))
      endwhere
   end do
   end do

   do k=1,nz
   do j=1,ny
   do i=1,nx
      u3(i,j,k) = ur3(i,j,k)*xx(i)/max( r2(i,j),.0001)
      v3(i,j,k) = ur3(i,j,k)*yy(j)/max( r2(i,j),.0001)
   end do
   end do
   end do


   do k=2,nz-1
   do j=2,ny-1
   do i=2,nx-1
      ztx3(i,j,k) =   ( w3(i,j+1,k) - w3(i,j-1,k) )/(yy(j+1)-yy(j-1)) & 
                    - ( v3(i,j,k+1) - v3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      zty3(i,j,k) = - ( w3(i+1,j,k) - w3(i-1,j,k) )/(xx(i+1)-xx(i-1)) & 
                    + ( u3(i,j,k+1) - u3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      ztz3(i,j,k) =   ( v3(i+1,j,k) - v3(i-1,j,k) )/(xx(i+1)-xx(i-1)) & 
                    - ( u3(i,j+1,k) - u3(i,j-1,k) )/(yy(j+1)-yy(j-1)) 
   end do
   end do
   end do






   do k=1,nz 
   do j=1,ny
   do i=2,nx-1 
      ux(i,j,k) =     ( u3(i+1,j,k) - u3(i-1,j,k) )/(xx(i+1)-xx(i-1))
      vx(i,j,k) =     ( v3(i+1,j,k) - v3(i-1,j,k) )/(xx(i+1)-xx(i-1))
      wx(i,j,k) =     ( w3(i+1,j,k) - w3(i-1,j,k) )/(xx(i+1)-xx(i-1))
   end do
   end do
   end do

   do k=1,nz
   do j=2,ny-1 
   do i=1,nx 
      uy(i,j,k) =      ( u3(i,j+1,k) - u3(i,j-1,k) )/(yy(j+1)-yy(j-1))
      vy(i,j,k) =      ( v3(i,j+1,k) - v3(i,j-1,k) )/(yy(j+1)-yy(j-1))
      wy(i,j,k) =      ( w3(i,j+1,k) - w3(i,j-1,k) )/(yy(j+1)-yy(j-1))
   end do
   end do
   end do

   do k=2,nz-1 
   do j=1,ny 
   do i=1,nx 
      uz(i,j,k) =      ( u3(i,j,k+1) - u3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      vz(i,j,k) =      ( v3(i,j,k+1) - v3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
      wz(i,j,k) =      ( w3(i,j,k+1) - w3(i,j,k-1) )/(zz(k+1)-zz(k-1)) 
   end do
   end do
   end do

   do k=2,nz-1 
   do j=2,ny-1 
   do i=2,nx-1 
      e2ij(i,j,k) =      4*(ux(i,j,k)**2)  +  2*((uy(i,j,k)+vx(i,j,k))**2)  +  2*((uz(i,j,k)+wx(i,j,k))**2)  &
                     +   4*(vy(i,j,k)**2)  +  2*((vz(i,j,k)+wy(i,j,k))**2)  &
                     +   4*(wz(i,j,k)**2) 
   end do
   end do
   end do

   e2ij=(1./4.)*e2ij


   ztmg = ztx3**2 + zty3**2 + ztz3**2



#ifdef SMOOTHDYNSRC
! Smooth 
   do n=1,5
     tmp = e2ij
     do k=2,nz-1 
     do j=2,ny-1 
     do i=2,nx-1 
        e2ij(i,j,k) =  (1./7.)* (tmp(i,j,k) + tmp(i+1,j,k)  + tmp(i-1,j,k)  + tmp(i,j+1,k)  + tmp(i,j-1,k)   + tmp(i,j,k+1)  + tmp(i,j,k-1) )
     end do
     end do
     end do
   end do

   do n=1,5
     tmp = ztmg
     do k=2,nz-1 
     do j=2,ny-1 
     do i=2,nx-1 
        ztmg(i,j,k) =  (1./7.)* (tmp(i,j,k) + tmp(i+1,j,k)  + tmp(i-1,j,k)  + tmp(i,j+1,k)  + tmp(i,j-1,k)   + tmp(i,j,k+1)  + tmp(i,j,k-1) )
     end do
     end do
     end do
   end do
#endif


     UR3 = e2ij - 0.5 * ztmg
     do k=1,nz
        UR3(:,:,k) = UR3(:,:,k) * Rho(k)  
     end do

#ifdef DEBUG
    if ( MOD( time, debug_wrt_intrv*DT ) .eq. 0. )  then
    !!if ( ( MOD( pces(1)%time, 100.) .eq. 0.) .or. (nshft>0) ) then
           write(*,*) " Writing from PR3d "
       write(313) nx,ny,nz,time
       write(313) xx,yy,zz
       write(313) u3
       write(313) v3
       write(313) w3
       write(313) e2ij
       write(313) ztmg
       write(313) ur3
   endif
#endif



return
end subroutine radial_wind_cart3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end MODULE CE1_PR3D
