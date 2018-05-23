!  $Id: ce1_diags.F90,v 1.24 2008/06/09 20:50:11 bacmj Exp $
!
!-----------------------------
MODULE CE1_DIAGS

use CE1_TYPES

 IMPLICIT NONE
 PRIVATE

 PUBLIC DIAGCTL
 !!PUBLIC T_PCE_DIAGS
 PUBLIC INIT_DIAGS
 PUBLIC DESTROY_DIAGS
 PUBLIC WRITE_CE_OUTPUTS
 PUBLIC SIMPLE_WRITE
 PUBLIC MORR_UPHYS_DIAGS

  PUBLIC PCEDIAGS
  
 type (T_DIAG_CTL), save :: DIAGCTL

 type (T_PCE_DIAGS), save     ::     PCEDIAGS
 character*72, parameter :: I_am_module = "ce1_diags" 

 real, parameter :: initval=0.

CONTAINS




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE INIT_DIAGS( LP,NP ) 
                   

 !!type (T_PCE_DIAGS), intent(inout)     ::     PCEDIAGS
 INTEGER,   intent(IN   )              ::     LP,NP

 allocate( PCEDIAGS%DHN(LP,NP) )
     PCEDIAGS%DHN=initval

 allocate( PCEDIAGS%PRACS(LP,NP) )
     PCEDIAGS%PRACS=initval
 allocate( PCEDIAGS%MUTT(LP,NP) )
     PCEDIAGS%MUTT=initval
 allocate( PCEDIAGS%TMUPH(LP,NP) )
     PCEDIAGS%TMUPH=initval
 allocate( PCEDIAGS%QMUPH(LP,NP) )
     PCEDIAGS%QMUPH=initval
 allocate( PCEDIAGS%QVSL(LP,NP) )
     PCEDIAGS%QVSL=initval
 allocate( PCEDIAGS%QVSI(LP,NP) )
     PCEDIAGS%QVSI=initval



 allocate( PCEDIAGS%E(LP,NP) )
     PCEDIAGS%E=initval
 allocate( PCEDIAGS%RHOA(LP,NP) )
     PCEDIAGS%RHOA=initval
 allocate( PCEDIAGS%RHOL(LP,NP) )
     PCEDIAGS%RHOL=initval
 allocate( PCEDIAGS%MASS(LP,NP) )
     PCEDIAGS%MASS=initval
 !allocate( PCEDIAGS%DZRW(LP,NP) )
 !    PCEDIAGS%DZRW=initval
 
allocate( PCEDIAGS%PP(LP+1,NP) )
     PCEDIAGS%MASS=initval
 allocate( PCEDIAGS%PHYD(LP+1,NP) )
     PCEDIAGS%MASS=initval
 allocate( PCEDIAGS%RCDOT(LP,NP) )
     PCEDIAGS%RCDOT=initval
 allocate( PCEDIAGS%PFZ(LP,NP) )
     PCEDIAGS%PFZ=initval


!========================================
! Surface, Integral, etc. ...,  quantities
!========================================
 allocate( PCEDIAGS%QTBKEN( NP) )
     PCEDIAGS%QTBKEN=initval
 allocate( PCEDIAGS%QTBKDE( NP) )
     PCEDIAGS%QTBKDE=initval
 allocate( PCEDIAGS%ABKEN( NP) )
     PCEDIAGS%ABKEN=initval
 allocate( PCEDIAGS%ABKDE( NP) )
     PCEDIAGS%ABKDE=initval
 allocate( PCEDIAGS%SHEXTOP( NP) )
     PCEDIAGS%SHEXTOP=initval
 allocate( PCEDIAGS%SHEXBOT( NP) )
     PCEDIAGS%SHEXBOT=initval
 allocate( PCEDIAGS%KEI( NP) )
     PCEDIAGS%KEI=initval
 allocate( PCEDIAGS%DKEI1( NP) )
     PCEDIAGS%DKEI1=initval

 end subroutine INIT_DIAGS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE DESTROY_DIAGS() 
                   

 !!type (T_PCE_DIAGS), intent(inout)     ::     PCEDIAGS


  deallocate( PCEDIAGS%PRACS )

 end subroutine DESTROY_DIAGS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE MORR_UPHYS_DIAGS ( PRACS , TTEND, TMUPH, QMUPH, QVSL, QVSI, IPC )
  real, intent(IN)    :: PRACS(:),TTEND(:),TMUPH(:),QMUPH(:),QVSL(:),QVSI(:)
  integer, intent(IN) :: IPC

  real, dimension(size(PRACS)):: tmp1,tmp2

  integer :: L,LP,I,J,N,LI

  LP = size(PRACS)

! invert height index
do L=1,LP
   Li      = LP+1-L
   PCEDIAGS%PRACS(Li,IPC) = PRACS(L)
   PCEDIAGS%MUTT(Li,IPC)  = TTEND(L)
   PCEDIAGS%TMUPH(Li,IPC) = TMUPH(L)
   PCEDIAGS%QMUPH(Li,IPC) = QMUPH(L)
   PCEDIAGS%QVSL(Li,IPC)  = QVSL(L)
   PCEDIAGS%QVSI(Li,IPC)  = QVSI(L)
end do



end SUBROUTINE MORR_UPHYS_DIAGS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE SIMPLE_WRITE( LUN   ,   &
                          PCES  , DT,  &
                          NP   )



  INTEGER, intent(in) :: LUN, NP
  type (T_ONE_PCE)   , intent(IN) :: PCES(NP)
  real               , intent(IN) :: DT

  type (T_ONE_PCE)        :: PCE
 !! type (T_DIAG_CTL), save :: DIAGCTL

  integer :: N , NPW, TKE_FLAG
 
  character*8 :: FIELD$

  real        :: TheTimeNow

  

     TKE_FLAG  =  0


     TheTimeNow = PCES(1)%TIME

!!if ( TheTimeNow >= 2000.) DIAGCTL%WRITE_INTERVAL=1.0

!!if( ( TheTimeNow - DIAGCTL%LAST_WRITE )  >= 10*DT ) then 
if( ( TheTimeNow - DIAGCTL%LAST_WRITE )  >= DIAGCTL%WRITE_INTERVAL ) then 

   DIAGCTL%LAST_WRITE = TheTimeNow

   write(*,*) " writing PCES at ",DIAGCTL%LAST_WRITE

     NPW=NP

   WRITE( LUN ) NPW

   do n=1,NPW

      PCE = PCES(n)

      WRITE( LUN )  PCE%IPC , PCE%L 
      WRITE( LUN)   PCE%I , PCE%J
      WRITE( LUN )  PCE%TYPE_OF_CE,TKE_FLAG
      WRITE( LUN )  PCE%MY_ENV
      WRITE( LUN )  PCE%MY_DNDRFT
      WRITE( LUN )  PCE%MY_PRSHFT
      WRITE( LUN )  PCE%MY_UPDRFT
      WRITE( LUN )  PCE%NPOP


      WRITE( LUN )  PCE%TIME 
      WRITE( LUN )  PCE%AGE 
      WRITE( LUN )  PCE%BIRTHTIME 
      WRITE( LUN )  PCE%DEATHTIME 
      WRITE( LUN )  PCE%CONFIG
      WRITE( LUN )  PCE%MAXW00
      WRITE( LUN )  PCE%MAXWTF
      WRITE( LUN )  PCE%STATUS
      WRITE( LUN )  PCE%RAIN_GAUGE
      WRITE( LUN )  PCE%SNOW_GAUGE
      WRITE( LUN )  PCE%HAIL_GAUGE
      WRITE( LUN )  PCE%XMEAN
      WRITE( LUN )  PCE%YMEAN
      WRITE( LUN )  PCE%A0
      WRITE( LUN )  PCE%ASH
      WRITE( LUN )  PCE%SHTOP


      WRITE( LUN )  PCE%X
      WRITE( LUN )  PCE%Y 
      WRITE( LUN )  PCE%ZGE 
                  ! Added 4/12/12
                  ! Dont know why not earlier
      WRITE( LUN )  PCE%PPE 

      WRITE( LUN )  PCE%U
      WRITE( LUN )  PCE%V
      WRITE( LUN )  PCE%W
      WRITE( LUN )  PCE%UR
      WRITE( LUN )  PCE%A
      WRITE( LUN )  PCE%Q
      WRITE( LUN )  PCE%QL
      WRITE( LUN )  PCE%QI
      WRITE( LUN )  PCE%QR
      WRITE( LUN )  PCE%QS
      WRITE( LUN )  PCE%QH
      WRITE( LUN )  PCE%THETA
      WRITE( LUN )  PCE%THBCK

          !---------------------------
          ! Additional Diag profiles 
          !---------------------------
      ! MuPHYS  
      if(allocated(PCEDIAGS%MUTT))        WRITE( LUN ) PCEDIAGS%MUTT(:,n)
      if(allocated(PCEDIAGS%PRACS))       WRITE( LUN ) PCEDIAGS%PRACS(:,n)
      if(allocated(PCEDIAGS%TMUPH))       WRITE( LUN ) PCEDIAGS%TMUPH(:,n)
      if(allocated(PCEDIAGS%QMUPH))       WRITE( LUN ) PCEDIAGS%QMUPH(:,n)
      if(allocated(PCEDIAGS%QVSL))        WRITE( LUN ) PCEDIAGS%QVSL(:,n)
      if(allocated(PCEDIAGS%QVSI))        WRITE( LUN ) PCEDIAGS%QVSI(:,n)
      !
      if(allocated(PCEDIAGS%DHN))         WRITE( LUN ) PCEDIAGS%DHN(:,n)
      if(allocated(PCEDIAGS%E))           WRITE( LUN ) PCEDIAGS%E(:,n)
      if(allocated(PCEDIAGS%MASS))        WRITE( LUN ) PCEDIAGS%MASS(:,n)
      if(allocated(PCEDIAGS%RHOA))        WRITE( LUN ) PCEDIAGS%RHOA(:,n)
      if(allocated(PCEDIAGS%RHOL))        WRITE( LUN ) PCEDIAGS%RHOL(:,n)
      if(allocated(PCEDIAGS%DZRW))        WRITE( LUN ) PCEDIAGS%DZRW(:,n)

      if(allocated(PCEDIAGS%PP))          WRITE( LUN ) PCEDIAGS%PP(:,n)
      if(allocated(PCEDIAGS%PHYD))        WRITE( LUN ) PCEDIAGS%PHYD(:,n)
      if(allocated(PCEDIAGS%RCDOT))       WRITE( LUN ) PCEDIAGS%RCDOT(:,n)


          !---------------------------------------------
          ! Additional Diag surface, integral quantities
          !---------------------------------------------
      if(allocated(PCEDIAGS%QTBKEN))      WRITE( LUN ) PCEDIAGS%QTBKEN(n)
      if(allocated(PCEDIAGS%QTBKDE))      WRITE( LUN ) PCEDIAGS%QTBKDE(n)
      if(allocated(PCEDIAGS%SHEXTOP))     WRITE( LUN ) PCEDIAGS%SHEXTOP(n)
      if(allocated(PCEDIAGS%SHEXBOT))     WRITE( LUN ) PCEDIAGS%SHEXBOT(n)



     if (allocated( PCEPHYS%RF ))      WRITE( LUN ) PCEPHYS%RF(:,n) 
     if (allocated( PCEPHYS%SF ))      WRITE( LUN ) PCEPHYS%SF(:,n) 
     if (allocated( PCEPHYS%HF ))      WRITE( LUN ) PCEPHYS%HF(:,n) 
  end do      
endif


 END SUBROUTINE SIMPLE_WRITE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE WRITE_CE_OUTPUTS( LUN      ,   &
                              PCE      ,   &
                              PCENV    ,   &
                              PCEDIAGS   )


  type (T_ONE_PCE)   , intent(IN) :: PCE
  type (T_ONE_PCE)   , intent(IN) :: PCENV
  type (T_PCE_DIAGS) , intent(IN) :: PCEDIAGS

  INTEGER, intent(in) :: LUN

  character*8 :: FIELD$



      WRITE( LUN )  PCE%IPC , PCE%L , PCE%TIME , PCE%AGE 
      WRITE( LUN )  PCE%CONFIG
      WRITE( LUN )  PCE%MAXW00
      WRITE( LUN )  PCE%MAXWTF
      WRITE( LUN )  PCE%STATUS
      WRITE( LUN )  PCE%X
      WRITE( LUN )  PCE%Y 
      !WRITE( LUN )  PCEDIAGS%ZGL 
      WRITE( LUN )  PCE%W
      WRITE( LUN )  PCE%UR
      WRITE( LUN )  PCE%A
      WRITE( LUN )  PCE%Q
      WRITE( LUN )  PCE%QL



 END SUBROUTINE WRITE_CE_OUTPUTS


END MODULE CE1_DIAGS
