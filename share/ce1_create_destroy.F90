!  $Id: ce1_create_destroy.F90,v 1.31 2008/06/09 20:50:11 bacmj Exp $
!
!-----------------------------
MODULE CE1_CREATE_DESTROY

use CE1_TYPES
use CE1_INFORM
use CE1_UTILS
use NUMERICAL_UTILITIES, only : QSAT, DQSAT
use CE1_CONSTS, only: PI 

 IMPLICIT NONE
 PRIVATE

 PUBLIC MANAGE_PCES

 character*72, parameter :: I_am_module = "ce1_create_destroy" 


 REAL   :: GLOM_INTERVAL  = 500.
 REAL   :: LAST_GLOM_TIME = 0.0


 CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MANAGE_PCES( PCES, DT )
  type (T_ONE_PCE), INTENT(INOUT) :: PCES(:)
  real,             intent(IN   ) :: DT

   !call CHECK_STATUS( PCES , DT )

   call ACT_PREEXIST_UPD( PCES )

   call ACT_PRSHAFTS( PCES )

   !call CHECK_STATUS( PCES , DT )

END SUBROUTINE MANAGE_PCES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CHECK_STATUS   ( PCES , DT )

   type (T_ONE_PCE)  ,  intent(inout), dimension(:)  :: PCES
   REAL,   INTENT(IN) :: DT

   type (T_ONE_PCE)  :: PCE
  
   REAL              :: TAU,WT,STATUS_NOW

   INTEGER           :: NP,I

   np = size( PCES%IPC)


   DO I=1,NP

     PCE = PCES(I)

     IF ( PCE%TYPE_OF_CE == UpdraftCode ) then
        IF(CONTROL%NEVER_KILL_UPDS) then
          PCE%STATUS = 1.0
        ENDIF
     ENDIF    

     IF ( PCE%TYPE_OF_CE == PrecipShaftCode ) then
        IF(CONTROL%NEVER_KILL_UPDS) then
          PCE%STATUS = 1.0
        ENDIF
     ENDIF    

     PCES(I) = PCE

   END DO



  END SUBROUTINE CHECK_STATUS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ACT_PREEXIST_UPD( PCES )

  type (T_ONE_PCE), INTENT(INOUT) :: PCES(:)
  INTEGER :: NP,I

  np = size( PCES%IPC)

  do i=1,np
     IF(                                                         &
       ( PCES(i)%TYPE_OF_CE == UpdraftCode )            .AND.  &
       ( PCES(i)%STATUS     <  DormantStatusCode )      .AND.  &
       ( PCES(i)%BIRTHTIME  >= 0.0 )                    .AND.  & 
       ( PCES(i)%TIME       >=  PCES(i)%BIRTHTIME ) )          then

          PCES(i)%STATUS  = 1.00

              write(*,*) "UPDRAFT ",PCES(i)%IPC," IS BORN AT ",PCES(i)%TIME 
      end if
  end do
          

end SUBROUTINE ACT_PREEXIST_UPD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ACT_PRSHAFTS( PCES )
  use CE1_CONSTS, only: grav
  type (T_ONE_PCE), INTENT(INOUT) :: PCES(:)
 



  INTEGER :: NP,IP, NEED, N ,L, NFREE, MAXIPC, NEWIPC, NPX,LP,K, NAVAIL,NPL,I, evin, NN
  LOGICAL,ALLOCATABLE :: WENEEDPRSHFTS (:), WEAREAVAILABLE(:)

  real :: MASS( PCES(1)%L ), tmass

  np = size( PCES%IPC)
  Lp = PCES(1)%L
  K  = 0 
  
#ifdef VERBY
 write(*,*) " LP, K" , LP,K
#endif


                      ! Allocate logicals for use by query functions
                      !--------------------------------------------
  ALLOCATE(WENEEDPRSHFTS(NP))
  ALLOCATE(WEAREAVAILABLE(NP))

                      ! Find and tally those PCEs that <<need>> precip shafts
                      ! and then those that are available, i.e. <<available for work>>
                      ! This routine is in this module. Compares PCE%QR with critical
                      ! value and sets a need_prshaft code in PCE
                      !------------------------------------------------------------
  CALL PRSHAFTS_NEEDED( PCEs, NP )

  WENEEDPRSHFTS   = WE_NEED_PRSHFTS(  PCEs , NP )
  WEAREAVAILABLE  = WE_ARE_AVAILABLE (  PCEs , NP )


  NEED    = COUNT( WENEEDPRSHFTS )
  NAVAIL  = COUNT( WEAREAVAILABLE )


  IF ( NEED > NAVAIL ) THEN
     WRITE(*,*) " ADMIN PROBLEM IN SETUP_SHAFTS, NEED \= AVAILABLE "
     STOP
  ELSE 
  ENDIF


               ! Now run through PCEs looking for those that need shafts
               ! Could prob use WENEEDPRSHFTS here, but dont ...
               !--------------------------------------------------------
  DO N=1,NP
                  ! if n-th PCE needs a shaft go into block
                  !----------------------------------------
     IF ( I_NEED_PRSHFT(PCES(N)) ) THEN
 
      !    Log info about shaft creation 
      !------------------------------------------------
      write(161,*) '                                            ==   '
      write(161,*) ' *-*_*_*_* ===================================   '
      write(161,*) ' A B O U T  T O  A C T I V A T E   P R E C I P  S H A F T '
      write(161,*) ' for PCE ',pces(N)%ipc,' at time=',pces(n)%time
            MASS   = C_LAYER_MASS(PCES(N))*100./grav
            tmass  = 0.
            do nn = 1,np 
            write( 161, *) "Mass ",NN, SUM( MASS * PCES(NN)%A , 1 )
            tmass = tmass + SUM( MASS * PCES(NN)%A , 1 )
            end do
            write( 161, *) "Total Mass ",tmass
      !-------------------------------------------------

                     ! Find first remaining available PCE.
                     ! NPL is index of this PCE in PCES(1:NP).
                     !----------------------------------------
         do i=1,np
            npl=i
            if(weareavailable(npl)) exit
         end do

         ! Find index for environment of updraft.
         ! Will be same for shaft
         !-----------------------------
         evin = indx_my_enviro( PCES(n) , PCES )
         !!
         !!write(*,*) " evin ", evin
         !!write(*,*) "pces(N)%my_env ",PCES(N  )%MY_ENV
         !!write(*,*) "pces(evin)%ipc ",PCES(evin)%IPC
         !!STOP
         !!

          WEAREAVAILABLE(NPL)     = .FALSE.  ! PCES(NPL) is no longer available!!

          PCES(NPL)%Type_Of_CE  = PrecipShaftCode  ! It is a shaft
          PCES(NPL)%STATUS      = 1.00             ! The shaft has STATUS=1.0
          PCES(NPL)%MAXWTF      = 1.00             ! Set backward time avg Q_prec_cond 
          PCES(NPL)%NEED_PRSHFT = FalseCode        ! shafts dont need shafts
          PCES(N  )%NEED_PRSHFT = FalseCode        ! its upd no longer needs one either

          PCES(N  )%MY_PRSHFT=PCES(NPL)%IPC        ! code in shaft ID to upd
          PCES(NPL)%MY_UPDRFT=PCES(N  )%IPC        ! code in upd ID to shaft
          PCES(NPL)%MY_ENV   =PCES(N  )%MY_ENV     ! give them the same env (obsolete?)
          PCES(NPL)%I        =PCES(N  )%I          ! give them the same gridbox i,j
          PCES(NPL)%J        =PCES(N  )%J          !               
          PCES(NPL)%NPOP     =PCES(N  )%NPOP       !               same population number

          PCES(NPL)%AGRID   =PCES(N)%AGRID
          PCES(NPL)%A0      =PCES(N)%A0

          PCES(NPL)%PPL     =PCES(N)%PPL
          PCES(NPL)%PPE     =PCES(N)%PPE
          PCES(NPL)%ZGE     =PCES(N)%ZGE
          PCES(NPL)%X       =PCES(N)%X
          PCES(NPL)%Y       =PCES(N)%Y

          PCES(NPL)%UR    =0.00
          PCES(NPL)%THETA =0.00
          !!++jtb 09/17
          PCES(NPL)%W     =0.00 ! PCES(evin)%W   
          PCES(NPL)%U     =PCES(evin)%U  
          PCES(NPL)%V     =PCES(evin)%V  
          PCES(NPL)%Q     =PCES(evin)%Q   
          PCES(NPL)%QL    =PCES(evin)%QL  
          PCES(NPL)%QI    =PCES(evin)%QI  
          PCES(NPL)%QR    =PCES(evin)%QR  
          PCES(NPL)%QS    =PCES(evin)%QS 
          PCES(NPL)%QH    =PCES(evin)%QH  

          PCES(NPL)%THBCK =PCES(evin)%THBCK  
          
          PCES(NPL)%A     = PCES(N)%A0*1.0e-5 !0.00
          PCES(NPL)%A_NM1 = PCES(N)%A0*1.0e-5
                    !--- remove area added to shaft
                    !--- from area of environment
          PCES(evin)%A    = PCES(evin)%A - PCES(NPL)%A 

          PCES(NPL)%RAIN_GAUGE= 0.0 
          PCES(NPL)%HAIL_GAUGE= 0.0 
          PCES(NPL)%SNOW_GAUGE= 0.0
          PCES(NPL)%BIRTHTIME = PCES(N  )%TIME  ! Record Time of Birth

          !call INIT_SHAFT_AREA( PCES(N) , PCES(NPL)  ) 


          NPL=NPL+1
          NEWIPC=NEWIPC+1

      !    Log info about shaft creation 
      !------------------------------------------------
      write(161,*) ' *-*_*_*_* '
      write(161,*) ' A C T I V A T E D   P R E C I P  S H A F T '
      write(161,*) '                  New Masses  '
            tmass = 0.
            do nn = 1,np 
            write( 161, *) "Mass ",NN, SUM( MASS * PCES(NN)%A , 1 )
            tmass = tmass + SUM( MASS * PCES(NN)%A , 1 )
            end do
            write( 161, *) "Total Mass ",tmass
      write(161,*) ' *-*_*_*_* '
      write(161,*) ' * '
      !------------------------------------------------

     END IF        
  END DO


 

  DEALLOCATE(WENEEDPRSHFTS)
  DEALLOCATE(WEAREAVAILABLE)


#ifdef VERBY
 write(*,*)" ip , ipc , prshft , updrft "
 do ip=1,npx
  write(*,*) ip,pces(ip)%ipc,pces(ip)%my_prshft,pces(ip)%my_updrft
 end do
#endif


END SUBROUTINE ACT_PRSHAFTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ADD_UPDRAFTS( PCES )

  type (T_ONE_PCE), POINTER :: PCES(:)
  type (T_ONE_PCE), pointer ::  PCES_COPY(:)

  type (T_ONE_PCE) ,    pointer  ::  PCENVS(:)
  type (T_ONE_PCE) ,    pointer  ::  PCECLST(:)



   integer :: NoOfClusters,NoInCluster,AddNUpdrafts 

   integer :: NP,LP,K,ic,ip,npx



  real :: ZeroBirthsFor, ADUM


  np = size( PCES%IPC)
  Lp = PCES(1)%L
  K  = 0
  
#ifdef VERBY
 write(*,*) " LP, K" , LP,K
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Identify clusters based on environments  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
              ! First find envs in PCES and
              ! pack into PCENVS     
   NoOfClusters = COUNT( FIND_ALL_ENVS(  PCES ,  NP) )
   allocate( PCENVS( NoOfClusters ) )

   PCENVS=PACK( PCES, FIND_ALL_ENVS(  PCES ,  NP) )

   AddNUpdrafts = 0


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Cluster Envs have been id d. Now start
   ! looping through clusters to update   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do ic=1,NoOfClusters

                    ! Env for current cluster

        NoInCluster = COUNT( FIND_MY_UPDS( PCES, PCENVS(ic), NP ) ) 

        write(*,*) " number of updrafts ", NoInCluster

        IF ( NoInCluster == 0 ) then

           ZeroBirthsFor =  PCENVS(ic)%AGE

        ELSE         
        
           allocate( PCECLST( NoInCluster ) )
           PCECLST = PACK( PCES , FIND_MY_UPDS( PCES, PCENVS(ic), NP ) ) 

           ZeroBirthsFor =  MINVAL( PCECLST%AGE )
           deallocate( PCECLST )

        ENDIF
              
        write(*,*) "No updrafts created for ",ZeroBirthsFor


        if ( ZeroBirthsFor > 500.0 ) AddNUpdrafts = 1

        NPX = NP + AddNUpdrafts

        IF ( AddNUpdrafts .GT. 0 ) THEN

        allocate(PCES_COPY(NP))

        do ip=1,np
           call CREATE_CE1(  LP , K , PCES_COPY(IP) )
           call COPY_PCE( PCES(IP) , PCES_COPY(IP) )
        enddo

        nullify(pces)
        allocate(PCES(NPX))

        do ip=1,npx
           call CREATE_CE1( LP , K ,  PCES(IP))
        end do

        do ip=1,np
           call COPY_PCE( PCES_COPY(IP) , PCES(IP) )
        enddo

        ADUM=PI*( 1000.**2)
        CALL STARTBUBBLE ( PCES(NPX)                 , & 
                           ADUM                      , & 
                           0.0                       , & 
                           0.0                       , & 
                           700.                      , & 
                           300.                        )


   
        end if
        
        end do




END SUBROUTINE ADD_UPDRAFTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INIT_SHAFT_AREA( PCE , PCESHFT )

  !INTEGER, intent(in   )               :: L
  type (T_ONE_PCE), INTENT(IN)         :: PCE
  type (T_ONE_PCE), INTENT(INOUT)      :: PCESHFT

  REAL, DIMENSION(PCE%L)        :: MASS
  REAL :: AQ, XQ, YQ

  INTEGER :: IP, NEED, N ,L, NFREE,  NEWIPC, NPL

  MASS = C_LAYER_MASS(PCE)
  AQ   =  SUM( MASS*PCE%A*(PCE%QR + PCE%QH + PCE%QS ) ) /SUM( MASS*(PCE%QR + PCE%QH + PCE%QS ) )
  XQ   =  SUM( MASS*PCE%X*(PCE%QR + PCE%QH + PCE%QS ) ) /SUM( MASS*(PCE%QR + PCE%QH + PCE%QS ) )
  YQ   =  SUM( MASS*PCE%Y*(PCE%QR + PCE%QH + PCE%QS ) ) /SUM( MASS*(PCE%QR + PCE%QH + PCE%QS ) )
  
  IF( MAXVAL(PCESHFT%A_NM1) == 0.0 ) THEN
     PCESHFT%A_NM1 = AQ
     PCESHFT%A0    = AQ
     PCESHFT%ASH   = AQ
     PCESHFT%A     = AQ
     PCESHFT%X     = XQ
     PCESHFT%Y     = YQ
  END IF


 END SUBROUTINE INIT_SHAFT_AREA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INIT_SHAFT_AREA_2( PCE , PCESHFT )

  !INTEGER, intent(in   )               :: L
  type (T_ONE_PCE), INTENT(IN)         :: PCE
  type (T_ONE_PCE), INTENT(INOUT)      :: PCESHFT

  REAL, DIMENSION(PCE%L)        :: MASS
  REAL :: AQ, XQ, YQ

  INTEGER :: IP, NEED, N ,L, NFREE,  NEWIPC, NPL, Lq(1)

  MASS = C_LAYER_MASS(PCE)
  XQ   =  SUM( MASS*PCE%X*(PCE%QR + PCE%QH + PCE%QS ) ) /SUM( MASS*(PCE%QR + PCE%QH + PCE%QS ) )
  YQ   =  SUM( MASS*PCE%Y*(PCE%QR + PCE%QH + PCE%QS ) ) /SUM( MASS*(PCE%QR + PCE%QH + PCE%QS ) )
  
                                  ! find loc of largest precip mass - LQ
  LQ = MAXLOC( PCE%A*MASS*(PCE%QR + PCE%QH + PCE%QS ) ) 
  AQ = PCE%A( LQ(1) )                ! Upd area at LQ

  IF( MAXVAL(PCESHFT%A_NM1) == 0.0 ) THEN
     PCESHFT%A0    = AQ
     PCESHFT%ASH   = AQ
     PCESHFT%X     = XQ
     PCESHFT%Y     = YQ

     PCESHFT%A(LQ(1):PCESHFT%L)     = AQ
     PCESHFT%A(1:LQ(1))             = 0.0

     PCESHFT%A_NM1 = AQ
  END IF


 END SUBROUTINE INIT_SHAFT_AREA_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION WE_NEED_DNDRFTS(  PCES , NP ) RESULT( NEED_DDRFTS )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(:)
   integer,                intent(in)    ::     NP
   
   logical, dimension(NP) :: NEED_DDRFTS


   NEED_DDRFTS =( PCES%NEED_DNDRFT.eq.TrueCode) 


 end FUNCTION WE_NEED_DNDRFTS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION I_NEED_PRSHFT(  PCE ) RESULT( NEED )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE
   
   logical  :: NEED


   NEED = (( PCE%NEED_PRSHFT.eq.TrueCode).and.(PCE%TYPE_OF_CE.ne.LatentCode) )


 end FUNCTION I_NEED_PRSHFT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION WE_NEED_PRSHFTS(  PCES , NP ) RESULT( NEED_DDRFTS )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(:)
   integer,                intent(in)    ::     NP
   
   logical, dimension(NP) :: NEED_DDRFTS


   NEED_DDRFTS = (( PCES%NEED_PRSHFT.eq.TrueCode).and.(PCES%TYPE_OF_CE.ne.LatentCode) )  


 end FUNCTION WE_NEED_PRSHFTS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION WE_ARE_AVAILABLE(  PCES , NP ) RESULT( WE_AVAIL )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(:)
   integer,                intent(in)    ::     NP
   
   logical, dimension(NP) :: WE_AVAIL

             ! PCEs are available if they are both latent and have a negative 
             ! birthtime

   WE_AVAIL = (( PCES%TYPE_OF_CE.eq.LatentCode) .AND.  (PCES%BIRTHTIME <  0.0 )) 


 end FUNCTION WE_ARE_AVAILABLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use to be in uPhys. More sensible to have here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PRSHAFTS_NEEDED( PCECLST , NCL )

   INTEGER,                INTENT(IN   )    ::     NCL
   type (T_ONE_PCE)   ,    intent(inout)    ::     PCECLST(NCL)

   INTEGER :: LP,I,L,N
                          !  
                          !  CreatePrecShaftValue moved to control% in ce1_types 
                          !  default value set in PceGridComp (4/14/12)
   DO N=1,NCL

        IF (  ( MAXVAL(PCECLST(N)%QR)   >  control%CreatePrecShaftValue ) .and. &
              ( PCECLST(N)%TYPE_OF_CE   == UpdraftCode)    .and. &   ! only updrafts
              ( PCECLST(N)%MY_PRSHFT    == MissingCode) )        &   ! only if you dont have one already
          THEN                 
          PCECLST(N)%NEED_PRSHFT=TrueCode          
        END IF
 
   END DO

           !!! write(*,*) I_am_module,"_PRSHAFTS ",PCECLST%NEED_PRSHFT

  end SUBROUTINE PRSHAFTS_NEEDED

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end MODULE CE1_CREATE_DESTROY


