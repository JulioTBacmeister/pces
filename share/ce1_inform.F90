!  $Id: ce1_inform.F90,v 1.17 2008/11/30 02:29:41 bacmj Exp $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Routines/functions that manage info about ce's, e.g., find relatives
!! locate environments, id types,.  Any code with logical
!! type(T_ONE_PCE), or integer result is a goo candidate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CE1_INFORM

use CE1_TYPES
use CE1_CONSTS, only: PI

 IMPLICIT NONE
 PRIVATE

 PUBLIC INDX_IS
 PUBLIC FIND_IDX_OF
 PUBLIC CLUSTER_ENV_IDX
 PUBLIC CLUSTER_ENV
 PUBLIC CE1_CLUSTER
 PUBLIC FIND_ALL_ENVS
 PUBLIC FIND_MY_ENV
 PUBLIC FIND_MY_KIDS
 PUBLIC FIND_US_ALL
 PUBLIC FIND_IN_I_J
 PUBLIC FIND_MY_UPDS
 PUBLIC FIND_MY_UNBORN
 PUBLIC FIND_MY_SISTERS
 PUBLIC IS_MY_ENV
 PUBLIC IS_SHAFT
 PUBLIC IS_UPDRAFT
 PUBLIC IS_DOWNDRAFT
 PUBLIC IS_ENVIRONMENT
 PUBLIC IS_DORMANT
 PUBLIC IS_ACTIVE_4_ENTR
 PUBLIC IS_PASSIVE_4_ENTR
 PUBLIC IS_ACTIVE_4_DYNMX
 PUBLIC IS_ACTIVE_4_UPHYS
 PUBLIC CAN_INDUCE_SUBSI
 PUBLIC WHO_TOUCH
 PUBLIC TYPE_OF
 PUBLIC INDX_MY_PRSHAFT 
 PUBLIC INDX_MY_ENVIRO
 PUBLIC WHERE_AM_I_IJ

 character*72, parameter :: I_am_module = "ce1_inform" 


 CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION TYPE_OF(  PCE  ) RESULT( TypeCode )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE
   integer     :: TypeCode

   TypeCode = PCE%TYPE_OF_CE
      
  end function TYPE_OF  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_SHAFT(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE    !! ME
   
   logical   :: IS_A_

   IS_A_ = ( PCE%TYPE_OF_CE .EQ. PrecipShaftCode )

 END FUNCTION IS_SHAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_UPDRAFT(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE    !! ME
   
   logical   :: IS_A_

   IS_A_ = ( PCE%TYPE_OF_CE .EQ. UpdraftCode )

 END FUNCTION IS_UPDRAFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_DORMANT(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE 
   
   logical   :: IS_A_

   IS_A_ = ( ( PCE%TYPE_OF_CE .EQ. LatentCode ) .OR. ( PCE%STATUS <= DormantStatusCode ) ) !++def

 END FUNCTION IS_DORMANT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_ACTIVE_4_ENTR(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE  
   
   logical   :: IS_A_

   IS_A_ = ( (.not.IS_DORMANT(PCE)) .and. (.not.IS_ENVIRONMENT(PCE))  ) !++def

 END FUNCTION IS_ACTIVE_4_ENTR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION CAN_INDUCE_SUBSI(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE   
   
   logical   :: IS_A_

   IS_A_ = ( (.not.IS_DORMANT(PCE)) .and. (.not.IS_ENVIRONMENT(PCE))  ) !++def

 END FUNCTION CAN_INDUCE_SUBSI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_PASSIVE_4_ENTR(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE 
   
   logical   :: IS_A_

   IS_A_ = ( (.not.IS_DORMANT(PCE)) ) !++def

  END FUNCTION IS_PASSIVE_4_ENTR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_ACTIVE_4_DYNMX(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE  
   
   logical   :: IS_A_

   IS_A_ = ( (.not.IS_DORMANT(PCE)) .and. (.not.IS_ENVIRONMENT(PCE))  ) !++def
   !IS_A_ = ( (.not.IS_DORMANT(PCE)) .and. (IS_UPDRAFT(PCE)) ) !++dbg

 END FUNCTION IS_ACTIVE_4_DYNMX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_ACTIVE_4_UPHYS(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE  
   
   logical   :: IS_A_

   IS_A_ = ( (.not.IS_DORMANT(PCE))  .and. (.not.IS_ENVIRONMENT(PCE)) ) !++def
   !!IS_A_ = ( (.not.IS_DORMANT(PCE))  .and. (IS_UPDRAFT(PCE)) ) !++dbg

 END FUNCTION IS_ACTIVE_4_UPHYS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_DOWNDRAFT(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE 
   
   logical   :: IS_A_

   IS_A_ = ( PCE%TYPE_OF_CE .EQ. DowndraftCode ) !++def

 END FUNCTION IS_DOWNDRAFT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_ENVIRONMENT(  PCE ) RESULT( IS_A_ )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE    !! ME
   
   logical   :: IS_A_

   IS_A_ = ( PCE%TYPE_OF_CE .EQ. EnvironmentCode )

 END FUNCTION IS_ENVIRONMENT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION FIND_MY_KIDS(  PCES , PCE, NP ) RESULT( IS_MY_KID )
   integer,                intent(in)    ::     NP
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(NP)
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE    !! ME
   
   integer                               ::     COMMON_ENV  
   logical, dimension(NP) :: IS_MY_KID

     IF ( PCE%TYPE_OF_CE.eq.EnvironmentCode ) THEN
        COMMON_ENV   = PCE%IPC   ! I can be the mommy
     ELSE
        COMMON_ENV   = -9999     ! I cant be a mommy
     ENDIF
                                                ! I cant be my own child
                                                ! AND kid has to have been
                                                ! born
     IS_MY_KID = (      (PCES%MY_ENV.eq.COMMON_ENV)  & 
                 .and.  (PCES%IPC.ne.PCE%IPC)        &
                 .and.  (PCES%TYPE_OF_CE.ne.LatentCode)  )

 end FUNCTION FIND_MY_KIDS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION FIND_IN_I_J(  PCES , NP, I, J ) RESULT( ONE_OF_US )
   integer,                intent(in)    ::     NP,I,J
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(NP)
   
   integer                               ::     COMMON_ENV  
   logical, dimension(NP) :: ONE_OF_US

                                                ! We share a geo. env. That is it.
     ONE_OF_US = ( ( PCES%I == I ) .and. (PCES%J == J  ) )   

 end FUNCTION FIND_IN_I_J

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION FIND_US_ALL(  PCES , PCE, NP ) RESULT( ONE_OF_US )
   integer,                intent(in)    ::     NP
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(NP)
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE    !! ME
   
   integer                               ::     COMMON_ENV  
   logical, dimension(NP) :: ONE_OF_US

                                                ! We share a geo. env. That is it.
     ONE_OF_US = ( PCES%MY_ENV == PCE%MY_ENV  )   

 end FUNCTION FIND_US_ALL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION FIND_MY_UPDS(  PCES , PCE, NP ) RESULT( IS_MY_UPD )
   integer,                intent(in)    ::     NP
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(NP)
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE    !! ME
   
   integer                               ::     COMMON_ENV  
   logical, dimension(NP) :: IS_MY_UPD

     IF ( PCE%TYPE_OF_CE.eq.EnvironmentCode ) THEN
        COMMON_ENV   = PCE%IPC   ! I can be the mommy
     ELSE
        COMMON_ENV   = -9999     ! I cant be a mommy
     ENDIF
                                                ! I cant be my own child
                                                ! AND kid has to have been
                                                ! born
     IS_MY_UPD = (      (PCES%MY_ENV.eq.COMMON_ENV)  & 
                 .and.  (PCES%IPC.ne.PCE%IPC)        &
                 .and.  (PCES%TYPE_OF_CE.eq.UpdraftCode)  )

 end FUNCTION FIND_MY_UPDS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION FIND_MY_UNBORN(  PCES , PCE, NP ) RESULT( IS_MY_KID )
   integer,                intent(in)    ::     NP
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(NP)
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE    !! ME
   
   integer                               ::     COMMON_ENV  
   logical, dimension(NP) :: IS_MY_KID

     IF ( PCE%TYPE_OF_CE.eq.EnvironmentCode ) THEN
        COMMON_ENV   = PCE%IPC   ! I can be the mommy
     ELSE
        COMMON_ENV   = -9999     ! I cant be a mommy
     ENDIF
                                                ! I cant be my own child
                                                ! AND kid has to be <<latent>>
     IS_MY_KID = (      (PCES%MY_ENV.eq.COMMON_ENV)  & 
                 .and.  (PCES%IPC.ne.PCE%IPC)        &
                 .and.  (PCES%TYPE_OF_CE.eq.LatentCode)  )

 end FUNCTION FIND_MY_UNBORN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION FIND_ALL_ENVS(  PCES , NP ) RESULT( I_AM_AN_ENV )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(:)
   integer,                intent(in)    ::     NP
   
   logical, dimension(NP) :: I_AM_AN_ENV


   I_AM_AN_ENV = (PCES%TYPE_OF_CE.eq.EnvironmentCode)


 end FUNCTION FIND_ALL_ENVS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION FIND_MY_ENV (  PCE , PCECLST ) RESULT( PCENV )   

   type (T_ONE_PCE)   ,    intent(in)    ::     PCECLST(:)
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE
   type (T_ONE_PCE)                      ::     PCENV

   integer :: n,ncl

   ncl=size(PCECLST%IPC)
   PCENV =PCE
   PCENV%TYPE_OF_CE = MissingCode

   do n=1,NCL
      if ( ( PCECLST(n)%IPC==PCE%MY_ENV ) .and. (PCECLST(n)%TYPE_OF_CE.eq.EnvironmentCode) ) then
         PCENV = PCECLST(n)
      endif
   end do

 end FUNCTION FIND_MY_ENV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION INDX_MY_PRSHAFT (  PCE , PCECLST ) RESULT( INDPRSH )   

   type (T_ONE_PCE)   ,    intent(in)    ::     PCECLST(:)
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE
   integer                               ::     INDPRSH

   integer :: n,ncl

   ncl=size(PCECLST%IPC)
   indprsh = -999

   do n=1,NCL
      if ( ( PCECLST(n)%IPC==PCE%MY_PRSHFT ) .and. (PCECLST(n)%TYPE_OF_CE.eq.PrecipShaftCode) ) then
         indprsh = n
      endif
   end do

 end FUNCTION INDX_MY_PRSHAFT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION INDX_MY_ENVIRO (  PCE , PCECLST ) RESULT( INDX )   

   type (T_ONE_PCE)   ,    intent(in)    ::     PCECLST(:)
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE
   integer                               ::     INDX

   integer :: n,ncl

   ncl=size(PCECLST%IPC)
   indx = -999

   do n=1,NCL
      if ( ( PCECLST(n)%IPC==PCE%MY_ENV ) .and. (PCECLST(n)%TYPE_OF_CE.eq.EnvironmentCode) ) then
         indx = n
      endif
   end do

 end FUNCTION INDX_MY_ENVIRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IS_MY_ENV (  PCEX , PCE ) RESULT( MY_ENV )   

   type (T_ONE_PCE)   ,    intent(in)    ::     PCEX
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE
   logical                               ::     MY_ENV

      MY_ENV = ( ( PCEX%IPC==PCE%MY_ENV ) .and. (PCEX%TYPE_OF_CE.eq.EnvironmentCode) )

 end FUNCTION IS_MY_ENV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION CLUSTER_ENV (  PCECLST ) RESULT( PCENV )   

   type (T_ONE_PCE)   ,    intent(in)    ::     PCECLST(:)
   type (T_ONE_PCE)                      ::     PCENV

   integer :: n,ncl,NoOfEnvs

   
   ncl=size(PCECLST%IPC)

   NoOfEnvs= COUNT( FIND_ALL_ENVS(  PCECLST ,  NCL) )

   PCENV =PCECLST(1)
   PCENV%TYPE_OF_CE = MissingCode

   do n=1,NCL
      if ( (PCECLST(n)%TYPE_OF_CE.eq.EnvironmentCode) ) then
         PCENV = PCECLST(n)
      endif
   end do

 end FUNCTION CLUSTER_ENV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION CLUSTER_ENV_IDX (  PCECLST ) RESULT( IDX )   

   type (T_ONE_PCE)   ,    intent(in)    ::     PCECLST(:)
   integer                               ::     IDX

   integer :: n,ncl,NoOfEnvs

   
   ncl=size(PCECLST%IPC)

   NoOfEnvs= COUNT( FIND_ALL_ENVS(  PCECLST ,  NCL) )

   idx   = MissingCode

   do n=1,NCL
      if ( (PCECLST(n)%TYPE_OF_CE.eq.EnvironmentCode) ) then
         idx = n
      endif
   end do

 end FUNCTION CLUSTER_ENV_IDX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION FIND_IDX_OF (  IPCS , IPC ) RESULT( IDX )   

   integer            ,    intent(in)    ::     IPCS(:), IPC
   integer                               ::     IDX

   integer :: n,ncl

   
   ncl=size(IPCS)

   idx   = MissingCode
    do N=1,ncl 
       if( IPCS(N) ==  IPC  ) IDX = N
    end do      

    !IF ( IDX == MissingCode ) write(*,*) "  Couldnt find a valid index "    

 end FUNCTION FIND_IDX_OF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION INDX_IS (  IS_True ) RESULT( INDS )   
   LOGICAL, dimension(:) ,  intent(in)    ::     IS_True
   integer, dimension(:) ,  pointer       ::     INDS
   integer                                ::     NumberT, NP, I,Icnt

   NP=size( IS_True )
   NumberT= COUNT( IS_True  )
   allocate(INDS(NumberT)) 
   icnt = 1
   do I = 1, NP
      if (IS_True(I) ) then
         INDS( icnt ) = I
         icnt         = icnt + 1
      endif
   end do

   


 end FUNCTION INDX_IS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION WHO_TOUCH(  PCES , PCE, NP, SELF ) RESULT( WE_TOUCH  )
   integer,                intent(in)    ::     NP
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(NP)
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE    !! ME
   logical, optional  ,    intent(in)    ::     SELF
   
   integer                               ::  L,LP,N
   logical, dimension(NP) :: WE_TOUCH
   logical :: SELF_
   real :: x0,y0,a0,r0
   real :: x1,y1,a1,r1
   integer :: i0 , i1
 
   real :: over_c, d01
 
   IF( PRESENT(SELF)) THEN
     SELF_ = SELF
   else
     SELF_ = .false.
   ENDIF

   over_c = 1.1
                              ! We overlap (mod threshhold) SOMEWHERE in the vertical
   LP=SIZE(PCE%A)

   DO N=1,NP
    WE_TOUCH(n)=.false.
    DO L=1,LP
      x0=PCE%X(l)
      y0=PCE%Y(l)
      a0=PCE%A(l)
      r0=SQRT(a0/PI)
      i0=PCE%IPC

      x1=PCES(n)%X(l)
      y1=PCES(n)%Y(l)
      a1=PCES(n)%A(l)
      r1=SQRT(a1/PI)
      i1=PCES(n)%IPC
    
      d01 = SQRT( (x0-x1)**2 + (y0-y1)**2 )

      if ( (r1+r0) > over_c*d01 ) WE_TOUCH(n)=.true.

      if ( i1 == i0 ) WE_TOUCH(n) = SELF_

    END DO
   END DO

 end FUNCTION WHO_TOUCH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION CE1_CLUSTER(  PCES , IENV, NP ) RESULT( MY_ENV_IS_IENV )
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(:)
   integer,                intent(in)    ::     IENV, NP
   !integer,   optional,    intent(in)    ::     SELFI
   
   logical, dimension(NP) :: MY_ENV_IS_IENV


     MY_ENV_IS_IENV = (   ( PCES%MY_ENV.eq.IENV )                  & 
                    .and. ( PCES%TYPE_OF_CE.ne.EnvironmentCode ) )

 end FUNCTION CE1_CLUSTER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION FIND_MY_SISTERS(  PCES , PCE, NP ) RESULT( IS_MY_SISTER )
   integer,                intent(in)    ::     NP
   type (T_ONE_PCE)   ,    intent(in)    ::     PCES(NP)
   type (T_ONE_PCE)   ,    intent(in)    ::     PCE    !! ME
   !integer,   optional,    intent(in)    ::     SELFI
   
   integer                               ::     COMMON_ENV  
   logical, dimension(NP) :: IS_MY_SISTER

     COMMON_ENV   = PCE%MY_ENV

                                                ! I am not my own sister
     IS_MY_SISTER = ((PCES%MY_ENV.eq.COMMON_ENV).and.(PCES%IPC.ne.PCE%IPC))

 end FUNCTION FIND_MY_SISTERS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE WHERE_AM_I_ij( BCKG, PCES, NP )   
  !------------------------------------------------                        
!+++++++++++++++++++++++++++++++++++++++++++++++++++
! This routine finds i's and j's of grid volumes
! containing PCES and puts into PCES(N)%{I,J}.
! Specific for logically rectangular XY grids with corners 
! indexed (1:IM+1) and (1:JM+1). Will need similar 
! routines for other grids.
!
! Note that this routine also assumes entire PCE is 
! contained within a volume - not neccesarily true
! in sheared env at high-res for example.
! (08/09/11)
!-----------------------------------------------------
   INTEGER, INTENT(IN) ::  NP
   type (T_BCKG), intent(in)                           ::     BCKG
   type (T_ONE_PCE)  , dimension(NP), intent(inout)    ::     PCES
   
    integer :: N,L,LP,KP,I,J,IM,JM,IXC,JYC
    logical,allocatable,dimension(:) :: gtxy
    real :: XWEST,XEAST,YSOUTH,YNORTH,X0,Y0

    IM = BCKG%IMS != IM
    JM = BCKG%JMS != JM
    
    allocate( GTXY( 1:MAXVAL((/IM,JM/)) ) )

    XWEST  = MINVAL( BCKG%XX )
    XEAST  = MAXVAL( BCKG%XX )
    YSOUTH = MINVAL( BCKG%YY )
    YNORTH = MAXVAL( BCKG%YY )
       
    !!write(*,*) " Domain bounds ",XWEST,XEAST,YSOUTH,YNORTH

    DO N=1,NP
       x0 = PCES(N)%XMEAN
       y0 = PCES(N)%YMEAN
       if( (X0 >= XWEST) .and. (X0 <= XEAST) .and. (Y0 >= YSOUTH) .and. ( Y0 <= YNORTH ) ) THEN 
       GTXY(:)=.FALSE.
       DO J=1,JM 
          gtxy(1:IM)=( BCKG%XX(1:IM,J) < X0 )
          IXC = COUNT( GTXY(1:IM) )   
       END DO 
       GTXY(:)=.FALSE.
       DO I=1,IM 
          gtxy(1:JM)=( BCKG%YY(I,1:JM) < Y0 )
          JYC = COUNT( GTXY(1:JM) )   
       END DO 
       
       PCES(N)%I = IXC
       PCES(N)%J = JYC

       !!write(*,*) " N, X Y ",N,PCES(N)%XMEAN,PCES(N)%YMEAN
       !!write(*,*) " N, I J ",N,PCES(N)%I,PCES(N)%J
       END IF
    end DO

    deallocate( GTXY )
    !!pause  

  end  SUBROUTINE WHERE_AM_I_ij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE CE1_INFORM
