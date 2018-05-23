!
!     file pois3d.f
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     SUBROUTINE POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,
!    +                   MDIMF,F,IERROR)
!
!
! DIMENSION OF           A(N), B(N), C(N), F(LDIMF,MDIMF,N)
! ARGUMENTS
!
! LATEST REVISION        June 2004
!
! PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
!                        FOR UNKNOWN X VALUES, WHERE I=1,2,...,L,
!                        J=1,2,...,M, AND K=1,2,...,N
!
!                        C1*(X(I-1,J,K) -2.*X(I,J,K) +X(I+1,J,K)) +
!                        C2*(X(I,J-1,K) -2.*X(I,J,K) +X(I,J+1,K)) +
!                        A(K)*X(I,J,K-1) +B(K)*X(I,J,K)+ C(K)*X(I,J,K+1)
!                        = F(I,J,K)
!
!                        THE INDICES K-1 AND K+1 ARE EVALUATED MODULO N,
!                        I.E. X(I,J,0)=X(I,J,N) AND X(I,J,N+1)=X(I,J,1).
!                        THE UNKNOWNS
!                        X(0,J,K), X(L+1,J,K), X(I,0,K), AND X(I,M+1,K)
!                        ARE ASSUMED TO TAKE ON CERTAIN PRESCRIBED
!                        VALUES DESCRIBED BELOW.
!
! USAGE                  CALL POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,
!                        N,A,B,C,LDIMF,MDIMF,F,IERROR)
!
! ARGUMENTS
!
! ON INPUT
!                        LPEROD
!                          INDICATES THE VALUES THAT X(0,J,K) AND
!                          X(L+1,J,K) ARE ASSUMED TO HAVE.
!                          = 0  X(0,J,K)=X(L,J,K), X(L+1,J,K)=X(1,J,K)
!                          = 1  X(0,J,K) = 0,      X(L+1,J,K) = 0
!                          = 2  X(0,J,K)=0,        X(L+1,J,K)=X(L-1,J,K)
!                          = 3  X(0,J,K)=X(2,J,K), X(L+1,J,K)=X(L-1,J,K)
!                          = 4  X(0,J,K)=X(2,J,K), X(L+1,J,K) = 0.
!
!                        L
!                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
!                          L MUST BE AT LEAST 3.
!
!                        C1
!                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
!                          OF EQUATIONS TO BE SOLVED.
!
!                        MPEROD
!                          INDICATES THE VALUES THAT X(I,0,K) AND
!                          X(I,M+1,K) ARE ASSUMED TO HAVE.
!                          = 0  X(I,0,K)=X(I,M,K), X(I,M+1,K)=X(I,1,K)
!                          = 1  X(I,0,K)=0,        X(I,M+1,K)=0
!                          = 2  X(I,0,K)=0,        X(I,M+1,K)=X(I,M-1,K)
!                          = 3  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=X(I,M-1,K)
!                          = 4  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=0
!
!                        M
!                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
!                          M MUST BE AT LEAST 3.
!
!                        C2
!                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
!                          OF EQUATIONS TO BE SOLVED.
!
!                        NPEROD
!                          = 0  IF A(1) AND C(N) ARE NOT ZERO.
!                          = 1  IF A(1) = C(N) = 0.
!
!                        N
!                          THE NUMBER OF UNKNOWNS IN THE K-DIRECTION.
!                          N MUST BE AT LEAST 3.
!
!                        A, B, C
!                          ONE-DIMENSIONAL ARRAYS OF LENGTH N THAT
!                          SPECIFY THE COEFFICIENTS IN THE LINEAR
!                          EQUATIONS GIVEN ABOVE.
!
!                          IF NPEROD = 0 THE ARRAY ELEMENTS MUST NOT
!                          DEPEND UPON INDEX K, BUT MUST BE CONSTANT.
!                          SPECIFICALLY,THE SUBROUTINE CHECKS THE
!                          FOLLOWING CONDITION
!                            A(K) = C(1)
!                            C(K) = C(1)
!                            B(K) = B(1)
!                          FOR K=1,2,...,N.
!
!                        LDIMF
!                          THE ROW (OR FIRST) DIMENSION OF THE THREE-
!                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
!                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
!                          USED TO SPECIFY THE VARIABLE DIMENSION
!                          OF F.  LDIMF MUST BE AT LEAST L.
!
!                        MDIMF
!                          THE COLUMN (OR SECOND) DIMENSION OF THE THREE
!                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
!                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
!                          USED TO SPECIFY THE VARIABLE DIMENSION
!                          OF F.  MDIMF MUST BE AT LEAST M.
!
!                        F
!                          A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE
!                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
!                          OF EQUATIONS GIVEN ABOVE.  F MUST BE
!                          DIMENSIONED AT LEAST L X M X N.
!
! ON OUTPUT
!
!                        F
!                          CONTAINS THE SOLUTION X.
!
!                        IERROR
!                          AN ERROR FLAG THAT INDICATES INVALID INPUT
!                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
!                          SOLUTION IS NOT ATTEMPTED.
!                          = 0  NO ERROR
!                          = 1  IF LPEROD .LT. 0 OR .GT. 4
!                          = 2  IF L .LT. 3
!                          = 3  IF MPEROD .LT. 0 OR .GT. 4
!                          = 4  IF M .LT. 3
!                          = 5  IF NPEROD .LT. 0 OR .GT. 1
!                          = 6  IF N .LT. 3
!                          = 7  IF LDIMF .LT. L
!                          = 8  IF MDIMF .LT. M
!                          = 9  IF A(K) .NE. C(1) OR C(K) .NE. C(1)
!                               OR B(I) .NE.B(1) FOR SOME K=1,2,...,N.
!                          = 10 IF NPEROD = 1 AND A(1) .NE. 0
!                               OR C(N) .NE. 0
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if N,M are too large
!                               for your computer)
!
!                          SINCE THIS IS THE ONLY MEANS OF INDICATING A
!                          POSSIBLY INCORRECT CALL TO POIS3D, THE USER
!                          SHOULD TEST IERROR AFTER THE CALL.
!
! SPECIAL CONDITIONS     NONE
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED files         fish.f,comf.f,fftpack.f
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
!                        1970'S.  RELEASED ON NCAR'S PUBLI! SOFTWARE
!                        LIBRARIES IN JANUARY, 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
!
! PORTABILITY            FORTRAN 90
!
! ALGORITHM              THIS SUBROUTINE SOLVES THREE-DIMENSIONAL BLOCK
!                        TRIDIAGONAL LINEAR SYSTEMS ARISING FROM FINITE
!                        DIFFERENCE APPROXIMATIONS TO THREE-DIMENSIONAL
!                        POISSON EQUATIONS USING THE FFT PACKAGE
!                        FFTPACK WRITTEN BY PAUL SWARZTRAUBER.
!
! TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
!                        IS ROUGHLY PROPORTIONAL TO
!                          L*M*N*(LOG2(L)+LOG2(M)+5)
!                        BUT ALSO DEPENDS ON INPUT PARAMETERS LPEROD
!                        AND MPEROD.
!
! ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
!                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
!                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
!                        IN THE 'PURPOSE' SECTION WITH
!                          A(K) = C(K) = -0.5*B(K) = 1,  K=1,2,...,N
!                        AND, WHEN NPEROD = 1
!                          A(1) = C(N) = 0
!                          A(N) = C(1) = 2.
!
!                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
!                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT
!                        SIDE Y WAS COMPUTED.  USING THIS ARRAY Y
!                        SUBROUTINE POIS3D WAS CALLED TO PRODUCE AN
!                        APPROXIMATE SOLUTION Z.  RELATIVE ERROR
!
!                        E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K
!
!                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
!                        OVER I=1,2,...,L, J=1,2,...,M AND K=1,2,...,N.
!                        VALUES OF E ARE GIVEN IN THE TABLE BELOW FOR
!                        SOME TYPICAL VALUES OF L,M AND N.
!
!                        L(=M=N)   LPEROD    MPEROD       E
!                        ------    ------    ------     ------
!
!                          16        0         0        1.E-13
!                          15        1         1        4.E-13
!                          17        3         3        2.E-13
!                          32        0         0        2.E-13
!                          31        1         1        2.E-12
!                          33        3         3        7.E-13
!
! REFERENCES              NONE
! ********************************************************************
      SUBROUTINE POIS3D(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, C &
      , LDIMF, MDIMF, F, IERROR)
      USE fish
      implicit none
      !TYPE (fishworkspace) :: w
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: LPEROD
      INTEGER  :: L
      INTEGER  :: MPEROD
      INTEGER  :: M
      INTEGER  :: NPEROD
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL  :: C1
      REAL  :: C2
#ifdef GFORTRAN
      REAL  :: A(N)
      REAL  :: B(N)
      REAL  :: C(N)
      REAL  :: F(LDIMF,MDIMF,N)
#else
      REAL  :: A(*)
      REAL  :: B(*)
      REAL  :: C(*)
      REAL  :: F(LDIMF,MDIMF,*)
#endif
      REAL  :: WSREW(1500)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LP, MP, NP, K, IRWK, ICWK,istatus
      REAL, DIMENSION(6) :: SAVE
!-----------------------------------------------
      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
!
!     CHECK FOR INVALID INPUT.
!
      IERROR = 0
      IF (LP<1 .OR. LP>5) IERROR = 1
      IF (L < 3) IERROR = 2
      IF (MP<1 .OR. MP>5) IERROR = 3
      IF (M < 3) IERROR = 4
      IF (NP<1 .OR. NP>2) IERROR = 5
      IF (N < 3) IERROR = 6
      IF (LDIMF < L) IERROR = 7
      IF (MDIMF < M) IERROR = 8
      IF (NP == 1) THEN
         DO K = 1, N
            IF (A(K) /= C(1)) GO TO 102
            IF (C(K) /= C(1)) GO TO 102
            IF (B(K) /= B(1)) GO TO 102
         END DO
         GO TO 104
  102    CONTINUE
         IERROR = 9
      ENDIF
      IF (NPEROD==1 .AND. (A(1)/=0. .OR. C(N)/=0.)) IERROR = 10
! 104 IF (IERROR .NE. 0) GO TO 122
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
!     allocate required work space length (generous estimate)
      IRWK=30+L+M+2*N+MAX0(L,M,N)+7*(INT((L+1)/2)+INT((M+1)/2))
      ICWK = 0
      !write(*,*) " ABOUT TO ALLOCATFISH "
!      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
!      IF (IERROR == 20) RETURN 
             write(*,*) "IRWK", IRWK
	!if (irwk > 0) then
	!     ALLOCATE(w%rew(irwk),STAT = istatus)
	!end if
!       !allocate icwk words of complex work space
	!if (icwk > 0) then
	!     ALLOCATE(w%cxw(icwk),STAT = istatus)
	!end if

      call pois3dd(LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,  &
#ifdef GFORTRAN
                   MDIMF,F,IERROR,wsrew,IRWK)
#else
                   MDIMF,F,IERROR,wsrew)
#endif
!
!     release work space
      !CALL FISHFIN (W)
	if (irwk > 0) then
	    ! DEALLOCATE(w%rew )
	end if
!       allocate icwk words of complex work space
	if (icwk > 0) then
	    ! DEALLOCATE(w%cxw)
	end if


                   
      RETURN 
      END SUBROUTINE POIS3D


 
      SUBROUTINE POIS3DD(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, &
#ifdef GFORTRAN
       C, LDIMF, MDIMF, F, IERROR, W,IWRK)
#else
       C, LDIMF, MDIMF, F, IERROR, W)
#endif
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LPEROD
      INTEGER  :: L
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER , INTENT(IN) :: NPEROD
#ifdef GFORTRAN
      INTEGER , INTENT(IN) :: IWRK
#endif
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL  :: C1
      REAL  :: C2
#ifdef GFORTRAN
      REAL  :: A(N)
      REAL  :: B(N)
      REAL  :: C(N)
      REAL  :: F(LDIMF,MDIMF,N)
      REAL  :: W(IWRK)
#else
      REAL  :: A(*)
      REAL  :: B(*)
      REAL  :: C(*)
      REAL  :: F(LDIMF,MDIMF,*)
      REAL  :: W(*)
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LP, MP, NP, IWYRT, IWT, IWD, IWBB, IWX, IWY, NH, NHM1, &
        NODD, I, J, K, NHPK, NHMK
      REAL, DIMENSION(6) :: SAVE
!-----------------------------------------------
      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
      IWYRT = L + 1
      IWT = IWYRT + M
      IWD = IWT + MAX0(L,M,N) + 1
      IWBB = IWD + N
      IWX = IWBB + N
      IWY = IWX + 7*((L + 1)/2) + 15
      GO TO (105,114) NP
!
!     REORDER UNKNOWNS WHEN NPEROD = 0.
!
  105 CONTINUE
      NH = (N + 1)/2
      NHM1 = NH - 1
      NODD = 1
      IF (2*NH == N) NODD = 2
      DO I = 1, L
         DO J = 1, M
            DO K = 1, NHM1
               W(K) = F(I,J,NH-K) - F(I,J,K+NH)
               W(K+NH) = F(I,J,NH-K) + F(I,J,K+NH)
            END DO
            W(NH) = 2.*F(I,J,NH)
            GO TO (108,107) NODD
  107       CONTINUE
            W(N) = 2.*F(I,J,N)
  108       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      SAVE(1) = C(NHM1)
      SAVE(2) = A(NH)
      SAVE(3) = C(NH)
      SAVE(4) = B(NHM1)
      SAVE(5) = B(N)
      SAVE(6) = A(N)
      C(NHM1) = 0.
      A(NH) = 0.
      C(NH) = 2.*C(NH)
      SELECT CASE (NODD) 
      CASE DEFAULT
         B(NHM1) = B(NHM1) - A(NH-1)
         B(N) = B(N) + A(N)
      CASE (2) 
         A(N) = C(NH)
      END SELECT
  114 CONTINUE
      CALL POS3D1 (LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, W, W(IWYRT) & 
        , W(IWT), W(IWD), W(IWX), W(IWY), C1, C2, W(IWBB))
      GO TO (115,122) NP
  115 CONTINUE
      DO I = 1, L
         DO J = 1, M
            W(NH-1:NH-NHM1:(-1))=0.5*(F(I,J,NH+1:NHM1+NH)+F(I,J,:NHM1))
            W(NH+1:NHM1+NH) = 0.5*(F(I,J,NH+1:NHM1+NH)-F(I,J,:NHM1))
            W(NH) = 0.5*F(I,J,NH)
            GO TO (118,117) NODD
  117       CONTINUE
            W(N) = 0.5*F(I,J,N)
  118       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      C(NHM1) = SAVE(1)
      A(NH) = SAVE(2)
      C(NH) = SAVE(3)
      B(NHM1) = SAVE(4)
      B(N) = SAVE(5)
      A(N) = SAVE(6)
  122 CONTINUE
      RETURN 
      END SUBROUTINE POIS3DD


      SUBROUTINE POS3D1(LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, XRT, &
      YRT, T, D, WX, WY, C1, C2, BB)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LP
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: MP
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: LDIMF
      INTEGER , INTENT(IN) :: MDIMF
      REAL , INTENT(IN) :: C1
      REAL , INTENT(IN) :: C2
#ifdef GFORTRAN
      REAL  :: A(N)
      REAL , INTENT(IN) :: B(N)
      REAL  :: C(N)
      REAL , INTENT(INOUT) :: F(LDIMF,MDIMF,N)
#else
      REAL  :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL  :: C(*)
      REAL , INTENT(INOUT) :: F(LDIMF,MDIMF,1)
#endif
      REAL , INTENT(INOUT) :: XRT(*)
      REAL , INTENT(INOUT) :: YRT(*)
      REAL  :: T(*)
      REAL  :: D(*)
      REAL  :: WX(*)
      REAL  :: WY(*)
      REAL  :: BB(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LR, MR, NR, LRDEL, I, MRDEL, J, IFWRD, IS, K
      REAL :: PI, DUM, SCALX, DX, DI, SCALY, DY, DJ
!-----------------------------------------------
      PI = 4.0*ATAN(1.0)
      LR = L
      MR = M
      NR = N
!
!     GENERATE TRANSFORM ROOTS
!
      LRDEL = ((LP - 1)*(LP - 3)*(LP - 5))/3
      SCALX = LR + LRDEL
      DX = PI/(2.*SCALX)
      GO TO (108,103,101,102,101) LP
  101 CONTINUE
      DI = 0.5
      SCALX = 2.*SCALX
      GO TO 104
  102 CONTINUE
      DI = 1.0
      GO TO 104
  103 CONTINUE
      DI = 0.0
  104 CONTINUE
      DO I = 1, LR
         XRT(I) = -4.*C1*SIN((FLOAT(I) - DI)*DX)**2
      END DO
      SCALX = 2.*SCALX
      GO TO (112,106,110,107,111) LP
  106 CONTINUE
      CALL SINTI (LR, WX)
      GO TO 112
  107 CONTINUE
      CALL COSTI (LR, WX)
      GO TO 112
  108 CONTINUE
      XRT(1) = 0.
      XRT(LR) = -4.*C1
      DO I = 3, LR, 2
         XRT(I-1) = -4.*C1*SIN(FLOAT(I - 1)*DX)**2
         XRT(I) = XRT(I-1)
      END DO
      CALL RFFTI (LR, WX)
      GO TO 112
  110 CONTINUE
      CALL SINQI (LR, WX)
      GO TO 112
  111 CONTINUE
      CALL COSQI (LR, WX)
  112 CONTINUE
      MRDEL = ((MP - 1)*(MP - 3)*(MP - 5))/3
      SCALY = MR + MRDEL
      DY = PI/(2.*SCALY)
      GO TO (120,115,113,114,113) MP
  113 CONTINUE
      DJ = 0.5
      SCALY = 2.*SCALY
      GO TO 116
  114 CONTINUE
      DJ = 1.0
      GO TO 116
  115 CONTINUE
      DJ = 0.0
  116 CONTINUE
      DO J = 1, MR
         YRT(J) = -4.*C2*SIN((FLOAT(J) - DJ)*DY)**2
      END DO
      SCALY = 2.*SCALY
      GO TO (124,118,122,119,123) MP
  118 CONTINUE
      CALL SINTI (MR, WY)
      GO TO 124
  119 CONTINUE
      CALL COSTI (MR, WY)
      GO TO 124
  120 CONTINUE
      YRT(1) = 0.
      YRT(MR) = -4.*C2
      DO J = 3, MR, 2
         YRT(J-1) = -4.*C2*SIN(FLOAT(J - 1)*DY)**2
         YRT(J) = YRT(J-1)
      END DO
      CALL RFFTI (MR, WY)
      GO TO 124
  122 CONTINUE
      CALL SINQI (MR, WY)
      GO TO 124
  123 CONTINUE
      CALL COSQI (MR, WY)
  124 CONTINUE
      IFWRD = 1
      IS = 1
  125 CONTINUE
!
!     TRANSFORM X
!
      DO J=1,MR
	 DO K=1,NR
	    DO I=1,LR
               T(I) = F(I,J,K)
	    END DO
            GO TO (127,130,131,134,135),LP
  127       GO TO (128,129),IFWRD
  128       CALL RFFTF (LR,T,WX)
            GO TO 138
  129       CALL RFFTB (LR,T,WX)
            GO TO 138
  130       CALL SINT (LR,T,WX)
            GO TO 138
  131       GO TO (132,133),IFWRD
  132       CALL SINQF (LR,T,WX)
            GO TO 138
  133       CALL SINQB (LR,T,WX)
            GO TO 138
  134       CALL COST (LR,T,WX)
            GO TO 138
  135       GO TO (136,137),IFWRD
  136       CALL COSQF (LR,T,WX)
            GO TO 138
  137       CALL COSQB (LR,T,WX)
  138       CONTINUE
	    DO I=1,LR
               F(I,J,K) = T(I)
	    END DO
	 END DO
      END DO
      GO TO (142,164) IFWRD
!
!     TRANSFORM Y
!
  142 CONTINUE
      DO I=1,LR
	 DO K=1,NR
	    DO J=1,MR
               T(J) = F(I,J,K)
	    END DO
            GO TO (144,147,148,151,152),MP
  144       GO TO (145,146),IFWRD
  145       CALL RFFTF (MR,T,WY)
            GO TO 155
  146       CALL RFFTB (MR,T,WY)
            GO TO 155
  147       CALL SINT (MR,T,WY)
            GO TO 155
  148       GO TO (149,150),IFWRD
  149       CALL SINQF (MR,T,WY)
            GO TO 155
  150       CALL SINQB (MR,T,WY)
            GO TO 155
  151       CALL COST (MR,T,WY)
            GO TO 155
  152       GO TO (153,154),IFWRD
  153       CALL COSQF (MR,T,WY)
            GO TO 155
  154       CALL COSQB (MR,T,WY)
  155       CONTINUE
	    DO J=1,MR
               F(I,J,K) = T(J)
	    END DO
	 END DO
      END DO
      GO TO (159,125) IFWRD
  159 CONTINUE
      DO I = 1, LR
         DO J = 1, MR
            BB(:NR) = B(:NR) + XRT(I) + YRT(J)
            T(:NR) = F(I,J,:NR)
            CALL TRID (NR, A, BB, C, T, D)
            F(I,J,:NR) = T(:NR)
         END DO
      END DO
      IFWRD = 2
      IS = -1
      GO TO 142
  164 CONTINUE
      F(:LR,:MR,:NR) = F(:LR,:MR,:NR)/(SCALX*SCALY)
      RETURN 
      END SUBROUTINE POS3D1


      SUBROUTINE TRID(MR, A, B, C, Y, D)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: MR
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL , INTENT(INOUT) :: Y(*)
      REAL , INTENT(INOUT) :: D(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: M, MM1, I, IP
      REAL :: Z
!-----------------------------------------------
      M = MR
      MM1 = M - 1
      Z = 1./B(1)
      D(1) = C(1)*Z
      Y(1) = Y(1)*Z
      DO I = 2, MM1
         Z = 1./(B(I)-A(I)*D(I-1))
         D(I) = C(I)*Z
         Y(I) = (Y(I)-A(I)*Y(I-1))*Z
      END DO
      Z = B(M) - A(M)*D(MM1)
      IF (Z == 0.) THEN
         Y(M) = 0.
      ELSE
         Y(M) = (Y(M)-A(M)*Y(MM1))/Z
      ENDIF
      DO IP = 1, MM1
         I = M - IP
         Y(I) = Y(I) - D(I)*Y(I+1)
      END DO
      RETURN 
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 Changes
!-----------------------------------------------------------------------
      END SUBROUTINE TRID
