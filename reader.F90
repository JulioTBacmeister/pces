SUBROUTINE ARM2(FILENAME, NT, NLEVEL, NLAYR,                 &
                PREF_MODEL_E,                                &
                TIME,                                        &
		YY,                                          &
		MO,                                          &
                DD,                                          &
                HH,                                          &
                MM,                                          &
                PCP_OBS,                                     &
		TS_AIR,                                      &
		TG_SOIL,                                     &
                TSKIN,                                       &
                QSFCAIR,                                     &
                QSKIN,                                       &
                PSFC,                                        &
		LHF,                                         & 
		SHF,                                         & 
                tt,                                          &
                qq,                                          &
                uu,                                          &
                vv,                                          &
                T_H_adv,                                     &
                T_V_adv,                                     &
                Q_H_adv,                                     &
                Q_V_adv,                                     &           
                Q1,                                          &
                Q2,                                          &
                P_MODEL_E                                    )

use NUMERICAL_UTILITIES, only : QSAT, DQSAT


        IMPLICIT NONE
! INPUT        XXXXX
! NT:  time slice number of ARM data
! NLEVEL: vertical levels of ARM data
! NLAYR:  vertical layer of Single column model
! FILENAME: file name of ARM data
!	
	INTEGER,  INTENT(IN) :: Nt, NLEVEL, NLAYR
	CHARACTER(LEN=*), INTENT(IN) :: FILENAME
	REAL,  DIMENSION(0:NLAYR),     INTENT(IN   ) :: PREF_MODEL_E

! OUTPUT
! single layer
	real, dimension(nt), intent(out) :: time, yy, mo, dd, hh, mm
	real, dimension(nt), intent(out) :: pcp_obs, tskin, qsfcair,&
	qskin, psfc, shf, lhf, ts_air, tg_soil
! multiple-layer
	real, dimension(nt,nlayr),intent(out):: tt, qq, uu, vv, &
	T_H_adv, T_V_adv,    &
        Q_H_adv, Q_V_adv, Q1, Q2
	real, dimension(nt,0:nlayr), intent(out):: P_MODEL_E
	real, dimension(nlayr) :: p_model
! temporiy array
	integer :: IVAR , IVAR_sfc, I, J, K 
	real, allocatable, dimension(:,:,:) :: TMP, DV
	real, allocatable, dimension(:,:)   :: TMP_sfc
	real, allocatable, dimension(:)     :: P
! working variables
	real, parameter :: ptop = 10.
	real :: pres, temd, PUPP, PDWN, PRGAT, PUPPK, PDWNK 
	real :: pmass, PRESK, TEM, TEMU
	integer :: IGD, IGU, IGTLEV, ITOP, UNIT
	
	allocate (P(NLEVEL))
!
        UNIT = 122
        OPEN( UNIT=UNIT, FILE=FILENAME, FORM =  "FORMATTED" )
	!!UNIT = GETFILE(FILENAME, form = "formatted")
!
	READ(UNIT,*)
	READ(UNIT,*)
	READ(UNIT,*)
	READ(UNIT,*)
	READ(UNIT,*)
	READ(UNIT,12) P
	READ(UNIT,*)
	READ(UNIT,12) TIME
	READ(UNIT,*)
	READ(UNIT,12) YY
	READ(UNIT,*)
	READ(UNIT,12) MO
	READ(UNIT,*)
	READ(UNIT,12) DD
	READ(UNIT,*)
	READ(UNIT,12) HH
	READ(UNIT,*)
	READ(UNIT,12) MM
	READ(UNIT,*)
	READ(UNIT,*)  IVAR
	allocate (TMP(IVAR,NLEVEL,NT))
	allocate (DV (IVAR,NLAYR, NT))
	do i = 1, IVAR
	READ(UNIT,*) 
	do j = 1, NLEVEL
	read(UNIT,12) (TMP(I,J,K),K=1,NT)
	END DO
	END DO
	READ(UNIT,*)
	READ(UNIT,*) IVAR_sfc
	READ(UNIT,*)
	allocate (TMP_sfc(IVAR_sfc,NT))
	DO I = 1, NT
	READ(UNIT,12) (TMP_sfc(K,i),K=1,IVAR_sfc)
	END DO
 12     format(5e15.7)
 	!!! CALL FREE_FILE(UNIT)
!  Make the fields for SCM
!  1.  SURFACE 
	
        do i = 1, NT
	PCP_OBS(i) = TMP_sfc(5,i)*24.0          ! mm/day
	PSFC(i)    = TMP_sfc(4,i)               ! mb
	TSKIN(i)   = TMP_sfc(3,i)               ! K
	end do
	TG_SOIL    = TSKIN
	TS_AIR     = TSKIN
	do i = 1, NT
	pmass = 100.*(PSFC(i) - ptop)/9.81
	SHF(i)    = TMP_sfc(1,i)
	LHF(I)    = TMP_sfc(2,i)
	end do
	TMP(6,:,:) = TMP(6,:,:)/3600.0          ! HTA k/s
	TMP(7,:,:) = TMP(7,:,:)/3600.0
	TMP(8,:,:) = TMP(8,:,:)/3600.0
	TMP(9,:,:) = TMP(9,:,:)/3600.0
	do i = 1, NT
	do j = 1, NLEVEL
	TMP(1,J,I) = TMP(1,J,I)*((1000.0/P(J))**0.286)
	END DO
	END DO

!       Interpolate vertically
	do i = 1, NT
	p_model_e(i,0:nlayr) = pref_model_e(0:nlayr)*psfc(i)*100./pref_model_e(nlayr)  ! output var in Pa
	p_model(1:nlayr) = (p_model_e(i,0:nlayr-1) + p_model_e(i,1:nlayr))*0.5/100.  ! Local var in hPa

! Note psfc is output directly in hPa, but is never used in SCM
	
	do k = 1, nlayr 
	PRES = P_MODEL(K)
	IF (PRES > P(1) ) THEN
	IGD  = 1
	IGU  = 1
	PUPP = 965.0
	PDWN = PSFC(I)
	GOTO 205
	END IF
	IF (PRES < P(NLEVEL)) THEN
	IGD = NLEVEL
	IGU = NLEVEL
	PUPP = PRES
	PDWN = P(NLEVEL)
	GOTO 205
	END IF
	PUPP = P(NLEVEL)
	PDWN = P(1)
	IGU  = NLEVEL
	IGD  = 2
	DO 204 IGTLEV = 2, NLEVEL
	PRGAT = P(IGTLEV)
	IF (PRGAT <= PRES ) THEN
	PUPP = PRGAT
	IGU  = IGTLEV
	GO TO 900
	END IF
	IF (PRGAT > PRES ) THEN
	PDWN = PRGAT
	IGD  = IGTLEV
	END IF
 204	CONTINUE
 900	CONTINUE
 205	CONTINUE
 	PUPPK = PUPP**0.286
	PDWNK = PDWN**0.286
	PRESK = PRES**0.286

	TEM   = PUPPK - PDWNK
	TEMU  = (PUPPK - PRESK)/TEM
	TEMD  = 1.0 - TEMU
	DO J = 1, IVAR
	DV(J,K,I) = TMP(J,IGU,I)*TEMD + TMP(J,IGD,I)*TEMU
	END DO
	END DO
 	END DO
	do K = 1, NLAYR
	   do I = 1, NT
	   tt(i,k) = dv(1,k,I)*((p_model(k)/1000.0)**0.286)
	   qq(i,k) = dv(2,K,I)*1000. ! kg/kg => g/kg
	   uu(i,k) = dv(3,K,I)
	   vv(i,k) = dv(4,K,I)
	   T_H_ADV(i,k) = -dv(6,K,I)
	   T_V_ADV(i,k) = -dv(7,K,I)
	   Q_H_ADV(i,k) = -dv(8,K,I)
	   Q_V_ADV(i,k) = -dv(9,K,I)
	   END DO
        END DO
	Q1 = -999.0
	Q2 = -999.0
	do i = 1, NT
	ITOP = 1
	do k = 1, nlayr-1
	if((p_model(k)<=minval(p)) .AND. (p_model(k+1)>minval(p))) ITOP = K
	END DO
	do k = 1, ITOP
	tt(i,k) = tt(i,itop)
	end do
	end do
	do i = 1, NT
	QSKIN(I) = 1000.*QSAT (TSKIN(i),PSFC(i))
	QSFCAIR(i) = QQ(i,nlayr)
	END DO
	END SUBROUTINE ARM2
