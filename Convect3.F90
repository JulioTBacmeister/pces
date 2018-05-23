!This stuff goes into cloude subroutine under "IF DO_TRACERS"

	   do ITR=1,ITRCR
	      XMASS0(ITR) = Dot_Product(XOI(ic:k,itr), &
			PRS(ic+1:k+1)-PRS(ic:k)) * 100.0/grav
	   end do

	   ! Calculate tracer mixing ratio in cloud after entrainment
	   ! at each level.  Weighted avg of XOI in each entrainment level.
	   XCD(K,1:itrcr) = XOI(K,1:itrcr)
	   do L=K-1,IC,-1
	      XCD(L,1:itrcr) = 1./ETA(L) * &
	     (XCD(L+1,1:itrcr)*ETA(L+1) + XOI(L,1:itrcr)*(ETA(L)-ETA(L+1)))
	   end do

	   !some unit weirdness means I have to set DT to 1.0
	   call convect(K,IC,ITRCR,ZET(IC-1:K+1),ZOL(IC-1:K), &
			PRS(IC-1:K+1),PRI(IC:K),GRAV,1.,WFN,ETA(IC:K), &
			XCD(IC:K,:), XOI(IC-1:K,:),XCU(IC:K,:))

	   do ITR=1,ITRCR
	      write(FNAME,'(I2.2)') ITR
	      FNAME = 'xcu'//trim(FNAME)//'.dat'
	      OPEN (UNIT=49, FILE=FNAME, POSITION='APPEND')
	      write(49,*) XCU(:,ITR)
	      CLOSE (UNIT=49)

	      XMASS(ITR) = Dot_Product(XOI(ic:k,itr), &
			PRS(ic+1:k+1)-PRS(ic:k)) * 100.0/grav
	      write(FNAME,'(I2.2)') ITR
	      FNAME = 'rasmass'//trim(FNAME)//'.dat'
	      OPEN (UNIT=51, FILE=FNAME, POSITION='APPEND')
	      write(51,*) XMASS0(ITR), XMASS(ITR)-XMASS0(ITR)
	      CLOSE (UNIT=51)
	   end do
           WFN     = WFN*0.5 *1.0           !*FRICFAC*0.5







   !RCooper's functions
   FUNCTION dLagrangeDz3(z1,z2,z3,z4,q1,q2,q3,q4,z)
   !INPUT:  4 z-coordinates and corresponding tracer values q, and
   !	    the value z at which we want to evaluate
   !OUTPUT: Fits a 3rd order Lagrangian polynomial to the 4 points
   !	    and calculates the derivative at the point in question

      REAL             :: dLagrangeDz3
      REAL, INTENT(IN) :: z1,z2,z3,z4, q1,q2,q3,q4, z

      dLagrangeDz3 = &
         q1 * (3*z**2. - 2*(z2+z3+z4)*z + (z2*z3+z2*z4+z3*z4)) / &
	 ((z1-z2)*(z1-z3)*(z1-z4)) + &
	 q2 * (3*z**2. - 2*(z1+z3+z4)*z + (z1*z3+z1*z4+z3*z4)) / &
	 ((z2-z1)*(z2-z3)*(z2-z4)) + &
	 q3 * (3*z**2. - 2*(z1+z2+z4)*z + (z1*z2+z1*z4+z2*z4)) / &
	 ((z3-z1)*(z3-z2)*(z3-z4)) + &
	 q4 * (3*z**2. - 2*(z1+z2+z3)*z + (z1*z2+z1*z3+z2*z3)) / &
	 ((z4-z1)*(z4-z2)*(z4-z3))

   RETURN
   END FUNCTION dLagrangeDz3


   FUNCTION dLagrangeDz2(z1,z2,z3,q1,q2,q3,z)
   ! INPUT:  3 z-coordinates and corresponding tracer values q, and
   !	     the value z at which we want to evaluate
   ! OUTPUT: Fits a 2nd order Lagrangian polynomial to the 3 points
   !	     and calculates the derivative at the point in question

      REAL		:: dLagrangeDz2
      REAL, INTENT(IN)	:: z1,z2,z3, q1,q2,q3, z

      dLagrangeDz2 = &
	 q1 * (2*z - (z2+z3)) / ((z1-z2)*(z1-z3)) + &
	 q2 * (2*z - (z1+z3)) / ((z2-z1)*(z2-z3)) + &
	 q3 * (2*z - (z1+z2)) / ((z3-z1)*(z3-z2))

      RETURN
   END FUNCTION dLagrangeDz2


   FUNCTION LPoly2(z1,z2,z3, q1,q2,q3, z)
   ! INPUT:  3 z-coordinates and corresponding tracer values q, and
   !	     the value z at which we want to evaluate
   ! OUTPUT: Fits a 2nd order Lagrangian polynomial to the 3 points
   !	     and calculates the value at the point in question
      REAL		 :: LPoly2
      REAL, INTENT(IN)	 :: z1,z2,z3, q1,q2,q3, z

      LPoly2 = &
         q1 * (z-z2)*(z-z3) / ((z1-z2)*(z1-z3)) + &
	 q2 * (z-z1)*(z-z3) / ((z2-z1)*(z2-z3)) + &
	 q3 * (z-z1)*(z-z2) / ((z3-z1)*(z3-z2))
      RETURN
   END FUNCTION LPoly2


   SUBROUTINE LPoly2co(z1,z2,z3, q1,q2,q3, a,b,c)
   ! INPUT:  3 z-coordinates and the corresponding values q, and 3
   !	     variables to hold thte coefficients
   ! RESULT: Fits a 2nd order Lagrangian polynomial to the 3 points
   !	     and stores the coefficients in the 3 output variables
   !	     i.e. the fit polynomial is q = a*z^2 + b*z + c

      REAL, INTENT(IN)	 :: z1,z2,z3, q1,q2,q3
      REAL, INTENT(OUT)	 :: a,b,c

      a = q1 / ((z1-z2)*(z1-z3)) + &
	  q2 / ((z2-z1)*(z2-z3)) + &
	  q3 / ((z3-z1)*(z3-z2))

      b = -1.*(q1*(z2+z3) / ((z1-z2)*(z1-z3)) + &
	       q2*(z1+z3) / ((z2-z1)*(z2-z3)) + &
	       q3*(z1+z2) / ((z3-z1)*(z3-z2)))

      c = q1*(z2*z3) / ((z1-z2)*(z1-z3)) + &
	  q2*(z1*z3) / ((z2-z1)*(z2-z3)) + &
	  q3*(z1*z2) / ((z3-z1)*(z3-z2))

      RETURN
   END SUBROUTINE LPoly2co


   FUNCTION antiDerivL2(z1,z2,z3, q1,q2,q3, a,b)
   ! INPUT:  3 z-coordinates and the corresponding values q, and
   !	     the 2 limits of integration
   ! OUTPUT: Fits a 2nd order Lagrangian polynomial to the 3 points
   !	     and returns the integral from a to b

      REAL		 :: antiDerivL2, evalA, evalB
      REAL, INTENT(IN)	 :: z1,z2,z3, q1,q2,q3, a,b

      evalA = &
	 q1 * ((a**3.)/3. - (a**2.)*(z2+z3)/2. + a*(z2*z3)) / &
	      ((z1-z2)*(z1-z3)) + &
	 q2 * ((a**3.)/3. - (a**2.)*(z1+z3)/2. + a*(z1*z3)) / &
	      ((z2-z1)*(z2-z3)) + &
	 q3 * ((a**3.)/3. - (a**2.)*(z1+z2)/2. + a*(z1*z2)) / &
	      ((z3-z1)*(z3-z2))

      evalB = &
	 q1 * ((b**3.)/3. - (b**2.)*(z2+z3)/2. + b*(z2*z3)) / &
	      ((z1-z2)*(z1-z3)) + &
	 q2 * ((b**3.)/3. - (b**2.)*(z1+z3)/2. + b*(z1*z3)) / &
	      ((z2-z1)*(z2-z3)) + &
	 q3 * ((b**3.)/3. - (b**2.)*(z1+z2)/2. + b*(z1*z2)) / &
	      ((z3-z1)*(z3-z2))

      antiDerivL2 = evalB - evalA
      RETURN
   END FUNCTION antiDerivL2

   FUNCTION antiDerivL1(z1,z2, q1,q2, a,b)
   ! INPUT:  2 z-coordinates and the corresponding values q, and
   !	     the limits of integration
   ! OUTPUT: Fits a 1st order Lagrangian polynomial (i.e. a line)
   !	     to the 2 points and returns the integral from a to b

      REAL		 :: antiDerivL1, m, intersect, evalA, evalB
      REAL, INTENT(IN)	 :: z1,z2, q1,q2, a,b

      m = (q1-q2)/(z1-z2)
      intersect   = (q1*z2 - q2*z1)/(z2-z1)
      evalA       = m/2.*a**2. + intersect*a
      evalB       = m/2.*b**2. + intersect*b
      antiderivL1 = evalB - evalA
      RETURN
   END FUNCTION antiDerivL1


   FUNCTION zprime(a,b,c, dt, h, incloud)
   ! INPUT:  The 3 coefficients of a 2nd order polynomial (parabola)
   !	     describing velocity as a function of height,
   !	     the timestep dt, the starting height, and whether we are
   !	     in a cloud (so air is moving up) or not (air subsiding)
   ! OUTPUT: The original height at time "now"-dt of the air which is
   !	     now at height h.  So the integral of air from h to zprime
   !	     (or zprime to h if we're in the cloud) is all the air/m^2
   !	     that passed through boundary at height h in one timestep.

      REAL		   :: zprime
      REAL, INTENT(IN)	   :: a,b,c, dt, h
      LOGICAL, INTENT(IN)  :: incloud
      REAL		   :: temp, temp2

      if (a.EQ.0.) then
	 if (b.EQ.0) then
	    if(incloud) then
	       zprime = h - c*dt
	    else
	       zprime = h + c*dt
	    end if
	 else
	   if(incloud) then
	      zprime = ((b*h+c)/b / EXP(b*dt)) - c/b
	   else
	      zprime = ((b*h+c)/b * EXP(b*dt)) - c/b
	   end if
	 end if

      else !a is not 0
         temp = 4.*a*c - b**2.
         if (temp .GT. 0.) then
            if (.NOT.incloud) then
	       zprime = (SQRT(temp) * &
	TAN(ATAN((2*a*h+b)/SQRT(temp)) + dt*SQRT(temp)/2.) - b)/(2*a) 
            else
               zprime = (SQRT(temp) * &
	TAN(ATAN((2*a*h+b)/SQRT(temp)) - dt*SQRT(temp)/2.) - b)/(2*a)
             end if

         else if (temp .LT. 0.) then
	    temp = -1.0*temp
	    if (.NOT.incloud) then
	       temp2 = EXP(dt*SQRT(temp)) * &
	          ((2*a*h+b-SQRT(temp)) / (2*a*h+b+SQRT(temp)))
	    else
	       temp2 = EXP(dt*SQRT(temp)*(-1.)) * &
	          ((2*a*h+b-SQRT(temp)) / (2*a*h+b+SQRT(temp)))
	    end if
	    zprime = ((b+SQRT(temp))*temp2 - b + SQRT(temp)) / &
		     (2*a*(1-temp2))

         else !4ac - b^2 == 0
            if(.NOT.incloud) then
	       zprime = 1./(a*(2./(2.*a*h+b) - dt)) - b/(2.*a)
	    else
	       zprime = 1./(a*(2./(2.*a*h+b) + dt)) - b/(2.*a)
	    end if
         end if
      end if
      RETURN
   END FUNCTION zprime


   SUBROUTINE convect(K,IC,ITRCR,ZET,ZOL,PRS,PRI,GRAV,dt, &
		      WFN,ETA,XCD,XOI,XCU)
   !This subroutine actually performs the 3rd order, mass-conserving
   !flux-based convection of tracer.  The steps are:
   ! 1) Fit parabolas to estimate the air density at layers IC-1:K and
   !     at edges IC-1:K+1, using     dens = -(1/g)*dP/dZ
   ! 2) Calculate the mass flux/m^2 at each edge IC:K by flux=wfn*eta,
   !	 and linearly extrapolate down to edge K+1
   ! 3) Fit a parabola to the wind field (u=flux/density) at each
   !	 edge IC+1:K, and use this to calculate from how far up in
   !     the atmosphere air traveled through this edge in a timestep
   !	 and from how far down in the cloud air traveled through
   !	 this edge in a timestep.
   ! 4) Use the above values as limits of integration to calculate
   !	 how much mass/m^2 traveled through each boundary in the
   !	 past timestep for edges IC+1:K.  The total mass through IC
   !	 and K+1 is set to 0.
   ! 5) Calculate the tendencies at layers IC:K as
   !	 delta(density*mixing ratio) = (mass in - mass out)/(dZ)
   !     delta(mixing ratio) = g/dP * (mass in - mass out)
   !
   ! CAUTION: In order to fit parabolas, etc. I had to extend ras's
   !	      definitions of variables beyond ICMIN to ICCALC, which
   !	      I define to be ICMIN-1.  I also extended the calculation
   !	      of ZET and ZOL.  Fortran doesn't always seem to warn you
   !	      if you call a variable beyond where it is explicitly
   !	      defined, so be careful!
   !
   ! ASSUMPTIONS: I am assuming here that WFN is given in units of
   !		  mass flux: kg/(s*m^2).  I'm not sure it is though.
   !		  Also, I'm not sure how to handle the timestep dt.
   !		  Should dt be 1 (as in timesteps) or 1800 (as in
   !		  seconds in each half-hour timestep)?

      INTEGER, INTENT(IN   )	:: K, IC, ITRCR
      REAL,    INTENT(IN   )	:: ZET(IC-1:K+1), ZOL(IC-1:K), PRS(IC-1:K+1)
      REAL,    INTENT(IN   )	:: PRI(IC:K), GRAV, dt
      REAL,    INTENT(IN   )	:: WFN, ETA(IC:K), XCD(IC:K,ITRCR)
      REAL,    INTENT(INOUT)	:: XOI(IC-1:K,ITRCR)
      REAL,    INTENT(  OUT)	:: XCU(IC:K,ITRCR)

      REAL			:: dens_E(IC-1:K+1), dens_L(IC-1:K)
      REAL			:: a, b, c, zp, zpc
      INTEGER			:: ITR
      REAL			:: flux(IC:K+1)
      REAL    :: totalflux(IC:K+1,ITRCR), totalfluxc(IC:K+1,ITRCR)



      do L=IC,K
         dens_E(L) = -100.0/GRAV*dLagrangeDz2(ZET(L+1),ZET(L),ZET(L-1), &
				PRS(L+1),PRS(L),PRS(L-1), ZET(L))
         dens_L(L) = -100.0/GRAV*dLagrangeDz2(ZET(L+1),ZET(L),ZET(L-1), &
				PRS(L+1),PRS(L),PRS(L-1), ZOL(L))
      end do
      !ZET is only declared down to K+1 and up to IC
 !   !  dens_E(IC-1) = -100.0/GRAV*dLagrangeDz2(ZET(IC+1),ZET(IC),ZET(IC-1), &
 !   !			       PRS(IC+1),PRS(IC),PRS(IC-1), ZET(IC-1))
      dens_L(IC-1) = -100.0/GRAV*dLagrangeDz2(ZET(IC+1),ZET(IC),ZET(IC-1), &
			       PRS(IC+1),PRS(IC),PRS(IC-1), ZOL(IC-1))
      dens_E(K+1) = -100.0/GRAV*dLagrangeDz2(ZET(K+1),ZET(K),ZET(K-1), &
			       PRS(K+1),PRS(K),PRS(K-1), ZET(K+1))

      OPEN(UNIT=77,FILE='density.dat',POSITION='APPEND')
      write(77,*) dens_E(IC:K+1)
      CLOSE(UNIT=77)

      ! How do I get flux in kg/(s*m^2) from eta, which seems
      ! to be unitless and wfn, which may be J/kg per timestep?
      flux(IC:K) = WFN*ETA(IC:K)
     ! flux(IC) = 0.
      flux(K+1) = flux(K) - &
	  (flux(K-1)-flux(K))*((ZET(K)-ZET(K+1))/(ZET(K-1)-ZET(K)))

write(*,*) 'wind:',flux(IC:K+1)/dens_E(IC:K+1)

      do L=IC+1,K
	 !Fit a parabola to wind field u=flux/density
	 call Lpoly2co(ZET(L+1),ZET(L),ZET(L-1), & 
	      flux(L+1)/dens_E(L+1), flux(L)/dens_E(L), & 
	      flux(L-1)/dens_E(L-1), a,b,c)

	 !Calculate how far to integrate tracer to get mass through
	 !a boundary in one timestep
	 zp = zprime(a,b,c, dt, ZET(L), .false.)
	 zpc = zprime(a,b,c, dt, ZET(L), .true.)

!	! OPEN (UNIT=77, FILE='z.dat',POSITION='APPEND')
!	! write(77,*) ZET(L),zp,zpc
!	! CLOSE(UNIT=77)

	 do itr=1,itrcr
	!   if(L.EQ.IC+1) then !We don't have XOI or ZOL above IC
	   !Be careful of undefined quantities, but don't compensate
	   !by looking downwind because this causes instability!
	!    totalflux(L,itr) = antiderivL2(ZOL(L  ),ZOL(L),ZOL(L-1), &
	!	       dens_L(L  ) * XOI(L  ,itr), &
	!	       dens_L(L  ) * XOI(L ,itr), &
	!	       dens_L(L-1) * XOI(L-1,itr), ZET(L), zp)
	!     totalflux(L,itr) = antiderivL1(ZOL(L),ZOL(L-1), &
	!	       dens_L(L  ) * XOI(L  ,itr), &
	!	       dens_L(L-1) * XOI(L-1,itr), ZET(L), zp)
	!      totalflux(L,itr) = dens_L(L-1)*XOI(L-1,itr) * (zp-ZET(L))
	!   else
	    totalflux(L,itr) = antiderivL2(ZOL(L),ZOL(L-1),ZOL(L-2), &
		       dens_L(L  ) * XOI(L  ,itr), &
		       dens_L(L-1) * XOI(L-1,itr), &
		       dens_L(L-2) * XOI(L-2,itr), ZET(L), zp)
	!   end if

	   if(L.EQ.K) then !We don't have XCD or ZOL below cloud base K
	    totalfluxc(L,itr) = antiderivL2(ZOL(L),ZOL(L-1),ZOL(L-2), &
		       dens_L(L  ) * XCD(L  ,itr), &
		       dens_L(L-1) * XCD(L-1,itr), &
		       dens_L(L-2) * XCD(L-2,itr), zpc, ZET(L))
	   else
	    totalfluxc(L,itr) = antiderivL2(ZOL(L+1),ZOL(L),ZOL(L-1), &
		       dens_L(L+1) * XCD(L+1,itr), &
		       dens_L(L  ) * XCD(L  ,itr), &
		       dens_L(L-1) * XCD(L-1,itr), zpc, ZET(L))
	   end if
	 end do
      end do

      !calculate totalflux(IC+1,itr) separately, since only have XOI
      !and dens_L up to IC.  Use lower values for parabola fit.
      !Hopefully the flux through IC+1 doesn't come from too high up
   !   call Lpoly2co(ZET(IC+2),ZET(IC+1),ZET(IC), & 
!	      flux(IC+2)/dens_E(IC+2), flux(IC+1)/dens_E(IC+1), & 
!	      flux(IC)/dens_E(IC), a,b,c)

!      zp = zprime(a,b,c, dt, ZET(IC+1), .false.)
!      zpc = zprime(a,b,c, dt, ZET(IC+1), .true.)

!      do itr=1,itrcr
!	 totalflux(IC+1,itr) = antiderivL2(ZOL(IC+2),ZOL(IC+1),ZOL(IC), &
!		       dens_L(IC+2) * XOI(IC+2,itr), &
!		       dens_L(IC+1) * XOI(IC+1,itr), &
!		       dens_L(IC  ) * XOI(IC  ,itr), ZET(IC+1), zp)

!	 totalfluxc(IC+1,itr) = antiderivL2(ZOL(IC+2),ZOL(IC+1),ZOL(IC), &
!		       dens_L(IC+2) * XCD(IC+2,itr), &
!		       dens_L(IC+1) * XCD(IC+1,itr), &
!		       dens_L(IC  ) * XCD(IC  ,itr), zpc, ZET(IC+1))
!      end do

      totalflux(K+1,1:itrcr)=0.
      totalflux(IC ,1:itrcr)=0.
      totalfluxc(K+1,1:itrcr)=0.
      totalfluxc(IC ,1:itrcr)=0.

      OPEN(UNIT=78,FILE='tflux.dat',POSITION='APPEND')
      write(78,*) totalflux
      CLOSE(UNIT=78)
      OPEN(UNIT=78,FILE='tcflux.dat',POSITION='APPEND')
      write(78,*) totalfluxc
      CLOSE(UNIT=78)

      do ITR=1,ITRCR
      do L=IC,K
	 XCU(L,itr) = grav*PRI(L) * &
	       (( totalflux(L,itr)- totalflux(L+1,itr)) - &
	        (totalfluxc(L,itr)-totalfluxc(L+1,itr)))
      end do
      end do

      ! Update tracer
      XOI(IC:K,1:ITRCR) = XOI(IC:K,1:ITRCR) + XCU(IC:K,1:ITRCR)
      RETURN
   END SUBROUTINE convect
