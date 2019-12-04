module functions
!
implicit none
!
contains
!
!   Function to calculate the number of simulation days
!
      function nodays(jtime,jy0)
      dimension ndmo(12)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
      jy=(jtime/10000)
      jrem=(jtime-jy*10000)
      jm=jrem/100
      jd=jrem-jm*100
      ny=jy-jy0
      nodays=365*ny+ndmo(jm)+jd
      return
      end
      function ndate(n,ny0)
      dimension ndmo(13)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334,365/
      ny=(n-1)/365
      njul=n-ny*365
      nm=1
 10   continue
      if(ndmo(nm+1).ge.njul) go to 50
      nm=nm+1
      go to 10
 50   continue
      nday=njul-ndmo(nm)
      nmon=100*nm
      nyear=10000*(ny0+ny)
      ndate=nyear+nmon+nday
!
      end function

!	Third-order polynomial interpolation using Lagrange
!     polynomials.  FUNCTION is SUBROUTINE POLINT from
!     Numerial Recipes
!
      FUNCTION tntrp(XA,YA,X,n)
!      PARAMETER (N=4)
!      DIMENSION XA(N),YA(N),C(N),D(N)
      DIMENSION XA(4),YA(4),C(4),D(4)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) DEN=0.001
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
	  tntrp=y
!
      END FUNCTION
!
      INTEGER FUNCTION Julian (YEAR,MONTH,DAY)
!
!---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!   DATE (YEAR,MONTH,DAY).
!
      INTEGER YEAR,MONTH,DAY,I,J,K
!

      I= YEAR
      J= MONTH
      K= DAY
!

      Julian=
     1   K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)
     2  /12-3*((I+4900+(J-14)/12)/100)/4
!

      END FUNCTION
!
      SUBROUTINE DENSE(DENSITY1,DENSITY2,T11,T22)
      REAL*8 DENSITY1,DENSITY2,T1,T2
         T1 = T11
         T2 = T22
!
          DENSITY1 = 6.793952E-02*T1 -9.095290E-03*(T1**2)
     &             + 1.001685E-04*(T1**3) - 1.120083E-06*(T1**4)
     &             + 6.536332E-09*(T1**5)
          DENSITY2 = 6.793952E-02*T2 -9.095290E-03*(T2**2)
     &             + 1.001685E-04*(T2**3) - 1.120083E-06*(T2**4)
     &             + 6.536332E-09*(T2**5)
!
      END SUBROUTINE  
!
      REAL FUNCTION KAPPA_Z(nd,res_nn,T_EPI,T_HYP)
      INTEGER res_nn
      REAL*8 RHO_EPI,RHO_HYP
!
!   
        PARAMETER (FACTOR = -1.96E-03,FACTOR3 = -9.8E-03,PI = 3.4159)
!
! FACTOR = -G/(RHO_H2O*DEPTH_APPROX)
!      G = 9.8 meters**2/sec; RHO_H2O = 1000.0 kg/meter**3; DEPTH_APPROX = 5.0 meters
!
        CALL DENSE (RHO_EPI,RHO_HYP,T_EPI,T_HYP)
        BRUNT = FACTOR*(RHO_EPI-RHO_HYP)
        if (res_nn .eq. 3) then 
          z_mod = 9.14*SIN(2.*PI*(DDAY-180.)/365.)
          BRUNT = FACTOR3*(RHO_EPI-RHO_HYP)*z_mod
        end if
        BRUNT = AMIN1(BRUNT,1.0E-04)
        BRUNT = AMAX1(BRUNT,1.0E-02)
!
!        KAPPA_Z = 5.6247E-10*(BRUNT**(-0.5028))
!        KAPPA_Z = 5.6247E-09/BRUNT
        KAPPA_Z = 5.6247E-09*(BRUNT**(-0.5028))
!
! Estimate the vertical eddy diffusivity using the method of Quay et al (1980)
!
      END FUNCTION
!
      INTEGER FUNCTION NLAYER(T11,T22)
      REAL*8 DENSITY1,DENSITY2,T1,T2
         T1 = T11
         T2 = T22
!
          DENSITY1 = 6.793952E-02*T1 -9.095290E-03*(T1**2)
     &             + 1.001685E-04*(T1**3) - 1.120083E-06*(T**4)
     &             + 6.536332E-09*(T1**5)
          DENSITY2 = 6.793952E-02*T2 -9.095290E-03*(T2**2)
     &             + 1.001685E-04*(T2**3) - 1.120083E-06*(T2**4)
     &             + 6.536332E-09*(T2**5)
         NLAYER = 1
         IF (DENSITY2 .GT. DENSITY1) NLAYER = 2
!
      END FUNCTION 
!
end module      
      








