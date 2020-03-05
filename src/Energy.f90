      SUBROUTINE ENERGY                        &
              (Tsurf,Qsurf,wind_fctr,A,B,ncell)
!
   use Block_Energy
   use Block_Network
!
   implicit none
!
   integer                     :: i,ncell
   real                        :: evap_rate,LVP, RB
   real                        :: A, B, e0
   real                        :: tsurf,qsurf,wind_fctr
   real                        :: qconv,qevap,qws,t_kelvin
   real,dimension(2),parameter :: evrte = (/1.5e-11,1.5e-11/)
   real,dimension(2)           :: q_fit, T_fit
!
      evap_rate=evrte(1)
      if (ind > 180) evap_rate=evrte(2)
      T_fit(1) = Tsurf-1.0
      T_fit(2) = Tsurf+1.0
      if (T_fit(1) .lt. 0.50) T_fit(1) = 0.50
      if (T_fit(2) .lt. 0.50) T_fit(2) = 1.00
      do i=1,2
         T_kelvin = T_fit(i) + 273.0
!
! Vapor pressure at water surface
         e0=2.1718E10*EXP(-4157.0/(T_kelvin-33.91))
!
! Bowen ratio
         RB=P_FCTR*(DBT(ncell)-T_fit(i))
!
! Latent heat of vaporization
         LVP=1.91846e06*(T_kelvin/(T_kelvin-33.91))**2
!
! Evaporative heat flux
         QEVAP=wind_fctr*RHO*LVP*evap_rate*WIND(ncell)
         if(qevap.lt.0.0) qevap=0.0
!
! Convective heat flux
         QCONV=RB*QEVAP
         QEVAP=QEVAP*(E0-EA(ncell))
! 
! Back radiation from the water surface
         QWS=280.23+6.1589*T_fit(i)
!
! Thermal energy budget for i = 1,2
         q_fit(i)=QNS(ncell)+0.97*QNA(ncell)-QWS-QEVAP+QCONV
      end do
!
! Linear relationship for equilibrium temperature and rate constant
!
      A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
      B=(T_fit(1)*q_fit(2)-T_fit(2)*q_fit(1))           &
       /(T_fit(1)-T_fit(2))
!
      qsurf=0.5*(q_fit(1)+q_fit(2))
      RETURN
      END SUBROUTINE ENERGY

