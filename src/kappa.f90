      REAL FUNCTION KAPPA_Z(nd,res_nn,T_EPI,T_HYP)
      integer, parameter          :: dp = selected_real_kind(15, 307)
      INTEGER                     :: nd,res_nn
      real                        :: dday      
      real                        :: brunt,t_epi,t_hyp,z_mod
      REAL(kind=dp)               :: RHO_EPI,RHO_HYP
!
     real(kind=dp),parameter      :: factor =-1.96e-03,factor3 = -9.8e-03
     real(kind=dp),parameter      :: PI = 3.14159265359
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
      END FUNCTION KAPPA_Z
