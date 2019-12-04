module Block_Energy
!
!
!   Energy budget variables
!
!   Incoming short wave radiation, Watts/m**2
!
    real, dimension(:), allocatable :: QNS
!
!   Incoming atmospheric radiation, Watts/m**2
!
    real, dimension(:), allocatable :: QNA
!
!   Air temperature at surface, deg. C
!
    real, dimension(:), allocatable :: DBT
!  
!   Wind speed, m/sec
!
    real, dimension(:), allocatable :: WIND
!
!   Vapor pressure of air at surface, MB
!
    real, dimension(:), allocatable :: EA
!
!   Air pressure at surface, mb
!
    real, dimension(:), allocatable :: PRESS 

!
    real, dimension (:), allocatable::mu,alphamu,beta,gmma,smooth_param
!
!   Some important constants
!
      integer          :: nwpd
      real,parameter   :: PI = 3.14159 
      real,parameter   :: EVRATE=1.5e-11        ! Lake Hefner coefficient, 1/meters
      real,parameter   :: P_FCTR=64.0,RHO=1000. ! Bowen ratio, water density
      real,parameter   :: RFAC=4.184e6          ! rho/Cp kg/meter**3/Joules/kg/Deg K    
end module Block_Energy  
