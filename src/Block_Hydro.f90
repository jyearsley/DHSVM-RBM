!
! Module for hydraulic characteristics and water quality constituents of the basin
!
module Block_Hydro
    integer, dimension(200):: no_dt,nstrt_elm
    real, dimension(200)   :: dt_part,x_part
!
! Hydraulic characteristics
!
    real, dimension(:),     allocatable  :: depth,D_a,D_b,D_min
    real, dimension(:),     allocatable  :: elev
    real, dimension(:),     allocatable  :: width
    real, dimension(:),     allocatable  :: u,U_a,U_b,U_min
    real, dimension(:),     allocatable  :: dt
    real, dimension(:),     allocatable  :: dx
    real, dimension(:,:,:), allocatable  :: x_dist
!
! Reservoir characteristics
!
    real, dimension(:), allocatable      :: Kappa
    real, dimension(:,:,:), allocatable  :: A_surf,Volume
!
! Flows
!
    real, dimension(:),     allocatable  :: Q_in
    real, dimension(:,:),   allocatable  :: Q_in_seg
    real, dimension(:),     allocatable  :: Q_trib
    real, dimension(:),     allocatable  :: Q_out
    real, dimension(:,:),   allocatable  :: Q_out_seg
    real, dimension(:),     allocatable  :: Q_diff
    real, dimension(:,:),   allocatable  :: Q_nps


!  Parameters
!
    real,parameter          :: cuft_cum = 0.028318 ! Cubic feet to cubic meters
!
end module Block_Hydro
