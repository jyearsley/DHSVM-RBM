Module Block_Network
!
! Module with stream topology variables
!

!


    integer, dimension(:), allocatable  :: last_seg,main_stem,ndelta
    integer, dimension(:), allocatable  :: ndelta,no_cells,node,no_tribs,no_wr_units
!
    integer, dimension(:,:), allocatable :: conflnce,head_cell,res_no,trib
    integer, dimension(:,:), allocatable :: no_celm,upstrm_cell,cells_wr,unit_type
    integer, dimension(:,:,:), allocatable :: segment_cell
!
    character (len=5) :: unit_type 

!
! Integer variables 
!
    integer           :: flow_cells,heat_cells
    integer           :: ndays,nyear,nreach,nwpd
    integer,parameter :: ns_max=500
    integer           :: start_year,start_month,start_day
    integer           :: end_year,end_month,end_day,end_hour
!
! Logical variables
!
    logical           :: IS_HEAD 
!
! Real variables
!
    real                           :: delta_n,dt_comp
    real                           :: day,day_fract,hr_fract,sin_incr,year,xd_year
    real                           :: n_default = 2.0
!
end module Block_Network
