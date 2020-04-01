Module Block_Network
!
! Module with stream topology variables
!
!
    integer, dimension(:), allocatable  :: first_seg
    integer, dimension(:), allocatable  :: ndelta,no_cells,no_tribs,no_wr_units
!    no_cells = 0
!    no_tribs = 0
!
    integer, dimension(:,:), allocatable    :: head_cell,res_no,trib
    integer, dimension(:,:), allocatable    :: no_celm,upstrm_cell,cells_wr
!    no_celm = 0
    integer, dimension(:,:,:), allocatable  :: segment_cell
!
    character(len=5),dimension(:,:), allocatable :: unit_type 
!
! Integer variables 
!
    integer           :: n1,n2,nhead,nrch,nreach,no_res
    integer           :: no_cells_ttl
    integer           :: ind,ipd,nd_year,ndays,nyear,nwpd
    integer           :: start_year,start_month,start_day,start_hour
    integer           :: end_year,end_month,end_day,end_hour
!
! Logical variables
!
    logical           :: IS_HEAD,Is_a_Riv,Is_a_Res 
!
! Real variables
!
    real                           :: delta_n,dt_comp
    real                           :: x_bndry,x_head
    real                           :: time,day,day_fract,hr_fract,sim_incr,year,xd_year,yr_days
    real,parameter                 :: n_default = 2.0
!
end module Block_Network
