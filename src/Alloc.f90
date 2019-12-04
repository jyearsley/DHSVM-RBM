Subroutine Alloc(nhead,nreach,no_res,res_seg_max,nseg_max,nwru_max)
!
!
!   Energy budget variables
!
!   Incoming short wave radiation, watts/m**2
!
    allocate (qns(ncell))
!
!   Incoming atmospheric radiation, kcal/m**2/sec
!
    allocate (qna(ncell))
!
!   Air temperature at surface, deg. C
!
    allocate (dbt(ncell))
!  
!   Wind speed, m/sec
!
    allocate (wind(ncell))
!
!   Vapor pressure of air at surface, mb
!
    allocate (ea(ncell))
!
!   Air pressure at surface, mb
!
    allocate (press(ncell)) 
!
!   Mohseni parameters
!
    allocate (mu(nhead),alpha_Mu(nhead),beta(nhead),gmma(nhead),smooth_param(nhead))
!
! Hydraulic characteristics
!
    allocate   (depth(ncell),D_a(ncell),D_b(ncell),D_min(ncell))
    allocate   (elev(ncell)_
    allocate   (width(ncell))
    allocate   (u(ncell),U_a(ncell),U_b(ncell),U_min(ncell))
    allocate   (dt(ncell))
    allocate   (dx(ncell))
    allocate   (x_dist(ncell,nwru_max,nseg_max))
!
! Reservoir characteristics
!
    allocate   (Kappa(no_res))
    allocate   (A_surf(no_res,res_seg_max,2),Volume(no_res,res_seg_max,2))
!
! Flows
!
    allocate   (Q_in(ncell))
    allocate   (Q_trib(ncell))
    allocate   (Q_out(ncell))
    allocate   (Q_diff(ncell))
!
    allocate   (Q_in_seg(ncell,nseg_max))
    allocate   (Q_out_seg(ncell,nseg_max))
    allocate   (Q_nps(ncell,nseg_max))
!
! Network variables
!
    allocate  (ndelta(ncell),no_cells,node,no_tribs,no_wr_units
!
    allocate (head_cell(nreach,no_wru_max),res_nonreach,no_wru_max))
    allocate (no_celm(nreach,no_wru_max),upstrm_cell(nreach,no_wru_max))
    allocate (cells_wr(nreach,no_wru_max),unit_type(nreach,no_wru_max))
    allocate (segment_cell(nreach,no_wru_max,nseg_max)
!
End Subroutine Alloc



