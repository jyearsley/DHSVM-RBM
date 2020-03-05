Subroutine Alloc(param_file)
!
use Block_Energy
use Block_Hydro
use Block_Network
use Block_WQ
!
implicit none
!
    character (len=15)  :: dummy
    character (len=5)   :: wr_type
    character (len=10)  :: dummy_a
    character (len=8)   :: dummy_b
    character (len=29)  :: dummy_c,basin_name
    character (len=50)  :: header
    character (len=200) :: param_file
!
! Integer variables
!
    integer,parameter                     :: nlayers = 2
    integer                               :: main_stem ! May not be needed
    integer                               :: n,nc,nnc,ncell,ncll,ndde,nseg
    integer                               :: nrch,no_wru,no_wru_0
    integer                               :: nseg_max=-9999,no_tribs_max=-9999
    integer                               :: no_res_max=-9999,no_wru_max=-9999,res_seg_max=-9999
    integer                               :: res_seg,trib_cell
    integer, dimension(:),allocatable     :: nn_tribs
    integer, dimension(:),allocatable     :: nn_wru_units
    integer, dimension(:),allocatable     :: nn_cells,nn_delta
    integer, dimension(:,:),allocatable   :: cclls_wr,nn_celm
    integer, dimension(:,:,:),allocatable :: ssgment_cell
!
! Real variables
!
    real                                 :: nndlta,source,x_0,x_1,z
!
!   allocate arrays
!
   allocate (nn_tribs(1000),nn_wru_units(1000))
   allocate (nn_cells(1000),nn_delta(1000))
   allocate (cclls_wr(500,50),nn_celm(200,200))
   allocate (ssgment_cell(200,50,100))
!
   nn_tribs     = 0
   nn_wru_units = 0
! 
!     Open network file
!
write(*,*) 'allocate1'
      OPEN(UNIT=90,FILE=param_file,STATUS='OLD')
!
!     Read header information from control file
      do n=1,5
        read(90,*) header
        write(*,*) header
      end do

!
!     Card Group I
!
      read(90,*)  nreach,dummy_a,ncell,dummy_c,basin_name
      write(*,*)  nreach,dummy_a,ncell,dummy_c,basin_name
      
!
!     Start reading the reach date and initialize the reach index, NR
!     and the cell index, NCELL
!
      ncell    = 0
      nrch     = 1
      no_res   = 0
      nseg     = 0
!
!     Card Group IIb. Reach characteristics
!
       do while (nrch .le. nreach) 
       write(*,*) 'nrch ',nreach,nrch
!
! Initialize index that counts the number of water resource units in this reach
!
         no_wru_0  = 0
!
!     Initialize NSEG, the total number of segments in this reach
      nseg=0
!
!     Read the number of cells in this reach, the headwater #,
!     the number of the cell where it enters the next higher order stream,
!     the headwater number of the next higher order stream it enters, and
!     the river mile of the headwaters.
!
!
        read(90,*) dummy_a,nn_cells(nrch)                          &
                  ,dummy_a,main_stem,dummy_b,trib_cell
        write(*,*) dummy_a,nn_cells(nrch)                          &
                  ,dummy_a,main_stem,dummy_b,trib_cell
!
!     If this is reach that is tributary to cell TRIB_CELL, give it the
!     pointer TRIB(TRIB_CELL) the index of this reach for further use.
!     Also keep track of the total number of tributaries for this cell
!
        if (trib_cell.gt.0) then
           nn_tribs(trib_cell)=nn_tribs(trib_cell)+1
           if (nn_tribs(trib_cell) .gt. no_tribs_max) no_tribs_max = nn_tribs(trib_cell)
        end if
!
!     Reading Reach Element information
!
!
! 
        nn_wru_units(nrch) = 0
        do nc=1,nn_cells(nrch)
          ncell=ncell+1
!
!     The headwaters index for each cell in this reach is given
!     in the order the cells are read*dt_calc/(z*rfac)
!
!     Card Type 3. Cell indexing #, Node # Row # Column Lat Long RM
!
          read(90,'(a3,i8,1x,a4,i6,1x,a3,f12.0,1x,a3,f12.0,1x,a9,f5.0,i3,1x,a4,i2,1x,a5)')         &
                   dummy,nnc,dummy,ndde                          &
                   ,dummy,x_0,dummy,x_1                          &
                   ,dummy,z,nn_delta(ncell)                      &
                   ,dummy,no_wru,wr_type
!
          if (no_wru .gt. no_wru_max) no_wru_max = no_wru
          if (no_wru .ne. no_wru_0) then
            nseg = 0
            nn_celm(nrch,no_wru) = 0
            cclls_wr(nrch,no_wru)  = 0
!
! If this is a reservoir, update reservoir index
!
            if (wr_type .eq. 'RSRVR') then
              no_res = no_res + 1
              if (no_res .gt. no_res_max) no_res_max = no_res
              res_seg = 0
            end if
!
! Update water resource unit index
!
            no_wru_0 = no_wru
          end if 
          if (wr_type .eq. 'RSRVR') then
             res_seg = res_seg + 1
             if (res_seg .gt. res_seg_max) res_seg_max = res_seg
             nn_celm(nrch,no_wru) = 1
             cclls_wr(nrch,no_wru) = cclls_wr(nrch,no_wru) + 1 
          else          
            do ncll = 1,nn_delta(ncell)
              nseg = nseg + 1
              nn_celm(nrch,no_wru) = nn_celm(nrch,no_wru)+1
              ssgment_cell(nrch,no_wru,nseg)=ncell
!
              if (nseg .gt. nseg_max) nseg_max = nseg
!
           end do
!
           cclls_wr(nrch,no_wru) = cclls_wr(nrch,no_wru) + 1 
          end if          
!   
        nn_wru_units(nrch) = no_wru
!
      end do
!
write(97,*) 'nrch ',nrch
      nrch = nrch + 1
      end do
!
! Close the network file
!
   nrch = nreach
!
   close(unit=90)
!
write(*,*) param_file
write(*,*) 'Results - ',ncell,nreach,nrch,nseg_max,no_res_max,res_seg_max,no_wru_max
!
!  Block_Energy
!
!   Incoming short wave radiation, Watts/m**2
!
    allocate  (QNS(ncell))
!
!   Incoming atmospheric radiation, Watts/m**2
!
    allocate  (QNA(ncell))
!
!   Air temperature at surface, deg. C
!
    allocate  (DBT(ncell))
!  
!   Wind speed, m/sec
!
    allocate  (WIND(ncell))
!
!   Vapor pressure of air at surface, MB
!
    allocate  (EA(ncell))
!
!   Air pressure at surface, mb
!
    allocate   (PRESS(ncell)) 
!
! Mohseni parameters
!
    allocate   (T_smth(nreach),alfa_Mu(nreach),beta(nreach),gmma(nreach),mu(nreach))

! Block Hydro
!
    allocate   (depth(ncell),D_a(ncell),D_b(ncell),D_min(ncell))
    allocate   (elev(ncell))
    allocate   (width(ncell))
    allocate   (u(ncell),U_a(ncell),U_b(ncell),U_min(ncell))
    allocate   (dt(ncell))
    allocate   (dx(ncell))
    allocate   (x_dist(ncell,no_wru_max,-2:nseg_max))
!
! Reservoir characteristics
!
    allocate   (Kappa(no_res_max))
    allocate   (A_surf(no_res_max,res_seg_max,nlayers)           &  
               ,Volume(no_res_max,res_seg_max,nlayers))
!
! Flows
!
    allocate   (Q_in(ncell))
    allocate   (Q_trib(no_tribs_max))
    allocate   (Q_out(ncell))
    allocate   (Q_diff(ncell))
!
    allocate   (Q_in_seg(ncell,nseg_max))
    allocate   (Q_out_seg(ncell,nseg_max))
    allocate   (Q_nps(ncell,nseg_max))
!
! Network variables
!
    allocate  (ndelta(ncell),no_tribs(ncell))
    allocate  (no_cells(nreach),no_wr_units(nreach))
!    allocate  (nstrt_elm(nseg_max),no_dt(nseg_max))
!
    allocate (head_cell(nreach,no_wru_max),res_no(nreach,no_wru_max))
    allocate (no_celm(nreach,no_wru_max),upstrm_cell(nreach,no_wru_max))
    allocate (cells_wr(nreach,no_wru_max),unit_type(nreach,no_wru_max))
    allocate (trib(ncell,no_tribs_max))
    allocate (segment_cell(nreach,no_wru_max,nseg_max))
!
! Water quality variables
!
    allocate (temp(nreach,no_wru_max,-2:nseg_max,2))
    allocate (temp_res(no_res_max,nlayers,2))
    allocate (temp_head(nreach,no_wru_max))
    allocate (temp_trib(no_tribs_max))
!   allocate (temp_source(----))
!   allocate (temp_nps----)) 
!
! Deallocate dummy variables
!
!
!   deallocate arrays
!
   deallocate (nn_tribs,nn_wru_units)
   deallocate (nn_cells,nn_delta)
   deallocate (cclls_wr,nn_celm)
   deallocate (ssgment_cell)
    write(*,*) ' Leaving Subroutine Alloc '
!
End Subroutine Alloc



