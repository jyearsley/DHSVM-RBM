!
Subroutine BEGIN (param_file)
!
use Block_Energy
use Block_Hydro
use Block_Network
use Block_WQ
!
!implicit none
!
    character (len=200) :: param_file,spatial_file
    integer             :: Julian
!
! Character variables
!
    character (len=5)   :: wr_type
    character (len=8)   :: end_date,start_date     
    character (len=8)   :: lat,dummy_b
    character (len=10)  :: long
    character (len=10)  :: dummy_a
    character (len=11)  :: end_time,start_time
    character (len=15)  :: dummy
    character (len=80)  :: header
!
! Integer variables
!
    integer           :: main_stem ! May not be needed
    integer           :: res_seg,trib_cell
    integer           :: jul_start,nyear1,nyear2,nc,ncell,nseg,seg_inp
    integer           :: n,nnc,ncll,ndde,nd_start,nr,nnr
    integer           :: nrch,ns_total,no_wru,no_wru_0,total_cells
    integer           :: nlayers = 2
! Logical variables
!
    logical           :: Test
!
! Real variables
!
  real            :: nndlta
  real            :: xwpd,x_0,x_1
  real            :: tds_inp,temp_inp
  real, parameter :: rho_Cp=1000. ! Units are kcal/m**3/deg K
  real, parameter :: miles_to_ft=5280. ! Convert miles to feet
!!
open(unit=90,file=param_file,status='old')
!
!     Read header information from control file
      do n=1,5
        read(90,*) header
        write(*,*) 'First file lines ',header
      end do
!
!     Card Group I
!
      read(30,*) start_time,end_time,nwpd,nd_start
!
      write(*,*) start_time,'  ',end_time,nwpd,nd_start
!
      read(start_time,'(i4,2i2,1x,i2)') start_year,start_month   &
                                       ,start_day,start_hour
      read(end_time,'(i4,2i2,1x,i2)') end_year,end_month         &
                                     ,end_day,end_hour
      nyear1=start_year
      nyear2=end_year
      write(*,*) start_year,start_month,start_day                &
                   ,end_year,end_month,end_day
!
!     Establish the Julian day for which simulations begin
! 
      jul_start=julian(start_year,start_month,start_day)
!
      read(90,*)  nreach,dummy_a,ncell

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
    allocate   (x_dist(ncell,10,0:100))
!
! Reservoir characteristics
!
    allocate   (Kappa(10))
    allocate   (A_surf(10,10,nlayers)           &  
               ,Volume(10,10,nlayers))
!
! Flows
!
    allocate   (Q_in(ncell))
    allocate   (Q_trib(nreach))
    allocate   (Q_out(ncell))
    allocate   (Q_diff(ncell))
!
    allocate   (Q_in_seg(ncell,100))
    allocate   (Q_out_seg(ncell,100))
    allocate   (Q_nps(ncell,100))
!
! Network variables
!
    allocate  (first_seg(ncell))
    allocate  (ndelta(ncell),no_tribs(ncell))
    allocate  (no_cells(nreach),no_wr_units(nreach))
!    allocate  (nstrt_elm(100),no_dt(100))
!
    allocate (head_cell(nreach,10),res_no(nreach,50))
    allocate (no_celm(nreach,50),upstrm_cell(nreach,50))
    allocate (cells_wr(nreach,50),unit_type(nreach,50))
    allocate (trib(ncell,50))
    allocate (segment_cell(nreach,50,100))
!
! Water quality variables
!
    allocate (temp(nreach,50,-2:100,2))
    allocate (temp_res(50,nlayers,2))
    allocate (temp_head(nreach,50))
    allocate (temp_trib(nreach))
    allocate (temp_nps(50,50))
    allocate (temp_source(50,50))
!
no_tribs=0
trib=0
temp_trib(:)=5.0
Temp_res(:,:,:) = 2.0
do n=1,nreach
 write(*,*) temp_trib(n)
end do
!
      write(*,*) "Number of stream reaches - ", nreach,ncell
      read(40,*) Test,a_smooth
      write(*,*) Test,a_smooth, nreach
      b_smooth=1.- a_smooth
      do nr=1, nreach
        read(40,*) nnr,alfa_Mu(nr),beta(nr),gmma(nr),mu(nr)
      end do
!
!     Start reading the reach date and initialize the reach index, NR
!     and the cell index, NCELL
!
      ncell    = 0
      nrch     = 1
      no_res   = 0
      ns_total = 0
      nseg     = 0
!
     upstrm_cell = -1
!
!  100 continue
!
!     Card Group IIb. Reach characteristics
!
       do while (nrch .le. nreach) 
!
! Initialize index that counts the number of water resource units in this reach
!
         no_wru_0  = 0
!
!     Initialize NSEG, the total number of segments in this reach
      nseg=0
!      write(*,*) ' Starting to read reach ', nrch
!
!     Read the number of cells in this reach, the headwater #,
!     the number of the cell where it enters the next higher order stream,
!     the headwater number of the next higher order stream it enters, and
!     the river mile of the headwaters.
!
!
        read(90,*) dummy_a,no_cells(nrch)                          &
                  ,dummy_a,main_stem,dummy_b,trib_cell
write(*,*) 'Trib cell ',trib_cell
!
!     If this is reach that is tributary to cell TRIB_CELL, give it the
!     pointer TRIB(TRIB_CELL) the index of this reach for further use.
!     Also keep track of the total number of tributaries for this cell
!
        if (trib_cell.gt.0) then
           no_tribs(trib_cell)=no_tribs(trib_cell)+1
           trib(trib_cell,no_tribs(trib_cell))= nrch
        end if
!
!     Reading Reach Element information
!
!
! 
        no_wr_units(nr) = 0
        Is_a_Riv = .FALSE.
        Is_a_Res = .FALSE.
        do nc=1,no_cells( nrch)
          ncell=ncell+1
          read(50,*) ncll,U_a(ncell),U_b(ncell),U_min(ncell)
          if (ncll .ne.ncell) then
            write(*,*) 'Mismatch in Leopold file',ncll,ncell
          end if  
          read(50,*) D_a(ncell),D_b(ncell),D_min(ncell)
!
!     The headwaters index for each cell in this reach is given
!     in the order the cells are read*dt_calc/(z*rfac)
!
!     Card Type 3. Cell indexing #, Node # Row # Column Lat Long RM
!

          read(90,'(a3,i8,1x,a4,i6,1x,a3,f12.0,1x,a3,f12.0,1x,a9,f5.0,i3,1x,a4,i2,1x,a5)')  &
                   dummy,nnc,dummy,ndde,                          &
                   dummy,x_0,dummy,x_1,                           &
                   dummy,z,ndelta(ncell),                         &
                   dummy,no_wru,wr_type
!
          unit_type(nrch,no_wru) = wr_type
!
!        
!
! Set reservior indices and geometry
!
          
          if (upstrm_cell(nrch,no_wru) .lt. 0) then
            upstrm_cell( nrch,no_wru) = ncell-1
          end if
          if (no_wru .ne. no_wru_0) then
            nseg = 0
            head_cell(nrch,no_wru) = ncell
            no_celm(nrch,no_wru) = 0
            cells_wr(nrch,no_wru)  = 0
            x_dist(nrch,no_wru,nseg) = x_1
!
! If this is a reservoir, update reservoir index
!
            if (wr_type .eq. 'RSRVR') then
              no_res = no_res + 1
              write(81,*) 'Reservoir _Number ',no_res
              res_seg = 0
            end if
!
! Update water resource unit index
!
            no_wru_0 = no_wru
          end if 
          if (wr_type .eq. 'RSRVR') then
             nseg = 1 !Temporary value while debugging JRY 02/28/2020
             first_seg(ncell) = nseg
             res_seg = res_seg + 1
             no_celm(nrch,no_wru) = 1
!            res_no(nrch,no_wru)   = no_res
             volume(no_res,res_seg,1) = U_a(ncell)
             volume(no_res,res_seg,2) = U_b(ncell)
             a_surf(no_res,res_seg,1) = D_a(ncell)
             a_surf(no_res,res_seg,2) = D_b(ncell)
             Kappa(no_res)            = D_min(ncell)
             cells_wr(nrch,no_wru) = cells_wr(nrch,no_wru) + 1 
            write(85,*) wr_type,nrch,no_wru,nseg                    &
                       ,upstrm_cell(nrch,no_wru)                   &
                       ,ncell,cells_wr(nrch,no_wru)                &
                       ,segment_cell(nrch,no_wru,nseg)             
          else     
write(85,*)        
            do ncll = 1,ndelta(ncell)
              nseg = nseg + 1
              if (ncll .eq. 1) first_seg(ncell) = nseg
              no_celm(nrch,no_wru) = no_celm(nrch,no_wru)+1
              segment_cell(nrch,no_wru,nseg)=ncell
              nndlta = ndelta(ncell)
              dx(ncell)=(x_1-x_0)/nndlta
              x_dist(nrch,no_wru,nseg) = x_dist(nrch,no_wru,nseg-1)-dx(ncell)
            write(85,*) wr_type,nrch,no_wru,first_seg(ncell),nseg  &
                       ,upstrm_cell(nrch,no_wru)                   &
                       ,ncell,cells_wr(nrch,no_wru)                &
                       ,segment_cell(nrch,no_wru,nseg)             &
                       ,x_dist(nrch,no_wru,nseg)                   &
                       ,x_dist(nrch,no_wru,nseg-1),x_1,x_0,dx(ncell)            
           end do
!
           cells_wr(nrch,no_wru) = cells_wr(nrch,no_wru) + 1 
          end if 
!
!  200 continue 
!   
        no_wr_units(nrch) = no_wru
!
      end do
!
!
!	go to 100

      nrch = nrch + 1
!
      end do
!
! Some final constants
!
      nrch= nreach
      xwpd=nwpd
      dt_comp=86400./xwpd
!
!     ******************************************************
!                         Return to RMAIN
!     ******************************************************
!
!
      return
end subroutine BEGIN
