!
Subroutine BEGIN(param_file,spatial_file)
!
use Block_Energy
use Block_Hydro
use Block_Network
use Block_WQ
!
implicit none
!
    character (len=200):: param_file,source_file,spatial_file
    integer:: Julian
!
! Character variables
!
    character (len=8) :: end_date,start_date     
    character (len=8) :: lat
    character (len=10):: long
!
! Integer variables
!
integer:: cell_check_tds,cell_check_temp,head_name,trib_cell
integer:: jul_start,main_stem,nyear1,nyear2,nc,ncell,nseg,seg_inp
integer:: ns_max_test,node,ncol,nrow,nr,cum_sgmnt
!
! Logical variables
!
logical:: first_cell,source
logical:: TRIBS_DONE
!
! Real variables
!
  real            :: nndlta
  real            :: rmile0,rmile1,xwpd
  real            :: tds_inp,temp_inp
  real, parameter :: rho_Cp=1000. ! Units are kcal/m**3/deg K
  real, parameter :: miles_to_ft=5280. ! Convert miles to feet
!
!
!
!   Mohseni parameters, if used
!
!
!
!     Card Group I
!
read(90,*) start_date,end_date
read(start_date,'(i4,2i2)') start_year,start_month,start_day
read(end_date,  '(i4,2i2)') end_year,end_month,end_day
nyear1=start_year
nyear2=end_year
write(*,'(2(2x,i4,2i2))')  &
 start_year,start_month,start_day,end_year,end_month,end_day
!
!     Establish the Julian day for which simulations begin
!
jul_start = Julian(start_year,start_month,start_day)
!
!
read(90,*) nreach,flow_cells,heat_cells,source
!
! Allocate dynamic arrays
!
 allocate(tds_head(nreach))
 allocate(temp_head(nreach))
 allocate(ndelta(heat_cells))
 allocate(mu(nreach))
 allocate(alphamu(nreach))
 allocate(beta(nreach))
 allocate(gmma(nreach))
 allocate (smooth_param(nreach))
 allocate(dx(heat_cells))
 allocate(no_celm(nreach))
 no_celm=0
 allocate(no_cells(nreach))
 no_cells=0
 allocate(no_tribs(heat_cells))
 no_tribs=0
 allocate(trib(heat_cells,20))
 trib=0
 allocate(head_cell(nreach))
 allocate (conflnce(heat_cells,20))
 conflnce=0
 allocate(reach_cell(nreach,ns_max))
 allocate(segment_cell(nreach,ns_max))
 allocate(x_dist(nreach,0:ns_max))
 allocate(tds_source(nreach,ns_max))
 allocate(temp_source(nreach,ns_max))
      SUBROUTINE BEGIN
      character*11 end_time,start_time
      character*5 Dummy_B
      character*5 wr_type
      character*10 Dummy_A
      integer head_name,res_seg,trib_cell,first_cell
      dimension ndmo(12)
      logical Is_a_Riv,Is_a_Res
      logical Test
      INCLUDE 'RBM.fi'
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
!      ndelta=2
!      delta_n=ndelta
!
!     Read the starting and ending times and the number of
!     periods per day of weather data from the forcing file
! 
! 
      read(30,*) start_time,end_time,nwpd,nd_start
!
      write(*,*) start_time,'  ',end_time,nwpd,nd_start
!
      read(start_time,'(i4,2i2,1x,i2)') start_year,start_month
     &                                 ,start_day,start_hour
      read(end_time,'(i4,2i2,1x,i2)') end_year,end_month
     &                               ,end_day,end_hour
      nyear1=start_year
      nyear2=end_year
      write(*,*) start_year,start_month,start_day
     &             ,end_year,end_month,end_day
!
!     Establish the Julian day for which simulations begin
! 
      jul_start=julian(start_year,start_month,start_day)
!
      read(90,*)  nreach
      write(*,*) "Number of stream reaches - ", nreach
      read(40,*) Test,a_smooth
      write(*,*) Test,a_smooth, nreach
      b_smooth=1.- a_smooth
      do nr=1, nreach
        read(40,*) nnr,alf_Mu(nr),beta(nr),gmma(nr),mu(nr)
        write(*,*) nnr,alf_Mu(nr),beta(nr),gmma(nr),mu(nr)
!
      end do
!
! Iniialize the upstream cell index
!
      do n = 1,500
        do nn = 1,50
          upstrm_cell(n,nn) = -1
        end do
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
  100 continue
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
      write(*,*) ' Starting to read reach ', nrch
!
!     Read the number of cells in this reach, the headwater #,
!     the number of the cell where it enters the next higher order stream,
!     the headwater number of the next higher order stream it enters, and
!     the river mile of the headwaters.
!
!
        read(90,*) Dummy_A,no_cells( nrch)
     &            ,Dummy_A,main_stem( nrch),Dummy_A,trib_cell
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

          read(90,*) Dummy_B,nnc,Dummy_B,node(nnc)
     &              ,Dummy_B,x_0,Dummy_B,x_1
     &              ,Dummy_B,elev( nrch),ndelta(ncell)
     &             ,Dummy_B,no_wru,wr_type
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
             res_seg = res_seg + 1
             no_celm(nrch,no_wru) = 1
!            res_no(nrch,no_wru)   = no_res
             volume(no_res,res_seg,1) = U_a(ncell)
             volume(no_res,res_seg,2) = U_b(ncell)
             a_surf(no_res,res_seg,1) = D_a(ncell)
             a_surf(no_res,res_seg,2) = D_b(ncell)
             Kappa(no_res)         = D_min(ncell)
            write(85,*) wr_type,no_wru,nseg
     &                 ,upstrm_cell(nrch,no_wru)
     &                 ,ncell,cells_wr(nrch,no_wru)
     &                 ,segment_cell(nrch,no_wru,nseg)
     &                 ,x_dist(nrch,no_wru,nseg)
     &                 ,x_dist(nrch,no_wru,nseg-1),x_1,x_0,dx(ncell)                        
                    cells_wr(nrch,no_wru) 
     &                   = cells_wr(nrch,no_wru) + 1 
          else             
            do ncll = 1,ndelta(ncell)
              nseg = nseg + 1
              no_celm(nrch,no_wru) = no_celm(nrch,no_wru)+1
              segment_cell(nrch,no_wru,nseg)=ncell
              dx(ncell)=(x_1-x_0)/ndelta(ncell)
              x_dist(nrch,no_wru,nseg) 
     &        = x_dist(nrch,no_wru,nseg-1)-dx(ncell)
           end do
!
           cells_wr(nrch,no_wru) 
     &                   = cells_wr(nrch,no_wru) + 1 
          end if
          
                   
!

  200 continue
! 
!   
        no_wr_units(nrch) = no_wru
!
      end do
!
!
!	go to 100

      nrch = nrch + 1

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
900 continue
!
!
end subroutine BEGIN
