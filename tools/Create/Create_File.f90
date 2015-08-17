Program Create_File
implicit none
!
! Integer variables
!
integer::delta_t,f_id,ierror,iostat,nc,nd_start,no_cycles,no_src
integer::narray,nf,nfile,n_head,n,nn,no_years,no_seg,nt,no_tribs
integer::no_dt,no_days,nobs_start,nobs_end
integer::start_day,start_mon,start_yr,end_day,end_mon,end_yr,start_hour,end_hour
integer::Julian,start_jul,end_jul
integer,dimension(1000)         :: file_id,net_ndx
integer,allocatable,dimension(:):: RBM_seg,seg_no,seg_seq,seg_net
integer,allocatable,dimension(:):: dummy
integer,allocatable,dimension(:):: trib_ndx_out
!
! Real variables
real::press=1013.
real,allocatable,dimension(:)::depth,out_flow,in_flow,lat_flow
real,allocatable,dimension(:)::trib_flow,trib_temp
real,allocatable,dimension(:,:)::forcing
!
! Character variables
!  
character (len=1)  :: colon=':'
character (len=12) :: Basin
character (len=19) :: start_date,end_date 
character (len=19) :: time_stamp0,time_stamp
character (len=4)  :: path,Date
character (len=6)  :: fluff
character (len=8)  :: sequence
character (len=200):: InDirectry,Project
!
!
integer iargc
integer numarg

!
! Command line input
!
      numarg = iargc ( )
      if (numarg .lt. 2) then
        write (*,*) 'Too few arguments were given'
        write (*,*) ' '
        write (*,*) 'First:  Directory with forcing files *(*.Only)'
        write (*,*) 'Second:  Project Name'
        write (*,*) 'eg: $ ./Create_File <Input directory> <Project Name>'
        write (*,*) ' '
        stop
      end if
      call getarg ( 1, InDirectry )
      call getarg ( 2, Project )
!
! Following are the files required for the Main Ste
!
!write(*,*) 'Name of Project.  Required files include:'
!write(*,*) 'ProjectName.map'
!write(*,*) 'ATP.Only     - Air temperature'
!write(*,*) 'NLW.Only     - Net longwave radiation'
!write(*,*) 'NSW.Only     - Net shortwave radiation'
!write(*,*) 'VP.Only      - Vapor pressure'
!write(*,*) 'WND.Only     - Wind speed'
!write(*,*) 'Inflow.Only  - Segment inflow'
!write(*,*) 'Outflow.Only - Segment outflow'
!
! Following are files required for the tributaries
!
!write(*,*) 'trib.Inflow.Only  - Tributary inflow'
!write(*,*) 'trib.Temp.Only    - Tributary temperature'
!
open(27,file='trib.Inflow.Only',status='old')
open(28,file='trib.Temp.Only',status='old')
!
! Forcings for Main Stem
!
! Inflow and stream temperature from the tributaries
!
open(35,file=TRIM(Project)//'.tribs',status='unknown')
!
! Allocate arrays after increasing expected size to account
! for differences in the number of segments and the size of
! the largest stream segment number in DHSVM.
!
narray=no_seg+no_seg/2
allocate (dummy(narray))
allocate (seg_no(narray))
allocate (seg_net(narray))
allocate (seg_seq(narray))
allocate (RBM_seg(narray))
allocate (in_flow(narray))
allocate (out_flow(narray))
!allocate (lat_flow(narray))
!allocate (depth(narray))
allocate (forcing(5,narray))
!
!
!
!
delta_t=no_dt
!
! Determine number of time steps per day
write(*,*) 'delta_t ',no_dt,delta_t
no_dt=24/delta_t
!
nobs_start=no_dt
nobs_end=no_dt
!
if (no_dt .gt. 1) then
  nobs_start = (24 - start_hour)/delta_t
  nobs_end   = 1+(end_hour/delta_t)
end if 
nd_start=1+(start_hour/delta_t)
if (nobs_start .lt. no_dt) then
  write(*,*)
  write(*,*) '!!!===WARNING: There is (are) only',nobs_start,'value(s) for the first simulation day'
  write(*,*) 'rather than ',no_dt, 'values. Daily averaging scripts will exclude the first' 
  write(*,*) 'day in the daily temperature output. ===!!!'
  write(*,*)
end if
if (nobs_end .lt. no_dt) then
  write(*,*)
  write(*,*) '!!!===WARNING: There is (are) only',nobs_end,'value(s) for the last simulation day'
  write(*,*) 'rather than ',no_dt, 'values. Daily averaging scripts will exclude the last'
  write(*,*) 'day in the daily temperature output. ===!!!'
  write(*,*)
end if
!

!
! Julian date of start and end of forcing records
!
end_jul=Julian(end_yr,end_mon,end_day)
no_days=end_jul - start_jul
!
write(*,*) 'Julian',start_jul,end_jul
!
no_cycles=nobs_start+nobs_end+no_dt*(no_days-1)

!
write(*,*) 'no_cycles ',no_cycles
!
!
! Read the headers of the .Source file
!
!
  read(15,*) no_tribs
!
! Allocate arrays
!
  allocate (trib_ndx_out(no_tribs+1))
  allocate (trib_flow(no_tribs+1))
  allocate (trib_temp(no_tribs+1))
!
! Read the *.Source file 
!
!   Upper_Conn     1   248  1 301 1 392
  do n=1,no_tribs+1
    read(15,*) Basin,nn,RBM_seg,no_src,net_ndx(n),nt,f_id
    file_id(f_id)=nt-1
  end do
!
! Set up the cross-references 
!
  do n=1,no_tribs
    trib_ndx_out(n)=file_id(net_ndx(n))
    write(*,*) n,trib_ndx_out(n),net_ndx(n)
  end do
!
! Cycle through reading and writing tributary advection
!
do nc=1,no_cycles
!
!  read(27,*) date,(trib_flow(nn),nn=1,no_tribs+1)
!  read(28,*) date,(trib_temp(nn),nn=1,no_tribs+1)
!
! Write the tributary file
!
!write(35,*) 
end do
end Program Create_File
