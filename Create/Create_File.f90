Program Create_File
implicit none
!
! Integer variables
!
integer::delta_t,hr_avg,ierror,iostat,n,nc,nd_start,no_steps
integer::narray,nf,nfile,n_head,nn,no_years,no_seg,nvg
integer::no_avg,no_dt,ndpnt,no_days,nobs,nobs_start,nobs_end
integer::start_day,start_mon,start_yr,end_day,end_mon,end_yr,start_hour,end_hour
integer::Julian,start_jul,end_jul
integer,allocatable,dimension(:):: seg_no,seg_indx,seg_seq,seg_net
integer,allocatable,dimension(:):: dummy
!
! Real variables
real                            :: dt_steps
real,parameter                  :: press=1013.
real,parameter,dimension(7) :: c_fctr = (/1.0,2.3884e-04,2.3884e-04,0.01,1.0   &
                                          ,35.315,35.315/)
real,allocatable,dimension(:)   :: depth,out_flow,in_flow,lat_flow
real,allocatable,dimension(:,:) :: DHSVM_out,forcing
!
! Character variables
!
character (len=1)  :: colon=':'
character (len=19) :: start_date,end_date 
character (len=19) :: time_stamp0,time_stamp
character (len=4)  :: path
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
!read(*,*) Project
!
write(*,*) TRIM(Project)//'.map'
open(10,file=TRIM(Project)//'.segmap',status='old')
open(20,file=TRIM(InDirectry)//'ATP.Only',status='old')
open(21,file=TRIM(InDirectry)//'NLW.Only',status='old')
open(22,file=TRIM(InDirectry)//'NSW.Only',status='old')
open(23,file=TRIM(InDirectry)//'VP.Only',status='old')
open(24,file=TRIM(InDirectry)//'WND.Only',status='old')
open(25,file=TRIM(InDirectry)//'Inflow.Only',status='old')
open(26,file=TRIM(InDirectry)//'Outflow.Only',status='old')
!
open(30,file=TRIM(Project)//'.forcing',status='unknown')
!
! Get some information about the network from the "<Project>.map" file
!
read(10,*) n_head,no_seg
!
! Allocate arrays after increasing expected size to account
! for differences in the number of segments and the size of
! the largest stream segment number in DHSVM.
!
narray=no_seg+no_seg/2
allocate (dummy(narray))
allocate (seg_no(narray))
allocate (seg_indx(narray))
allocate (seg_net(narray))
allocate (seg_seq(narray))
allocate (in_flow(narray))
allocate (out_flow(narray))
!allocate (lat_flow(narray))
!allocate (depth(narray))
allocate (DHSVM_out(7,narray))
allocate (forcing(7,narray))
!
do n=1,no_seg
  read(10,*) sequence,nn,path,seg_no(n)
end do
!
nfile=20
do nf=1,7
  read(nfile,'(A19,1x,A19,1x,i2)') start_date,end_date, no_dt
  write(*,*) 'start ',start_date
  write(*,*) 'end ',end_date
  read(start_date,'(i2,1x,i2,1x,i4,1x,i2,a6)') start_mon,start_day,start_yr  &
                                              ,start_hour,fluff
  start_jul=Julian(start_yr,start_mon,start_day)
  read(end_date,'(i2,1x,i2,1x,i4,1x,i2,a6)') end_mon,end_day,end_yr          &
                                              ,end_hour,fluff
  nfile=nfile+1
end do
!
write(*,*) 'Number of hours in averaging period,hr_avg. Note: hr_avg/no_dt must be integer'
read(*,*) hr_avg
no_avg = hr_avg/no_dt
dt_steps = no_avg
!
! Number of daily values in a day in DHSVM file
ndpnt = 24/no_dt
!
! Number of daily values written to the forcing file
nforce = 24/hr_avg
!
! Determine number of time steps per day
write(*,*) 'no_dt,hr_avg,no_avg ',no_dt,hr_avg,no_avg,dt_steps
!
delta_t = no_dt
!
nd_start=1+(start_hour/delta_t)
!
! Julian date of start and end of forcing records
!
end_jul=Julian(end_yr,end_mon,end_day)
!
! Total number of data points in each file
!
nobs=ndpnt*(end_jul - start_jul)
!
write(*,*) 'Julian',start_jul,end_jul
!
write(*,*) 'no_cycles ',no_days
!
write(30,'(2(i4.4,i2.2,i2.2,a1,i2.2,1x),2i4)')          &
     start_yr,start_mon,start_day,colon,start_hour          &
    ,end_yr,end_mon,end_day,colon,end_hour                  &
    ,no_dt,nd_start
!
! Read segment mapping
!
nfile=19
!
! Read the segment sequencing from the header of the ATP.Only file
! and establish the relationship between DHSVM segment numbers and
! the indexed location in the forcing files
!
  nfile=nfile+1
  read(nfile,*) (seg_seq(n),n=1,no_seg)
  do nf=1,no_seg
    seg_net(seg_seq(nf))=nf
  end do
do nf=2,7
  nfile=nfile+1
  read(nfile,*) (dummy(n),n=1,no_seg)
end do!
! Read the forcings from the DHSVM file
!
do nc=1,nobs
!
  forcing = 0.0
!
  do nvg = 1,no_avg
    nfile=19
     do nf=1,5
      nfile=nfile+1
      read(nfile,*) time_stamp,(DHSVM_out(nf,n),n=1,no_seg)
      if (nf .eq. 1) write(35,*) time_stamp
      do n=1,no_seg
        forcing(nf,n)=forcing(nf,n)+c_fctr(nf)*DHSVM_out(nf,n)/dt_steps
      end do
    end do
    do nf=6,7
      nfile=nfile+1
      read(nfile,*) time_stamp,(DHSVM_out(nf,n),n=1,no_seg)
!
! Minimum flow is 1.0 (m**3/sec) 
!
      if (DHSVM_out(nf,n) .lt. 1.0) DHSVM_out(nf,n) = 1.0
      do n=1,no_seg
        forcing(nf,n)=forcing(nf,n)+c_fctr(nf)*DHSVM_out(nf,n)/dt_steps
      end do
    end do
  end do
!
! Write the output file
!
  do n=1,no_seg
    nf=seg_no(n)
    nn=seg_net(nf)
    write(30,*) n,nn,press,(forcing(nf,nn),nf=1,7)
  end do
end do
end Program Create_File
