      Subroutine systmm
!
      use Block_Hydro
      use Block_Network
      use Block_WQ
      use Block_Energy
!
      integer                           :: nd_start = 1
      integer                           :: nseg,nseg_trib,nseq,res_nn
      logical                           :: Leap_Year
      integer,dimension(12,2),parameter :: ndmo = reshape((/0,31,59,90,120,151,181,212,243,273,304,334, &
                                                  0,31,60,91,121,152,182,213,244,274,305,335/),shape(ndmo))
!
!
      hour_inc=1./nwpd
!
! Initialize the headwaters temperature and smoothed air temperature 
!
      Temp_head = 0.0
      T_smth = 0.0
!
      n1=1
      n2=2
      nobs=0
!
!     Initialize the day counter used for calculating the
!     record position for direct access files
!
!
      write(*,*) 'Start_year ',start_year,start_month,start_day
      write(*,*) 'End_year   ',end_year,end_month,end_day
      ndays=-Julian(start_year,start_month,start_day)            &
           +Julian(end_year,end_month,end_day)+1
      write(*,*) 'Number of days ',ndays
!      
!     Setup the timing of the simulation
! 
      lp_year = 1
      if (leap_year(start_year)) lp_year = 2
!
      day_fract=ndmo(start_month,lp_year)+start_day-1
      day_fract=day_fract/yr_days
      hr_fract=start_hour/(24.*yr_days)
      sim_incr=dt_comp/(365.*86400.)
      year=start_year
      time=year+day_fract+hr_fract
!
!     Year loop starts
!
!
end_year =1996
      do nyear=start_year,end_year
         write(*,*) ' Simulation Year - ',nyear
         nd_year=365
         yr_days=365
         xd_year=nd_year
!
         if (leap_year(nyear)) then
            nd_year = 366
            yr_days = 366
         end if
!
!        Day loop starts
         if (nyear.eq.start_year) then
           nd1=ndmo(start_month,lp_year)+start_day
         else
           nd1=1
         end if
         if (nyear.eq.end_year) then
           nd2=ndmo(end_month,lp_year)+end_day
         else
           nd2=nd_year
         end if
!
!
         DO ND=nd1,nd2
!
!     Start the numbers of days-to-date counter
           ndays=ndays+1
!
!     Daily period loop starts
           DO ndd=nd_start,nwpd 
!
!     Begin reach computations
!      
             ind=nd
             ipd=ndd
             day=nd
             period=ndd
             hour_inc = 3.*(period - 0.5)/24.
             time = (day - 1.)/xd_year
             time = time + hour_inc/xd_year
             time = year + time
!
! Initialize number of reservoirs (index is basin-wide)
!
             res_nn = 0
             nseq = 0
!
             do nr=1,nreach
               IS_HEAD = .TRUE.
               nrch = nr
               nseg = 0
               do no_wr = 1,no_wr_units(nr)
!
! Select either the RIVER or RSRVR case
!
!
               Select Case(unit_type(nrch,no_wr))
                 case('RIVER')
                   call RIVER(nseg,nseq,nr,no_wr)
!
                 case('RSRVR')      
                   res_nn = res_nn + 1
                   call RSRVR(nseg,nseq,nrch,no_wr,res_nn)
!
               End Select
               IS_HEAD = .FALSE. 
!
! End SELECT
!

               end do
!
! Establish tributary temperatures at reach end
!
!
               nseg_trib = no_celm(nr,no_wr_units(nr))
               Temp_trib(nr)=temp(nr,no_wr_units(nr),nseg_trib,n2)
!
             end do
             ntmp=n1
             n1=n2
             n2=ntmp
!
!     End of weather period loop (NDD=1,NWPD)
!
           end do
!
! Reset daily loop counter
!
          nd_start=1
! 
!
!     End of main loop (ND=1,365/366)
!

         end do
!
!    Update initial time for new year
!
      year=year+1
!
!     End of year loop
!
      end do
!
! Finish
!
!
!
!     ******************************************************
!                        return to rmain
!     ******************************************************
!
END SUBROUTINE Systmm
