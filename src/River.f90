!
!    Subroutine that simulates water temperature in advective system
!
      SUBROUTINE RIVER(nseg,nseq,nr,no_wr)
!
      use Block_Hydro
      use Block_Network
      use Block_WQ
!
      integer                            :: nr,no_wr
      integer, dimension(:), allocatable :: no_dt,nstrt_elm
      integer, dimension(:), allocatable :: ndltp,nterp
      logical                            :: DONE,pp_T
!
! Real variables
!
      real                               :: Kappa_bed,lat_flow
      real, dimension(:)                 :: dt_part(2000),x_part(2000)
      real, dimension(:), allocatable    :: xa,ta
      real, dimension(4), parameter      :: ndltp=(/-2,-1,-2,-2/),nterp=(/4,3,2,3/)
!
      real, parameter                    :: ft_to_km = 1./(3.2808*1000.)
      real, parameter                    :: (PI=3.14159,RFAC=4.184e6,wind_fctr = 1.0)
      real, parameter                    :: rho_sed = 2200., Cps = 0.210*4184.
      real, parameter                    :: alpha_s = 6.45e-07,H_sed = 0.10, T_sed = 15.0
!
      l_cell = 0
!
      do nc=1,cells_wr(nr,no_wr)
        l_cell = l_cell+1
        nseq  = nseq + 1
        read(30,*) l1
     &                    ,press(nseq),dbt(nseq)
     &                    ,qna(nseq),qns(nseq)
     &                    ,ea(nseq),wind(nseq)
     &                    ,qin(nseq),qout(nseq)
      if (l1 .ne. nseq) then
         write(*,*) 'Input file error at cell - ',l1,nseq
         stop
      end if
!
! Find the headwaters temperature
!
      nc_head=upstrm_cell(nr,no_wr)
      if(nseq .eq. nc_head+1) then
        call HEAD_TEMP (ncell,nr,no_wr)
!      if (nr.eq.95) write(72,*) nd,nr,no_wr,ncell,nseq
!     &              ,dbt(ncell),TEMP_head(nr,no_wr) 
!      if (nr .eq. 95) write(72,*)             
      end if
!
! Average flow in segment for estimating hydraulic conditions

        qavg=0.5*(qin(nseq)+qout(nseq))
!
!    Stream speed estimated with Leopold coefficients
! 
                 u(nseq)=U_a(nseq)*(qavg**U_b(nseq)) 
                 u(nseq) = amax1(u_min(nseq),u(nseq))
! 
!                 qdiff(nseq)=(qout(nseq)-qin(nseq))/delta_n
! 
! 
        dt(nseq)=dx(nseq)/u(nseq)
!
!    Depth estimated with Leopold coefficients
! 
        depth(nseq)=D_a(nseq)*(qavg**D_b(nseq))
        depth(nseq)=amax1(D_min(nseq),depth(nseq))
!
!      end do
!
!  Set the value of the tributary flow at the downstream segment
!  of the RIVER water resource segment
!
        q_trib(nr)=qout(nseq)
!
!     Main stem inflows and outflows for each reach first
!     Flows are cumulative and do not include tributaries if
!     tributaries are downstream of the inflow junction
!
      end do
!
! Temporary fix to the shortwave radiation in Segment 236  JRY 7/3/2018
!
      qns(236) = qns(237)
!
! Temporary method for fixing headwaters location - best done in Begin.f90
! 
      x_head=x_dist(nr,no_wr,0)
      x_bndry=x_head-1.0


!   Do the reverse particle tracking
!
      do ns=no_celm(nr,no_wr),1,-1
!
!     Segment is in cell SEGMENT_CELL(NC)
!

        ncell=segment_cell(nr,no_wr,ns)
        nx_s=1
        nx_part=ns
        dt_part(ns)=dt(ncell)
        dt_total=dt_part(ns)
        x_part(ns)=x_dist(nr,no_wr,ns)

 100    continue
!
!     Determine if the total elapsed travel time is equal to the
!     computational interval
!

        if(dt_total.lt.dt_comp) then
           x_part(ns)=x_part(ns)+dx(segment_cell(nr,no_wr,nx_part))
!     
!     If the particle has started upstream from the boundary point, give it
!     the value of the boundary
!

          if(x_part(ns).ge.x_bndry) then
            x_part(ns)=x_head
            dt_part(ns)=dt(segment_cell(nr,no_wr,nx_part))
            dt_total=dt_total+dt_part(ns)
            go to 200
          end if
!
!     Increment the segment counter if the total time is less than the
!     computational interval
!
          nx_s=nx_s+1
          nx_part=nx_part-1
          dt_part(ns)=dt(segment_cell(nr,no_wr,nx_part))
          dt_total=dt_total+dt_part(ns)
          go to 100
        else
!
!     For the last segment of particle travel, adjust the particle location
!     such that the total particle travel time is equal to the computational
!     interval.
!

          dt_before=dt_part(ns)
          xpart_before =x_part(ns)
          dt_part(ns)=dt_comp-dt_total+dt_part(ns)
          x_part(ns)=x_part(ns)+u(segment_cell(nr,no_wr,nx_part))
     &              *dt_part(ns)
          if(x_part(ns).ge.x_head) then
            x_part(ns)=x_head
            nx_s=nx_s-1
            dt_part(ns)=dt(head_cell(nr,no_wr))
          end if
        end if
 200    continue
        if(nx_part.lt.1) nx_part=1
        nstrt_elm(ns)=nx_part
        no_dt(ns)=nx_s
      end do
      DONE=.FALSE.
!
! Begin simulations of temperature in this water resource unit
!
      do ns=1,no_celm(nr,no_wr)
        ncell=segment_cell(nr,no_wr,ns)
!
!     Net solar radiation (kcal/meter^2/second)
!
!     qns
!
!     Net atmospheric radiation (kcal/meter^2/second)
!
!     qna
!
!     Dry bulb temperature (deg C)
!
!     dbt
!
!     Wind speed (meters/second)
!
!     wind
!      time=day_fract+hr_fract
!     Factor for Bowen ratio ((deg C)^-1)
!
!     pf
!
!     Vapor pressure at given air temperature (mb)
!
!     ea
!
!     Photo period (fraction of a day.  Not used in the energy budget)
!
!     phper
!


 250  continue
!
!     Now do the third-order interpolation to
!     establish the starting temperature values
!     for each parcel
!
        nseg=nstrt_elm(ns)
        npndx=1
!
!     If starting element is the first one, then set
!     the initial temperature to the boundary value
!
        if (nseg.eq.1) then
          t0=TEMP_head(nr,no_wr)
        else
!
!     Perform polynomial interpolation
!
          do ntrp=1,nterp(npndx)
            npart=nseg+ntrp+ndltp(npndx)-1
            xa(ntrp)=x_dist(nr,no_wr,npart)
            ta(ntrp)=temp(nr,no_wr,npart,n1)
          end do
          x=x_part(ns)
  280   continue
!
!     Call the interpolation function
!

          T0=tntrp(xa,ta,x,nterp(npndx))
          ttrp=T0
!
! End of headwaters or interpolation block
!
        end if
!
        dt_calc=dt_part(ns)
        nncell=segment_cell(nr,no_wr,nstrt_elm(ns))
!
!    Set NCELL0 for purposes of tributary input
!
        ncell0=nncell
        dt_total=0.0
        do nm=no_dt(ns),1,-1         
          u_river=u(nncell)/3.2808
          z=depth(nncell)/3.2808
!
! Calculate the transfer of energy across the air-water interface
!
          call energy(t0,QSURF,wind_fctr,A,B,nncell)
!
! Bed conduction
!
!          q_sed = rho_sed*Cps*alpha_s*(T_sed - T0)/H_sed
          q_sed = 0.0
          dt_total=dt_total+dt_calc
          dvsr = 1.0/(z*rfac)
          alpha_1 = Rate_eq*dvsr
          alpha_2 = Kappa_bed*dvsr
          qdot=(qsurf + q_sed)/(z*rfac) 
          T0 = T0+qdot*dt_calc
          t00 = t0
          qd_calc = qdot*dt_calc 
          if(t0.lt.0.0) t0=0.0
 400      continue
!
!     Look for a tributary.
! 
          q1      = qin(nncell)
          q2      = qin(nncell)
          Q_advct = 0.0    
          ntribs = no_tribs(nncell)
!
          do while (ntribs.gt.0.and..not.DONE)
            do ntrb=1,ntribs
              nr_trib = trib(nncell,ntrb)
              q2      = q2+q_trib(nr_trib)
              Q_advct = Q_advct + q_trib(nr_trib)*T_trib(nr_trib)
            end do
            DONE=.TRUE.
          end do
          t0 = (q1*t0 + Q_advct)/q2
          lat_flow = qout(nncell) - q2
          if (nyear .eq. 2000) write(60,*) 'q1,q2,qout(nncell),lat_flow'
     &                   ,nr,ncell,nncell,q1,q2,qout(nncell),lat_flow 
!          if (lat_flow.gt.0) then
!
!  Modified nonpoint source temperature so as to be the same
!  as the instream simulated temperature for Connecticut River 7/2015
!            T_dist=t0
!            Q_advct = Q_advct + lat_flow*T_dist
!          end if
 500      continue
!

!
          nseg=nseg+1
          nncell=segment_cell(nr,no_wr,nseg)
!
!     Reset tributary flag is this is a new cell
! 
          if (ncell0.ne.nncell) then
            ncell0=nncell
            DONE=.FALSE.
          end if
          dt_calc=dt(nncell)
        end do
        if (t0.lt.0.5) t0=0.5
        temp(nr,no_wr,ns,n2)=t0
!
!   Write file 20 with all temperature output 11/19/2008
! 
        ndout=nd
        nyear_out=nyear
!        time=year+(day-1.+hour_inc*period)/xd_year
        if (ndd .eq. nwpd)then
          time=year + day/xd_year
          ndout = nd
          if (nd .eq.nd_year) then
            ndout = 1
            nyear_out=nyear+1
          end if
        end if
        qsw_out=qns(ncell)
        xkm_plot=AMAX1(0.001,x_dist(nr,no_wr,ns)*3.048e-04)
!
! Write output to unit 20
!
        write(20,'(f11.5,i5,1x,2i4,1x,5i5,1x,5f7.2,4f9.2,2f9.1)') 
     &             time,nyear_out,ndout,ndd,nr,no_wr,ncell,ns,nseg
     &            ,t0,TEMP_head(nr,no_wr),dbt(ncell)
     &            ,depth(ncell),u(ncell),qin(ncell),q1,q2
     &            ,lat_flow,xkm_plot,qsw_out
!
!  write output for longitudinal plots, if requested
!
!        if (ndout .eq. jjl_day .and. nr .eq. nnl_rch) then
!          xll_dist = x_dist(nr,no_wr,ns)*ft_to_km
!          nll_unit = 70
!          write(nll_unit,'(4i5,f10.3,f10.2)') 
!     &          nyear,nnl_day,ncell,nseg,xll_dist,T0
!        end if
!
!     End of computational element loop
!
      end do
!
! Reset NSEG back one for the next water resource unit
!
      nseg = ns
!
      RETURN
      END
