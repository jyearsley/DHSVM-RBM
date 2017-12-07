c
c      PROGRAM RMAIN
C
C     Dynamic river basin model for simulating water quality in
C     branching river systems with freely-flowing river segments. 
c
c     This version uses Reverse Particle Tracking in the Lagrangian
c     mode and Lagrangian interpolation in the Eulerian mode.
c
c     Topology and routing is set up to be consistent with output
c     from the Distributed Hydrologic Soil and Vegetation Model (DHSVM)
c     model developed by the Land Surface Hydrology Group at 
c     the University of Washington.
c
C     For additional information visit:
c
c     http://www.hydro.washington.edu/Lettenmaier/Models/DHSVM/
c
c     or contact:
c
C     John Yearsley
C     Land Surface Hydrology Group
C     Dept. of Civil and Environmental Engineering
C     Box 352700
C     University of Washington
C     Seattle, Washington
C     98195-2700
C     yearsley@hydro.washington.edu
C
      character*8  start_data,end_data     
      character*200 Prefix1, Prefix2
      integer iargc
      integer numarg
 
c     Command line input
c
      numarg = iargc ( )
      if (numarg .lt. 2) then
        write (*,*) 'Too few arguments were given'
        write (*,*) ' '
        write (*,*) 'First:  Location and prefix of input file'
        write (*,*) '        (networkfile)'
        write (*,*) 'Second:  Location and prefix of output file'
        write (*,*) '        (networkfile)'
        write (*,*) 'eg: $ <program-name> <Project Name> <Output File>'
        write (*,*) ' '
        stop
      end if
      call getarg ( 1, Prefix1 )
      call getarg ( 2, Prefix2 )
c
c     Identify and open necessary files
c
c     open the output file 
      open(unit=20,file=TRIM(Prefix2)//'.temp',status='unknown')
c
c     Open file with weather and infPlow data
      write(*,*) 'Forcing file -  ', TRIM(Prefix1)//'.forcing'
      open(unit=30,file=TRIM(Prefix1)//'.forcing',STATUS='old')
C
c     open Mohseni file 
      open(40,file=TRIM(Prefix1)//'.Mohseni',STATUS='old')    
c
c     open Leopold file 
      open(50,file=TRIM(Prefix1)//'.Leopold',STATUS='old')    
c
c     Open network file
      OPEN(UNIT=90,FILE=TRIM(Prefix1)//'.net',STATUS='OLD')
c
c     Read header information from control file
      do n=1,5
        read(90,*)
      end do
c
C     Call systems programs to get started
C
C     SUBROUTINE BEGIN reads control file, sets up topology and
C     important properties of reaches
      write(*,*) 'Calling BEGIN'
      CALL BEGIN
C
C     SUBROUTINE SYSTMM performs the simulations
C
      CALL SYSTMM
C
C     Close files after simulation is complete
C
      write(*,*) ' Closing files after simulation'
      CLOSE(30)
      CLOSE(90)
      STOP
      END
      SUBROUTINE BEGIN
      character*11 end_time,start_time
      character*5 Dummy_B
      character*5 wr_type
      character*10 Dummy_A
      integer head_name,trib_cell,first_cell
      dimension ndmo(12)
      logical Is_a_Riv,Is_a_Res
      logical Test
      INCLUDE 'RBM.fi'
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
c      ndelta=2
c      delta_n=ndelta
c
c     Read the starting and ending times and the number of
c     periods per day of weather data from the forcing file
c 
c 
      read(30,*) start_time,end_time,nwpd,nd_start
c
c Kludge here for averaged forcing 
C
      nwpd = 1
      write(*,*) start_time,'  ',end_time,nwpd,nd_start
c 
      write(*,*) 'Number of simulations per day - ',nwpd
c
      read(start_time,'(i4,2i2,1x,i2)') start_year,start_month
     &                                 ,start_day,start_hour
      read(end_time,'(i4,2i2,1x,i2)') end_year,end_month
     &                               ,end_day,end_hour
      nyear1=start_year
      nyear2=end_year
      write(*,*) start_year,start_month,start_day
     &             ,end_year,end_month,end_day
c
c     Establish the Julian day for which simulations begin
c 
      jul_start=julian(start_year,start_month,start_day)
c
      read(90,*)  nreach
      write(*,*) "Number of stream reaches - ", nreach
      read(40,*) Test,a_smooth
      write(*,*) Test,a_smooth, nreach
      b_smooth=1.- a_smooth
      do nr=1, nreach
        read(40,*) nnr,alf_Mu(nr),beta(nr),gmma(nr),mu(nr)
        write(*,*) nnr,alf_Mu(nr),beta(nr),gmma(nr),mu(nr)
c
      end do
c
c Iniialize the upstream cell index
c
      do n = 1,500
        do nn = 1,50
          upstrm_cell(n,nn) = -1
        end do
      end do
c
c     Start reading the reach date and initialize the reach index, NR
c     and the cell index, NCELL
c
      ncell    = 0
      nrch     = 1
      no_res   = 0
      ns_total = 0
      nseg     = 0
  100 continue
C
C     Card Group IIb. Reach characteristics
C
       do while (nrch .le. nreach) 
c
c Initialize index that counts the number of water resource units in this reach
c
         no_wru_0  = 0
c
c     Initialize NSEG, the total number of segments in this reach
      nseg=0
      write(*,*) ' Starting to read reach ', nrch
c
c     Read the number of cells in this reach, the headwater #,
c     the number of the cell where it enters the next higher order stream,
c     the headwater number of the next higher order stream it enters, and
c     the river mile of the headwaters.
c
c
        read(90,*) Dummy_A,no_cells( nrch)
     &            ,Dummy_A,main_stem( nrch),Dummy_A,trib_cell
c
c     If this is reach that is tributary to cell TRIB_CELL, give it the
c     pointer TRIB(TRIB_CELL) the index of this reach for further use.
c     Also keep track of the total number of tributaries for this cell
c
        if (trib_cell.gt.0) then
           no_tribs(trib_cell)=no_tribs(trib_cell)+1
           trib(trib_cell,no_tribs(trib_cell))= nrch
        end if
c
c     Reading Reach Element information
c
c
c 
        no_wr_units(nr) = 0
        Is_a_Riv = .FALSE.
        Is_a_Res = .FALSE.
        do nc=1,no_cells( nrch)
          ncell=ncell+1
          read(50,*) ncll,U_a(ncell),U_b(ncell),U_min(ncell)
          if (ncll .ne.ncell) then
            write(*,*) 'Mismatch in Leopold file'
          end if  
          read(50,*) D_a(ncell),D_b(ncell),D_min(ncell)


c
c     The headwaters index for each cell in this reach is given
c     in the order the cells are read*dt_calc/(z*rfac)
c
C     Card Type 3. Cell indexing #, Node # Row # Column Lat Long RM
C

          read(90,*) Dummy_B,nnc,Dummy_B,node(nnc)
     &              ,Dummy_B,x_0,Dummy_B,x_1
     &              ,Dummy_B,elev( nrch),ndelta(ncell)
     &             ,Dummy_B,no_wru,wr_type
c
          unit_type(nrch,no_wru) = wr_type
c
        
c
c Set reservior indices and geometry
c
          if (upstrm_cell(nrch,no_wru) .lt. 0) then
            upstrm_cell( nrch,no_wru) = ncell-1
          end if
          if (no_wru .ne. no_wru_0) then
            nseg = 0
            head_cell(nrch,no_wru) = ncell
            no_celm(nrch,no_wru) = 0
            cells_wr(nrch,no_wru)  = 0
            x_dist(nrch,no_wru,nseg) = x_1
c
c If this is a reservoir, update reservoir index
c
            if (wr_type .eq. 'RSRVR') then
              no_res = no_res + 1
            end if
c
c Update water resource unit index
c
            no_wru_0 = no_wru
          end if 
          do ncll = 1,ndelta(ncell)
            nseg = nseg + 1
            if (wr_type .eq. 'RSRVR') then
              no_res = no_res + 1
c            res_no(nrch,no_wru)   = no_res
              volume(no_res,nseg,1) = U_a(ncell)
              volume(no_res,nseg,2) = U_b(ncell)
              a_surf(no_res,nseg,1) = D_a(ncell)
              a_surf(no_res,nseg,2) = D_b(ncell)
              Kappa(no_res)         = D_min(ncell)
            end if             
            no_celm(nrch,no_wru) = no_celm(nrch,no_wru)+1
            segment_cell(nrch,no_wru,nseg)=ncell
            dx(ncell)=(x_1-x_0)/ndelta(ncell)
            x_dist(nrch,no_wru,nseg) 
     &      = x_dist(nrch,no_wru,nseg-1)-dx(ncell)
            write(85,*) wr_type,no_wru,nseg
     &                 ,upstrm_cell(nrch,no_wru)
     &                 ,ncell,cells_wr(nrch,no_wru)
     &                 ,segment_cell(nrch,no_wru,nseg)
     &                 ,x_dist(nrch,no_wru,nseg)
     &                 ,x_dist(nrch,no_wru,nseg-1),x_1,x_0,dx(ncell)                        
          end do
                    cells_wr(nrch,no_wru) 
     &                   = cells_wr(nrch,no_wru) + 1
          write(85,*) 
          write(85,*) 'WRU cells ',nrch,no_wru,cells_wr(nrch,no_wru)
c                      cells_wr(nrch,no_wru) = ncell
          write(85,*)                          
c

  200 continue
c 
c   
        no_wr_units(nrch) = no_wru
        write(85,*) 'Total WRU in NR ',no_wru
c
      end do
c
c
c	go to 100

      nrch = nrch + 1

      end do
c  500	continue
      nrch= nreach
      xwpd=nwpd
      dt_comp=86400./xwpd
C
C     ******************************************************
C                         Return to RMAIN
C     ******************************************************
C
c
      RETURN
  900 END
      SUBROUTINE SYSTMM
c
      integer res_nn
      integer ndmo(12,2)
c 

      INCLUDE 'RBM.fi'
      data lat/47.6/,pi/3.14159/,rfac/304.8/
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334
     &         ,0,31,60,91,121,152,182,213,244,274,305,335/
c
c
      hour_inc=1./nwpd
      do nr=1,500
        T_smth(nr)=mu(nr)
        do no_wr = 1,10
          T_head(nr,no_wr)=mu(nr)
        end do
      end do
c
      n1=1
      n2=2
      nobs=0
c
c     Initialize the day counter used for calculating the
c     record position for direct access files
c
c
      write(*,*) 'Start_year ',start_year,start_month,start_day
      write(*,*) 'End_year   ',end_year,end_month,end_day
      ndays=-Julian(start_year,start_month,start_day)
     &     +Julian(end_year,end_month,end_day)+1
      write(*,*) 'Number of days ',ndays
c      
c     Setup the timing of the simulation
c 
      lp_year=1
      yr_days=365.
      if (mod(start_year,4).eq.0) then
        lp_year=2
        yr_days=366.
      end if
      day_fract=ndmo(start_month,lp_year)+start_day-1
      day_fract=day_fract/yr_days
      hr_fract=start_hour/(24.*yr_days)
      sim_incr=dt_comp/(365.*86400.)
      year=start_year
      time=year+day_fract+hr_fract
c Temporary for debugging
      end_year =start_year+1
c
c     Year loop starts
      do nyear=start_year,end_year
         write(*,*) ' Simulation Year - ',nyear
         nd_year=365
         lp_year=1
         if (mod(nyear,4).eq.0) then 
           nd_year=366
           lp_year=2
         end if
         xd_year=nd_year
c
c        Day loop starts
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
c
c
         DO ND=nd1,nd2
c
c     Start the numbers of days-to-date counter
           ndays=ndays+1
c
c     Daily period loop starts
           DO ndd=nd_start,nwpd 
c
c     Begin reach computations
c      
             ind=nd
             ipd=ndd
             day=nd
             period=ndd
             time=year+(day-1.+hour_inc*period)/xd_year
c
c Initialize number of reservoirs (index is basin-wide)
c
             res_nn = 0
             ncell = 0
c
             do nr=1,nreach
             write(29,*) 'reach ',nr,no_wr_units(nr)
               IS_HEAD = .TRUE.
               nrch = nr
               nseg = 0
               do no_wr = 1,no_wr_units(nr)
                write(29,*) 'Systemm ',no_wr,unit_type(nrch,no_wr)
c
c Select either the RIVER or RSRVR case
c
               Select Case(unit_type(nrch,no_wr))
                 case('RIVER')
                   call RIVER(nseg,ncell,nr,no_wr)
                 case('RSRVR')      
                   res_nn = res_nn + 1
                   call RSRVR(nseg,ncell,nrch,no_wr,res_nn)
               End Select
               IS_HEAD = .FALSE. 
c
c End SELECT
c

               end do
c
c Establish tributary temperatures at reach end
c
c
               nseg_trib = no_celm(nr,no_wr_units(nr))
	       T_trib(nr)=temp(nr,no_wr_units(nr),nseg_trib,n2)
               write(29,*) 'end ',nseg_trib,nrch,no_wru,T_trib(nr)
c
             end do
             ntmp=n1
             n1=n2
             n2=ntmp
c
c     End of weather period loop (NDD=1,NWPD)
c
           end do
c
c Reset daily loop counter
c
          nd_start=1
c 
C
c     End of main loop (ND=1,365/366)
c

         end do
c
c    Update initial time for new year
c
      year=year+1
c
c     End of year loop
c
      end do
c
c Finish
c
  900 Continue
c
c
c     ******************************************************
c                        return to rmain
c     ******************************************************
c

  950 return
      end
      SUBROUTINE ENERGY
     &           (TSURF,QSURF,A,B,ncell)
      REAL*4 Ksw,LVP
      real*4 q_fit(2),T_fit(2),evrate
      INCLUDE 'RBM.fi'
c      data evrate/1.5e-9/,pf/0.640/
      parameter (pi=3.14159)
      parameter (P_FCTR=64.0,RHO=1000.,EVRATE=1.5e-11)
c     
      td=nd
      evap_rate=EVRATE
      T_fit(1)=tsurf-0.5
      T_fit(2)=tsurf+0.5
      do i=1,2
         T_kelvin = T_fit(i) + 273.0
c
c Vapor pressure at water surface
         e0=2.1718E10*EXP(-4157.0/(T_kelvin-33.91))
c
c Bowen ratio
         rb=P_FCTR*(DBT(ncell)-T_fit(i))
c
c Latent heat of vaporization
         lvp=1.91846e06*(T_kelvin/(T_kelvin-33.91))**2
c
c Evaporative heat flux
         QEVAP=rho*lvp*evap_rate*WIND(ncell)
         if(qevap.lt.0.0) qevap=0.0
c
c Convective heat flux
         QCONV=rb*QEVAP
         QEVAP=QEVAP*(E0-EA(ncell))
c 
c Back radiation from the water surface
         QWS=280.23+6.1589*T_fit(i)
c
c Thermal energy budget for i = 1,2
         q_fit(i)=QNS(ncell)+0.97*QNA(ncell)-QWS-QEVAP+QCONV
      end do

      A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
      B=(T_fit(1)*q_fit(2)-T_fit(2)*q_fit(1))
     .     /(T_fit(1)-T_fit(2))
      qsurf=0.5*(q_fit(1)+q_fit(2))
      RETURN
      END
C
c
c    Subroutine that simulates water temperature in advective system
c
      SUBROUTINE RIVER(nseg,ncell,nr,no_wr)
      integer nr,no_wr
      integer no_dt(2000),nstrt_elm(2000)
     .     ,ndltp(4),nterp(4)
      logical DONE,pp_T
      real*4 xa(4),ta(4)
      real*4 lat_flow
      real*4 dt_part(2000),x_part(2000)
      INCLUDE 'RBM.fi'
      data ndltp/-2,-1,-2,-2/,nterp/4,3,2,3/
c      data pi/3.14159/,rfac/304.8/
      parameter (PI=3.14159,RFAC=4.184e6)
c
      l_cell = 0
c
      do nc=1,cells_wr(nr,no_wr)
        l_cell = l_cell+1
        ncell  = ncell + 1
        read(30,*) ntemp,l1
     &                    ,press(ncell),dbt(ncell)
     &                    ,qna(ncell),qns(ncell)
     &                    ,ea(ncell),wind(ncell)
     &                    ,qin(ncell),qout(ncell)
c
c Find the headwaters temperature
c
      nc_head=upstrm_cell(nr,no_wr)
      write(29,*) 'nchead ',nr,no_wr,nc_head,ncell
      if(ncell .eq. nc_head+1) then
        call HEAD_TEMP (ncell,nr,no_wr)
        write(28,*) 'head ',IS_HEAD,nr,no_wr,ncell,nseg
      end if
c
c Average flow in segment for estimating hydraulic conditions

        qavg=0.5*(qin(ncell)+qout(ncell))
c
c    Stream speed estimated with Leopold coefficients
c 
                 u(ncell)=U_a(ncell)*(qavg**U_b(ncell)) 
                 u(ncell) = amax1(u_min(ncell),u(ncell))
c 
c                 qdiff(ncell)=(qout(ncell)-qin(ncell))/delta_n
c 
c 
        dt(ncell)=dx(ncell)/u(ncell)
c
c    Depth estimated with Leopold coefficients
c 
        depth(ncell)=D_a(ncell)*(qavg**D_b(ncell))
        depth(ncell)=amax1(D_min(ncell),depth(ncell))
c
c      end do
c
c  Set the value of the tributary flow at the downstream segment
c  of the RIVER water resource segment
c
        q_trib(nr)=qout(ncell)
c
c     Main stem inflows and outflows for each reach first
c     Flows are cumulative and do not include tributaries if
c     tributaries are downstream of the inflow junction
c
      end do
c
c Temporary method for fixing headwaters location - best done in Begin.f90
c 
      x_head=x_dist(nr,no_wr,0)
      x_bndry=x_head-1.0


c   Do the reverse particle tracking
c
c      write(28,*) 'nr,no_celm ',nr,no_wr,no_celm(nr,no_wr)
      do ns=no_celm(nr,no_wr),1,-1
c
c     Segment is in cell SEGMENT_CELL(NC)
c

        ncell=segment_cell(nr,no_wr,ns)
        nx_s=1
        nx_part=ns
        dt_part(ns)=dt(ncell)
        dt_total=dt_part(ns)
        x_part(ns)=x_dist(nr,no_wr,ns)
 100    continue
c
c     Determine if the total elapsed travel time is equal to the
c     computational interval
c

        if(dt_total.lt.dt_comp) then
           x_part(ns)=x_part(ns)+dx(segment_cell(nr,no_wr,nx_part))
c     
c     If the particle has started upstream from the boundary point, give it
c     the value of the boundary
c

          if(x_part(ns).ge.x_bndry) then
            x_part(ns)=x_head
            dt_part(ns)=dt(segment_cell(nr,no_wr,nx_part))
            dt_total=dt_total+dt_part(ns)
            go to 200
          end if
c
c     Increment the segment counter if the total time is less than the
c     computational interval
c
          nx_s=nx_s+1
          nx_part=nx_part-1
          dt_part(ns)=dt(segment_cell(nr,no_wr,nx_part))
          dt_total=dt_total+dt_part(ns)
          go to 100
        else
c
c     For the last segment of particle travel, adjust the particle location
c     such that the total particle travel time is equal to the computational
c     interval.
c

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
c
c Begin simulations of temperature in this water resource unit
c
      do ns=1,no_celm(nr,no_wr)
        ncell=segment_cell(nr,no_wr,ns)
c
c     Net solar radiation (kcal/meter^2/second)
c
c     qns
c
c     Net atmospheric radiation (kcal/meter^2/second)
c
c     qna
c
c     Dry bulb temperature (deg C)
c
c     dbt
c
c     Wind speed (meters/second)
c
c     wind
c      time=day_fract+hr_fract
c     Factor for Bowen ratio ((deg C)^-1)
c
c     pf
c
c     Vapor pressure at given air temperature (mb)
c
c     ea
c
c     Photo period (fraction of a day.  Not used in the energy budget)
c
c     phper
c


 250  continue
c
c     Now do the third-order interpolation to
c     establish the starting temperature values
c     for each parcel
c
        nseg=nstrt_elm(ns)
        npndx=1
c
c     If starting element is the first one, then set
c     the initial temperature to the boundary value
c
        if (nseg.eq.1) then
          t0=T_head(nr,no_wr)
        else
c
c     Perform polynomial interpolation
c
          do ntrp=1,nterp(npndx)
            npart=nseg+ntrp+ndltp(npndx)-1
            xa(ntrp)=x_dist(nr,no_wr,npart)
            ta(ntrp)=temp(nr,no_wr,npart,n1)
          end do
          x=x_part(ns)
  280   continue
c
c     Call the interpolation function
c

          t0=tntrp(xa,ta,x,nterp(npndx))
          ttrp=t0
c
c End of headwaters or interpolation block
c
        end if
c
        dt_calc=dt_part(ns)
        nncell=segment_cell(nr,no_wr,nstrt_elm(ns))
c
c    Set NCELL0 for purposes of tributary input
c
        ncell0=nncell
        dt_total=0.0
        do nm=no_dt(ns),1,-1         
          u_river=u(nncell)/3.2808
          z=depth(nncell)/3.2808
c
c Calculate the transfer of energy across the air-water interface
c
          call energy(t0,QSURF,A,B,nncell)
          t_eq=-B/A
          qdot=qsurf/(z*rfac)
c          write(28,*) 'nd qdot dt_calc ',nd,nr,no_wr,nncell
          dt_calc=dt_part(nm)
          dt_total=dt_total+dt_calc
          too = t0
          t0=t0+qdot*dt_calc 
          if(t0.lt.0.0) t0=0.0
 400      continue
c
c     Look for a tributary.
c 
          q1      = qin(nncell)
          q2      = qin(nncell)
          Q_advct = 0.0    
          ntribs = no_tribs(nncell)
          do while (ntribs.gt.0.and..not.DONE)
            do ntrb=1,ntribs
              nr_trib = trib(nncell,ntrb)
              q2      = q2+q_trib(nr_trib)
              Q_advct = Q_advct + q_trib(nr_trib)*T_trib(nr_trib)
            end do
            DONE=.TRUE.
          end do
c
          lat_flow = qout(nncell) - q2
          if (lat_flow.gt.0) then
c
c  Modified nonpoint source temperature so as to be the same
c  as the instream simulated temperature for Connecticut River 7/2015
            T_dist=t0
            Q_advct = Q_advct + lat_flow*T_dist
          end if
 500      continue
c

c
          nseg=nseg+1
          nncell=segment_cell(nr,no_wr,nseg)
c
c     Reset tributary flag is this is a new cell
c 
          if (ncell0.ne.nncell) then
            ncell0=nncell
            DONE=.FALSE.
          end if
          dt_calc=dt(nncell)
        end do
        if (t0.lt.0.5) t0=0.5
        temp(nr,no_wr,ns,n2)=t0
c
c   Write file 20 with all temperature output 11/19/2008
c 
        ndout=nd
        nyear_out=nyear
        time=year+(day-1.+hour_inc*period)/xd_year
        if (ndd .eq. nwpd)then
          time=year + day/xd_year
          ndout = nd + 1
          if (nd .eq.nd_year) then
            ndout = 1
            nyear_out=nyear+1
          end if
        end if
        qsw_out=qns(ncell)
        rmile_plot=x_dist(nr,no_wr,ns)/5280.
c
c Write output to unit 20
c
        write(20,'(f11.5,i5,1x,2i4,1x,4i5,1x,5f7.2,f9.2,f9.1)') 
     &                  time,nyear_out,ndout,ndd,no_wr,ncell,ns,nseg
     &                 ,t0,T_head(nr,no_wr),dbt(ncell)
     &                 ,depth(ncell),u(ncell),qin(ncell),qsw_out
c                     end if
c
c 
c Update NSEG
c

c     End of computational element loop
c
      end do
c
c Reset NSEG back one for the next water resource unit
c
      nseg = ns
c      ntmp=n1
c      n1=n2
c      n2=ntmp
      RETURN
      END
c
c Reservoir subroutine
c
      SUBROUTINE RSRVR(nseg,ncell,nr,no_wr,res_nn)
c
c Reservoir inlows and advected thermal energy
c
      real*4 q_in_res(2),q_out_res(2)
      real*4 T_in(2),T_res(500,2,2),Q_advect(2)

!
! Reservoir temperatures are saved as T_res(m1,m2,m3,m4)
! m1 = reach #, m2 = water resource unit #, m3 = res #, m4 = layer #, m5 = time index
!
      integer res_nn
      logical DONE
      real*4 Area_sum(2),Vol_sum(2)
c      real*4 T_res(500,2,2)
      INCLUDE 'RBM.fi'
      SAVE T_res
      data Pi/3.1415927/,rho_cp/4.186e06/,cuft_cum/0.028318/
c
c Find the input temperature
c
      call HEAD_TEMP (ncell+1,nr,no_wr)

c
c Initialize  reservoir properties
c
      Area_sum(:) = 0.0
      Vol_sum(:)  = 0.0

c
c Initialize epilimnion and hypolimnion temperatures
c
      T_epi = T_res(res_nn,1,n1)
      T_hyp = T_res(res_nn,2,n1)
    
c*************************************************************
c Convert upstream inflow to m**3/second
c
      q_inflow  = cuft_cum*qout(ncell)
c
c*************************************************************
c
c
c Use inflow temperature to determine placement of upstream input
c
      T_inflow    = T_head(nr,no_wr)
      q_in_res(:) = 0.0
      T_in(:) = 0.0
c
      layer = NLAYER(T_inflow,T_epi)
      T_in(layer)     = T_inflow
      q_in_res(layer) = q_inflow
c
c
      Q_advect(1) = q_in_res(1)*T_in(1)
      Q_advect(2) = q_in_res(2)*T_in(2)
c       
      nseg = 0
c 
c Initialize important reservoir properties
c
      Q_netsurf    = 0.0

c
c Cycle through the cells (NCELL) and segments (NSEG)
c
      do nc = 1,cells_wr(nr,no_wr)
c
c Update segment and cell #'s
c
        ncell  = ncell + 1  
        nseg   = nseg + 1
        DONE   = .FALSE.
        write(29,*) 'Reservoir ',ncell,nseg,res_nn
c
c Accumulate surface areas and volumes of each reservoir cell
c
        do nl = 1,2
          Area_sum(nl) = Area_sum(nl) + a_surf(nr,res_nn,nl)
          Vol_sum(nl)  = Vol_sum(nl) + volume(nr,res_nn,nl)
        end do
c
c Read the forcing file
c
        read(30,*) ntemp,l1
     &                    ,press(ncell),dbt(ncell)
     &                    ,qna(ncell),qns(ncell)
     &                    ,ea(ncell),wind(ncell)
     &                    ,qin(ncell),qout(ncell)
c
c Call energy budget routine
c
        Call ENERGY(T_epi,Q_surface,A,B,ncell)
        Q_netsurf = Q_netsurf + a_surf(res_nn,nseg,1)*Q_surface/rho_cp
c
c     Look for a tributary.
c 
          q1      = qin(ncell)
          q2      = qin(ncell)
          q_tmp   = 0.0
c
          Q_advct_trb = 0.0  
          T_tmp       = 0.0  
          ntribs = no_tribs(ncell)
          do while (ntribs.gt.0.and..not.DONE)
            do ntrb=1,ntribs
              nr_trib = trib(ncell,ntrb)
              q2      = q2+q_trib(nr_trib)
              q_tmp   = q_tmp + q_trib(nr_trib)
              Q_advct_trb = Q_advct_trb 
     &                    + q_trib(nr_trib)*T_trib(nr_trib)
              T_tmp = T_tmp + q_trib(nr_trib)*T_trib(nr_trib)
            end do
            DONE=.TRUE.
          end do
c
c Add the advected thermal energy from the tributaries in this cell
c
          T_tmp = T_tmp/q_tmp
          layer = NLAYER(T_tmp,T_epi)
          Q_advect(layer) = Q_advect(layer) + Q_advct_trb
c
c End of the loop that accumlates area, volume, surface transfer and advected sources
c from each reservoir cell
c      
      end do
c
c Reservoir outflows are based on the results from the last reservoir cell
c
        q_out_res(1) = 0.0
        q_out_res(2) = qout(ncell)
        q_total =  q_out_res(1) + q_out_res(2)

c
c Hypolimnion
c
        q_vert = q_in_res(1)
c
        T_res(res_nn,2,n2) = T_res(res_nn,2,n1)
     &        + ((q_vert*T_res(res_nn,1,n1) + Q_advect(2)
     &        - q_out_res(2)*T_res(res_nn,2,n1))/ Vol_sum(2))
     &      * dt_comp 
        if (T_res(1,2,n2) .lt. 0.0) T_res(1,2,n2) = 0.0
!
! Epilimnion
!
        T_res(res_nn,1,n2) = T_res(res_nn,1,n1) 
     &        + ((Q_netsurf 
     &        +  q_in_res(1)*T_in(1)- q_out_res(1)*T_res(res_nn,1,n1))
     &        / Vol_sum(1))*dt_comp 
c         a1 = q_surface*a_surf(res_nn,1)*dt_comp/V_epi
c         a2 = Q_in_res(1)*T_in(1)*dt_comp/V_epi
c         a3 = Q_out_res(1)*T_res(nres_nn,1,n1)*dt_comp/V_epi
c         a4 = Q_in_res(1)*T_in(1)*dt_comp/V_epi
c
        if (T_res(res_nn,1,n2) .lt. 0.5) T_res(res_nn,1,n2) = 0.5
c
c
c
        T_out = (q_out_res(1)*T_res(res_nn,1,n2) 
     &        +  q_out_res(2)*T_res(res_nn,2,n2)) / q_total
        temp(nr,no_wr,nseg,n2) = T_out
        write(20,'(f11.5,i5,1x,2i4,1x,4i5,1x,5f7.2,f9.1,f9.1)') 
     &                       time,nyear,nd,ndd,no_wr,ncell,nseg,nseg
     &                      ,T_out,T_out,dbt(ncell)
     &                      ,depth(ncell),u(ncell),qin(ncell)
     &                      ,q_surface
c
c      ntmp=n1
c      n1=n2
c      n2=ntmp
c
      return
      end
c
c
c Headwaters subroutine
c
      SUBROUTINE HEAD_TEMP (ncell,nr,no_wr)
!
!
      INCLUDE 'RBM.fi'
c
c     Headwaters flow and temperature
c
      write(28,*) IS_HEAD,ncell,nr,no_wr
      if(IS_HEAD) then
        T_smth(nr)=b_smooth*T_smth(nr)
     &            +a_smooth*dbt(ncell)
        T_head(nr,no_wr)=mu(nr)
     &            +(alf_Mu(nr)
     &            /(1.+exp(gmma(nr)*(beta(nr)-T_smth(nr))))) 
      else
        nwr_1 = no_wr-1
        T_head(nr,no_wr) = temp(nr,nwr_1,no_celm(nr,nwr_1),n1)
      end if
c        
      temp(nr,no_wr,0,n1)=T_head(nr,no_wr)
      temp(nr,no_wr,-1,n1)=T_head(nr,no_wr)
      temp(nr,no_wr,-2,n1)=T_head(nr,no_wr)
      temp(nr,no_wr,no_celm(nr,no_wr)+1,n1)
     &     =temp(nr,no_wr,no_celm(nr,no_wr),n1)
      return
      end
c
c   Function to calculate the number of simulation days
c
      function nodays(jtime,jy0)
      dimension ndmo(12)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
      jy=(jtime/10000)
      jrem=(jtime-jy*10000)
      jm=jrem/100
      jd=jrem-jm*100
      ny=jy-jy0
      nodays=365*ny+ndmo(jm)+jd
      return
      end
      function ndate(n,ny0)
      dimension ndmo(13)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334,365/
      ny=(n-1)/365
      njul=n-ny*365
      nm=1
 10   continue
      if(ndmo(nm+1).ge.njul) go to 50
      nm=nm+1
      go to 10
 50   continue
      nday=njul-ndmo(nm)
      nmon=100*nm
      nyear=10000*(ny0+ny)
      ndate=nyear+nmon+nday
      return
      end

c	Third-order polynomial interpolation using Lagrange
c     polynomials.  FUNCTION is SUBROUTINE POLINT from
c     Numerial Recipes
c
      FUNCTION tntrp(XA,YA,X,n)
c      PARAMETER (N=4)
c      DIMENSION XA(N),YA(N),C(N),D(N)
      DIMENSION XA(4),YA(4),C(4),D(4)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) DEN=0.001
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
	  tntrp=y
      RETURN
      END
c
      INTEGER FUNCTION Julian (YEAR,MONTH,DAY)
C
C---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
C   DATE (YEAR,MONTH,DAY).
C
      INTEGER YEAR,MONTH,DAY,I,J,K
C

      I= YEAR
      J= MONTH
      K= DAY
C

      Julian=
     1   K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)
     2  /12-3*((I+4900+(J-14)/12)/100)/4
C

      RETURN
      END
c
      INTEGER FUNCTION NLAYER(T1,T2)
         DENSITY1 = ((((6.536332E-9*T1-1.120083E-6)*T1+1.001685E-4)
     &           * T1-909529E-3)*T1+6.793952E-2)*T1+0.842594
         DENSITY1 = DENSITY1 + 999.0
         DENSITY2 = ((((6.536332E-9*T2-1.120083E-6)*T2+1.001685E-4)
     &           * T2-909529E-3)*T2+6.793952E-2)*T2+0.842594
         DENSITY2 = DENSITY2 + 999.0
         NLAYER = 1
         IF (DENSITY2 .GT. DENSITY1) NLAYER = 2
      RETURN
      END   








