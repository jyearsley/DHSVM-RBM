!
! Reservoir subroutine
!
      SUBROUTINE RSRVR(nseg,nseq,nr,no_wr,res_nn)
      use Block_Hydro
      use Block_Network
      use Block_WQ
!
      implicit none
!
! Reservoir inlows and advected thermal energy
!
      integer, parameter                  :: dp=kind(0.d0)
      integer                             :: res_nn,res_seg
!
      logical                             :: DONE
!
      real*4                              :: kappa_z
      real, dimension(:), allocatable     :: q_advct_trb,q_in_res,q_out_res,q_in_trb,q_vert
      real, dimension(:), allocatable     :: area_sum,vol_sum,T_in,Q_advect
      real, dimension(:,:,:), allocatable :: temp_res
!
      real(dp)                            :: density1,density2

!
! Reservoir temperatures are saved as TEMP_res(m1,m2,m3,m4)
! m1 = reach #, m2 = water resource unit #, m3 = res #, m4 = layer #, m5 = time index
!
      SAVE TEMP_res
!
! Find the input temperature
!
      call HEAD_TEMP (nseq+1,nr,no_wr)
!
! Initialize  reservoir properties
!
      Area_sum(:) = 0.0
      Vol_sum(:)  = 0.0

!
! Initialize epilimnion and hypolimnion temperatures
!
      T_epi = TEMP_res(res_nn,1,n1)
      T_hyp = TEMP_res(res_nn,2,n1)
    
!*************************************************************
!
! Convert upstream inflow to m**3/second
!
      q_inflow  = cuft_cum*qout(nseq)
!
!*************************************************************
!
!
! Use inflow temperature to determine placement of upstream input
!
      T_inflow    = TEMP_head(nr,no_wr)
      q_in_res(:) = 0.0
      T_in(:) = 0.0
!
      nsq_1 = nseq + 1
!
!      layer = 1
!      call DENSE(DENSITY1,DENSITY2,T_inflow,T_epi)
!      if (DENSITY1 .gt. DENSITY2) layer = 2
!
      layer = NLAYER(T_inflow,T_epi)
      T_in(layer)     = T_inflow
      q_in_res(layer) = q_inflow
!
!
      Q_advect(1) = q_in_res(1)*T_in(1)
      Q_advect(2) = q_in_res(2)*T_in(2)
!       
      res_seg = 0
! 
! Initialize important reservoir properties
!
      Q_netsurf    = 0.0
!
! Vertical eddy diffusivity
!
      Q_dff_epi    = 0.0
      Q_dff_hyp    = 0.0
!
! Cycle through the cells (nseq) and segments (NSEG)
!
      do nc = 1,cells_wr(nr,no_wr)
!
! Update segment and cell #'s
!
        nseq     = nseq + 1  
        res_seg   = res_seg + 1
        DONE   = .FALSE.
!
! Accumulate surface areas and volumes of each reservoir cell
!
        do nl = 1,2
          Area_sum(nl) = Area_sum(nl) + a_surf(res_nn,res_seg,nl)
          Vol_sum(nl)  = Vol_sum(nl) + volume(res_nn,res_seg,nl)
        end do
!
! Read the forcing file
!
        read(30,*) l1
     &                    ,press(nseq),dbt(nseq)
     &                    ,qna(nseq),qns(nseq)
     &                    ,ea(nseq),wind(nseq)
     &                    ,qin(nseq),qout(nseq)
!
!
! Call energy budget routine
!
        Call ENERGY(T_epi,Q_surface,wind_fctr,A,B,nseq)
        Q_netsurf = Q_netsurf
     &            + a_surf(res_nn,res_seg,1)*Q_surface/rho_cp
!
! Vertical eddy diffusivity
!
!        eddy_dff = KAPPA_Z(nd,res_nn,T_epi,T_hyp)
!
        eddy_dff = 0.0
!
        Q_dff_epi = Q_dff_epi 
     &            + eddy_dff*(T_hyp-T_epi)*a_surf(res_nn,res_seg,1)
!
        Q_dff_hyp = Q_dff_hyp 
     &            + eddy_dff*(T_epi-T_hyp)*a_surf(res_nn,res_seg,2)
!     Look for a tributary.
! 
          q1      = cuft_cum*qin(nseq)
          q2      = cuft_cum*qin(nseq)
          q_tmp   = 0.0
!
          Q_advct_sum = 0.0  
          T_tmp       = 0.0  
          ntribs = no_tribs(nseq)
          do while (ntribs.gt.0.and..not.DONE)
            do ntrb=1,ntribs
              nr_trib = trib(nseq,ntrb)
              qq_trib = cuft_cum*q_trib(nr_trib)
              q2      = q2+qq_trib
              q_tmp   = q_tmp + qq_trib
              Q_advct_sum = Q_advct_sum 
     &                    + qq_trib*T_trib(nr_trib)
              T_tmp = T_tmp + qq_trib*T_trib(nr_trib)
            end do
            DONE=.TRUE.
          end do
!
! Add the advected flow and thermal energy from the tributaries in this cell
!
          T_tmp = T_tmp/q_tmp
          q_in_trb(:)    = 0.0
          Q_advct_trb(:) = 0.0
!
! Place reservoir inflow
!
!          layer = 1
!          call DENSE(DENSITY1,DENSITY2,T_tmp,T_epi)
!          if (DENSITY1 .gt. DENSITY2) layer = 2
            layer = NLAYER(T_tmp,T_epi)
!
! Place flow and thermal energy in composited layer
!
          q_in_trb(layer)    = q_tmp
          Q_advct_trb(layer) = Q_advct_sum
!
! Calculate residence time
!
          res_time = Vol_sum(1)/(86400.*q2)
     
!
! End of the loop that accumlates area, volume, surface transfer and advected sources
! from each reservoir cell
!      
      end do
!
! Reservoir outflows are based on the results from the last reservoir cell
! At present (Farmington River paper) Reservoirs Nepaug (1), Barkhamstead (2) and 
! Lake McDonough (3) are run-of-the-river and assumed to release water from the
! surface layer (epilimnion). Only Colebrook Lake and West Branch. which generate hydropower are 
! assumed to release water from the hypolimnion
!
      if (res_nn .le. 3) then
        q_out_res(1) = q2
        q_out_res(2) = 0.0
      else
        q_out_res(1) = 0.0
        q_out_res(2) = q2
      end if 
!
!	q_vert(:) = 0.0
!        q_vert(1) = q_out_res(1) - (q_in_res(1) + q_in_trb(1))
!        q_vert(2) = q_out_res(2) - (q_in_res(2) + q_in_trb(2))
!
! Special case for the Farmington River paper
!
         Q_vert(:) = 0.
         q_v = q_in_res(1) + q_in_trb(1) - q_out_res(1)
         if (q_v .gt.0.0001) then
           Q_vert(2) = q_v
         else
           Q_vert(1) = -q_v
         end if

!
        q_total =  q_out_res(1) + q_out_res(2)
!
        TTEMP_head = TEMP_head(nr,no_wr)

!
! Calculate residence time
!
        res_time = Vol_sum(1)/(86400.*q_total)
          q_vert_epi = q_in_res(1) - q_out
!        if (res_time .gt. 2.0000) then     
!
! Hypolimnion
!

!
          TEMP_res(res_nn,2,n2) = TEMP_res(res_nn,2,n1)
     &          + Q_dff_hyp
     &          + ((q_vert(2)*TEMP_res(res_nn,1,n1) 
     &          + Q_advect(2) + Q_advct_trb(2)
     &          - (q_vert(1) + q_out_res(2))*TEMP_res(res_nn,2,n1))
     &          / Vol_sum(2))* dt_comp 
         if (TEMP_res(1,2,n2) .lt. 0.0) TEMP_res(1,2,n2) = 0.0
!
! Epilimnion
!
          TEMP_res(res_nn,1,n2) = TEMP_res(res_nn,1,n1) 
     &         + ((Q_netsurf
     &         +  Q_dff_epi
     &         +  q_vert(1)*TEMP_res(res_nn,2,n1) 
     &         +  Q_advect(1) + Q_advct_trb(1)
     &         - (q_vert(2) + q_out_res(1))*TEMP_res(res_nn,1,n1))
     &        / Vol_sum(1))*dt_comp 
          if (TEMP_res(res_nn,1,n2) .lt. 0.5) TEMP_res(res_nn,1,n2) = 0.5
          if (TEMP_res(res_nn,2,n2) .gt. TEMP_res(res_nn,1,n2)) 
     &        TEMP_res(res_nn,2,n2) = TEMP_res(res_nn,1,n2)
!
!
!
          T_out = (q_out_res(1)*TEMP_res(res_nn,1,n2) 
     &          +  q_out_res(2)*TEMP_res(res_nn,2,n2)) / q_total
!    
! 
          nseg = no_celm(nr,no_wr)
!
        nnseg = 1
        temp(nr,no_wr,nnseg,n2) = T_out
        xkm_plot=x_dist(nr,no_wr,nseg)*3.048e-04
        write(20,'(f11.5,i5,1x,2i4,1x,5i5,1x,5f7.2,f9.1,2f9.1,3E12.3)') 
     &             time,nyear,nd,ndd,nr,no_wr,nseq,res_nn,nseg
     &            ,T_out,TEMP_head(nr,no_wr),dbt(nseq)
     &            ,TEMP_res(res_nn,1,n2),TEMP_res(res_nn,2,n2)
     &            ,qin(nseq),xkm_plot,res_time,Q_netsurf
     &            ,Q_dff_epi,Q_dff_hyp
!      ntmp=n1
!      n1=n2
!      n2=ntmp
!
      return
      end
