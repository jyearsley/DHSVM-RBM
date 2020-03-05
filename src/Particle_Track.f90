SUBROUTINE Particle_Track(no_wr,nr,ns,nx_s)
USE Block_Hydro
USE Block_Network
IMPLICIT NONE
integer, intent(IN):: no_wr,nr,ns
integer            :: ncell, nx_part, nx_s
real               :: dt_before,dt_total,dt_xcess
real               :: xpart_before
!
!     First do the reverse particle tracking
!
!     Segment is in cell SEGMENT_CELL(NC)
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
!          x_part(ns)=x_part(ns)+u(segment_cell(nr,no_wr,nx_part))            &
!                    *dt_part(ns)
          x_part(ns) = xpart_before
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
END SUBROUTINE Particle_Track
