!
!   Function to calculate the number of simulation days
!
      integer function nodays(jtime,jy0)
      integer,dimension(12),parameter :: ndmo = (/0,31,59,90,120,151,181,212,243,273,304,334/)
      integer                         :: jd,jm,jrem,jtime,jy0,jy,ny
!
      jy=(jtime/10000)
      jrem=(jtime-jy*10000)
      jm=jrem/100
      jd=jrem-jm*100
      ny=jy-jy0
      nodays=365*ny+ndmo(jm)+jd
!
      return
      end function nodays



