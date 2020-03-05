      integer function ndate(n,ny0)
      integer,dimension(13),parameter :: ndmo = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer                         :: n,nday,njul,nm,nmon,ny0,ny,nyear
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
!
      end function ndate

