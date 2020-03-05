      SUBROUTINE DENSE(DENSITY1,DENSITY2,T11,T22)
      integer, parameter          :: dp = selected_real_kind(15, 307)
      real                        :: t11,t22
      REAL(kind=dp)               :: DENSITY1,DENSITY2,T1,T2
         T1 = T11
         T2 = T22
!
          DENSITY1 = 6.793952E-02*T1 -9.095290E-03*(T1**2)       &
                  + 1.001685E-04*(T1**3) - 1.120083E-06*(T1**4)  &
                  + 6.536332E-09*(T1**5)
          DENSITY2 = 6.793952E-02*T2 -9.095290E-03*(T2**2)       &
                  + 1.001685E-04*(T2**3) - 1.120083E-06*(T2**4)  &
                  + 6.536332E-09*(T2**5)
!
      END SUBROUTINE DENSE 
