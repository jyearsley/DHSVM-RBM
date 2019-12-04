!
!
! Headwaters subroutine
!
      SUBROUTINE HEAD_TEMP (ncell,nr,no_wr)
!
      use Block_Hydro
      use Block_Network
      use Block_WQ
!
      implicit none
!
!     Headwaters flow and temperature
!
      if(IS_HEAD) then
        T_smth=b_smooth*T_smth
     &            +a_smooth*dbt(ncell)
        TEMP_head(nr,no_wr)=mu(nr)
     &            +(alf_Mu(nr)
     &            /(1.+exp(gmma(nr)*(beta(nr)-T_smth)))) 
      else
        nwr_1 = no_wr-1
        TEMP_head(nr,no_wr) = temp(nr,nwr_1,no_celm(nr,nwr_1),n1)
      end if
!        
      temp(nr,no_wr,0,n1)=TEMP_head(nr,no_wr)
      temp(nr,no_wr,-1,n1)=TEMP_head(nr,no_wr)
      temp(nr,no_wr,-2,n1)=TEMP_head(nr,no_wr)
      temp(nr,no_wr,no_celm(nr,no_wr)+1,n1)
     &     =temp(nr,no_wr,no_celm(nr,no_wr),n1)
      return
      end

