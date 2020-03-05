      SUBROUTINE RK4(Y,DYDT,N,X,H,YOUT)
      INTEGER,PARAMETER             :: NMAX=10
      REAL,DIMENSION(:),ALLOCATABLE :: Y(N),DYDT(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      ALLOCATE (Y(N),DYDT(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX))
!                     
      HH=H*0.5
      H6=H/6.
      XH=X+HH
!
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDT(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDT(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END SUBFROUTINE RK4
