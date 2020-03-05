      PROGRAM RBM10_VIC

!      PROGRAM RMAIN
!
!     Dynamic river basin model for simulating water quality in
!     branching river systems with freely-flowing river sfegments. 
!
!     This version uses Reverse Particle Tracking in the Lagrangian
!     mode and Lagrangian interpolation in the Eulerian mode.
!
!     Topology and routing is set up to be consistent with output
!     from the Distributed Hydrologic Soil and Vegetation Model (DHSVM)
!     model developed by the Land Surface Hydrology Group at 
!     the University of Washington.
!
!     For additional information visit:
!
!     http://github.com/jyearsley/DHSVM-RBM/tree/reservoir_mod/
!
!     or contact:
!
!     John Yearsley
!     UW Hydro/Computational Hydrology
!     Dept. of Civil and Environmental Engineering
!     Box 352700
!     University of Washington
!     Seattle, Washington
!     98195-2700
!     jyearsley@uw.edu
!
  implicit none
      character (len=200)    :: param_file,source_file,spatial_file
      character (len=200)    :: Prefix1, Prefix2
      integer                :: iargc
      integer                :: numarg
 
!     Command line input
!
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
!
!     Identify and open necessary files
!
!     open the output file 
      open(unit=20,file=TRIM(Prefix2)//'.temp',status='unknown')
!
!     Open file with weather and inflow data
      write(*,*) 'Forcing file -  ', TRIM(Prefix1)//'.forcing'
      open(unit=30,file=TRIM(Prefix1)//'.forcing',STATUS='old')
!
!     open Mohseni file 
      open(40,file=TRIM(Prefix1)//'.Mohseni',STATUS='old')    
!
!     open Leopold file 
      open(50,file=TRIM(Prefix1)//'.Leopold',STATUS='old')    
!
!     Open network file
      param_file = TRIM(Prefix1)//'.net'
      write(*,*) 'Network file ',param_file
      OPEN(UNIT=90,FILE=param_file,STATUS='OLD')

!
!     Call systems programs to get started
!
!     SUBROUTINE BEGIN reads control file, sets up topology and
!     important properties of reaches
!
      CALL BEGIN (param_file)
      write(*,*) 'Calling SYSTMM'
!
!     SUBROUTINE SYSTMM performs the simulations
!
      CALL SYSTMM
!
!     Close files after simulation is complete
!
      write(*,*) ' Closing files after simulation'
      CLOSE(30)
      CLOSE(90)
      STOP
!
      END PROGRAM RBM10_VIC
