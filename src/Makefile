#######################################################################
###  RBM_DHSVM.f makefile #############################################
#######################################################################
#
# Semi-Lagrangian stream temperature model with reservoir module
#
# The compiler
FC = gfortran

FFFLAGS = -03 -C

FCFLAGS = -03

HFILES = RBM.fi

OBJECTS = RBM.o

exe:   $(OBJECTS)
  $(FC) $(FFLAGS) $(OBJECTS) -o RBM
  
  RBM.o: RBM.f
  
  clean:
    /bin/rm *.o
    
    %.o: %.f
      $(FC) $(FCFLAGS) -c $<
