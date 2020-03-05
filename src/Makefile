#
# Makefile for the grid-based semi-Lagrangian water temperature model, RBM_X
#
# Start of the makefile
#
# Defining variables
#
objects = RBM_X.o Begin.o Energy.o \
    Particle_Track.o River.o Reservoir.o Systmm.o\
    Block_Energy.o Block_Hydro.o Block_Network.o Block_WQ.o\
    density.o head_temp.o julian.o kappa.o leap_year.o \
    ndate.o nlayer.o nodays.o tntrp.o
#
# Use the open source fortran compiler
#	  
f90comp = gfortran
# Makefile
RBM_X: $(objects)
	$(f90comp) -fbounds-check -o RBM_X $(objects) 
#
Block_WQ.o: Block_WQ.f90
	$(f90comp) -c Block_WQ.f90 -Wall
block_wq.mod: Block_WQ.o Block_WQ.f90
	$(f90comp) -c Block_WQ.f90 -Wall
#
Block_Energy.o: Block_Energy.f90
	$(f90comp) -c Block_Energy.f90 -Wall
block_energy.mod: Block_Energy.o Block_Energy.f90
	$(f90comp) -c Block_Energy.f90 -Wall
#
Block_Hydro.o: Block_Hydro.f90
	$(f90comp) -c Block_Hydro.f90 -Wall
block_hydro.mod: Block_Hydro.o Block_Hydro.f90
	$(f90comp) -c Block_Hydro.f90 -Wall
#
Block_Network.o: Block_Network.f90
	$(f90comp) -c Block_Network.f90 -Wall
block_network.mod: Block_Network.o Block_Network.f90
	$(f90comp) -c Block_Network.f90 -Wall
#
#Alloc.o: block_energy.mod block_network.mod block_hydro.mod block_wq.mod Alloc.f90
#	$(f90comp) -c Alloc.f90 -Wall
#
Begin.o: block_energy.mod block_network.mod block_hydro.mod block_wq.mod Begin.f90
	$(f90comp) -c Begin.f90 -Wall
#
Read_Forcing.o: block_energy.mod block_hydro.mod block_network.mod Read_Forcing.f90
	$(f90comp) -c Read_Forcing.f90 -Wall
#
Reservoir.o: block_energy.mod block_hydro.mod block_network.mod Reservoir.f90
	$(f90comp) -c Reservoir.f90 -Wall
#
River.o: block_energy.mod block_hydro.mod block_network.mod River.f90
	$(f90comp) -c River.f90 -Wall
#
Systmm.o: block_network.mod block_energy.mod block_hydro.mod Systmm.f90
	$(f90comp) -c Systmm.f90 -Wall
#
Energy.o: block_energy.mod block_network.mod Energy.f90
	$(f90comp) -c Energy.f90 -Wall
#
Particle_Track.o: block_hydro.mod block_network.mod Particle_Track.f90
	$(f90comp) -c Particle_Track.f90 -Wall
#
Write.o: Write.f90
	$(f90comp) -c Write.f90 -Wall
#
julian.o: julian.f90
	$(f90comp) -c julian.f90 -Wall
density.o: density.f90
	$(f90comp) -c density.f90 -Wall
head_temp.o: head_temp.f90
	$(f90comp) -c head_temp.f90 -Wall
kappa.o: kappa.f90
	$(f90comp) -c kappa.f90 -Wall
leap_year.o: leap_year.f90
	$(f90comp) -c leap_year.f90 -Wall
#
ndate.o: ndate.f90
	$(f90comp) -c ndate.f90 -Wall
nlayer.o: nlayer.f90
	$(f90comp) -c nlayer.f90 -Wall
nodays.o: nodays.f90
	$(f90comp) -c nodays.f90 -Wall
tntrp.o: tntrp.f90
	$(f90comp) -c tntrp.f90 -Wall
#
RBM_X.o: RBM_X.f90
	$(f90comp) -c RBM_X.f90 -Wall

# Cleaning everything
clean:
	rm block_energy.mod block_hydro.mod block_network.mod block_wq.mod\
           RBM_X
	rm $(objects)
