include ../Make.inc

##############################################################################

all: distclean lib mpilib

clean:
	@echo "Cleaning build files..."
	@rm -rf *.o;

distclean: clean
	@echo "Cleaning libraries..."
	@rm -rf libpso.a;
	@rm -rf libpso_mpi.a;

##############################################################################

Swarm.o: Swarm.cpp Swarm.h
	$(CPP) $(CPPFLAGS) -c Swarm.cpp -o Swarm.o

Particle.o: Particle.cpp Particle.h
	$(CPP) $(CPPFLAGS) -c Particle.cpp -o Particle.o

Point.o: Point.cpp Point.h
	$(CPP) $(CPPFLAGS) -c Point.cpp -o Point.o

##############################################################################

Swarm_mpi.o: Swarm.cpp Swarm.h
	$(MPICPP) $(MPICPPFLAGS) -DUSE_MPI -c Swarm.cpp -o Swarm_mpi.o

Particle_mpi.o: Particle.cpp Particle.h
	$(MPICPP) $(MPICPPFLAGS) -DUSE_MPI -c Particle.cpp -o Particle_mpi.o

Point_mpi.o: Point.cpp Point.h
	$(MPICPP) $(MPICPPFLAGS) -DUSE_MPI -c Point.cpp -o Point_mpi.o

##############################################################################

lib: Swarm.o Particle.o Point.o
	ar crs libpso.a Swarm.o Particle.o Point.o

mpilib: Swarm_mpi.o Particle_mpi.o Point_mpi.o
	ar crs libpso_mpi.a Swarm_mpi.o Particle_mpi.o Point_mpi.o

##############################################################################

test: test.o $(OBJ)
	$(CPP) -o test test.o $(OBJ) $(LDFLAGS)

