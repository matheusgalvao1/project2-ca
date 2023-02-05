CC=gcc
MFLAGS=-lm
GFLAGS=-lGL -lGLU  -lglut

CC_MPI=mpicc
PTHREAD=-pthread

#use the 1st line in alternative to the 2nd (1st is for valgrind, 2nd is for benchmarking)
#CFLAGS=-Wall -O0 -g
CFLAGS=-Wall -O2

all: mandelbrot-gui-serial.exe

mandelbrot-gui-serial.exe: mandelbrot-gui-serial.c
	$(CC) $(CFLAGS) mandelbrot-gui-serial.c $(MFLAGS) $(GFLAGS) -o mandelbrot-gui-serial.exe 

#uncomment one (only if you are running in the cluster)
#MPIFLAGS=-I/share/apps/openmpi-4.1.0/include
#MPIFLAGS=-I/share/apps/mpich-3.3.2/include

mpi: mandelbrot-gui-mpi.c
	$(CC_MPI) $(CFLAGS) mandelbrot-gui-mpi.c $(MFLAGS) $(GFLAGS) $(PTHREAD) -o mandelbrot-gui-mpi.exe

clean:
	rm -f *.exe *.o a.out

