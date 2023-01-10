CC=gcc
MFLAGS=-lm
GFLAGS=-lGL -lGLU  -lglut

#use the 1st line in alternative to the 2nd (1st is for valgrind, 2nd is for benchmarking)
CFLAGS=-Wall -O0 -g
#CFLAGS=-Wall -O2

all: mandelbrot-gui-serial.exe

mandelbrot-gui-serial.exe: mandelbrot-gui-serial.c
	$(CC) $(CFLAGS) mandelbrot-gui-serial.c $(MFLAGS) $(GFLAGS) -o mandelbrot-gui-serial.exe 

clean:
	rm -f *.exe *.o a.out
