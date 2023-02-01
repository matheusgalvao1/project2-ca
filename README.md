# project2-ca

Project for the course of Advanced Computing from Instituto Politecnico de Braganca, this project applies parallel computing using MPI to a serial code for the generation of Mandelbrot set. 

1- Create the executable files
- To makefile has two different commands
- "make" will generate the .exe file for serial version
- "make mpi" will generate the .exe file for the MPI version

2- Execute the program
- Serial version command: ./mandelbrot-gui-serial.exe
- MPI version command: mpirun -np 2 --hostfile localhost.OPENMPI ./mandelbrot-gui-mpi.exe
- For the MPI version you can change the number of tasks used by changing the number that comes after "-np"
