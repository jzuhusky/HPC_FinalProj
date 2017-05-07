# HPC_FinalProj

Multigrid Code for Solutions to 2D Laplace / Poisson Eq. on Unit Square:

Makefile to be pushed at some point

As of now:

If makefile avail:

# make 
# mpirun -n [# MPI TASKS] ./Executable ['N' Points Per Dim, solves (N+1)x(N+1) Mesh, N must be power of 2] (Opt # Grids)

Serial Run
# ./Executable (Npoints, power of 2) 

If no make file:
mpicc -O3 -g MPI_CFILE -lrt -lm -o EXECNAME
gcc -O3 SerialCFile -lrt -lm -o EXECNAME
