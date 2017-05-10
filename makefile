
CC=gcc 
MPICC=mpicc
flags=-O3 -lrt -lm

METHODS=serial2DMultigrid MPI_Multigrid gs2D_serial

all:${METHODS}

serial2DMultigrid: serial2DMultigrid.c
	${CC} ${flags} $^ -o serial.run

MPI_Multigrid: MPI_Multigrid.c
	${MPICC} ${flags} $^ -o parallel.run

gs2D_serial: gs2D_serial.c
	${CC} ${flags} $^ -o gsSerial.run

clean:
	rm -f *.run *.dat *.out
