
serialCC=gcc
FLAGS=-O3 -lm -lrt
EXEC=serialGSPoisson.run

METHODS=serialGS

all:${METHODS}

serialGS: gs2D_serial.c
	${serialCC} ${FLAGS} -o ${EXEC} gs2D_serial.c

clean:
	rm -f *.run *.dat
