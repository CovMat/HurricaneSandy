CC = gcc
CFLAGS = -Wall -O1
BIN = .
all:cal_omp
z_sacio.o:z_sacio.c
	$(CC) -c $^ ${CFLAGS}
distaz.o:distaz.c
	$(CC) -c $^ ${CFLAGS}
cal_omp.o:cal_omp.c
	$(CC) -c -fopenmp $^ ${CFLAGS}
cal_omp:z_sacio.o distaz.o cal_omp.o
	$(CC) -o $(BIN)/$@ -fopenmp $^ ${CFLAGS} -L. -lm -lsac -lsacio 
clean:
	rm -f *.o
