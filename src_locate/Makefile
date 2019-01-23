CC = gcc
CFLAGS = -Wall -O1
BIN = .
all:project
distaz.o:distaz.c
	$(CC) -c $^ ${CFLAGS}
project.o:project.c
	$(CC) -c -fopenmp $^ ${CFLAGS}
project:project.o distaz.o
	$(CC) -o $(BIN)/$@ $^ ${CFLAGS} -L. -lm -fopenmp 
clean:
	rm -f *.o
