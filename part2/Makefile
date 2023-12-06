# Makefile
CC=mpicxx
CFLAGS=-Wall
TARGET=matmul
OBJS=main.o matrix_operations.o
HOSTFILE=hosts_matrix.txt
DATAFILE=/home/.grade/HW4/data-set/data2_7
# DATAFILE=data0_1

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

main.o: main.cc src/matrix_operations.h
	$(CC) $(CFLAGS) -c main.cc

matrix_operations.o: src/matrix_operations.cpp src/matrix_operations.h
	$(CC) $(CFLAGS) -c src/matrix_operations.cpp

clean:
	rm -f $(OBJS) $(TARGET)

run: $(TARGET)
	parallel-scp -A -h $(HOSTFILE) $(TARGET) ~
	mpirun -np 4 --hostfile $(HOSTFILE) ./$(TARGET) < $(DATAFILE) > output.txt
