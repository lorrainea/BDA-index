SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -pthread -std=c++0x -DNDEBUG -O3 -march=native

all: construct_sa

construct_sa:
	$(CC) $(CFLAGS) -I ./libsdsl/include/ -o construct_sa ./src/main.cpp ./src/utils.cpp -L ./libsdsl/lib/ -lsdsl -ldivsufsort -ldivsufsort64 -fopenmp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f construct_sa *.o
