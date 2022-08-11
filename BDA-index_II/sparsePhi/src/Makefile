SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -std=c++0x -O3 -DNDEBUG -pthread
AUX_PAR_FLAGS = -fopenmp
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++0x -g2 -O3 -pthread

all: construct_lcp_sequential construct_lcp_parallel

construct_lcp_sequential:
	$(CC) $(CFLAGS) -o construct_lcp_sequential em_sparse_phi_src/utils.cpp main.cpp

construct_lcp_parallel:
	$(CC) $(CFLAGS) $(AUX_PAR_FLAGS) -o construct_lcp_parallel em_sparse_phi_src/utils.cpp main.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f construct_lcp_sequential construct_lcp_parallel *.o
