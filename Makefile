CC=gcc
OPTS=-g -Wall -O3 -march=native 

simulation: simulation.c
	$(CC) $(OPTS) simulation.c -o simulation -fopenmp

.PHONY: clean
clean:
	rm -rf *.o simulation