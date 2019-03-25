all: tiling

tiling: tiling.c
	icc -O3 -g -Wall -fopenmp tiling.c -lm
clean:
	rm -f tiling *.o
