all: tiling

CC = icc
DEPS = tile_auxiliary.h


tiling: tiling.c
	$(CC) -O3 -g -ffast-math -w -fopenmp tiling.c -lm
clean:
	rm -f tiling *.o
