#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include <stdio.h>


void print_array(int nrows, int ncols, float **B)
{
    int xi,yi;
    printf("\n");
    for (xi = 0; xi < nrows; xi++) {
      printf("\n");
      for (yi = 0; yi < ncols; yi++) {
        printf("%1.1f ", B[xi][yi]);
      }
    }
}

void jacobi_omp(int timesteps, int nrows, int ncols, float **A, float **B)
{
  int t, xi, yi;
  for (t = 0; t < timesteps; t++) {
    //#pragma omp parallel for
    for (xi = 1; xi < nrows-1; xi++) {
      for (yi = 1; yi < ncols-1; yi++) {
        B[xi][yi] = (A[xi][yi] + A[xi-1][yi] + A[xi+1][yi] + A[xi][yi-1] + A[xi][yi+1])/6;
      }
    }
    for (xi = 1; xi < nrows-1; xi++) {
      //#pragma omp simd
      for (yi = 1; yi < ncols-1; yi++) {
        A[xi][yi] = B[xi][yi];
      }
    }
  }
}




  int main(int argc, const char * argv[]) {
    int timesteps = 90000;
    int nrows = 1000;
    int ncols = 1000;
    int i = 0;

    struct timeval t1, t2;
    double elapsedTime;

        float ** A;
    A = malloc(nrows * sizeof(float * ));
    if (A == NULL) {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
    for (i = 0; i < nrows; i++) {
      A[i] = malloc(ncols * sizeof(float));
      if (A[i] == NULL) {
        fprintf(stderr, "out of memory\n");
        return 0;
      }
    }

    float ** B;
    B = malloc(nrows * sizeof(float * ));
    if (B == NULL) {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
    for (i = 0; i < nrows; i++) {
      B[i] = malloc(ncols * sizeof(float));
      if (B[i] == NULL) {
        fprintf(stderr, "out of memory\n");
        return 0;
      }
    }

// Zero the matrix
    int t = 0;
    int xi = 0;
    int yi = 0;
    for (xi = 0; xi < nrows; xi++) {
      for (yi = 0; yi < ncols; yi++) {
        A[xi][yi] = 1;
        B[xi][yi] = 1;
      }
    }
// Print the matrix
//print_array(nrows, ncols, A);

gettimeofday(&t1, NULL);
jacobi_omp(timesteps, nrows, ncols, A, B);
gettimeofday(&t2, NULL);
elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
printf("%f ms.\n ", elapsedTime);

//print_array(nrows, ncols, B);


    printf("\n");
    return 0;
  }
