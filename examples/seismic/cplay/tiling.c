#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include <time.h>
#include <stdio.h>
#include <cilk/cilk.h>

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

void jacobi_omp_par(int timesteps, int nrows, int ncols, float **A, float **B)
{
  int t, xi, yi;
  for (t = 0; t < timesteps; t++) {
    #pragma omp parallel
    #pragma omp for
    for (xi = 1; xi < nrows-1; xi++) {
      for (yi = 1; yi < ncols-1; yi++) {
        B[xi][yi] = (A[xi][yi] + A[xi-1][yi] + A[xi+1][yi] + A[xi][yi-1] + A[xi][yi+1])/6;
      }
    }
    #pragma omp for
    for (xi = 1; xi < nrows-1; xi++) {
      for (yi = 1; yi < ncols-1; yi++) {
        A[xi][yi] = B[xi][yi];
      }
    }
  }
}

void jacobi_omp_parsimd(int timesteps, int nrows, int ncols, float **A, float **B)
{
  int t, xi, yi;
  for (t = 0; t < timesteps; t++) {
    #pragma omp parallel for
    //#pragma omp simd
    for (xi = 1; xi < nrows-1; xi++) {
      #pragma omp simd
      for (yi = 1; yi < ncols-1; yi++) {
        B[xi][yi] = (A[xi][yi] + A[xi-1][yi] + A[xi+1][yi] + A[xi][yi-1] + A[xi][yi+1])/6;
      }
    }
    //#pragma omp simd
    for (xi = 1; xi < nrows-1; xi++) {
      #pragma omp simd
      for(yi = 1; yi < ncols-1; yi++) {
        A[xi][yi] = B[xi][yi];
      }
    }
  }
}

void jacobi_omp_simd(int timesteps, int nrows, int ncols, float **A, float **B)
{
  int t, xi, yi;
  for (t = 0; t < timesteps; t+=1) {
    //#pragma omp simd
    for (xi = 1; xi < nrows-1; xi+=1) {
      #pragma omp simd
      for (yi = 1; yi < ncols-1; yi+=1) {
        B[xi][yi] = (A[xi][yi] + A[xi-1][yi] + A[xi+1][yi] + A[xi][yi-1] + A[xi][yi+1])/6;
      }
    }
    for (xi = 1; xi < nrows-1; xi+=1) {
      #pragma omp simd
      for (yi = 1; yi < ncols-1; yi+=1) {
        A[xi][yi] = B[xi][yi];
      }
    }
  }
}

void jacobi_cilk(int timesteps, int nrows, int ncols, float **A, float **B)
{
  int t, xi, yi;
  for (t = 0; t < timesteps; t+=1) {
    //#pragma omp simd
    cilk_for (xi = 1; xi < nrows-1; xi+=1) {
      for (yi = 1; yi < ncols-1; yi+=1) {
        B[xi][yi] = (A[xi][yi] + A[xi-1][yi] + A[xi+1][yi] + A[xi][yi-1] + A[xi][yi+1])/6;
      }
    }
    cilk_for (xi = 1; xi < nrows-1; xi+=1) {
      for (yi = 1; yi < ncols-1; yi+=1) {
        A[xi][yi] = B[xi][yi];
      }
    }
  }
}






void initialize(int nrows, int ncols, float **A_init)
{
  // Initialize the matrix
    int t = 0;
    int xi = 0;
    int yi = 0;
    for (xi = 0; xi < nrows; xi+=1) {
      for (yi = 0; yi < ncols; yi+=1) {
        A_init[xi][yi] = 1;
      }
    }
}


  int main(int argc, const char * argv[]) {
    int timesteps = 40;
    int nrows = 20000;
    int ncols = 20000;
    int i = 0;

    struct timeval t1, t2;
    double elapsedTime;

    float ** A;
    A = malloc(nrows * sizeof(float * ));
    if (A == NULL) {
      fprintf(stderr, "out of memory\n");
      return 0;
    }
    for (i = 0; i < nrows; i+=1) {
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
    for (i = 0; i < nrows; i+=1) {
      B[i] = malloc(ncols * sizeof(float));
      if (B[i] == NULL) {
        fprintf(stderr, "out of memory\n");
        return 0;
      }
    }


// Jacobi  #pragma omp parallel for
initialize(nrows, ncols, A); initialize(nrows, ncols, B);
gettimeofday(&t1, NULL);
jacobi_omp_par(timesteps, nrows, ncols, A, B);
gettimeofday(&t2, NULL);
elapsedTime = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
printf("Jacobi OpenMP par-for, Time taken by program is : %f\n",elapsedTime);


if (nrows < 15) {
  print_array(nrows, ncols, B);
}

// Jacobi  #pragma omp simd
initialize(nrows, ncols, A); initialize(nrows, ncols, B);
gettimeofday(&t1, NULL);
jacobi_omp_simd(timesteps, nrows, ncols, A, B);
gettimeofday(&t2, NULL);
elapsedTime = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
printf("Jacobi OpenMP SIMD, Time taken by program is : %f\n",elapsedTime);

if (nrows < 15) {
  print_array(nrows, ncols, B);
}
// Jacobi  #pragma omp parfor + simd
initialize(nrows, ncols, A); initialize(nrows, ncols, B);
gettimeofday(&t1, NULL);
jacobi_omp_parsimd(timesteps, nrows, ncols, A, B);
gettimeofday(&t2, NULL);
elapsedTime = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
printf("Jacobi OpenMP parfor SIMD, Time taken by program is : %f\n",elapsedTime);

if (nrows < 15) {
  print_array(nrows, ncols, B);
}


// Jacobi CILK
initialize(nrows, ncols, A); initialize(nrows, ncols, B);
gettimeofday(&t1, NULL);
jacobi_omp_parsimd(timesteps, nrows, ncols, A, B);
gettimeofday(&t2, NULL);
elapsedTime = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
printf("Jacobi CILK, Time taken by program is : %f\n",elapsedTime);

if (nrows < 15) {
  print_array(nrows, ncols, B);
}





















    printf("\n");
    return 0;
  }
