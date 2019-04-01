#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include <time.h>
#include <stdio.h>
#include <cilk/cilk.h>
#include "xmmintrin.h"
#include "pmmintrin.h"
#include "tile_auxiliary.h"

int main(int argc, const char * argv[]) {
  int timesteps = 40000;
  int nrows = 1000;
  int ncols = 1000;
  int i = 0;
  struct timeval t1, t2;
  double elapsedTime;

  //double ** A = mem_allocate(nrows, ncols, A);

  double ** A;
  malloc2d(&A, nrows, ncols);
  double ** B;
  malloc2d(&B, nrows, ncols);

  /* Flush denormal numbers to zero in hardware */
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);


// Jacobi  #pragma omp parallel for
initialize(nrows, ncols, A);
initialize(nrows, ncols, B);

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
