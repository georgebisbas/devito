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

double*** createMatrix(int nrows, int ncols, int timestamps) {

  double ***matrix;

  matrix = (double***)malloc(nrows*sizeof(double**));
  for (int row = 0; row < nrows; row++) {
    matrix[row] = (double**)malloc(ncols*sizeof(double*));
    for (int col = 0; col < ncols; col++) {
      matrix[row][col] = (double*)malloc(timestamps*sizeof(double));
    }
  }

  //printf("Insert the elements of your matrix:\n");
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      for (int k = 0; k < timestamps; k++) {
        //printf("Insert element [d][d][%d]: ", i, j, k);
        matrix[i][j][k] = 0;
  //      printf("matrix[%d][%d][%d]: %lf\n", i, j, k, matrix[i][j][k]);
      }
    }
  }
  printf("allocate3d...");

  return matrix;
}


int main(int argc, char **argv) {
    int nrows = atoi(argv[1]);
    int ncols = atoi(argv[2]);
    int timestamps = 2;
    int timesteps = atoi(argv[3]);

    int nsrc = 1;
    int nrecs = 10;

    struct timeval t1, t2;
    double elapsedTime;
    printf("Problem setup is \nnrows: %d\n ncols: %d\nTimesteps: %d\nTimestamps: %d\n",nrows, ncols, timesteps, timestamps);

    //float ** A = mem_allocate(nrows, ncols, A);

    float ** A;
    malloc2d(&A, nrows, ncols);
    float ** B; // B is like A(t-1)
    malloc2d(&B, nrows, ncols);
    float ** src_coords; // B is like A(t-1)
    malloc2d(&src_coords, nsrc, 2);
    float ** rec_coords; // B is like A(t-1)
    malloc2d(&rec_coords, nsrc, 2);
    float ** src; // B is like A(t-1)
    malloc2d(&src, nsrc, 2);
    float *** u;
    u = createMatrix(nrows, ncols, timestamps);


    /* Flush denormal numbers to zero in hardware */
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

    // Jacobi  #pragma omp parallel for
    //initialize2(nrows, ncols, A);
    //initialize2(nrows, ncols, B);

    //for (int i = 0; i < timestamps; i++) {
    initialize3(nrows, ncols,timestamps, u);
    //}
    //for (int k = 0; k < timestamps; k++) {
    if(0){

      printf("\n ------------------------\n");
        for (int i = 0; i < nrows; i++) {
          printf("\n");
          for (int j = 0; j < ncols; j++) {
          //printf(" u[%d][%d][%d]: %1.1f", i, j, k, u[i][j][k]);
          printf(" %1.1f", u[i][j][1]);

        }
      }
    //}
}
    printf("Starting Jacobi..."); gettimeofday(&t1, NULL);
    u = jacobi_omp_par_src_rcv(timesteps, nrows, ncols, timestamps, nsrc, nrecs, u, src_coords, rec_coords);
    gettimeofday(&t2, NULL);    printf("... Finished \n");
    elapsedTime = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
    printf("Jacobi OpenMP par-for, Time taken by program is : %3.3f\n",elapsedTime);

    //for (int i = 0; i < timestamps; i++) {
      //print_array_2d(nrows, ncols, timestamps, u[][][i]);
     //}

    // printMatrix3d(nrows,ncols,timestamps,u);
    if(0){
    for (int k = 0; k < timestamps; k++) {
      printf("\n ------------------------\n");
        for (int i = 0; i < nrows; i++) {
          printf("\n");
          for (int j = 0; j < ncols; j++) {
          //printf(" u[%d][%d][%d]: %1.1f", i, j, k, u[i][j][k]);
          printf(" %1.1f", u[i][j][k]);

        }
      }
    }
    }
    //print_array_3d(nrows, ncols, timestamps, u[0]);

    //print_array_3d(nrows, ncols, u );

    //validate3d(nrows, ncols, u);
    printf(" Success ");
		free(u);

    printf("\n");
    return 0;
}
