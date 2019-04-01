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
#include "tile_cilk.h"
#include "tile_blk.h"



double*** createMatrix(int nrows, int ncols, int timestamps) {

  double ***matrix;

  matrix = (double***)malloc(nrows*sizeof(double**));
  for (int row = 0; row < nrows; row++) {
    matrix[row] = (double**)malloc(ncols*sizeof(double*));
    for (int col = 0; col < ncols; col++) {
      matrix[row][col] = (double*)malloc(timestamps*sizeof(double));
    }
  }
  printf("allocate3d...");
  return matrix;
}



int main(int argc, char **argv) {
    int nrows = atoi(argv[1]);
    int ncols = atoi(argv[2]);
    int timesteps = atoi(argv[3]);

    int nsrc = 1;
    int nrecs = 10;
    int timestamps = 2;
    int omp_opt=0;



    struct timeval t1, t2;
    double elapsedTime;
    printf("Problem setup is \nnrows: %d\n ncols: %d\nTimesteps: %d\nTimestamps: %d\n",nrows, ncols, timesteps, timestamps);

    double ** A;
    malloc2d(&A, nrows, ncols);
    double ** B; // B is like A(t-1)
    malloc2d(&B, nrows, ncols);
    double ** src_coords; // B is like A(t-1)
    malloc2d(&src_coords, nsrc, 2);
    double ** rec_coords; // B is like A(t-1)
    malloc2d(&rec_coords, nsrc, 2);
    double ** src; // B is like A(t-1)
    malloc2d(&src, nsrc, 2);
    double *** u;
    u = createMatrix(nrows, ncols, timestamps);
    double *** u2;
    u2 = createMatrix(nrows, ncols, timestamps);


    /* Flush denormal numbers to zero in hardware */
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

    initialize3(nrows, ncols,timestamps, u);
    initialize3(nrows, ncols,timestamps, u2);


    printf("\n Starting Jacobi..."); gettimeofday(&t1, NULL);

    u = jacobi_3d_all(timesteps, nrows, ncols, u, omp_opt=0);

    gettimeofday(&t2, NULL);
    printf("... Finished \n");
    elapsedTime = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
    printf("Jacobi OpenMP par-for, Time taken by program is : %3.3f\n",elapsedTime);

    printf("Starting Jacobi..."); gettimeofday(&t1, NULL);
    u2 = jacobi_3d_all_SKEW(timesteps, nrows, ncols, u2, omp_opt=0);
        //u2 = jacobi_omp_par_src_rcv(timesteps, nrows, ncols, timestamps, nsrc, nrecs, u2, src_coords, rec_coords);
    gettimeofday(&t2, NULL);    printf("... Finished \n");
    elapsedTime = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
    printf("Jacobi NAIVE, Time taken by program is : %3.3f\n",elapsedTime);



int validate_flag=1;
if(validate_flag){
    for (int k = 0; k < timestamps; k++) {
        for (int i = 0; i < nrows; i++) {
          for (int j = 0; j < ncols; j++) {
          if ((u[i][j][k]-u2[i][j][k])> 0.01) {
            printf(" Failed \n");
          }
        }
      }
    }
}


int si = 40; int ei = 50; int sc = 40; int ec = 50;
if(1){
  printf("\n ------------------------\n");
    for (int i = si; i < ei; i++) {
      printf("\n");
      for (int j = sc; j < ec; j++) {
      printf(" %3.3f", u[i][j][1]);
    }
  }
}

if(1){
  printf("\n ------------------------\n");
    for (int i = si; i < ei; i++) {
      printf("\n");
      for (int j = sc; j < ec; j++) {
      printf(" %3.3f", u2[i][j][1]);
    }
  }
}


		free(u);
    free(u2);
    printf("\n");
    return 0;
}
