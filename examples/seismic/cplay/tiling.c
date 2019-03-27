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

int main(int argc, char **argv) {
    int nrows = atoi(argv[1]);
    int ncols = atoi(argv[2]);
    int timesteps = atoi(argv[3]);

    int nsrc = 1;
    int nrecs = 10;

    struct timeval t1, t2;
    double elapsedTime;
    printf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\n",nrows, ncols, timesteps);

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

    /* Flush denormal numbers to zero in hardware */
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);


    // Jacobi  #pragma omp parallel for
    initialize(nrows, ncols, A);
    initialize(nrows, ncols, B);

    printf("Starting Jacobi..."); gettimeofday(&t1, NULL);
    jacobi_omp_par_src_rcv(timesteps, nrows, ncols, nsrc, nrecs, A, B, src_coords, rec_coords);
    gettimeofday(&t2, NULL);    printf("... Finished \n");
    elapsedTime = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
    printf("Jacobi OpenMP par-for, Time taken by program is : %f\n",elapsedTime);

    //print_array(nrows, ncols, A);
    //print_array(nrows, ncols, B);

    validate(nrows, ncols, A, B);

    printf("\n");
    return 0;
}
