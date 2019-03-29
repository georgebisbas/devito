

float malloc2d(float *** C, int nrows, int ncols) {
    int i;
    *C = malloc( sizeof(float *) * nrows);

    if (*C == NULL) {
        printf("ERROR: out of memory\n");
        return 1;
    }

    for (i=0; i<nrows; i++) {
        (*C)[i] = malloc( sizeof(float) * ncols);
        if ((*C)[i] == NULL) {
            printf("ERROR: out of memory\n");
            return 1;
        }
    }
    //printf("Allocated!\n");
    return 0;
}

float malloc3d(float *** C, int timestamps, int nrows, int ncols) {
    int i,t_i;
    C = (float ***) malloc (sizeof(float **) * timestamps);

    for (t_i = 0; t_i < timestamps; t_i++) {
      C[t_i] = (float **) malloc( sizeof(float *) * nrows);

      if (*C == NULL) {
          printf("ERROR: out of memory\n");
          return 1;
        }

        for (i=0; i<nrows; i++) {
          C[t_i][i] = (float *)malloc( sizeof(float) * ncols);
          if (C[t_i][i] == NULL) {
              printf("ERROR: out of memory\n");
              return 1;
        }
    }
  }     printf("Allocated 3D!\n");

  printf("\n%1.1f ", C[0][5][9]);
  printf("\n%1.1f ", C[0][10][7]);
  printf("\n%1.1f ", C[0][20][5]);

  return 0;

}


void mem_allocate(int nrows, int ncols, float **C)
{
    C = malloc(nrows * sizeof(float * ));
    if (C == NULL) {
        fprintf(stderr, "out of memory\n");
    }
    for (int i = 0; i < nrows; i+=1) {
        C[i] = malloc(ncols * sizeof(float));
        if (C[i] == NULL) {
            fprintf(stderr, "out of memory\n");
        }
    }
    //return C;
}

void initialize2(int nrows, int ncols, float **A_init)
{
    // Initialize the matrix
    int ti = 0;
    int xi = 0;
    int yi = 0;
    for (xi = 0; xi < nrows; xi+=1) {
        for (yi = 0; yi < ncols; yi+=1) {
            A_init[xi][yi] = 10;
        }
  }
}

void initialize3(int nrows, int ncols, int timestamps, float ***A_init)
{
    // Initialize the matrix
    int ti = 0;
    int xi = 0;
    int yi = 0;
    for (xi = 0; xi < nrows; xi+=1) {
        for (yi = 0; yi < ncols; yi+=1) {
          for (ti = 0; ti < timestamps; ti+=1) {
            A_init[xi][yi][ti] = fabs(xi-yi);
        }
    }
  }
}

void printMatrix3d(int nrows, int ncols, int timestamps, float ***matrix) {

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      for (int k = 0; k < timestamps; k++) {
        printf("matrix[%d][%d][%d]: %lf\n", i, j, k, matrix[i][j][k]);
      }
    }
  }

  //return matrix;

}



void print_array_2d(int nrows, int ncols,int timestamps, float B[nrows][ncols][timestamps])
{
    int xi = 0;
    int yi = 0;
    int ti = 0;
    printf("\n");

    for (xi = 0; xi < nrows; xi++) {
        printf("\n");
        for (yi = 0; yi < ncols; yi++) {
            printf("%3.3f ", B[xi][yi][0]);
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



void validate3d(int nrows, int ncols, float ***u, float ***u2)
{
if (u[4:6][5:10] == u2[4:6][5:10]) {
   //printf("Validated \n");
 }
 else {printf("Not equal \n");}
}
