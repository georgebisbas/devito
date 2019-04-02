#define min(a, b) (((a) < (b)) ? (a) : (b))

double ** * jacobi_omp_par_src_rcv(int timesteps, int nrows, int ncols, int timestamps, int nsrc, int nrcv, double ** * u, double ** src_coords, double ** rec_coords) {
  int t, xi, yi;
  for (int titer = 0; titer < timesteps - 1; titer++) {
    //for (int time = 0, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3); time <= timestamps; time += 1, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3))
    //{
    t = 0;

    u[1: nrows - 2][1: ncols - 2][t + 1] = (u[1: nrows - 2][1: ncols - 2][t] + u[0: nrows - 3][1: ncols - 2][t] + u[2: nrows - 1][1: ncols - 2][t] + u[1: nrows - 2][0: ncols - 3][t] + u[1: nrows - 2][2: ncols - 1][t]) / 6;
    u[0: nrows - 1][0: ncols - 1][t] = u[0: nrows - 1][0: ncols - 1][t + 1]; //Cilk array notation

    int x_m = 1;
    int x_M = nrows;
    int y_m = 1;
    int y_M = ncols;
    cilk_for(int p_src = 0; p_src < nsrc; p_src++) {
      double r1 = (int)(floor(5.0e-1F * src_coords[p_src][0]));
      int ii_src_0 = r1 + 10;
      double r2 = (int)(floor(5.0e-1F * src_coords[p_src][1]));
      int ii_src_1 = r2 + 10;
      int ii_src_2 = r2 + 11;
      int ii_src_3 = r1 + 11;
      double px = (double)(-2.0F * r1 + src_coords[p_src][0]);
      double py = (double)(-2.0F * r2 + src_coords[p_src][1]);
      if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1) {
        int r3 = ii_src_0 + 2;
        int r4 = ii_src_1 + 2;
        int r5 = ii_src_0 + 1;
        int r6 = ii_src_1 + 1;
        double r7 = 5.7121e-2F * (2.5e-1F * px * py - 5.0e-1F * px - 5.0e-1F * py + 1);
        //*src[time][p_src]/m[r5][r6];
        u[r3][r4][t + 1] += r7;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_2 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_2 <= y_M + 1) {
        int r8 = ii_src_0 + 2;
        int r9 = ii_src_2 + 2;
        int r10 = ii_src_0 + 1;
        int r11 = ii_src_2 + 1;
        double r12 = 5.7121e-2F * (-2.5e-1F * px * py + 5.0e-1F * py);
        //*src[time][p_src]/m[r10][r11];
        u[r8][r9][t + 1] += r12;
      }
      if (ii_src_1 >= y_m - 1 && ii_src_3 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= x_M + 1) {
        int r13 = ii_src_3 + 2;
        int r14 = ii_src_1 + 2;
        int r15 = ii_src_3 + 1;
        int r16 = ii_src_1 + 1;
        double r17 = 5.7121e-2F * (-2.5e-1F * px * py + 5.0e-1F * px);
        //*src[time][p_src]/m[r15][r16];
        u[r13][r14][t + 1] += r17;
      }
      if (ii_src_2 >= y_m - 1 && ii_src_3 >= x_m - 1 && ii_src_2 <= y_M + 1 && ii_src_3 <= x_M + 1) {
        int r18 = ii_src_3 + 2;
        int r19 = ii_src_2 + 2;
        int r20 = ii_src_3 + 1;
        int r21 = ii_src_2 + 1;
        double r22 = 1.428025e-2F * px * py;
        //*src[time][p_src]/m[r20][r21];
        u[r18][r19][t + 1] += r22;
      }
    }

    cilk_for(int p_rec = 0; p_rec <= nrcv; p_rec += 1) {
      double r23 = (int)(floor(-5.0e-1F + 5.0e-1F * rec_coords[p_rec][0]));
      int ii_rec_0 = r23 + 10;
      double r24 = (int)(floor(-5.0e-1F + 5.0e-1F * rec_coords[p_rec][1]));
      int ii_rec_1 = r24 + 10;
      int ii_rec_2 = r24 + 11;
      int ii_rec_3 = r23 + 11;
      double px = (double)(-2.0F * r23 + rec_coords[p_rec][0]);
      double py = (double)(-2.0F * r24 + rec_coords[p_rec][1]);
      double sum = 0.0F;
      if (ii_rec_0 >= x_m - 1 && ii_rec_1 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_1 <= y_M + 1) {
        int r25 = ii_rec_0 + 2;
        int r26 = ii_rec_1 + 2;
        sum += (2.5e-1F * px * py - 5.0e-1F * px - 5.0e-1F * py + 1) * u[r25][r26][t + 1];
      }
      if (ii_rec_0 >= x_m - 1 && ii_rec_2 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_2 <= y_M + 1) {
        int r27 = ii_rec_0 + 2;
        int r28 = ii_rec_2 + 2;
        sum += (-2.5e-1F * px * py + 5.0e-1F * py) * u[r27][r28][t + 1];
      }
      if (ii_rec_1 >= y_m - 1 && ii_rec_3 >= x_m - 1 && ii_rec_1 <= y_M + 1 && ii_rec_3 <= x_M + 1) {
        int r29 = ii_rec_3 + 2;
        int r30 = ii_rec_1 + 2;
        sum += (-2.5e-1F * px * py + 5.0e-1F * px) * u[r29][r30][t + 1];
      }
      if (ii_rec_2 >= y_m - 1 && ii_rec_3 >= x_m - 1 && ii_rec_2 <= y_M + 1 && ii_rec_3 <= x_M + 1) {
        int r31 = ii_rec_3 + 2;
        int r32 = ii_rec_2 + 2;
        sum += 2.5e-1F * px * py * u[r31][r32][t + 1];
      }
      //rec[time][p_rec] = sum;
    }

  }

  return u;

}

double ** * jacobi_3d_all(int timesteps, int nrows, int ncols, double ** * grid, int omp_opt) {
  int t = 0;
  int xi = 0;
  int yi = 0;
  int tile_size = 50;

  for (int titer = 0; titer < timesteps; titer++) {

    t = 0;
    int x_m = 1;
    int x_M = nrows - 1;
    int y_m = 1;
    int y_M = ncols - 1;

    // Update core
    if (omp_opt == 1) {
      /* code */
      //printf("\n PARALLEL \n");

      for (int xi = x_m; xi < x_M; xi++) {

        #pragma omp simd

        for (int yi = y_m; yi < y_M; yi += 1) {
          grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi - 1][yi][t] + grid[xi + 1][yi][t] + grid[xi][yi - 1][t] + grid[xi][yi + 1][t]) / 5;
        }
      }
    } else {

    for (int xblk = x_m; xblk < x_M; xblk+=tile_size) {
      for (int xi = xblk; xi < min(x_M, (xblk + tile_size)); xi++) {
        for (int yblk = y_m; yblk < y_M; yblk += tile_size) {
        for (int yi = yblk; yi < min(y_M, (yblk + tile_size)); yi++) {
          grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi - 1][yi][t] + grid[xi + 1][yi][t] + grid[xi][yi - 1][t] + grid[xi][yi + 1][t]) / 5;
        }
      }
    }
   }
  }
    // Update boundary
    for (int yi = 1; yi < y_M - 1; yi++) {
      xi = 0;
      grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi][yi - 1][t] + grid[xi][yi + 1][t] + grid[xi + 1][yi][t]) / 4;
      xi = x_M;
      grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi][yi - 1][t] + grid[xi][yi + 1][t] + grid[xi - 1][yi][t]) / 4;
    }
    for (int xi = 1; xi < x_M - 1; xi++) {
      yi = 0;
      grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi - 1][yi][t] + grid[xi + 1][yi][t] + grid[xi][yi + 1][t]) / 4;
      yi = y_M;
      grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi - 1][yi][t] + grid[xi + 1][yi][t] + grid[xi][yi - 1][t]) / 4;
    }

    grid[0][0][t + 1] = (grid[0][0][t] + grid[0][1][t] + grid[1][0][t]) / 3;
    grid[x_M][y_M][t + 1] = (grid[x_M][y_M][t] + grid[x_M - 1][y_M][t] + grid[x_M][y_M - 1][t]) / 3;

    grid[0][y_M][t + 1] = (grid[0][y_M][t] + grid[0][y_M - 1][t] + grid[1][y_M][t]) / 3;
    grid[x_M][0][t + 1] = (grid[x_M][0][t] + grid[x_M - 1][1][t] + grid[x_M][1][t]) / 3;

    // Update old grid
    for (int xi = 0; xi < nrows; xi++) {
      for (int yi = 0; yi < ncols; yi++) {
        grid[xi][yi][t] = grid[xi][yi][t + 1];
      }
    }

  }
  return grid;

}




double *** jacobi_3d_all_SKEW(int timesteps, int nrows, int ncols, double ** * grid, int omp_opt) {
  int t = 0;
  int xi = 0;
  int yi = 0;
  int tile_size = 50;

  for (int titer = 0; titer < timesteps; titer++) {

    t = 0;
    int x_m = 1;
    int x_M = nrows - 1;
    int y_m = 1;
    int y_M = ncols - 1;

    int tx_skewed = 4;
    int ty_skewed = 4;
    // Update core
    if (omp_opt == 1) {
      /* code */
      //printf("\n PARALLEL \n");


      // TO ADD blocksss


      for (int xi = x_m; xi < x_M; xi++) {
        #pragma omp simd
        for (int yi = y_m; yi < y_M; yi += 1) {
          grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi - 1][yi][t] + grid[xi][yi - 1][t]);
          grid[xi][yi][t + 1] = (grid[xi][yi][t + 1]  + grid[xi + 1][yi][t]  + grid[xi][yi + 1][t]);
          grid[xi][yi][t + 1] = grid[xi][yi][t + 1]/ 5;
        }
      }
    } else {
      //printf("\n Skewed version...");


      for (int xblk = x_m; xblk < x_M; xblk+=tile_size) {
        for (int xi = (xblk + tx_skewed); xi < min( (x_M + tx_skewed), (xblk + tile_size + tx_skewed)); xi++) {
          for (int yblk = y_m; yblk < y_M; yblk += tile_size) {
          for (int yi = (yblk + ty_skewed); yi < min((y_M + ty_skewed), (yblk + tile_size + ty_skewed)); yi++) {



      //for (int xi = (x_m + tx_skewed) ; xi < (x_M +  tx_skewed) ; xi++) {
      //  for (int yi = (y_m + ty_skewed); yi < (y_M + ty_skewed); yi++) {




          grid[xi - tx_skewed][yi - ty_skewed][t + 1] = (grid[xi - tx_skewed][yi - ty_skewed][t] + grid[(xi - tx_skewed) - 1][yi - ty_skewed][t] + grid[(xi - tx_skewed)][(yi - ty_skewed) - 1][t]);
          grid[xi - tx_skewed][yi - ty_skewed][t + 1] = (grid[xi - tx_skewed][yi - ty_skewed][t + 1]  + grid[(xi - tx_skewed) + 1][yi - ty_skewed][t]  + grid[(xi - tx_skewed)][(yi - ty_skewed) + 1][t]);
          grid[xi - tx_skewed][yi - ty_skewed][t + 1] =  grid[xi - tx_skewed][yi - ty_skewed][t + 1]/ 5;        }

        }
      }
      }
    }
    // Update boundary
    for (int yi = 1; yi < y_M - 1; yi++) {
      xi = 0;
      grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi][yi - 1][t] + grid[xi][yi + 1][t] + grid[xi + 1][yi][t]) / 4;
      xi = x_M;
      grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi][yi - 1][t] + grid[xi][yi + 1][t] + grid[xi - 1][yi][t]) / 4;
    }
    for (int xi = 1; xi < x_M - 1; xi++) {
      yi = 0;
      grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi - 1][yi][t] + grid[xi + 1][yi][t] + grid[xi][yi + 1][t]) / 4;
      yi = y_M;
      grid[xi][yi][t + 1] = (grid[xi][yi][t] + grid[xi - 1][yi][t] + grid[xi + 1][yi][t] + grid[xi][yi - 1][t]) / 4;
    }

    grid[0][0][t + 1] = (grid[0][0][t] + grid[0][1][t] + grid[1][0][t]) / 3;
    grid[x_M][y_M][t + 1] = (grid[x_M][y_M][t] + grid[x_M - 1][y_M][t] + grid[x_M][y_M - 1][t]) / 3;

    grid[0][y_M][t + 1] = (grid[0][y_M][t] + grid[0][y_M - 1][t] + grid[1][y_M][t]) / 3;
    grid[x_M][0][t + 1] = (grid[x_M][0][t] + grid[x_M - 1][1][t] + grid[x_M][1][t]) / 3;

    // Update old grid
    for (int xi = 0; xi < nrows; xi++) {
      for (int yi = 0; yi < ncols; yi++) {
        grid[xi][yi][t] = grid[xi][yi][t + 1];
      }
    }

  }
  return grid;

}
