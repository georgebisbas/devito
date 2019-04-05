#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

double ** * buffer_jacobi_3d(int timesteps, int nrows, int ncols, double ** * grid, int omp_opt, int tile_size) {
  int t = 0;
  int xi = 0;
  int yi = 0;
  int t1=1;int t0=0;
  int x_m = 2;
  int x_M = nrows - 2;
  int y_m = 2;
  int y_M = ncols - 2;

  for (int titer = 0; titer < timesteps-1; titer++) {
    //#pragma omp parallel
    // Update core
  //#pragma omp for
  t1 = titer+1 ;
  t0 = titer ;


  for (int xi = x_m; xi < x_M; xi++) {
    //#pragma omp parallel for
        for (int yi = y_m; yi < y_M; yi++) {
        grid[xi][yi][t1] = (grid[xi][yi][t0] + grid[xi - 2][yi][t0] + grid[xi + 2][yi][t0] + grid[xi][yi - 2][t0] + grid[xi][yi + 2][t0]) / 5;
        }
      }
  }
  return grid;

}

double ** * tile_buffer_jacobi_3d(int timesteps, int nrows, int ncols, double ** * grid, int omp_opt, int tile_size) {
  int t = 0;
  int xi = 0;
  int yi = 0;
  int t1=1;int t0=0;


  int x_m = 2;
  int x_M = nrows - 2;
  int y_m = 2;
  int y_M = ncols - 2;

  int t_blk_size = 5;
  int x_blk_size = 5;

  int tx_skewed = 2;
  int ty_skewed = 2;

//timesteps should be more than tile_size

  printf(" Valid loop interchange!!! \n");


  for (int t_blk = 0; t_blk  < timesteps-1; t_blk  +=  t_blk_size){
    for (int titer = max ( 0, t_blk ); titer < min(timesteps-1,t_blk+t_blk_size); titer++) {
    t1 = titer+1 ;
    t0 = titer;
    //for (int xi = max(x_m + titer, xblk); xi < min((x_M + titer), (xblk + tile_size)); xi++) {
    for (int xblk = x_m; xblk < x_M + timesteps-1; xblk+= x_blk_size) {
    for (int xi = max(x_m + titer,xblk); xi < min( (x_M + titer), (xblk + x_blk_size)); xi++) {
    //#pragma omp simd
    for (int yi = y_m + ty_skewed; yi < y_M + ty_skewed; yi++) {

    //#pragma omp parallel for
        //grid[xi][yi][titer] = (grid[xi][yi][titer-2] + grid[xi - 2][yi][titer-1] + grid[xi + 1][yi][titer-1] + grid[xi][yi - 2][titer-1] + grid[xi][yi + 2][titer-1]) / 5;
        grid[xi - tx_skewed][yi - ty_skewed][t1] = (grid[xi - tx_skewed][yi - ty_skewed][t0] + grid[(xi - tx_skewed) - 2][yi - ty_skewed][t0] + grid[(xi - tx_skewed) + 2][yi - ty_skewed][t0] + grid[xi - tx_skewed][yi - ty_skewed - 2][t0] + grid[xi - tx_skewed][yi - ty_skewed + 2][t0]) / 5;

        }
      }
  }
}
}
  return grid;

}
