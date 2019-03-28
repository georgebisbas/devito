

float *** jacobi_cilk_rcv_src(int timesteps, int nrows, int ncols, int timestamps, int nsrc, int nrcv, float ***u, float **src_coords, float **rec_coords)
{
    int t, xi, yi;
    for (int titer = 0; titer < timesteps-1; titer++) {
    //for (int time = 0, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3); time <= timestamps; time += 1, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3))
//{
        t = 0;

        u[1:nrows-2][1:ncols-2][t+1] = (u[1:nrows-2][1:ncols-2][t] + u[0:nrows-3][1:ncols-2][t] + u[2:nrows-1][1:ncols-2][t] + u[1:nrows-2][0:ncols-3][t] + u[1:nrows-2][2:ncols-1][t] )/6;
        u[0:nrows-1][0:ncols-1][t] = u[0:nrows-1][0:ncols-1][t+1]; //Cilk array notation

        int x_m = 1;
        int x_M = nrows;
        int y_m = 1;
        int y_M = ncols;
        cilk_for (int p_src = 0; p_src < nsrc; p_src++) {
            float r1 = (int)(floor(5.0e-1F*src_coords[p_src][0]));
            int ii_src_0 = r1 + 10;
            float r2 = (int)(floor(5.0e-1F*src_coords[p_src][1]));
            int ii_src_1 = r2 + 10;
            int ii_src_2 = r2 + 11;
            int ii_src_3 = r1 + 11;
            float px = (float)(- 2.0F*r1 + src_coords[p_src][0]);
            float py = (float)(- 2.0F*r2 + src_coords[p_src][1]);
            if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1)
            {
                int r3 = ii_src_0 + 2;
                int r4 = ii_src_1 + 2;
                int r5 = ii_src_0 + 1;
                int r6 = ii_src_1 + 1;
                float r7 = 5.7121e-2F*(2.5e-1F*px*py - 5.0e-1F*px - 5.0e-1F*py + 1);
                //*src[time][p_src]/m[r5][r6];
                u[r3][r4][t+1] += r7;
            }
            if (ii_src_0 >= x_m - 1 && ii_src_2 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_2 <= y_M + 1)
            {
                int r8 = ii_src_0 + 2;
                int r9 = ii_src_2 + 2;
                int r10 = ii_src_0 + 1;
                int r11 = ii_src_2 + 1;
                float r12 = 5.7121e-2F*(-2.5e-1F*px*py + 5.0e-1F*py);
                //*src[time][p_src]/m[r10][r11];
                u[r8][r9][t+1] += r12;
            }
            if (ii_src_1 >= y_m - 1 && ii_src_3 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= x_M + 1)
            {
                int r13 = ii_src_3 + 2;
                int r14 = ii_src_1 + 2;
                int r15 = ii_src_3 + 1;
                int r16 = ii_src_1 + 1;
                float r17 = 5.7121e-2F*(-2.5e-1F*px*py + 5.0e-1F*px);
                //*src[time][p_src]/m[r15][r16];
                u[r13][r14][t+1] += r17;
            }
            if (ii_src_2 >= y_m - 1 && ii_src_3 >= x_m - 1 && ii_src_2 <= y_M + 1 && ii_src_3 <= x_M + 1)
            {
                int r18 = ii_src_3 + 2;
                int r19 = ii_src_2 + 2;
                int r20 = ii_src_3 + 1;
                int r21 = ii_src_2 + 1;
                float r22 = 1.428025e-2F*px*py;
                //*src[time][p_src]/m[r20][r21];
                u[r18][r19][t+1] += r22;
            }
        }

        cilk_for (int p_rec = 0; p_rec <= nrcv; p_rec += 1)
            {
              float r23 = (int)(floor(-5.0e-1F + 5.0e-1F*rec_coords[p_rec][0]));
              int ii_rec_0 = r23 + 10;
              float r24 = (int)(floor(-5.0e-1F + 5.0e-1F*rec_coords[p_rec][1]));
              int ii_rec_1 = r24 + 10;
              int ii_rec_2 = r24 + 11;
              int ii_rec_3 = r23 + 11;
              float px = (float)( - 2.0F*r23 + rec_coords[p_rec][0]);
              float py = (float)( - 2.0F*r24 + rec_coords[p_rec][1]);
              float sum = 0.0F;
              if (ii_rec_0 >= x_m - 1 && ii_rec_1 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_1 <= y_M + 1)
              {
                int r25 = ii_rec_0 + 2;
                int r26 = ii_rec_1 + 2;
                sum += (2.5e-1F*px*py - 5.0e-1F*px - 5.0e-1F*py + 1)*u[r25][r26][t+1];
              }
              if (ii_rec_0 >= x_m - 1 && ii_rec_2 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_2 <= y_M + 1)
              {
                int r27 = ii_rec_0 + 2;
                int r28 = ii_rec_2 + 2;
                sum += (-2.5e-1F*px*py + 5.0e-1F*py)*u[r27][r28][t+1];
              }
              if (ii_rec_1 >= y_m - 1 && ii_rec_3 >= x_m - 1 && ii_rec_1 <= y_M + 1 && ii_rec_3 <= x_M + 1)
              {
                int r29 = ii_rec_3 + 2;
                int r30 = ii_rec_1 + 2;
                sum += (-2.5e-1F*px*py + 5.0e-1F*px)*u[r29][r30][t+1];
              }
              if (ii_rec_2 >= y_m - 1 && ii_rec_3 >= x_m - 1 && ii_rec_2 <= y_M + 1 && ii_rec_3 <= x_M + 1)
              {
                int r31 = ii_rec_3 + 2;
                int r32 = ii_rec_2 + 2;
                sum += 2.5e-1F*px*py*u[r31][r32][t+1];
              }
              //rec[time][p_rec] = sum;
            }

    }



return u;

}
