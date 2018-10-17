/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include <do_thread.h>

void compute_dest_mac(GV gv, int dir, int my_rank_x, int my_rank_y, int my_rank_z, int Px, int Py, int Pz){
  int dest_x, dest_y, dest_z;
  switch(dir){
    case 1: {
      //dir 1
      dest_x = (my_rank_x+1) % Px;
      dest_y = my_rank_y;
      dest_z = my_rank_z;
      printf("mac%d: dest_x=%d, dest_y=%d, dest_z=%d\n",
        gv->my_rank, dest_x, dest_y, dest_z);
      fflush(stdout);
      // if(dest_x >= Px)
      //   dest_x = 0;
      break;
    }

    case 2: {//dir 2
      dest_x = (my_rank_x-1) % Px;
      dest_y = my_rank_y;
      dest_z = my_rank_z;
      if(dest_x < 0)
        dest_x += Px;
      break;
    }

    case 3: {
      dest_x = my_rank_x;
      dest_y = my_rank_y;
      dest_z = (my_rank_z+1) % Pz;
      // if(dest_z >= Pz)
      //   dest_z = 0;
      break;
    }

    case 4: {
      dest_x = my_rank_x;
      dest_y = my_rank_y;
      dest_z = (my_rank_z-1) % Pz;
      if(dest_z < 0)
        dest_z += Pz;
      break;
    }

    case 5: { //dir 5
      dest_x = my_rank_x;
      dest_y = (my_rank_y-1) % Py;
      dest_z = my_rank_z;
      if(dest_y < 0)
        dest_y += Py;
      break;
    }

    case 6: { //dir 6
      dest_x = my_rank_x;
      dest_y = (my_rank_y+1) % Py;
      dest_z = my_rank_z;
      break;
    }

    case 7: {
      dest_x = (my_rank_x+1) % Px;
      dest_y = my_rank_y;
      dest_z = (my_rank_z+1) % Pz;

      break;
    }

    case 8: {
      dest_x = (my_rank_x-1) % Px;
      dest_y = my_rank_y;
      dest_z = (my_rank_z-1) % Pz;
      if(dest_x < 0)
        dest_x += Px;
      if(dest_z < 0)
        dest_z += Pz;
      break;
    }

    case 9: {
      dest_x = (my_rank_x+1) % Px;
      dest_y = my_rank_y;
      dest_z = (my_rank_z-1) % Pz;
      if(dest_z < 0)
        dest_z += Pz;
      break;
    }

    case 10: {
      dest_x = (my_rank_x-1) % Px;
      dest_y = my_rank_y;
      dest_z = (my_rank_z+1) % Pz;
      if(dest_x < 0)
        dest_x += Px;
      break;
    }

    case 11: {
      dest_x = my_rank_x;
      dest_y = (my_rank_y-1) % Py;
      dest_z = (my_rank_z+1) % Pz;
      if(dest_y < 0)
        dest_y += Py;
      break;
    }

    case 12: {
      dest_x = my_rank_x;
      dest_y = (my_rank_y+1) % Py;
      dest_z = (my_rank_z-1) % Pz;
      if(dest_z < 0)
        dest_z += Pz;
      break;
    }

    case 13: {
      dest_x = my_rank_x;
      dest_y = (my_rank_y+1) % Py;
      dest_z = (my_rank_z+1) % Pz;
      break;
    }

    case 14: {
      dest_x = my_rank_x;
      dest_y = (my_rank_y-1) % Py;
      dest_z = (my_rank_z-1) % Pz;
      if(dest_y < 0)
        dest_y += Py;
      if(dest_z < 0)
        dest_z += Pz;
      break;
    }

    case 15: {
      dest_x = (my_rank_x+1) % Px;
      dest_y = (my_rank_y-1) % Py;
      dest_z = my_rank_z;
      if(dest_y < 0)
        dest_y += Py;
      break;
    }

    case 16: {
      dest_x = (my_rank_x-1) % Px;
      dest_y = (my_rank_y+1) % Py;
      dest_z = my_rank_z;
      if(dest_x < 0)
        dest_x += Px;
      break;
    }

    case 17: {
      dest_x = (my_rank_x+1) % Px;
      dest_y = (my_rank_y+1) % Py;
      dest_z = my_rank_z;
      break;
    }

    case 18: {
      dest_x = (my_rank_x-1) % Px;
      dest_y = (my_rank_y-1) % Py;
      dest_z = my_rank_z;
      if(dest_x < 0)
        dest_x += Px;
      if(dest_y < 0)
        dest_y += Py;
      break;
    }
  }

  gv->dest_mac[dir] = dest_x * Py * Pz + dest_y * Pz + dest_z;
}

void init_stream_msg_surface(GV gv, int start_dir, int end_dir, int dim_1, int dim_2){
  int cube_size =  gv->cube_size;
  int i;
  for (i = start_dir; i <= end_dir; i++){
    gv->stream_last_pos[i] = 0;     // Initialize gv->ifd_last_pos
    // Initialize mutex lock_stream_msg
    if (pthread_mutex_init(&gv->lock_stream_msg[i], NULL)){
      fprintf(stderr, "Unable to initialize lock_stream_msg mutex\n");
      exit(1);
    }
    //1d, along x axis
    gv->stream_msg[i] = (char*)malloc((sizeof(int)* 7 + sizeof(double))* cube_size * dim_1 * cube_size * dim_2 * 5);
  }
}

// void init_stream_msg_surface(GV gv, int dir, int dim_1, int dim_2){

//   gv->stream_last_pos[dir] = 0;     // Initialize gv->ifd_last_pos
//   // Initialize mutex lock_stream_msg
//   if (pthread_mutex_init(&gv->lock_stream_msg[dir], NULL)){
//     fprintf(stderr, "Unable to initialize lock_stream_msg mutex\n");
//     exit(1);
//   }
//   gv->stream_msg[dir] = (char*)malloc((sizeof(int)* 7 + sizeof(double))* cube_size * dim_1 * cube_size * dim_2 * 5);

// }

void init_stream_msg_line(GV gv, int start_dir, int end_dir, int length){
  int cube_size =  gv->cube_size;
  int i;
  for (i = start_dir; i <= end_dir; i++){
    gv->stream_last_pos[i] = 0;     // Initialize gv->ifd_last_pos
    // Initialize mutex lock_stream_msg
    if (pthread_mutex_init(&gv->lock_stream_msg[i], NULL)){
      fprintf(stderr, "Unable to initialize lock_stream_msg mutex\n");
      exit(1);
    }
    //1d, along x axis
    gv->stream_msg[i] = (char*)malloc((sizeof(int)* 7 + sizeof(double))* cube_size * length);
  }
}

// void init_stream_msg_line(GV gv, int dir, int length){
//   gv->stream_last_pos[dir] = 0;     // Initialize gv->ifd_last_pos
//   // Initialize mutex lock_stream_msg
//   if (pthread_mutex_init(&gv->lock_stream_msg[dir], NULL)){
//     fprintf(stderr, "Unable to initialize lock_stream_msg mutex\n");
//     exit(1);
//   }

//   gv->stream_msg[dir] = (char*)malloc((sizeof(int)* 7 + sizeof(double))* cube_size * length);
// }

void init_gv(GV gv, Fibershape* fiber_shape, Fluidgrid* fluid_grid, int cube_size) {
  Fluidnode     *nodes;

  int          i, j;
  int          BI, BJ, BK; //to identify the Sub grids
  int          cube_idx, node_idx;
  int          dim_x, dim_y, dim_z;
  int          total_sub_grids;
  Sub_Fluidgrid *sub_fluid_grid;
  int          num_cubes_x, num_cubes_y, num_cubes_z;
  int          temp_mac_rank;
  int          li, lj, lk;//local access point inside cube
  int          fiber_mac_rank = gv->num_macs - 1;
  int          P,Q,R;
  int          Px, Py, Pz;
  int          my_rank, my_rank_x, my_rank_y, my_rank_z;
  // if (gv->my_rank == gv->Px * gv->Py * gv->Pz){
  // printf("Fiber Machine %d Enter init_gv\n", gv->my_rank);
  if (fiber_shape == NULL) {
    fprintf(stderr, "Error(%s): Null fiber shape!\n", __func__); exit(1);
  }

  gv->fiber_shape = fiber_shape;

  // }

  // if (gv->my_rank != gv->Px * gv->Py * gv->Pz){
  if (fluid_grid == NULL) {
    fprintf(stderr, "Error(%s): Null fluid grid!\n", __func__); exit(1);
  }
  // }



  gv->fluid_grid = fluid_grid;
  sub_fluid_grid = fluid_grid->sub_fluid_grid;
  gv->cube_size = cube_size;
  nodes = gv->fluid_grid->sub_fluid_grid->nodes;


  dim_x = gv->fluid_grid->x_dim;
  dim_y = gv->fluid_grid->y_dim;
  dim_z = gv->fluid_grid->z_dim;

  total_sub_grids = (dim_x*dim_y*dim_z) / pow(cube_size, 3);

  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  P=gv->P;
  Q=gv->Q;
  R=gv->R;

  Px=gv->Px;
  Py=gv->Py;
  Pz=gv->Pz;

  my_rank   = gv->my_rank;
  my_rank_x = gv->my_rank_x;
  my_rank_y = gv->my_rank_y;
  my_rank_z = gv->my_rank_z;

  gv->ib = 2;       //0 and 1 are used for buffer zone
  gv->ie = dim_x - 3; //dim_x-1 (i.e., the last node), dim_x-2 are used for buffer zone
  gv->jb = 2;
  gv->je = dim_y - 3;
  gv->kb = 2;
  gv->ke = dim_z - 3;

  gv->Re = 1.5e2;
  gv->rho_l = 1.0e0;
  gv->u_l = 0.001;       /* choice of u_l should make Ma <0.1 */

  //Todo: move it to fiber sheet!!
  gv->L_l = 2.0e1;       /*the dimensionless char LENGTH, should shortest fiber or width of the sheet should be fibersheet_w*/

  gv->nu_l = gv->u_l * gv->L_l / gv->Re; //viscosity
  gv->tau = 0.5 * (6.0 * gv->nu_l + 1.0);
  gv->Fr = 1.0e10;        /* huge number of Fr means gra is not significant */
  //NO IB if Kbhat and Kshat = 0.0
  gv->Kbhat = 0.005; /* Dimensionless flexure Modulus Kb(.0001-.05) Prateek */
  gv->Kshat = 2.0e1;  /*Streching Compression Coefficient of fiber Kst(20) Prateek */
  gv->Ksthat = 0.0e0;     //TODO:add Ksthat value    /*relative stiffness of the virtual spring tethered to the fixed point*/
  gv->cs_l = 1.0 / (sqrt(3));    /* non-dimensional */ /*cs is the speed of the sound of the model and c is taken to be 1....Why other values of j(from paper not considered) Prateek */

  /* calculate the values of M, Kb, Ks, and g in LB world */
  gv->g_l = gv->u_l * gv->u_l / (gv->Fr * gv->L_l);  /*How...? PN*/
  gv->g_l = 0.0; // gravity is equal to zero.

  gv->Kb_l = gv->Kbhat * gv->rho_l * gv->u_l * gv->u_l * pow(gv->L_l, 4); // For bending Force
  //Not used any more.gv->Kst_l = gv->Ksthat * gv->rho_l * gv->u_l * gv->u_l * gv->L_l;       // For stretching Force may be used fo tethered points only
  gv->Ks_l = gv->Kshat * gv->rho_l * gv->u_l * gv->u_l * gv->L_l;
  //printf("GV->Ks_l:::::%f\n",gv->Ks_l);

  /* discrete particle velocity, 19 velocities: see Figure 1 */
  gv->c[0][0] = 0.0;
  gv->c[0][1] = 0.0;
  gv->c[0][2] = 0;   //direction 0
  gv->c[1][0] = 1.0;
  gv->c[1][1] = 0.0;
  gv->c[1][2] = 0.0; //direction 1
  gv->c[2][0] = -1.0;
  gv->c[2][1] = 0.0;
  gv->c[2][2] = 0.0;
  gv->c[3][0] = 0.0;
  gv->c[3][1] = 0.0;
  gv->c[3][2] = 1;
  gv->c[4][0] = 0.0;
  gv->c[4][1] = 0.0;
  gv->c[4][2] = -1;
  gv->c[5][0] = 0.0;
  gv->c[5][1] = -1.0;
  gv->c[5][2] = 0.0;
  gv->c[6][0] = 0.0;
  gv->c[6][1] = 1.0;
  gv->c[6][2] = 0.0;
  gv->c[7][0] = 1.0;
  gv->c[7][1] = 0.0;
  gv->c[7][2] = 1;
  gv->c[8][0] = -1.0;
  gv->c[8][1] = 0.0;
  gv->c[8][2] = -1;
  gv->c[9][0] = 1.0;
  gv->c[9][1] = 0.0;
  gv->c[9][2] = -1.0;
  gv->c[10][0] = -1.0;
  gv->c[10][1] = 0.0;
  gv->c[10][2] = 1.0;
  gv->c[11][0] = 0.0;
  gv->c[11][1] = -1.0;
  gv->c[11][2] = 1.0;
  gv->c[12][0] = 0.0;
  gv->c[12][1] = 1.0;
  gv->c[12][2] = -1.0;
  gv->c[13][0] = 0.0;
  gv->c[13][1] = 1.0;
  gv->c[13][2] = 1.0;
  gv->c[14][0] = 0.0;
  gv->c[14][1] = -1.0;
  gv->c[14][2] = -1.0;
  gv->c[15][0] = 1.0;
  gv->c[15][1] = -1.0;
  gv->c[15][2] = 0.0;
  gv->c[16][0] = -1.0;
  gv->c[16][1] = 1.0;
  gv->c[16][2] = 0.0;
  gv->c[17][0] = 1.0;
  gv->c[17][1] = 1.0;
  gv->c[17][2] = 0.0;
  gv->c[18][0] = -1.0;
  gv->c[18][1] = -1.0;
  gv->c[18][2] = 0.0; //velocity 18

  // max size may occupy too much memory
  gv->ifd_max_bufsize = 64 * (sizeof(int)* 6 + sizeof(double)* 3) * (gv->fiber_shape->sheets[0].num_rows) * (gv->fiber_shape->sheets[0].num_cols);


  if (gv->my_rank != fiber_mac_rank){
    for (BI = 0; BI < num_cubes_x; ++BI)
    for (BJ = 0; BJ < num_cubes_y; ++BJ)
    for (BK = 0; BK < num_cubes_z; ++BK){
      cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank);
      if (gv->my_rank == temp_mac_rank){ //MPI changes
        cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        nodes = fluid_grid->sub_fluid_grid[cube_idx].nodes;
        for (li = 0; li < cube_size; ++li)
        for (lj = 0; lj < cube_size; ++lj)
        for (lk = 0; lk < cube_size; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk;
          nodes[node_idx].rho = gv->rho_l;//same as rho_l
          nodes[node_idx].vel_x = gv->u_l;
          nodes[node_idx].vel_y = 0.0;
          nodes[node_idx].vel_z = 0.0;
        }
      }
    }

    // printf("row:%d, col:%d, ifd_max_bufsize = %d \n", gv->fiber_shape->sheets[0].num_rows, gv->fiber_shape->sheets[0].num_cols, gv->ifd_max_bufsize);
    // fflush(stdout);
    gv->ifd_recv_count = 0;
    gv->ifd_recv_buf = (char*)malloc(sizeof(char)*gv->ifd_max_bufsize);

    // Initilize stream_msg
    int stream_max_dim = max(max(P*Q, Q*R), P*R);
    // printf("stream_max_dim=%d\n", stream_max_dim);
    gv->stream_recv_max_bufsize = (sizeof(int) * 7 + sizeof(double)) * cube_size * cube_size * stream_max_dim * 5;

    int direction, dest;
    if(Px==1 && Py==1 && Pz==1){
      printf("Fluid mac %d doesn't need to stream\n", my_rank);
      fflush(stdout);
    }
    else if (check_1d(Px, Py, Pz, &direction)){ //1d, any direction=1, function return 0, at most 2 msg
      printf("Fluid mac %d Enter 1D Init configuration! direction=%d\n", my_rank, direction);
      fflush(stdout);

      gv->stream_msg = (char**)malloc(sizeof(char*) * 7);
      // gv->dest_mac = (int *)malloc(sizeof(int) * 7);

      gv->lock_stream_msg = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)* 7);
      gv->stream_last_pos = (int*)malloc(sizeof(int)*7);

      gv->stream_msg[0] = (char*)malloc(gv->stream_recv_max_bufsize);
      gv->stream_recv_buf = gv->stream_msg[0];
      // printf("Here! stream_recv_max_bufsize=%d\n", gv->stream_recv_max_bufsize);
      // fflush(stdout);

      if(direction == X_transfer_1D){
        printf("Fluid mac %d Enter 1D X_transfer_1D configuration!\n", my_rank);
        fflush(stdout);
        //1d, along x axis
        init_stream_msg_surface(gv, 1, 2, Q, R);
        //dir 1
        compute_dest_mac(gv, 1, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);

        //dir 2
        compute_dest_mac(gv, 2, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);

        for(i=1; i<=18; i++){
          printf("Fluid mac%d: dest_mac[%d]=%d\n", my_rank, i, gv->dest_mac[i]);
          fflush(stdout);
        }

        // find_stream_dest(gv, 1);

        //   init_stream_msg_surface(gv, 1, Q, R);

        //   init_stream_msg_surface(gv, 2, Q, R);
      }
      else if(direction == Z_transfer_1D){
        //1d, along z axis
        init_stream_msg_surface(gv, 3, 4, P, Q);
        //dir 3
        compute_dest_mac(gv, 3, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);
        //   init_stream_msg_surface(gv, 3, P, Q);

        //dir 4
        compute_dest_mac(gv, 4, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);
        //   init_stream_msg_surface(gv, 4, P, Q);

      }
      else if(direction == Y_transfer_1D){
        //1d, along y axis
        init_stream_msg_surface(gv, 5, 6, P, R);
        //dir 5
        compute_dest_mac(gv, 5, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);
        //   init_stream_msg_surface(gv, 5, P, R);

        //dir 6
        compute_dest_mac(gv, 6, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);
        // dest = my_rank + 1;
        // if(dest >= Py)
        //   dest = 0;
        // gv->dest_mac[6] = dest;
        //   init_stream_msg_surface(gv, 6, P, R);
      }
      else{
        printf("check_1d init Error!\n");
        fflush(stdout);
      }
    }
    else if (check_2d(Px, Py, Pz, &direction)){//2d, at most 8
      printf("Fluid mac Enter 2D Init configuration!\n");
      fflush(stdout);
      gv->stream_msg = (char**)malloc(sizeof(char*)* 19);

      gv->lock_stream_msg = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)* 19);
      gv->stream_last_pos = (int*)malloc(sizeof(int) * 19);

      gv->stream_msg[0] = (char*)malloc(gv->stream_recv_max_bufsize);
      gv->stream_recv_buf = gv->stream_msg[0];

      if(direction==X_transfer_2D){ //1, 2, 5, 6, 15~18

        init_stream_msg_surface(gv, 1, 2, Q, R);
        init_stream_msg_surface(gv, 5, 6, P, R);
        init_stream_msg_line(gv,15, 18, R);

        //dir 1
        compute_dest_mac(gv, 1, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);

        //dir 2
        compute_dest_mac(gv, 2, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);

        //dir 5
        compute_dest_mac(gv, 5, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);

        //dir 6
        compute_dest_mac(gv, 6, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);


        //dir 15
        for(i=15; i<=18; i++)
          compute_dest_mac(gv, i, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);

      }
      else if(direction==Y_transfer_2D){
        init_stream_msg_surface(gv, 3, 4, P, Q);
        init_stream_msg_surface(gv, 5, 6, P, R);
        init_stream_msg_line(gv,11, 14, P);

        for(i=3; i<=6; i++)
          compute_dest_mac(gv, i, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);
        for(i=11; i<=14; i++)
          compute_dest_mac(gv, i, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);
      }
      else if(direction==Z_transfer_2D){
        init_stream_msg_surface(gv, 1, 2, Q, R);
        init_stream_msg_surface(gv, 3, 4, P, Q);
        init_stream_msg_line(gv,7, 10, Q);

        for(i=1; i<=4; i++)
          compute_dest_mac(gv, i, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);
        for(i=7; i<=10; i++)
          compute_dest_mac(gv, i, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);
      }
      else{
        printf("check_2d init Error!\n");
        fflush(stdout);
      }
    }
    else{ //3d
      printf("Fluid mac Enter 3D Init configuration!\n");
      fflush(stdout);
      gv->stream_msg = (char**)malloc(sizeof(char*) * 19);

      gv->lock_stream_msg = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)* 19);
      gv->stream_last_pos = (int*)malloc(sizeof(int)* 19);

      gv->stream_msg[0] = (char*)malloc(gv->stream_recv_max_bufsize);
      gv->stream_recv_buf = gv->stream_msg[0];

      init_stream_msg_surface(gv, 1, 2, Q, R);
      init_stream_msg_surface(gv, 3, 4, P, Q);
      init_stream_msg_surface(gv, 5, 6, P, R);
      init_stream_msg_line(gv,7, 10, Q);
      init_stream_msg_line(gv,11, 14, P);
      init_stream_msg_line(gv,15, 18, R);

      for(i=1; i<=18; i++)
          compute_dest_mac(gv, i, my_rank_x, my_rank_y, my_rank_z, Px, Py, Pz);
    }


  }

  /*PTHREAD_Change*/

  gv->timesteps= 5;

  //gv->N_WR = gv->TIME_STOP / gv->TIME_WR + 1;//param for Cd need to add later....

  gv->dt = 1;      // Time step size, always equal 1!
  gv->time = gv->dt; // Time step#, starting from 1.


  //Fiber machine initialize
  /******* Shift fiber sheet for one time step ************/
  //TODO: for all the fiber sheets, shift.
  if (gv->my_rank == fiber_mac_rank){
    for (i = 0; i < gv->fiber_shape->sheets[0].num_rows; i++) {
      for (j = 0; j < gv->fiber_shape->sheets[0].num_cols; j++){
        gv->fiber_shape->sheets[0].fibers[i].nodes[j].x += gv->u_l*gv->dt;
        //TODO: if v_l, w_l !=0, then ADD .y, .z = v_l, w_l times gv->dt.
      }
    }

    /*MPI changes*/
    int num_fluid_macs = gv->num_macs - 1;
    gv->ifd_bufpool = (char **) malloc(sizeof(char*) * num_fluid_macs);
    gv->ifd_bufpool_msg_size = (int *) malloc(sizeof(int) * num_fluid_macs);
    gv->ifd_last_pos = (int*) malloc(sizeof(int) * num_fluid_macs);
    gv->lock_ifd_fluid_mac_msg = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t) * num_fluid_macs);
    gv->num_influenced_macs = 0;
    gv->influenced_macs = (int*) malloc(sizeof(int) * num_fluid_macs);

    for (i = 0; i<num_fluid_macs; i++){

      // Initialize ifd_bufpool
      gv->ifd_bufpool[i] = (char*)malloc(sizeof(char)* gv->ifd_max_bufsize);
      gv->ifd_bufpool_msg_size[i] = 0;

      gv->ifd_last_pos[i] = 0;     // Initialize gv->ifd_last_pos

      // Initialize mutex lock_ifd_fluid_mac_msg
      if (pthread_mutex_init(&gv->lock_ifd_fluid_mac_msg[i], NULL)){
        fprintf(stderr, "Unable to initialize lock_ifd_fluid_mac_msg mutex\n");
        exit(1);
      }

      gv->influenced_macs[i] = 0;     // Initialize gv->influenced_macs

    }
  } //endif fiber machine init

  // printf("***********Gv Init exit*****\n");

}
