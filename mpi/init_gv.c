/*  -- Distributed-LB-IB --
 * Copyright 2018 Indiana University Purdue University Indianapolis 
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *   
 * @author: Yuankun Fu (Purdue University, fu121@purdue.edu)
 *
 * @file:
 *
 * @date:
 */

#include <do_thread.h>

extern std::vector<IFDMap> vecOfIfdmap;

int c[19][3] = { //{x, y, z}
            { 0, 0, 0},

            { 1, 0, 0}, {-1, 0, 0}, { 0, 0, 1}, //1, 2, 3
            { 0, 0,-1}, { 0,-1, 0}, { 0, 1, 0}, //4, 5, 6
            { 1, 0, 1}, {-1, 0,-1}, { 1, 0,-1}, //7, 8, 9

            {-1, 0, 1}, { 0,-1, 1}, { 0, 1,-1}, //10, 11, 12
            { 0, 1, 1}, { 0,-1,-1}, { 1,-1, 0}, //13, 14, 15
            {-1, 1, 0}, { 1, 1, 0}, {-1,-1, 0}  //16, 17, 18
}; // stores 19 different velocity directions for ksi

double t[19] = {
            1./3.,

            1./18., 1./18., 1./18., //1~6
            1./18., 1./18., 1./18.,

            1./36., 1./36., 1./36., //7~18
            1./36., 1./36., 1./36.,
            1./36., 1./36., 1./36.,
            1./36., 1./36., 1./36.
};

void init_gv_constant(GV gv){
  int dim_x, dim_y, dim_z;

  // gv->timesteps= 5; // add default value
  //gv->N_WR = gv->TIME_STOP / gv->TIME_WR + 1;//param for Cd need to add later....
  gv->dt = 1;      // Time step size, always equal 1!
  gv->time = gv->dt; // Time step#, starting from 1.

  dim_x = gv->fluid_grid->x_dim;
  dim_y = gv->fluid_grid->y_dim;
  dim_z = gv->fluid_grid->z_dim;

  // fiber also needs fluid grid dimension infomation
  int cube_size = gv->cube_size;
  int num_cubes_x = gv->fluid_grid->num_cubes_x = dim_x / cube_size;
  int num_cubes_y = gv->fluid_grid->num_cubes_y = dim_y / cube_size;
  int num_cubes_z = gv->fluid_grid->num_cubes_z = dim_z / cube_size;

  gv->ib = 2;         //0 and 1 are used for buffer zone
  gv->ie = dim_x - 3; //dim_x-1 (i.e., the last node), dim_x-2 are used for buffer zone
  gv->jb = 2;
  gv->je = dim_y - 3;
  gv->kb = 2;
  gv->ke = dim_z - 3;

  gv->Re = 1.5e2;
  gv->rho_l = 1.0e0;
  gv->u_l = 0.001;       /* choice of u_l should make Ma <0.1 */
  // gv->u_l = 8;

  //Todo: move it to fiber sheet!!
  gv->L_l = 2.0e1;       /*the dimensionless char LENGTH, should shortest fiber or width of the sheet should be fibersheet_w*/

  gv->nu_l = gv->u_l * gv->L_l / gv->Re; //viscosity
  gv->tau = 0.5 * (6.0 * gv->nu_l + 1.0);
  gv->Fr = 1.0e10;        /* huge number of Fr means gra is not significant */
  //NO IB if Kbhat and Kshat = 0.0
  gv->Kbhat = 0.005;  /* Dimensionless flexure Modulus Kb(.0001-.05) Prateek */
  gv->Kshat = 2.0e1;  /*Streching Compression Coefficient of fiber Kst(20) Prateek */
  gv->Ksthat = 0.0e0; //TODO:add Ksthat value /*relative stiffness of the virtual spring tethered to the fixed point*/
  gv->cs_l = 1.0 / (sqrt(3));    /* non-dimensional */ /*cs is the speed of the sound of the model and c is taken to be 1....Why other values of j(from paper not considered) Prateek */

  /* calculate the values of M, Kb, Ks, and g in LB world */
  gv->g_l = gv->u_l * gv->u_l / (gv->Fr * gv->L_l);  /*How...? PN*/
  gv->g_l = 0.0; // gravity is equal to zero.

  gv->Kb_l = gv->Kbhat * gv->rho_l * gv->u_l * gv->u_l * pow(gv->L_l, 4); // For bending Force
  //Not used any more.gv->Kst_l = gv->Ksthat * gv->rho_l * gv->u_l * gv->u_l * gv->L_l;       // For stretching Force may be used fo tethered points only
  gv->Ks_l = gv->Kshat * gv->rho_l * gv->u_l * gv->u_l * gv->L_l;
  //printf("GV->Ks_l:::::%f\n",gv->Ks_l);
}

void init_stream_msg(GV gv, int dir, int points){
  int cube_size =  gv->cube_size;
  int i;

  gv->stream_last_pos[dir] = 0;     // Initialize gv->stream_last_pos
  // Initialize mutex lock_stream_msg
  if (pthread_mutex_init(&gv->lock_stream_msg[dir], NULL)){
    fprintf(stderr, "Unable to initialize lock_stream_msg mutex\n");
    exit(1);
  }
  
  gv->stream_msg[dir] = (char*) malloc ((sizeof(int)* 4 + sizeof(double))* points);
  // printf("Fluid%d: Init stream_msg[%d] = %d\n", gv->taskid, dir, (sizeof(int)* 4 + sizeof(double))* points);

}

void init_gv(GV gv) {
  int i, j;
  int BI, BJ, BK; //to identify the Sub grids
  int cube_idx, node_idx;
  int temp_taskid;
  int li, lj, lk; //local access point inside cube
  int P, Q, R;
  int Px, Py, Pz;
  int my_rank; //my_rank_x, my_rank_y, my_rank_z;
  int cube_size = gv->cube_size;
  int num_fluid_tasks = gv->num_fluid_tasks;
  int num_fiber_tasks = gv->fiber_shape->num_sheets;
  size_t tmp;
  Fluidgrid* fluid_grid = gv->fluid_grid;

  P = gv->tx;
  Q = gv->ty;
  R = gv->tz;

  Px = gv->num_fluid_task_x;
  Py = gv->num_fluid_task_y;
  Pz = gv->num_fluid_task_z;

  int dim_x = fluid_grid->x_dim;
  int dim_y = fluid_grid->y_dim;
  int dim_z = fluid_grid->z_dim;

  my_rank   = gv->taskid;
  // my_rank_x = gv->my_rank_x;
  // my_rank_y = gv->my_rank_y;
  // my_rank_z = gv->my_rank_z;

  // printf("Task%d: Enter init_gv\n", my_rank);

  // determine the ifd_max_bufsize received by a fluid task
  // int ifd_size = 64; //4*4*4, x,y,z: (-2,2)
  gv->ifd_max_bufsize = 0;
  for(i = 0; i < num_fiber_tasks; i++){
    tmp = (sizeof(int)*3 + sizeof(double)*3) * IFD_SIZE * IFD_SIZE *
                        (gv->fiber_shape->sheets[i].width + IFD_SIZE) * 
                        (gv->fiber_shape->sheets[i].height + IFD_SIZE);
    printf("tmp=%d B, width=%f, height=%f\n", 
      tmp, gv->fiber_shape->sheets[i].width,gv->fiber_shape->sheets[i].height);
    if(tmp > gv->ifd_max_bufsize)
      gv->ifd_max_bufsize = tmp;
  }
  // printf("ifd_max_bufsize=%d\n", gv->ifd_max_bufsize);

  // Fluid task
  if (gv->taskid < gv->num_fluid_tasks){

    // printf("row:%d, col:%d, ifd_max_bufsize = %d \n", gv->fiber_shape->sheets[0].num_rows, gv->fiber_shape->sheets[0].num_cols, gv->ifd_max_bufsize);
    // fflush(stdout);
    gv->ifd_recv_count = 0;
    gv->ifd_recv_buf = (char*) malloc(sizeof(char) * gv->ifd_max_bufsize);

    // Initilize stream_msg
    int max_stream_msg_points = max(max(dim_x*dim_y/(Px*Py), dim_x*dim_z/(Px*Pz)), dim_z*dim_y/(Pz*Py)) * 5; //TOO BIG: max_stream_msg_points: stream a 2D surface to neighbour

    gv->stream_recv_max_bufsize = (sizeof(int) * 4 + sizeof(double)) * max_stream_msg_points;

    // use msg[0] as recv buffer
    gv->stream_msg[0] = (char*)malloc(gv->stream_recv_max_bufsize);
    gv->stream_recv_buf = gv->stream_msg[0];
    gv->stream_msg_recv_cnt = 0;

    if (my_rank == 0){
      printf("Fluid%d: stream_recv_max_bufsize=%d\n", my_rank, gv->stream_recv_max_bufsize);
      fflush(stdout);
    }
    
    int destCoord[3], srcCoord[3];
    for (int streamdir = 0; streamdir < 19; ++streamdir) {
      destCoord[0] = gv->rankCoord[0] + c[streamdir][0];
      destCoord[1] = gv->rankCoord[1] + c[streamdir][1];
      destCoord[2] = gv->rankCoord[2] + c[streamdir][2];

      srcCoord[0] = gv->rankCoord[0] - c[streamdir][0];
      srcCoord[1] = gv->rankCoord[1] - c[streamdir][1];
      srcCoord[2] = gv->rankCoord[2] - c[streamdir][2];

      int dest, src;
      MPI_Cart_rank(gv->cartcomm, destCoord, &dest);
      MPI_Cart_rank(gv->cartcomm, srcCoord, &src);

      gv->streamDest[streamdir] = dest;
      gv->streamSrc[streamdir] = src;

#if 1
      printf("Fluid%2d: streamdir=%2d, dest(x,y,z)=(%2d, %2d, %2d), dest=%2d || src(x,y,z)=(%2d, %2d, %2d), src=%2d\n", 
        gv->rank[0], streamdir, destCoord[0], destCoord[1], destCoord[2], dest,
        srcCoord[0], srcCoord[1], srcCoord[2], src);
      fflush(stdout);
#endif
      // if(gv->rank[0] != dest){
        
      //   gv->streamdir[streamdir] = dest;

      //   // if(streamdir < 7){
      //   //   init_stream_msg(gv, streamdir, * 5);
      //   // }
      //   // else{
      //   //   init_stream_msg(gv, streamdir, Q);
      //   // }
      // }
      // else{
      //   gv->streamdir[streamdir] = -1;
      // }

      init_stream_msg(gv, streamdir, max_stream_msg_points); //need to optimize
    }
  }

  //Fiber task initialization
  /******* Shift fiber sheet for one time step ************/
  //TODO: for all the fiber sheets, shift.
  if (gv->taskid >= gv->num_fluid_tasks){

    // fiber taskid in group = gv->rank[1]
    for (i = 0; i < gv->fiber_shape->sheets[gv->rank[1]].num_rows; i++) {
      for (j = 0; j < gv->fiber_shape->sheets[gv->rank[1]].num_cols; j++){
        gv->fiber_shape->sheets[gv->rank[1]].fibers[i].nodes[j].x += gv->u_l*gv->dt;
        //TODO: if v_l, w_l !=0, then ADD .y, .z = v_l, w_l times gv->dt.
      }
    }

    /*MPI changes*/
    vecOfIfdmap.resize(num_fluid_tasks);
    gv->ifd_bufpool = (char **) malloc(sizeof(char*) * num_fluid_tasks);
    // gv->ifd_bufpool_msg_size = (int *) malloc(sizeof(int) * num_fluid_tasks);
    gv->ifd_last_pos = (int*) malloc(sizeof(int) * num_fluid_tasks);
    gv->lock_ifd_fluid_task_msg = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t) * num_fluid_tasks);
    gv->num_influenced_macs = 0;
    gv->influenced_macs = (int*) malloc(sizeof(int) * num_fluid_tasks);

    int max_msg_size = (sizeof(int) * 3 + sizeof(double) * 3) * IFD_SIZE * IFD_SIZE *
                        (gv->fiber_shape->sheets[gv->rank[1]].width + IFD_SIZE) * 
                        (gv->fiber_shape->sheets[gv->rank[1]].height + IFD_SIZE);

    printf("Fiber%d of %d: Init width+IFD_SIZE=%f, height+IFD_SIZE=%f, max_msg_size=%d B\n", 
      gv->rank[1], gv->size[1], 
      gv->fiber_shape->sheets[gv->rank[1]].width + IFD_SIZE, 
      gv->fiber_shape->sheets[gv->rank[1]].height + IFD_SIZE, 
      max_msg_size);                    

    for (i = 0; i < num_fluid_tasks; i++){

      // Initialize ifd_bufpool
      gv->ifd_bufpool[i] = (char*) malloc(sizeof(char) * max_msg_size);
      // gv->ifd_bufpool_msg_size[i] = 0;

      gv->ifd_last_pos[i] = 0;     // Initialize gv->ifd_last_pos

      // Initialize mutex lock_ifd_fluid_task_msg
      if (pthread_mutex_init(&gv->lock_ifd_fluid_task_msg[i], NULL)){
        fprintf(stderr, "Unable to initialize lock_ifd_fluid_task_msg mutex\n");
        exit(1);
      }

      gv->influenced_macs[i] = 0;     // Initialize gv->influenced_macs

    }
  } //endif fiber machine init

  // printf("***********Gv Init exit*****\n");

}
