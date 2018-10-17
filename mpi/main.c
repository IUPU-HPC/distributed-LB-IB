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


#include "do_thread.h"

void needs_argument(int i, int argc, const char *flag) {
  if (i+1 >= argc) {
    fprintf(stderr, "error: Flag \"%s\" requires an argument\n", flag);
    abort();
  }
}

int main(int argc, char* argv[]) {

/*  if (argc != 18) {*/
    /*fprintf(stderr,*/
      /*"Usage: %s -> fibersheet_width, fibersheet_height, \*/
                      /*total_fibers_row, total_fibers_clmn, \*/
                        /*FluidgridZDim, FluidgridYDim, FluidGridXDim, \*/
                          /*fibersheet_xo,  fibersheet_yo,  fibersheet_zo \*/
                            /*cubedimension \*/
                              /*P, Q, R, \*/
                                /*Px, Py, Pz \n",*/
                                  /*argv[0]); //Enter #fiber and #gripoints*/
    /*exit(1);*/
  /*}*/

  // Parse command line
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-steps")) {
      needs_argument(i, argc, "-steps");
      long value = atol(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-steps %ld\" must be > 0\n", value);
        abort();
      }
      gv->timesteps = value;
    }

    if (!strcmp(argv[i], "-numfibers")) {
      needs_argument(i, argc, "-numfibers");
      long value = atol(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-numfibers%ld\" must be > 0\n", value);
        abort();
      }
      gv->num_fiber_task = value;
    }

    if (!strcmp(argv[i], "-fiber_x0")) {
      needs_argument(i, argc, "-fiber_x0");
      long value = atol(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fiber_x0 %ld\" must be > 0\n", value);
        abort();
      }
      gv->timesteps = value;
    }

    if (!strcmp(argv[i], "-setps")) {
      needs_argument(i, argc, "-steps");
      long value = atol(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-steps %ld\" must be > 0\n", value);
        abort();
      }
      gv->timesteps = value;
    }

    if (!strcmp(argv[i], "-setps")) {
      needs_argument(i, argc, "-steps");
      long value = atol(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-steps %ld\" must be > 0\n", value);
        abort();
      }
      gv->timesteps = value;
    }

  }
  double fibersheet_w = atof(argv[1]);
  double fibersheet_h = atof(argv[2]);
  int total_fibers_row = atoi(argv[3]);   /* no of fibres along height */
  int total_fibers_clmn = atoi(argv[4]);  /* no of fibres along width  or column should be 512 for 1024:*/

  int fluidgrid_z = atoi(argv[5]);      // fluid dim z in whole domain
  int fluidgrid_y = atoi(argv[6]);      // fluid dim y in whole domain
  int fluidgrid_x = atoi(argv[7]);      // fluid dim x in whole domain
  double fibersheet_xo = atof(argv[8]); // initial position of the fiber sheet chosen somewher in the middle of grid
  double fibersheet_yo = atof(argv[9]);
  double fibersheet_zo = atof(argv[10]);

  int k_cubedim = atoi(argv[11]); // cubesize

  int P = atoi(argv[12]); //along x num of fluid thread per fluid process
  int Q = atoi(argv[13]); //along y
  int R = atoi(argv[14]); //along z

  int num_fluid_task_x = atoi(argv[15]); //along x num of fluid process
  int num_fluid_task_y = atoi(argv[16]); //along y
  int num_fluid_task_z = atoi(argv[17]); //along z

  int numtasks_= Px*Py*Pz + 1;
  int num_threads_per_process = P*Q*R;

  if (fluidgrid_z%R != 0){
    fprintf(stderr, "Check Fluid_elem_z== :%d and  R==:%d :Should be multiple of 2 \n", fluidgrid_z, R);
    exit(1);
  }

  if (fluidgrid_y%Q != 0){
    fprintf(stderr, "Check Fluid_elem_y== :%d and  Q==:%d :Should be multiple of 2 \n", fluidgrid_y, Q);
    exit(1);
  }

  if (fluidgrid_x%P != 0){
    fprintf(stderr, "Check Fluid_elem_x== :%d and  P==:%d :Should be multiple of 2 \n", fluidgrid_x, P);
    exit(1);
  }

  //New Error Check
  if (fluidgrid_x%Px != 0 || fluidgrid_y%Py != 0 || fluidgrid_z%Pz != 0){
    fprintf(stderr, "Check Fluid_elem_x:%dor elem_y:%d or elem_z== :%d and #machines==:%d :\
                     Lateral Distribution with one machine reserved for Fiber structure\n", 
                    fluidgrid_x, fluidgrid_y, fluidgrid_z, num_macs);
    exit(1);
  }

  //Error check for Fluid nodes and cubes to be divisible
  if (fluidgrid_x % (k_cubedim) != 0 || fluidgrid_y % (k_cubedim) != 0 || fluidgrid_z % (k_cubedim) != 0){
    fprintf(stderr, "Check Fluid_elem_x:%dor elem_y:%d or elem_z== :%d and CubeSIZE:K==:%d : \
                     The cubes should be distributed equally on the entire fluid grid\n", 
                    fluidgrid_x, fluidgrid_y, fluidgrid_z, k_cubedim);
    exit(1);
  }

  /*MPI code starts*/
  int ierr, my_rank;
  int provided;

  ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_macs); // num_macs == num of process
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  LV lvs;
  GV gv;
  gv = (GV) calloc (1, sizeof(*gv));
  gv->num_cubes_x = fluidgrid_x / k_cubedim;
  gv->num_cubes_y = fluidgrid_y / k_cubedim;
  gv->num_cubes_z = fluidgrid_z / k_cubedim;
  gv->P = P;
  gv->Q = Q;
  gv->R = R;
  gv->total_threads = num_threads_per_process;
  gv->num_macs = num_macs;
  gv->my_rank = my_rank;
  gv->Px = Px;
  gv->Py = Py;
  gv->Pz = Pz;
  int fiber_mac_rank = gv->num_macs - 1;

  if (my_rank == fiber_mac_rank){
    printf("***********distributed-LB-IB Simulation using Pthreads cube starts************\n");
    printf("*****INPUT*********\n");
    printf("provided=%d\n", provided);
    printf("    Fibersheet: width %f, height %f, row %d, clmn %d, x0 %f, y0 %f, z0 %f\n ", 
           fibersheet_w, fibersheet_w, total_fibers_row, total_fibers_row, 
           fibersheet_xo, fibersheet_yo, fibersheet_zo);
    printf("    Fluid: elem_z:%d , elem_y :%d, elem_x: %d\n ", fluidgrid_z, fluidgrid_y, fluidgrid_x);
    fflush(stdout);
  }
  else{
    // fluid mac init my_rank_x,y,z
    int direction;
    if (check_1d(Px, Py, Pz, &direction)){
      if (direction == X_transfer_1D){
        gv->my_rank_z = 0;
        gv->my_rank_y = 0;
        gv->my_rank_x = my_rank;
      }
      else if (direction == Z_transfer_1D){
        gv->my_rank_z = my_rank;
        gv->my_rank_y = 0;
        gv->my_rank_x = 0;
      }
      else{//Y_transfer_1D
        gv->my_rank_z = my_rank;
        gv->my_rank_y = 0;
        gv->my_rank_x = 0;
      }
    }
    else if (check_2d(Px, Py, Pz, &direction)){
      if (direction == X_transfer_2D){
        gv->my_rank_z = 0;
        gv->my_rank_y = my_rank % Py;
        gv->my_rank_x = my_rank / Py;
      }
      else if (direction == Y_transfer_2D){
        gv->my_rank_z = my_rank % Pz;
        gv->my_rank_y = my_rank / Pz;
        gv->my_rank_x = 0;
      }
      else{//direction==Z_transfer_2D
        gv->my_rank_z = my_rank % Pz;
        gv->my_rank_y = 0;
        gv->my_rank_x = my_rank / Pz;
      }
    }
    else{
      gv->my_rank_z = my_rank % Pz;
      gv->my_rank_y = my_rank / Pz;
      gv->my_rank_x = my_rank / (Py*Pz);
    }
  }

  // Last machine is the fiber machine
  Fibershape* fiber_shape;
  Fluidgrid*  fluid_grid;

  // if (my_rank == num_macs - 1){
  fiber_shape = gen_fiber_shape(fibersheet_w, fibersheet_h,
                                total_fibers_clmn, total_fibers_row,
                                fibersheet_xo, fibersheet_yo, fibersheet_zo, gv);
  // }
  // else{
  if(my_rank == 0){
    printf("Fluid Machine Rank %d print info of fluid machine\n", my_rank);
    printf("thread dimension per fluid process -- P:%d, Q:%d, R:%d \n TuningFactor:cubedim::%d \n", P, Q, R, k_cubedim);
    printf("***********num_threads_per_process -- P*Q*R =%d  ************\n\n", num_threads_per_process);
    printf("fluid machines dimension -- Px:%d, Py:%d, Pz:%d\n", Px, Py, Pz);
    fflush(stdout);
  }
  fluid_grid = gen_fluid_grid(fluidgrid_x, fluidgrid_y, fluidgrid_z, k_cubedim, gv);
  // }

  init_gv(gv, fiber_shape, fluid_grid, k_cubedim);
  // printf("%d Pass init gv!\n", my_rank);
  // fflush(stdout);
  //allocating memory for mutex
  // gv->lock_Fluid = (pthread_mutex_t*)malloc(sizeof(*gv->lock_Fluid)*num_threads_per_process);


  //Error Chek#3
  //initiallise mutex
  // for (i = 0; i < num_threads_per_process; i++){
  //   if (pthread_mutex_init(&gv->lock_Fluid[i], NULL)){
  //     fprintf(stderr, "Unable to initialize a mutex\n");
  //     exit(1);
  //   }
  // }

  pthread_barrier_t barr;
  //Error Chek#4

  if (pthread_barrier_init(&barr, NULL, num_threads_per_process)){
    fprintf(stderr, "Could not create a barrier\n");
    exit(1);
  }

  gv->barr = barr;

  if (gv->my_rank == fiber_mac_rank){
    printf("After generating the fiber shape:\n");
    printf("Printing for Corner Points(z,y) : 0,0 \n");
    print_fiber_sub_grid(gv, 0, 0, 0, 0);
    printf("Printing for Corner Points(z,y) : 51,0 \n");
    print_fiber_sub_grid(gv, 0, 51, 0, 51);
    printf("Printing for Corner Points(z,y) : 0,51 \n");
    print_fiber_sub_grid(gv, 51, 0, 51, 0);
    printf("Printing for Corner Points (z,y): 51,51 \n");
    print_fiber_sub_grid(gv, 51, 51, 51, 51);
    fflush(stdout);
  }
  else{
    init_eqlbrmdistrfuncDF0(gv);
    init_df1(gv);
    init_df_inout(gv);
  }

  // check_point : fluid grid info
  // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z,
  //   lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
  // print_fluid_sub_grid(gv, 127, 127, 127, 127, 127, 127, gv->cube_size);
  // print_fluid_cube(gv, 4, 0, 0, 4, 0, 0, gv->cube_size);


  double startTime = get_cur_time();
  //Shared Distribution to Thread in each machine
  pthread_t *threads;
  pthread_attr_t *attrs;
  void           *retval;
  threads = (pthread_t*)malloc(sizeof(pthread_t)*num_threads_per_process);
  attrs = (pthread_attr_t*)malloc(sizeof(pthread_attr_t)*num_threads_per_process);

  lvs = (LV)malloc(sizeof(*lvs) * num_threads_per_process);

  gv->threads = threads;
  for (i = 0; i < num_threads_per_process; ++i){
    lvs[i].tid = i;
    lvs[i].gv = gv;
    if (pthread_create(threads + i, NULL, do_thread, lvs + i)){
      //if(pthread_create(threads +i, &attr, do_thread,lvs+i )){
      perror("Thread creation do_thread\n");
      printf("Error in Thread creation LBM_IB simulation********\n");
      exit(1);
    }
  }

  for (i = 0; i < num_threads_per_process; i++) {
    pthread_join(threads[i], &retval); //to check why retval not allocated
#ifdef DEBUG_PRINT
    printf("Thread %d is finished\n", i);
#endif //DEBUG_PRINT
  }

  double endTime = get_cur_time();
  if (gv->my_rank == fiber_mac_rank){
    printf("after Running for Timesteps:%d:\n", gv->time - 1);
    printf("Printing for Corner Points(z,y) : 0,0 \n");
    print_fiber_sub_grid(gv, 0, 0, 0, 0);
    printf("Printing for Corner Points(z,y) : 51,0 \n");
    print_fiber_sub_grid(gv, 0, 51, 0, 51);
    printf("Printing for Corner Points(z,y) : 0,51 \n");
    print_fiber_sub_grid(gv, 51, 0, 51, 0);
    printf("Printing for Corner Points (z,y): 51,51 \n");
    print_fiber_sub_grid(gv, 51, 51, 51, 51);
    fflush(stdout);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  printf("Proc%d: TOTAL TIME TAKEN IN Seconds:%f by machine =%d\n", gv->my_rank, endTime - startTime);
  fflush(stdout);


  MPI_Finalize();

  // pthread_mutex_destroy(&gv->lock_Fluid[num_threads_per_process]);
  // Need to do
  pthread_mutex_destroy(&gv->lock_ifd_fluid_mac_msg[num_threads_per_process]);
  free(gv->fluid_grid->sub_fluid_grid);
  free(gv->fluid_grid);
  free(gv);

  return 0;
}

