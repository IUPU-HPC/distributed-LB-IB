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
#include <cstring>
#include "timer.h"

void needs_argument(int i, int argc, const char *flag) {
  if (i+1 >= argc) {
    fprintf(stderr, "error: Flag \"%s\" requires an argument\n", flag);
    abort();
  }
}

// IFDMap ifdmap;
std::vector<std::vector<IFDMap> > Ifdmap_proc_thd;
MsgMap msg_pos;

int main(int argc, char* argv[]) {

  LV lvs;
  GV gv;
  gv = (GV) calloc (1, sizeof(*gv));

#if 0
  int X=0, Y=0, Z=0;
  std::array<int, 3> a;
  a[0]=X; a[1]=Y; a[2]=Z;

  for(int i=0; i<3; i++){
    std::vector<IFDMap> threadmap; //each proc has one threadmap, each thread has a IFDMap
    threadmap.resize(4);
    threadmap[0].insert(std::pair<std::array<int, 3>, int>(a, 0));
    std::cout << "thread[0] size = "<< threadmap[0].size() << "\n";
    std::cout << "Max size is " << threadmap[0].max_size() << "\n";

    Ifdmap_proc_thd.push_back(threadmap);
    printf("Proc[%d] size = %d\n", i, Ifdmap_proc_thd[i].size());
    std::cout << "Ifdmap_proc_thd[0][0] size = "<< Ifdmap_proc_thd[i][0].size() << "\n";
  }
#endif

  // Parse command line
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-steps")) {
      needs_argument(i, argc, "-steps");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-steps %d\" must be > 0\n", value);
        abort();
      }
      gv->timesteps = value;
    }

    if (!strcmp(argv[i], "-dump")) {
      needs_argument(i, argc, "-dump");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-dump %d\" must be > 0\n", value);
        abort();
      }
      gv->dump = value;
    }

    if (!strcmp(argv[i], "-fluid_grid_xyz")) {
      gv->fluid_grid = (Fluidgrid*) calloc(1, sizeof(Fluidgrid));
      // gv->fluid_grid->sub_fluid_grid = (Sub_Fluidgrid*) calloc (1, sizeof(Sub_Fluidgrid));

      needs_argument(i, argc, "-fluid_grid_x");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fluid_grid_x %d\" must be > 0\n", value);
        abort();
      }
      gv->fluid_grid->x_dim = value;

      needs_argument(i, argc, "-fluid_grid_y");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fluid_grid_y %d\" must be > 0\n", value);
        abort();
      }
      gv->fluid_grid->y_dim = value;

      needs_argument(i, argc, "-fluid_grid_z");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fluid_grid_z %d\" must be > 0\n", value);
        abort();
      }
      gv->fluid_grid->z_dim = value;
    }

    if (!strcmp(argv[i], "-cube_size")) {
      needs_argument(i, argc, "-cube_size");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-cube_size %d\" must be > 0\n", value);
        abort();
      }
      gv->cube_size = value;
    }

    if (!strcmp(argv[i], "-fluid_task_xyz")) {
      needs_argument(i, argc, "-fluid_task_x");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fluid_task_x %d\" must be > 0\n", value);
        abort();
      }
      gv->num_fluid_task_x = value;

      needs_argument(i, argc, "-fluid_task_y");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fluid_task_y %d\" must be > 0\n", value);
        abort();
      }
      gv->num_fluid_task_y = value;

      needs_argument(i, argc, "-fluid_task_z");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fluid_task_z %d\" must be > 0\n", value);
        abort();
      }
      gv->num_fluid_task_z = value;
    }

    if (!strcmp(argv[i], "-thread_per_task_xyz")) {
      needs_argument(i, argc, "-thread_per_task_x");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-thread_per_task_x %d\" must be > 0\n", value);
        abort();
      }
      gv->tx = value;

      needs_argument(i, argc, "-thread_per_task_y");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-thread_per_task_y %d\" must be > 0\n", value);
        abort();
      }
      gv->ty = value;

      needs_argument(i, argc, "-thread_per_task_z");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-thread_per_task_z %d\" must be > 0\n", value);
        abort();
      }
      gv->tz = value;
    }

    if (!strcmp(argv[i], "-num_fibersht")) {
      needs_argument(i, argc, "-num_fibersht");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-num_fibersht %d\" must be > 0\n", value);
        abort();
      }
      gv->fiber_shape = (Fibershape*) calloc (1, sizeof(Fibershape));
      gv->fiber_shape->num_sheets = value;
      gv->fiber_shape->sheets = (Fibersheet*) calloc (value, sizeof(Fibersheet));
    }

    if (!strcmp(argv[i], "-fibersht_width_height")) {
      for(int j = 0; j < gv->fiber_shape->num_sheets; j++){
        needs_argument(i, argc, "-fibersht_width");
        double value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_width %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].width = value;

        needs_argument(i, argc, "-fibersht_height");
        value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_height %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].height = value;
      }
    }

    if (!strcmp(argv[i], "-fibersht_row_clmn")) {
      for(int j = 0; j < gv->fiber_shape->num_sheets; j++){
        needs_argument(i, argc, "-fibersht_row");
        int value = atoi(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_row %d\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].num_rows = value;

        needs_argument(i, argc, "-fibersht_clmn");
        value = atoi(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_clmn %d\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].num_cols = value;
      }
    }

    if (!strcmp(argv[i], "-fibersht_xyz_0")) {
      for(int j = 0; j < gv->fiber_shape->num_sheets; j++){
        needs_argument(i, argc, "-fibersht_x0");
        double value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_x0 %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].x_orig = value;

        needs_argument(i, argc, "-fibersht_y0");
        value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_y0 %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].y_orig = value;

        needs_argument(i, argc, "-fibersht_z0");
        value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_z0 %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].z_orig = value;
      }
    }
  }

  //Error Check for Fluid grid can be divided by fluid_task
  if ( (gv->fluid_grid->z_dim % gv->num_fluid_task_z != 0) ||
       (gv->fluid_grid->y_dim % gv->num_fluid_task_y != 0) ||
       (gv->fluid_grid->x_dim % gv->num_fluid_task_x != 0) ){
    fprintf(stderr, "Check Fluid_task x: %d y: %d z: %d\n",
                    gv->num_fluid_task_x, gv->num_fluid_task_y, gv->num_fluid_task_z);
    exit(1);
  }

  //Error Check for Fluid grid can be divided by thread_per_task_xyz
  if ( (gv->fluid_grid->z_dim % gv->tz != 0) ||
       (gv->fluid_grid->y_dim % gv->ty != 0) ||
       (gv->fluid_grid->x_dim % gv->tx != 0) ){
    fprintf(stderr, "Check Fluid_grid_z: %ld and tz: %d Should be multiple of 2 \n", gv->fluid_grid->z_dim, gv->tz);
    fprintf(stderr, "Check Fluid_grid_y: %ld and ty: %d Should be multiple of 2 \n", gv->fluid_grid->y_dim, gv->ty);
    fprintf(stderr, "Check Fluid_grid_x: %ld and tx: %d Should be multiple of 2 \n", gv->fluid_grid->z_dim, gv->tx);
    exit(1);
  }

  //Error check for Fluid grid can be divided by cube_size
  if ( (gv->fluid_grid->x_dim % (gv->cube_size) != 0) ||
       (gv->fluid_grid->y_dim % (gv->cube_size) != 0) ||
       (gv->fluid_grid->z_dim % (gv->cube_size) != 0) ){
    fprintf(stderr, "cube_size: %d \
                     The cubes should be distributed equally on the entire fluid grid\n",
                    gv->cube_size);
    exit(1);
  }

  /*MPI code starts*/
  init_gv_constant(gv);
  // gv->total_tasks = gv->num_fluid_task_x * gv->num_fluid_task_y * gv->num_fluid_task_z + gv->num_fibersht;
  gv->threads_per_task = gv->tx * gv->ty * gv->tz;

  int ierr;
  int provided;

  ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &gv->total_tasks);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &gv->taskid);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &gv->rank[0]);
  MPI_Get_processor_name(gv->processor_name, &gv->namelen);
  printf("Hello world! I’m rank %d of %d on %s\n", gv->rank[0], gv->size[0], gv->processor_name);
  fflush(stdout);

  gv->num_fluid_tasks = gv->total_tasks - gv->fiber_shape->num_sheets;

  printf("%d: width=%d, height=%d\n",
      gv->taskid, gv->fiber_shape->sheets[0].width,
      gv->fiber_shape->sheets[0].height);

  int color;
  if(gv->taskid < gv->num_fluid_tasks)
    //Fuild Group
    color = 0;
  else
    //Fiber Group
    color = 1;
  MPI_Comm mygroup;
  MPI_Comm_split(MPI_COMM_WORLD, color, gv->taskid, &mygroup);
  MPI_Comm_rank(mygroup, &gv->rank[1]);
  MPI_Comm_size(mygroup, &gv->size[1]);

  // print info
  if (gv->taskid == 0){
    printf("***********distributed-LB-IB Simulation using Pthreads cube starts************\n");
    printf("provided=%d\n", provided);
    printf("    Fluidgrid: z %ld, y %ld, x %ld\n",
      gv->fluid_grid->z_dim, gv->fluid_grid->y_dim, gv->fluid_grid->x_dim);
    printf("    Fluid task dimension: Px %d, Py %d, Pz %d\n",
      gv->num_fluid_task_x, gv->num_fluid_task_y, gv->num_fluid_task_z);
    printf("    Fluid threads_per_task: x %d, y %d, z %d ; TuningFactor: cube_size %d\n",
      gv->tx, gv->ty, gv->tz, gv->cube_size);
    fflush(stdout);
  }

  // Fiber tasks: generate fibersheet
  if (gv->taskid >= gv->num_fluid_tasks){
    // ifdmap.insert(std::make_pair(0, 0));
    // ifdmap.clear();
    printf("Friber task%d start gen_fiber_sheet!\n", gv->taskid);
    for(int i = 0; i < gv->fiber_shape->num_sheets; ++i)
      if (gv->taskid == (gv->num_fluid_tasks + i)){
        gen_fiber_sheet(gv->fiber_shape->sheets + i);
        printf("    Fibersheet %d : width %f, height %f, row %ld, clmn %ld, x0 %f, y0 %f, z0 %f\n ",
                     i, gv->fiber_shape->sheets[i].width, gv->fiber_shape->sheets[i].height,
                        gv->fiber_shape->sheets[i].num_rows, gv->fiber_shape->sheets[i].num_cols,
                        gv->fiber_shape->sheets[i].x_orig, gv->fiber_shape->sheets[i].y_orig,
                        gv->fiber_shape->sheets[i].z_orig);
      }
  }
  else{ //Fluid tasks: generate fluid
    // fluid task init my_rank_x,y,z
    int direction;
    int Px, Py, Pz;
    Px = gv->num_fluid_task_x;
    Py = gv->num_fluid_task_y;
    Pz = gv->num_fluid_task_z;

    int ndims, reorder;
    ndims = 3;
    int dim_size[ndims], periods[ndims];

    // old_comm = mygroup;
    int nprocs = gv->num_fluid_tasks;
    dim_size[0] = Px;
    dim_size[1] = Py;
    dim_size[2] = Pz;
    periods[0] = 1; //periodic = false
    periods[1] = 1;
    periods[2] = 1;
    reorder = 1;

    // ierr = MPI_Dims_create (nprocs, ndims, dim_size);
    ierr = MPI_Cart_create(mygroup, ndims, dim_size,
               periods, reorder, &gv->cartcomm);
    if(ierr != 0) printf("ERROR[%d] creating CART\n", ierr);

    /* find my coordinates in the cartesian communicator group */
    ierr = MPI_Cart_coords(gv->cartcomm, gv->taskid, ndims, gv->rankCoord);
    if(ierr != 0) printf("ERROR[%d] MPI_Cart_rank\n", ierr);

    /* use my coordinates to find my rank in cartesian group*/
    int my_cart_rank;
    ierr = MPI_Cart_rank(gv->cartcomm, gv->rankCoord, &my_cart_rank);
    MPI_Comm_size(gv->cartcomm, &gv->size[2]);
    if(ierr != 0) printf("ERROR[%d] MPI_Cart_rank\n", ierr);
    printf("Fluid PW[%d]: my_cart_rank PCM[%d], my coords(x,y,z) = (%d, %d, %d)\n",
      gv->rank[0], my_cart_rank, gv->rankCoord[0], gv->rankCoord[1], gv->rankCoord[2]);

    gen_fluid_grid(gv->fluid_grid, gv->cube_size, gv->taskid, gv);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  init_gv(gv);
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef INIT
  printf("Task%d Pass init gv!\n", gv->taskid);
  fflush(stdout);
#endif

#if 0
  //allocating memory for mutex
  gv->lock_Fluid = (pthread_mutex_t*)malloc(sizeof(*gv->lock_Fluid) * gv->threads_per_task);

  //Error Chek#3
  //initiallise mutex
  for (int i = 0; i < gv->threads_per_task; i++){
    if (pthread_mutex_init(&gv->lock_Fluid[i], NULL)){
      fprintf(stderr, "Unable to initialize a mutex\n");
      exit(1);
    }
  }
#endif

  // Init barrier
  pthread_barrier_t barr;
  if (pthread_barrier_init(&barr, NULL, gv->threads_per_task)){
    fprintf(stderr, "Could not create a barrier\n");
    exit(1);
  }
  gv->barr = barr;

  // Fiber print corner points
  if (gv->taskid >= gv->num_fluid_tasks){
    printf("Friber task%d gen_fiber_sheet complete!\n", gv->taskid);
#if 0
    printf("Printing for Corner Points(z,y) : 0,0 \n");
    print_fiber_sub_grid(gv, 0, 0, 0, 0);
    printf("Printing for Corner Points(z,y) : 51,0 \n");
    print_fiber_sub_grid(gv, 0, 51, 0, 51);
    printf("Printing for Corner Points(z,y) : 0,51 \n");
    print_fiber_sub_grid(gv, 51, 0, 51, 0);
    printf("Printing for Corner Points (z,y): 51,51 \n");
    print_fiber_sub_grid(gv, 51, 51, 51, 51);
    fflush(stdout);
#endif

// #ifdef SAVE //PASS
//     char filename[80];
//     sprintf(filename, "Fiber%d_init.dat", gv->taskid);
//     save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
// #endif
  }
  else{
    init_eqlbrmdistrfuncDF0(gv);
    init_df1(gv);
    init_df_inout(gv);
#ifdef INIT
    printf("Fluid task%d init_eqlbrmdistrfuncDF0, init_df1, init_df_inout complete!\n", gv->taskid);
#endif

#if 0 //PASS init
    char filename[80];
    sprintf(filename, "Fluid%d_init.dat", gv->taskid);
    save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
#endif
  }

  // check_point : fluid grid info
  // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z,
  //   lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
  // print_fluid_sub_grid(gv, 127, 127, 127, 127, 127, 127, gv->cube_size);
  // print_fluid_cube(gv, 4, 0, 0, 4, 0, 0, gv->cube_size);

  double t_start = Timer::get_cur_time();

  //Shared Distribution to Thread in each machine
  pthread_t *threads;
  pthread_attr_t *attrs;
  void *retval;
  threads = (pthread_t*) malloc(sizeof(pthread_t) * gv->threads_per_task);
  attrs = (pthread_attr_t*) malloc(sizeof(pthread_attr_t) * gv->threads_per_task);

  lvs = (LV) malloc(sizeof(*lvs) * gv->threads_per_task);

  gv->threads = threads;
  for (int i = 0; i < gv->threads_per_task; ++i){
    lvs[i].tid = i;
    lvs[i].gv = gv;
    if (pthread_create(threads + i, NULL, do_thread, lvs + i)){
      //if(pthread_create(threads +i, &attr, do_thread,lvs+i )){
      perror("Thread creation do_thread\n");
      printf("Error in Thread creation LBM_IB simulation********\n");
      exit(1);
    }
  }

  for (int i = 0; i < gv->threads_per_task; i++) {
    pthread_join(threads[i], &retval); //to check why retval not allocated
#ifdef DEBUG_PRINT
    printf("Task%d Thread%d is finished\n", gv->taskid, i);
#endif //DEBUG_PRINT
  }

#if 0
  // Fiber print corner points
  if (gv->taskid >= gv->num_fluid_tasks){
    printf("after Running for timesteps: %ld\n", gv->timesteps);
    printf("Printing for Corner Points(z,y) : 0 , 0 \n");
    print_fiber_sub_grid(gv, 0, 0, 0, 0);
    printf("Printing for Corner Points(z,y) : 51, 0 \n");
    print_fiber_sub_grid(gv, 0, 51, 0, 51);
    printf("Printing for Corner Points(z,y) : 0 ,51 \n");
    print_fiber_sub_grid(gv, 51, 0, 51, 0);
    printf("Printing for Corner Points (z,y): 51,51 \n");
    print_fiber_sub_grid(gv, 51, 51, 51, 51);
    fflush(stdout);
  }
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  double time_elapsed = Timer::get_cur_time();

  printf("Task%d: TOTAL TIME TAKEN IN Seconds: %f\n", gv->taskid, time_elapsed - t_start);
  fflush(stdout);

  if(gv->taskid < gv->num_fluid_tasks) MPI_Comm_free(&gv->cartcomm);
    MPI_Comm_free(&mygroup);

  MPI_Finalize();

#if 0
  pthread_mutex_destroy(&gv->lock_Fluid[gv->threads_per_task]);
#endif
  // Need to do
  // pthread_mutex_destroy(&gv->lock_ifd_fluid_task_msg[gv->threads_per_task]);
  free(gv->fluid_grid->sub_fluid_grid);
  free(gv->fluid_grid);
  free(gv);

  return 0;
}

