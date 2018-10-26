/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

/* Step1: Influential domain for force-spreading and velocity interpolation. */
/* Step2: Do actual force spreading. */
// eqn 19 left assign
void fluid_get_SpreadForce(LV lv){//Fiber influences fluid
#ifdef DEBUG_PRINT
  // printf("****Inside fluid_get_SpreadForce******\n");
  fflush(stdout);
#endif //DEBUG_PRINT
  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

  int i = 0, j = 0;//inner i corresponding to fluid index

  int total_fibers_row, total_fibers_clmn;
  Fiber* fiberarray;

  double tmp_dist;

  Fluidnode* nodes;
  long BI, BJ, BK;
  long li, lj, lk;
  long  total_sub_grids, dim_x, dim_y, dim_z;
  long cube_idx, node_idx;
  int  P, Q, R, total_threads;

  /*MPI changes*/
  int temp_mac_rank, fluid_owner_mac;
  int fiber_mac_rank = gv->num_fluid_tasks;
  int my_rank = gv->taskid;
  double elastic_force_x, elastic_force_y, elastic_force_z;
  MPI_Status status;

  total_fibers_row = gv->fiber_shape->sheets[0].num_rows;
  total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  fiberarray = gv->fiber_shape->sheets[0].fibers;

  /*Pthread chages*/
  int owner_tid;
  dim_x = gv->fluid_grid->x_dim;
  dim_y = gv->fluid_grid->y_dim;
  dim_z = gv->fluid_grid->z_dim;

  int cube_size = gv->cube_size;
  long num_cubes_x = gv->fluid_grid->num_cubes_x;
  long num_cubes_y = gv->fluid_grid->num_cubes_y;
  long num_cubes_z = gv->fluid_grid->num_cubes_z;
  total_sub_grids = (dim_x*dim_y*dim_z) / pow(cube_size, 3);

  //Annuling Forces on Fluid grid ::TOO Expensive
  int fluid_mac_rank;

  for (BI = 0; BI <num_cubes_x; BI++)
  for (BJ = 0; BJ <num_cubes_y; BJ++)
  for (BK = 0; BK <num_cubes_z; BK++){
    if (cube2thread_and_task(BI, BJ, BK, gv, &temp_mac_rank) == tid){
      if (temp_mac_rank == my_rank){
        cube_idx = BI*num_cubes_y*num_cubes_z + BJ*num_cubes_z + BK;
        for (li = 0; li< cube_size; li++)
        for (lj = 0; lj< cube_size; lj++)
        for (lk = 0; lk< cube_size; lk++){
          node_idx = li* cube_size * cube_size + lj * cube_size + lk;
          gv->fluid_grid->sub_fluid_grid[cube_idx].nodes[node_idx].elastic_force_x = 0.0e0;
          gv->fluid_grid->sub_fluid_grid[cube_idx].nodes[node_idx].elastic_force_y = 0.0e0;
          gv->fluid_grid->sub_fluid_grid[cube_idx].nodes[node_idx].elastic_force_z = 0.0e0;
        }
      }//if machine check ends
    }//if cube2thread ends
  }


  // fluid machine starts
  // fluid machine thread 0 do receive
  if (tid == 0){
#ifdef DEBUG_PRINT
    printf("fluid_mac%d gv->ifd_max_bufsize=%d\n", my_rank, gv->ifd_max_bufsize);
#endif //DEBUG_PRINT
    MPI_Recv(gv->ifd_recv_buf, gv->ifd_max_bufsize, MPI_CHAR, fiber_mac_rank, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &gv->ifd_recv_count);
    // printf("fluid_mac%d receive a message with ifd_recv_count=%d\n", my_rank, gv->ifd_recv_count);
    // fflush(stdout);
  }

  // All other thread wait until thread 0 finish receiving msg
  pthread_barrier_wait(&(gv->barr));


  if (gv->ifd_recv_count == 1){
    // receive stop message
// #ifdef DEBUG_PRINT
//   printf("**** Fluid_mac%d_tid%d: fluid_get_SpreadForce recv STOP and Exit******\n", my_rank, tid);
// #endif //DEBUG_PRINT
    return;
  }
  else{
// #ifdef DEBUG_PRINT
//   printf("**** Fluid_mac%d_tid%d: fluid_get_SpreadForce recv MSG******\n", my_rank, tid);
// #endif //DEBUG_PRINT
    int position = 0;
    while (position<gv->ifd_recv_count){

      BI = *((int*)(gv->ifd_recv_buf + position));
      BJ = *((int*)(gv->ifd_recv_buf + position + sizeof(int)));
      BK = *((int*)(gv->ifd_recv_buf + position + sizeof(int)* 2));
      li = *((int*)(gv->ifd_recv_buf + position + sizeof(int)* 3));
      lj = *((int*)(gv->ifd_recv_buf + position + sizeof(int)* 4));
      lk = *((int*)(gv->ifd_recv_buf + position + sizeof(int)* 5));
      elastic_force_x = *((double*)(gv->ifd_recv_buf + position + sizeof(int)* 6));
      elastic_force_y = *((double*)(gv->ifd_recv_buf + position + sizeof(int)* 6 + sizeof(double)));
      elastic_force_z = *((double*)(gv->ifd_recv_buf + position + sizeof(int)* 6 + sizeof(double)* 2));

      position += sizeof(int)* 6 + sizeof(double)* 3;

      /*Spreading force*/
      cube_idx = BI *num_cubes_y *num_cubes_z + BJ *num_cubes_z + BK;
      node_idx = (li)*cube_size*cube_size + (lj)*cube_size + lk;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      int fluid_owner_mac;
      owner_tid = cube2thread_and_task(BI, BJ, BK, gv, &fluid_owner_mac);//owner_tid is thread id in the fluid machine
      // Need check my_rank == fluid_owner_mac?
      if (tid == owner_tid){// since stop message alonhg with data is sent to all fluid machines, so N-1 machine will recv wrong cube
        // Don't need lock here
        // pthread_mutex_lock(&(gv->lock_Fluid[owner_tid]));
        nodes[node_idx].elastic_force_x += elastic_force_x;
        nodes[node_idx].elastic_force_y += elastic_force_y;
        nodes[node_idx].elastic_force_z += elastic_force_z;
        // pthread_mutex_unlock(&(gv->lock_Fluid[owner_tid]));
      }

    }

  }


#ifdef DEBUG_PRINT
  // printf("**** Fluid_mac%d: fluid_get_SpreadForce recv MSG and Exit******\n", my_rank);
#endif //DEBUG_PRINT
}
