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

  int total_fibers_row, total_fibers_clmn;
  Fiber* fiberarray;

  Fluidnode* nodes;
  int BI, BJ, BK;
  int li, lj, lk;
  int total_sub_grids, dim_x, dim_y, dim_z;
  int cube_idx, node_idx;
  int  P, Q, R, total_threads;

  /*MPI changes*/
  int temp_mac_rank, fluid_owner_mac;
  int fiber_mac_rank = gv->num_fluid_tasks;
  int my_rank = gv->taskid;
  double elastic_force_x, elastic_force_y, elastic_force_z;
  MPI_Status status;
  int toProc, X, Y, Z;

  total_fibers_row = gv->fiber_shape->sheets[0].num_rows;
  total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  fiberarray = gv->fiber_shape->sheets[0].fibers;

  /*Pthread chages*/
  int owner_tid;
  dim_x = gv->fluid_grid->x_dim;
  dim_y = gv->fluid_grid->y_dim;
  dim_z = gv->fluid_grid->z_dim;

  int cube_size = gv->cube_size;
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;
  total_sub_grids = (dim_x*dim_y*dim_z) / pow(cube_size, 3);

  //Annuling Forces on Fluid grid :: Necessary to annul in every time step, because it is between fluid and fiber
  int fluid_mac_rank;

  for (BI = 0; BI < num_cubes_x; BI++)
  for (BJ = 0; BJ < num_cubes_y; BJ++)
  for (BK = 0; BK < num_cubes_z; BK++){
    owner_tid = cube2thread_and_task(BI, BJ, BK, gv, &toProc);

    if (tid == owner_tid && my_rank == toProc){

      cube_idx = BI*num_cubes_y*num_cubes_z + BJ*num_cubes_z + BK;
      for (li = 0; li < cube_size; li++)
      for (lj = 0; lj < cube_size; lj++)
      for (lk = 0; lk < cube_size; lk++){
#if 0
        X = BI * cube_size + li;
        Y = BJ * cube_size + lj;
        Z = BK * cube_size + lk; 
        printf("Tid%d write (%d, %d, %d)\n", tid, X, Y, Z);
#endif
        node_idx = li* cube_size * cube_size + lj * cube_size + lk;
        gv->fluid_grid->sub_fluid_grid[cube_idx].nodes[node_idx].elastic_force_x = 0.0e0;
        gv->fluid_grid->sub_fluid_grid[cube_idx].nodes[node_idx].elastic_force_y = 0.0e0;
        gv->fluid_grid->sub_fluid_grid[cube_idx].nodes[node_idx].elastic_force_z = 0.0e0;
      }
    }//if cube2thread ends
  }


  // fluid task starts
  // fluid task thread 0 do receive
  if (tid == 0){
#ifdef DEBUG_PRINT
    printf("Fluid task%d gv->ifd_max_bufsize=%d\n", my_rank, gv->ifd_max_bufsize);
#endif //DEBUG_PRINT
    MPI_Recv(gv->ifd_recv_buf, gv->ifd_max_bufsize, MPI_CHAR, fiber_mac_rank, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &gv->ifd_recv_count);
    printf("Fluid task%d tid%d receive a message with ifd_recv_count=%d\n", my_rank, tid, gv->ifd_recv_count);
    fflush(stdout);
  }

  // All other thread wait until thread 0 finish receiving msg
  pthread_barrier_wait(&(gv->barr));


  if (gv->ifd_recv_count == 1){
    // receive stop message
// #ifdef DEBUG_PRINT
//   printf("**** Fluid task%d_tid%d: fluid_get_SpreadForce recv STOP and Exit******\n", my_rank, tid);
// #endif //DEBUG_PRINT
    return;
  }
  else{
// #ifdef DEBUG_PRINT
//   printf("**** Fluid task%d_tid%d: fluid_get_SpreadForce recv MSG******\n", my_rank, tid);
// #endif //DEBUG_PRINT
    int position = 0;
    while (position < gv->ifd_recv_count){

      X = *((int*)(gv->ifd_recv_buf + position));
      Y = *((int*)(gv->ifd_recv_buf + position + sizeof(int)));
      Z = *((int*)(gv->ifd_recv_buf + position + sizeof(int) * 2));
      elastic_force_x = *((double*)(gv->ifd_recv_buf + position + sizeof(int) * 3));
      elastic_force_y = *((double*)(gv->ifd_recv_buf + position + sizeof(int) * 3 + sizeof(double)));
      elastic_force_z = *((double*)(gv->ifd_recv_buf + position + sizeof(int) * 3 + sizeof(double)* 2));

      BI = X / cube_size;
      BJ = Y / cube_size;
      BK = Z / cube_size;

      li = X % cube_size;
      lj = Y % cube_size;
      lk = Z % cube_size;

      // printf("(X,Y,Z)=(%ld, %ld, %ld), (BI,BJ,BK)=(%ld, %ld, %ld), (li, lj, lk)=(%ld, %ld, %ld)\n", 
      //   X, Y, Z, BI, BJ, BK, li, lj, lk);
      // fflush(stdout);

      /*Spreading force*/
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      node_idx = li * cube_size * cube_size + lj * cube_size + lk;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      owner_tid = cube2thread_and_task(BI, BJ, BK, gv, &toProc);//owner_tid is thread id in the fluid task

      // printf("(my_rank,toProc)=(%d, %d), (X,Y,Z)=(%ld, %ld, %ld), (BI,BJ,BK)=(%ld, %ld, %ld), (li, lj, lk)=(%ld, %ld, %ld)\n", 
      //   my_rank, toProc, X, Y, Z, BI, BJ, BK, li, lj, lk);
      // fflush(stdout);

      assert(my_rank == toProc);

      if (tid == owner_tid){// since stop message alonhg with data is sent to all fluid machines, so N-1 task will recv wrong cube
        // Don't need lock here
        // pthread_mutex_lock(&(gv->lock_Fluid[owner_tid]));
        nodes[node_idx].elastic_force_x += elastic_force_x;
        nodes[node_idx].elastic_force_y += elastic_force_y;
        nodes[node_idx].elastic_force_z += elastic_force_z;
        // pthread_mutex_unlock(&(gv->lock_Fluid[owner_tid]));
        // printf("Tid%d update (%d,%d,%d) || elastic_force (%.6f,%.24f,%.24f)\n", 
        //   tid, X, Y, Z, elastic_force_x, elastic_force_y, elastic_force_z);
      }

      position += sizeof(int) * 3 + sizeof(double) * 3;

    }

  }


#ifdef DEBUG_PRINT
  printf("**** Fluid task%d: fluid_get_SpreadForce recv MSG and Exit******\n", my_rank);
#endif //DEBUG_PRINT
}
