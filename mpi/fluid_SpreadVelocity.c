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

void fluid_SpreadVelocity(LV lv){ //Fluid spread velocity to fiber

  int tid;
  GV gv = lv->gv;
  tid = lv->tid;
  int my_rank = gv->taskid;

#ifdef DEBUG_PRINT
  // printf("Fluid%dtid%d: ****Inside fluid_SpreadVelocity******\n", my_rank, tid);
  fflush(stdout);
#endif //DEBUG_PRINT  

  int total_fibers_row, total_fibers_clmn;
  Fiber* fiberarray;

  Fluidnode* nodes;
  int BI, BJ, BK;
  int li, lj, lk;
  int total_sub_grids, dim_x, dim_y, dim_z;
  int cube_idx, node_idx;

  /*MPI changes*/
  int temp_mac_rank, fluid_owner_mac;
  int fiber_mac_rank = gv->num_fluid_tasks;

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
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;
  total_sub_grids = (dim_x * dim_y * dim_z) / pow(cube_size, 3);

  int fluid_mac_rank;

  // fluid task starts
  if (gv->ifd_recv_count == 1){
    // receive stop message
  	char stop = 100;

#ifdef DEBUG_PRINT
        printf("Fluid%d: send EMPTY msg to %d ifd_max_bufsize=%d\n",
          my_rank, fiber_mac_rank, gv->ifd_max_bufsize);
        fflush(stdout);
#endif //DEBUG_PRINT

    MPI_Send(&stop, 1, MPI_CHAR, fiber_mac_rank, 0, MPI_COMM_WORLD);

    return;
  }
  else{
// #ifdef DEBUG_PRINT
//   printf("**** Fluid task%d_tid%d: fluid_SpreadVelocity prepare MSG******\n", my_rank, tid);
// #endif //DEBUG_PRINT
    int position = 0;
    while (position < gv->ifd_recv_count){

      int X = *((int*)(gv->ifd_recv_buf + position));
      int Y = *((int*)(gv->ifd_recv_buf + position + sizeof(int)));
      int Z = *((int*)(gv->ifd_recv_buf + position + sizeof(int) * 2));

#if 0      
      elastic_force_x = *((double*)(gv->ifd_recv_buf + position + sizeof(int) * 3));
      elastic_force_y = *((double*)(gv->ifd_recv_buf + position + sizeof(int) * 3 + sizeof(double)));
      elastic_force_z = *((double*)(gv->ifd_recv_buf + position + sizeof(int) * 3 + sizeof(double) * 2));
#endif

      position += sizeof(int)* 3 + sizeof(double)* 3;

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
      int toProc;
      owner_tid = cube2thread_and_task(BI, BJ, BK, gv, &toProc);//owner_tid is thread id in the fluid task

      // printf("(my_rank,toProc)=(%d, %d), (X,Y,Z)=(%ld, %ld, %ld), (BI,BJ,BK)=(%ld, %ld, %ld), (li, lj, lk)=(%ld, %ld, %ld)\n", 
      //   my_rank, toProc, X, Y, Z, BI, BJ, BK, li, lj, lk);
      // fflush(stdout);

      assert(my_rank == toProc);

      if (tid == owner_tid){// since stop message alonhg with data is sent to all fluid machines, so N-1 task will recv wrong cube
        // Don't need lock here
        // pthread_mutex_lock(&(gv->lock_Fluid[owner_tid]));
        *((double*)(gv->ifd_recv_buf + position + sizeof(int) * 3)) = nodes[node_idx].vel_x;
        *((double*)(gv->ifd_recv_buf + position + sizeof(int) * 3 + sizeof(double))) = nodes[node_idx].vel_y;
        *((double*)(gv->ifd_recv_buf + position + sizeof(int) * 3 + sizeof(double) * 2)) = nodes[node_idx].vel_z;
        // pthread_mutex_unlock(&(gv->lock_Fluid[owner_tid]));
      }

    }

    pthread_barrier_wait(&(gv->barr));

    if(tid == 0){
    	printf("Fluid%d: send velocity to Fiber task%d, sendcnt=%d, ifd_max_bufsize=%d\n",
          my_rank, fiber_mac_rank, gv->ifd_recv_count, gv->ifd_max_bufsize);
        fflush(stdout);

        assert(gv->ifd_recv_count <= gv->ifd_max_bufsize);

        MPI_Send(gv->ifd_recv_buf, gv->ifd_recv_count, MPI_CHAR, fiber_mac_rank, 0, MPI_COMM_WORLD);
    }

    pthread_barrier_wait(&(gv->barr));

  }


#ifdef DEBUG_PRINT
  printf("**** Fluid%dtid%d: fluid_SpreadVelocity recv MSG and Exit******\n", my_rank, tid);
#endif //DEBUG_PRINT
}