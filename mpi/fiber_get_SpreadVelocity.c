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
#include "timer.h"

void search_velocity(char* msg, int recv_cnt, int inneri, int innerj, int innerk, double* vel_x, double* vel_y, double* vel_z){
	int position = 0;
    while (position < recv_cnt){

      int X = *((int*)(msg + position));
      int Y = *((int*)(msg + position + sizeof(int)));
      int Z = *((int*)(msg + position + sizeof(int) * 2));

      if ( inneri==X && innerj==Y && innerk==Z){
      	*vel_x = *((double*)(msg + position + sizeof(int) * 3));
      	*vel_y = *((double*)(msg + position + sizeof(int) * 3 + sizeof(double)));
      	*vel_z = *((double*)(msg + position + sizeof(int) * 3 + sizeof(double) * 2));
      	return;
      }

      position += sizeof(int)* 3 + sizeof(double)* 3;
    }

}

void fiber_get_SpreadVelocity(LV lv){ //Fiber recv spread velocity from Fluid

#ifdef DEBUG_PRINT
  // printf("****Inside fiber_SpreadForce******\n");
  fflush(stdout);
#endif //DEBUG_PRINT

  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

  int    i = 0, j = 0, inneri = 0, innerj = 0, innerk = 0; //inner i corresponding to fluid index
  int    istart, istop, jstart, jstop, kstart, kstop;
  double dx = 1.0, dy = 1.0, dz = 1.0; //fluid distance
  double rx = 0.0, ry = 0.0, rz = 0.0; //local temporary variable

  int    total_fibers_row, total_fibers_clmn;
  Fiber  *fiberarray;

  double tmp; //spreading_velocity coefficient

  Fluidnode *nodes;
  int BI, BJ, BK;
  int li, lj, lk;
  int total_sub_grids, dim_x, dim_y, dim_z;
  int cube_idx, node_idx;
  int P, Q, R, total_threads;
  /* Todo: Move the following variables to GV to save time */
  double PI = 3.14159265358979;
  double c_x = PI / (2.0 * dx);
  double c_y = PI / (2.0 * dy);
  double c_z = PI / (2.0 * dz);
  double c_64 = 1.0 / 64.0;

  /*MPI changes*/
  int ifd2FluidProc;
  int num_fluid_tasks = gv->num_fluid_tasks;
  int my_rank = gv->taskid;
  double elastic_force_x, elastic_force_y, elastic_force_z;
  MPI_Status status;

  double t0 = 0, t1 = 0, t_search = 0;

  total_fibers_row = gv->fiber_shape->sheets[0].num_rows;
  total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  fiberarray = gv->fiber_shape->sheets[0].fibers;

  /*Pthread chages*/
  dim_x = gv->fluid_grid->x_dim;
  dim_y = gv->fluid_grid->y_dim;
  dim_z = gv->fluid_grid->z_dim;
  
  int cube_size = gv->cube_size;
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;

  // printf("fiber_SpreadForce: num_cubes_x %ld\n", num_cubes_x);

  total_sub_grids = (dim_x * dim_y * dim_z) / pow(cube_size, 3);

  P = gv->tx;
  Q = gv->ty;
  R = gv->tz;
  total_threads = P*Q*R; //gv->total_threads

  int recv_cnt; // SpreadVelocity_recv_count
  double s1 = 0, s2 = 0, s3 = 0;

  Timer::time_start();

#if 1
  //fiber task thread 0 send the message out
  if (tid == 0){
    Timer::time_start();
    for (int fromProc = 0; fromProc < num_fluid_tasks; fromProc++){
      MPI_Recv(gv->ifd_bufpool[fromProc], gv->ifd_max_bufsize, MPI_CHAR, fromProc, 0, MPI_COMM_WORLD, &status);	
      MPI_Get_count(&status, MPI_CHAR, &recv_cnt);

      if (recv_cnt > 1){

      	assert(recv_cnt == gv->ifd_last_pos[fromProc]);
        printf("Fiber%d: recv velocity from Fluid%d, recv_cnt=%d, ifd_max_bufsize=%d\n",
          my_rank, fromProc, recv_cnt, gv->ifd_max_bufsize);
        fflush(stdout);

      }
      else{

#ifdef DEBUG_PRINT
        printf("Fiber%d: get EMPTY msg from Fluid%d ifd_max_bufsize=%d\n",
          my_rank, fromProc, gv->ifd_max_bufsize);
        fflush(stdout);
#endif //DEBUG_PRINT

      }
    }
    double time_elapsed = Timer::time_end();
    printf("Fiber task%d: recv_velocity_msg_time=%f,\n", my_rank, time_elapsed);
    fflush(stdout);
  }

#endif

  // other fiber threads wait until tid=0 complete receive
  pthread_barrier_wait(&(gv->barr));

  for (i = 0; i < total_fibers_row; ++i){
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
      for (j = 0; j < total_fibers_clmn; ++j){
        // for fibre machine
        //Find influential domain
        istart = floor(fiberarray[i].nodes[j].x / dx - 2) + 1; //x dimension
        istop = istart + 3;
        jstart = floor(fiberarray[i].nodes[j].y / dy - 2) + 1; //y dimension
        jstop = jstart + 3;
        kstart = floor(fiberarray[i].nodes[j].z / dz - 2) + 1; //z dimension
        kstop = kstart + 3;

        for (inneri = istart; inneri <= istop; inneri++) //x direction
        for (innerj = jstart; innerj <= jstop; innerj++) //y direction
        for (innerk = kstart; innerk <= kstop; innerk++){//z direction

          /*Used for calculating eqn 21 ....PN*/
          /*distance between the fiber node and all the fluid nodes within the infulential doman.*/
          rx = dx*inneri - fiberarray[i].nodes[j].x;
          ry = dy*innerj - fiberarray[i].nodes[j].y;
          rz = dz*innerk - fiberarray[i].nodes[j].z;

          /*notice the difference for spreading and for interpolation here, cf is used for spreading*/
          /*factors in nondimensionalizing h^(-3) and dxdydz cancel each other*/
          tmp = c_64 * (1.0e0 + cos(c_x*rx)) * (1.0e0 + cos(c_y*ry))*(1.0e0 + cos(c_z*rz));

          BI = inneri / cube_size; //BI should be inneri/cube_size but BI_end seems correct
          BJ = innerj / cube_size; 
          BK = innerk / cube_size; 
          li = inneri % cube_size;
          lj = innerj % cube_size;
          lk = innerk % cube_size;
            
          cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          //node_idx = (inneri%cube_size) *cube_size*cube_size +(innerj%cube_size) *cube_size +innerk%cube_size;
          node_idx = (li) * cube_size * cube_size +(lj) * cube_size + lk;

          //find ifd Fluid toProc
          int ifd2FluidProc = global2task(inneri, innerj, innerk, gv);

          double vel_x = 0, vel_y = 0, vel_z = 0;

          t0 = Timer::get_cur_time();
          search_velocity(gv->ifd_bufpool[ifd2FluidProc], gv->ifd_last_pos[ifd2FluidProc], inneri, innerj, innerk, &vel_x, &vel_y, &vel_z); //need to optimize
          t1 = Timer::get_cur_time();
          t_search += t1 - t0;

          s1 += vel_x * tmp; 
      	  s2 += vel_y * tmp;
          s3 += vel_z * tmp;

        }//for innerk ends
      } //for fiberclmn ends
    }//if fiber2thread ends
  }//for fiber row ends

  // wait until all fiber threads complete computation
  pthread_barrier_wait(&(gv->barr));

  // if(tid == 0){
  	printf("Fiber%d prepare to final update location, t_search=%f\n", my_rank, t_search);
  // }
  for (i = 0; i < total_fibers_row; ++i){
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
      for (j = 0; j < total_fibers_clmn; ++j){
        fiberarray[i].nodes[j].x += gv->dt * s1;
      	fiberarray[i].nodes[j].y += gv->dt * s2;
      	fiberarray[i].nodes[j].z += gv->dt * s3;
      } //for fiberclmn ends
    }//if fiber2thread ends
  }//for fiber row ends

  // wait until all fiber threads complete computation
  pthread_barrier_wait(&(gv->barr));

  double time_elapsed = Timer::time_end();
  printf("Fiber task%d tid%d: T_prepare_ifd_send_msg=%f\n", gv->taskid, tid, time_elapsed);
  fflush(stdout);

  //reset
  gv->num_influenced_macs = 0;
  for (int toProc = 0; toProc < num_fluid_tasks; toProc++){  
    // set ifd_last_pos to default value 0
    gv->ifd_last_pos[toProc] = 0;
  }

#ifdef DEBUG_PRINT
  // printf("----- Fiber task:%d fiber_SpreadForce Exit! -----\n", my_rank);
#endif //DEBUG_PRINT
}