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

/* Step1: Influential domain for force-spreading and velocity interpolation. */
/* Step2: Do actual force spreading. */
// eqn 19 21, ..
void fiber_SpreadForce(LV lv){//Fiber influences fluid
#ifdef DEBUG_PRINT
  // printf("****Inside fiber_SpreadForce******\n");
  fflush(stdout);
#endif //DEBUG_PRINT
  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

  int     i = 0, j = 0, inneri = 0, innerj = 0, innerk = 0;//inner i corresponding to fluid index
  int     istart, istop, jstart, jstop, kstart, kstop;
  double  dx = 1.0, dy = 1.0, dz = 1.0; //fluid distance
  double  rx = 0.0, ry = 0.0, rz = 0.0; //local temporary variable

  int     total_fibers_row, total_fibers_clmn;
  Fiber   *fiberarray;

  double  tmp_dist; //distance

  Fluidnode     *nodes;
  long BI, BJ, BK;
  long li, lj, lk;
  long total_sub_grids, dim_x, dim_y, dim_z;
  long cube_idx, node_idx;
  int P, Q, R, total_threads;
  /* Todo: Move the following variables to GV to save time */
  double  PI = 3.14159265358979;
  double  c_x = PI / (2.0*dx);
  double  c_y = PI / (2.0*dy);
  double  c_z = PI / (2.0*dz);
  double  cf = 1.0 / (64.0*dx*dy*dz);

  /*MPI changes*/
  int fluid_owner_mac;
  int num_fluid_tasks = gv->num_fluid_tasks;
  int my_rank = gv->taskid;
  double elastic_force_x, elastic_force_y, elastic_force_z;
  MPI_Status status;

  total_fibers_row = gv->fiber_shape->sheets[0].num_rows;
  total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  fiberarray = gv->fiber_shape->sheets[0].fibers;
  /*Pthread chages*/

  dim_x = gv->fluid_grid->x_dim;
  dim_y = gv->fluid_grid->y_dim;
  dim_z = gv->fluid_grid->z_dim;
  
  int cube_size = gv->cube_size;
  long num_cubes_x = gv->fluid_grid->num_cubes_x;
  long num_cubes_y = gv->fluid_grid->num_cubes_y;
  long num_cubes_z = gv->fluid_grid->num_cubes_z;

  total_sub_grids = (dim_x*dim_y*dim_z) / pow(cube_size, 3);

  P = gv->tx;
  Q = gv->ty;
  R = gv->tz;
  total_threads = P*Q*R;//gv->total_threads

  Timer::time_start();
  for (i = 0; i < total_fibers_row; ++i){
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
      for (j = 0; j < total_fibers_clmn; ++j){
        // for fibre machine
        //Find influential domain
        istart = floor(fiberarray[i].nodes[j].x / dx - 2) + 1; //x dimension
        istop = istart + 3;
        jstart = floor(fiberarray[i].nodes[j].y / dy - 2) + 1; //y dimension
        jstop = jstart + 3;
        kstart = floor(fiberarray[i].nodes[j].z / dz - 2) + 1; // z dimension
        kstop = kstart + 3;

        for (inneri = istart; inneri <= istop; inneri++)//x direction
        for (innerj = jstart; innerj <= jstop; innerj++)//y direction
        for (innerk = kstart; innerk <= kstop; innerk++){//z direction

          /*Used for calculating eqn 21 ....PN*/
          /*distance between the fiber node and all the fluid nodes within the infulential doman.*/
          rx = dx*inneri - fiberarray[i].nodes[j].x;
          ry = dy*innerj - fiberarray[i].nodes[j].y;
          rz = dz*innerk - fiberarray[i].nodes[j].z;

          BI = inneri / cube_size;
          BJ = innerj / cube_size;
          BK = innerk / cube_size;
          li = inneri % cube_size;
          lj = innerj % cube_size;
          lk = innerk % cube_size;

          tmp_dist = cf * (1.0e0 + cos(c_x*rx)) * (1.0e0 + cos(c_y*ry))*(1.0e0 + cos(c_z*rz));
          elastic_force_x = fiberarray[i].nodes[j].elastic_force_x * tmp_dist;
          elastic_force_y = fiberarray[i].nodes[j].elastic_force_y * tmp_dist;
          elastic_force_z = fiberarray[i].nodes[j].elastic_force_z * tmp_dist;

          //find fluid_owner_mac for this fluid node
          int fluid_owner_mac = cube2task(BI, BJ, BK, gv);

          pthread_mutex_lock(&gv->lock_ifd_fluid_task_msg[fluid_owner_mac]);

          *((int*)(gv->ifd_bufpool[fluid_owner_mac] + gv->ifd_last_pos[fluid_owner_mac]))                   = BI;
          *((int*)(gv->ifd_bufpool[fluid_owner_mac] + gv->ifd_last_pos[fluid_owner_mac] + sizeof(int)))     = BJ;
          *((int*)(gv->ifd_bufpool[fluid_owner_mac] + gv->ifd_last_pos[fluid_owner_mac] + sizeof(int)* 2))  = BK;
          *((int*)(gv->ifd_bufpool[fluid_owner_mac] + gv->ifd_last_pos[fluid_owner_mac] + sizeof(int)* 3))  = li;
          *((int*)(gv->ifd_bufpool[fluid_owner_mac] + gv->ifd_last_pos[fluid_owner_mac] + sizeof(int)* 4))  = lj;
          *((int*)(gv->ifd_bufpool[fluid_owner_mac] + gv->ifd_last_pos[fluid_owner_mac] + sizeof(int)* 5))  = lk;
          *((double*)(gv->ifd_bufpool[fluid_owner_mac] + gv->ifd_last_pos[fluid_owner_mac] + sizeof(int)* 6))    = elastic_force_x;
          *((double*)(gv->ifd_bufpool[fluid_owner_mac] + gv->ifd_last_pos[fluid_owner_mac] + sizeof(int)* 6 + sizeof(double))) = elastic_force_y;
          *((double*)(gv->ifd_bufpool[fluid_owner_mac] + gv->ifd_last_pos[fluid_owner_mac] + sizeof(int)* 6 + sizeof(double)* 2)) = elastic_force_z;

          gv->ifd_last_pos[fluid_owner_mac] += sizeof(int) * 6 + sizeof(double) * 3;

          pthread_mutex_unlock(&gv->lock_ifd_fluid_task_msg[fluid_owner_mac]);


        }//for innerk ends
      } //for fiberclmn ends
    }//if fiber2thread ends
  }//for fiber row ends

  // wait until all fiber threads complete computation
  pthread_barrier_wait(&(gv->barr));
  double time_elapsed = Timer::time_end();
  printf("Fiber task %d tid %d: lock_send_buffer=%f\n", gv->taskid, tid, time_elapsed);
  fflush(stdout);

  //fiber task thread 0 send the message out
  if (tid == 0){
    Timer::time_start();
    for (int tmp_fluid_mac = 0; tmp_fluid_mac < num_fluid_tasks; tmp_fluid_mac++){
      if (gv->ifd_last_pos[tmp_fluid_mac] > 0){
        if(gv->ifd_last_pos[tmp_fluid_mac] > gv->ifd_max_bufsize){
          printf("Error: send message cannot be bigger than ifd_max_bufsize! ifd_last_pos[%d]=%d, ifd_max_bufsize=%d\n",
            tmp_fluid_mac, gv->ifd_last_pos[tmp_fluid_mac], gv->ifd_max_bufsize);
          exit(1);
        }
        // gv->bufferpool_msg_size[tmp_fluid_mac] = gv->ifd_last_pos[tmp_fluid_mac];
        gv->influenced_macs[gv->num_influenced_macs]=tmp_fluid_mac;
        printf("Fiber_mac:%d send msg to %d ifd_last_pos[%d]=%d, ifd_max_bufsize=%d, gv->influenced_macs[%d]=%d,\n",
          my_rank, tmp_fluid_mac, tmp_fluid_mac, gv->ifd_last_pos[tmp_fluid_mac], gv->ifd_max_bufsize,
          gv->num_influenced_macs, tmp_fluid_mac);
        fflush(stdout);
        MPI_Send(gv->ifd_bufpool[tmp_fluid_mac], gv->ifd_last_pos[tmp_fluid_mac], MPI_CHAR, tmp_fluid_mac, 0, MPI_COMM_WORLD);


        gv->num_influenced_macs++;
        // set ifd_last_pos to default value 0
        gv->ifd_last_pos[tmp_fluid_mac] = 0;
      }
      else{
        char stop = 100;

#ifdef DEBUG_PRINT
        printf("Fiber_mac:%d send STOP msg to %d ifd_max_bufsize=%d\n",
          my_rank, tmp_fluid_mac, gv->ifd_max_bufsize);
        fflush(stdout);
#endif //DEBUG_PRINT

        MPI_Send(&stop, 1, MPI_CHAR, tmp_fluid_mac, 0, MPI_COMM_WORLD);
      }
    }
    double time_elapsed = Timer::time_end();
    printf("Fiber_mac:%d send_msg_time=%f,\n", my_rank, time_elapsed);
    fflush(stdout);
  }

#ifdef DEBUG_PRINT
  // printf("----- Fiber_mac:%d fiber_SpreadForce Exit! -----\n", my_rank);
#endif //DEBUG_PRINT
}
