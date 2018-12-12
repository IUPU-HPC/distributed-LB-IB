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

#if 0
void annul_ifd_msg(char* msg, int ifd_size, int point_size){
  for(i = 0; i < ; i){

  }

}

void insert_ifd_msg(char* msg, int position, int point_size, int inneri, int innerj, int innerk, double elastic_force_x, double elastic_force_y, double elastic_force_z){

  for(int i = 0; i < 64 * point_size; i += point_size){
    int X = *((int*)(msg + position + i));
    int Y = *((int*)(msg + position + i + sizeof(int)));
    int Z = *((int*)(msg + position + i + sizeof(int) * 2));

    if ( X == -1 && Y == -1 && Z == -1){ // init state
      *((int*)(msg + position + i)) = inneri;
      *((int*)(msg + position + i + sizeof(int))) = innerj;
      *((int*)(msg + position + i + sizeof(int) * 2)) = innerk;
      *((double*)(msg + position + i + sizeof(int) * 3)) += elastic_force_x;
      *((double*)(msg + position + i + sizeof(int) * 3) + sizeof(double)) += elastic_force_y;
      *((double*)(msg + position + i + sizeof(int) * 3) + sizeof(double) * 2) += elastic_force_z;
      return;
    }
    
    if( X==inneri && Y==innerj && Z==innerk ){ //same fluid point add up elastic force
      *((double*)(msg + position + i + sizeof(int) * 3)) += elastic_force_x;
      *((double*)(msg + position + i + sizeof(int) * 3) + sizeof(double)) += elastic_force_y;
      *((double*)(msg + position + i + sizeof(int) * 3) + sizeof(double) * 2) += elastic_force_z;
      return;
    }
  }
}
#endif

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
  int BI, BJ, BK;
  int li, lj, lk;
  int total_sub_grids, dim_x, dim_y, dim_z;
  int cube_idx, node_idx;
  int P, Q, R, total_threads;
  /* Todo: Move the following variables to GV to save time */
  double  c_x = PI / (2.0*dx);
  double  c_y = PI / (2.0*dy);
  double  c_z = PI / (2.0*dz);
  double  cf = 1.0 / (64.0*dx*dy*dz);
  int ifd_size = 64;
  int point_size = sizeof(int) * 3 + sizeof(double) * 3;

  /*MPI changes*/
  int ifd2FluidProc;
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
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;

  // printf("fiber_SpreadForce: num_cubes_x %ld\n", num_cubes_x);

  total_sub_grids = (dim_x * dim_y * dim_z) / pow(cube_size, 3);

  P = gv->tx;
  Q = gv->ty;
  R = gv->tz;
  total_threads = P*Q*R; //gv->total_threads

#if 0
  //Annul coordinate(x,y,z) in ifd_bufpool to -1, elastic(x,y,z) to 0
  for (int toProc = 0; toProc < num_fluid_tasks; toProc++){
    if(tid == (toProc%total_threads))
      annul_ifd_msg(gv->ifd_bufpool[toProc], point_size);
  }
#endif

  // wait until all fiber threads complete computation
  pthread_barrier_wait(&(gv->barr));

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

          tmp_dist = cf * (1.0e0 + cos(c_x*rx)) * (1.0e0 + cos(c_y*ry)) * (1.0e0 + cos(c_z*rz) );
          elastic_force_x = fiberarray[i].nodes[j].elastic_force_x * tmp_dist;
          elastic_force_y = fiberarray[i].nodes[j].elastic_force_y * tmp_dist;
          elastic_force_z = fiberarray[i].nodes[j].elastic_force_z * tmp_dist;

          //find ifd Fluid toProc
          int ifd2FluidProc = global2task(inneri, innerj, innerk, gv);
#if 0         
          int position = (j + i * total_fibers_clmn) * ifd_size * (sizeof(int) * 3 + sizeof(double) * 3);
#endif 
          pthread_mutex_lock(&gv->lock_ifd_fluid_task_msg[ifd2FluidProc]); //Have repeated (inneri, innerj, innerk) inserted here; Can be updated!

          // printf("Tid%d: Fluid(%d, %d, %d) ifd_last_pos[%d]=%ld\n", 
          //   tid, inneri, innerj, innerk, ifd2FluidProc, gv->ifd_last_pos[ifd2FluidProc]);
          // fflush(stdout);
#if 0
          // map ifd_fluid_coordinate to a specific location in the message
          insert_ifd_msg(gv->ifd_bufpool[ifd2FluidProc], position, point_size, inneri, innerj, innerk, elastic_force_x, elastic_force_y, elastic_force_z);
#endif          
          // naive implementation
          *((int*)(gv->ifd_bufpool[ifd2FluidProc] + gv->ifd_last_pos[ifd2FluidProc]))                  = inneri;
          *((int*)(gv->ifd_bufpool[ifd2FluidProc] + gv->ifd_last_pos[ifd2FluidProc] + sizeof(int)))    = innerj;
          *((int*)(gv->ifd_bufpool[ifd2FluidProc] + gv->ifd_last_pos[ifd2FluidProc] + sizeof(int)* 2)) = innerk;
          *((double*)(gv->ifd_bufpool[ifd2FluidProc] + gv->ifd_last_pos[ifd2FluidProc] + sizeof(int) * 3))                      = elastic_force_x;
          *((double*)(gv->ifd_bufpool[ifd2FluidProc] + gv->ifd_last_pos[ifd2FluidProc] + sizeof(int) * 3 + sizeof(double)))     = elastic_force_y;
          *((double*)(gv->ifd_bufpool[ifd2FluidProc] + gv->ifd_last_pos[ifd2FluidProc] + sizeof(int) * 3 + sizeof(double) * 2)) = elastic_force_z;

          gv->ifd_last_pos[ifd2FluidProc] += sizeof(int) * 3 + sizeof(double) * 3;

          pthread_mutex_unlock(&gv->lock_ifd_fluid_task_msg[ifd2FluidProc]);


        }//for innerk ends
      } //for fiberclmn ends
    }//if fiber2thread ends
  }//for fiber row ends

  // wait until all fiber threads complete computation
  pthread_barrier_wait(&(gv->barr));
  double time_elapsed = Timer::time_end();
  printf("Fiber task%d tid%d: T_prepare_ifd_send_msg=%f\n", gv->taskid, tid, time_elapsed);
  fflush(stdout);

  //fiber task thread 0 send the message out
  // 2018-11-06
  if (tid == 0){
    Timer::time_start();
    for (int toProc = 0; toProc < num_fluid_tasks; toProc++){
      if (gv->ifd_last_pos[toProc] > 0){

        // gv->bufferpool_msg_size[toProc] = gv->ifd_last_pos[toProc];

        gv->influenced_macs[gv->num_influenced_macs] = toProc;
        printf("Fiber task%d: send msg to Fluid task%d, ifd_last_pos[%d]=%d, ifd_max_bufsize=%d, gv->influenced_macs[%d]=%d\n",
          my_rank, toProc, toProc, gv->ifd_last_pos[toProc], gv->ifd_max_bufsize,
          gv->num_influenced_macs, toProc);
        fflush(stdout);

        assert(gv->ifd_last_pos[toProc] <= gv->ifd_max_bufsize);

        MPI_Send(gv->ifd_bufpool[toProc], gv->ifd_last_pos[toProc], MPI_CHAR, toProc, 0, MPI_COMM_WORLD);

        gv->num_influenced_macs++;
        // set ifd_last_pos to default value 0
        // gv->ifd_last_pos[toProc] = 0;
      }
      else{
        char stop = 100;

#ifdef DEBUG_PRINT
        printf("Fiber task%d: send STOP msg to %d ifd_max_bufsize=%d\n",
          my_rank, toProc, gv->ifd_max_bufsize);
        fflush(stdout);
#endif //DEBUG_PRINT

        MPI_Send(&stop, 1, MPI_CHAR, toProc, 0, MPI_COMM_WORLD);
      }
    }
    double time_elapsed = Timer::time_end();
    printf("Fiber task%d: SpreadForce_send_msg_time=%f\n", my_rank, time_elapsed);
    fflush(stdout);
  }
  // 2018-10-06

  //reset
  gv->num_influenced_macs = 0;

#ifdef DEBUG_PRINT
  // printf("----- Fiber task:%d fiber_SpreadForce Exit! -----\n", my_rank);
#endif //DEBUG_PRINT
}
