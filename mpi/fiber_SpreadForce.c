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

extern IFDMap ifdmap;

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

  int     i = 0, j = 0, ii = 0, jj = 0, kk = 0;//inner i corresponding to fluid index
  int     istart, jstart, kstart;
  double  dx = 1.0, dy = 1.0, dz = 1.0; //fluid distance
  double  rx = 0.0, ry = 0.0, rz = 0.0; //local temporary variable

  int     total_fibers_row, total_fibers_clmn;
  Fiber   *fiberarray;
  Fibernode *fibernode;

  double  tmp_dist; //distance

  Fluidnode     *nodes;
  int BI, BJ, BK;
  int li, lj, lk;
  int total_sub_grids, dim_x, dim_y, dim_z;
  int cube_idx, node_idx;
  int total_threads;
  /* Todo: Move the following variables to GV to save time */
  double  c_x = PI / (2.0*dx);
  double  c_y = PI / (2.0*dy);
  double  c_z = PI / (2.0*dz);
  double  cf = 1.0 / (64.0*dx*dy*dz);
  int point_size = sizeof(int) * 3 + sizeof(double) * 3;

  /*MPI changes*/
  int num_fluid_tasks = gv->num_fluid_tasks;
  int my_rank = gv->taskid;
  double elastic_force_x, elastic_force_y, elastic_force_z;
  MPI_Status status;
  double t0, t1, t2, t3, t4, t5;
  double t_lock=0, t_insert=0, t_find=0;

  // ifd
  int ifd_fld_proc;
  int X, Y, Z; //influenced fluid point coordintate
  int ifd_fld_thd;
  long global_index;
  int ifd_msg_pos;
  IFDMap ifdmap;
  std::pair<IFDMap::iterator, bool> ret;
  char* msg;
  int fl_tid;

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

  total_sub_grids = num_cubes_x * num_cubes_y * num_cubes_z;

  total_threads = gv->tx * gv->ty * gv->tz; //gv->total_threads


#ifdef VERIFY
  char fName[80];
  sprintf(fName, "Fiber%d_spread_elastic_force.dat", tid);
  FILE* oFile = fopen(fName, "w");
#endif

  t0 = Timer::get_cur_time();
  for (i = 0; i < total_fibers_row; ++i){
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
      for (j = 0; j < total_fibers_clmn; ++j){

        fibernode = &(fiberarray[i].nodes[j]);

        //Find influential domain
        istart = floor(fibernode->x / dx - 2) + 1; //x dimension
        jstart = floor(fibernode->y / dy - 2) + 1; //y dimension
        kstart = floor(fibernode->z / dz - 2) + 1; //z dimension

        for (ii = 0; ii < IFD_SIZE; ii++) //x direction
        for (jj = 0; jj < IFD_SIZE; jj++) //y direction
        for (kk = 0; kk < IFD_SIZE; kk++){//z direction

          X = istart + ii; 
          Y = jstart + jj;
          Z = kstart + kk;
          global_index = X*dim_y*dim_z + Y*dim_z + Z;

          /*Used for calculating eqn 21 ....PN*/
          /*distance between the fiber node and all the fluid nodes within the infulential doman.*/
          rx = dx*X - fibernode->x;
          ry = dy*Y - fibernode->y;
          rz = dz*Z - fibernode->z;

          tmp_dist = cf * (1.0e0 + cos(c_x*rx)) * (1.0e0 + cos(c_y*ry)) * (1.0e0 + cos(c_z*rz) );

          // NEED this check in case of fiber moving out of bound
          if (X < 0 || X >= dim_x) {
            fprintf(stderr, "X out of bound: %d, x bound=%d!\n", X, gv->fluid_grid->x_dim); exit(1);
          }
          if (Y < 0 || Y >= dim_y) {
            fprintf(stderr, "Y out of bound: %d, y bound=%d!\n", Y, gv->fluid_grid->y_dim); exit(1);
          }
          if (Z < 0 || Z >= dim_z) {
            fprintf(stderr, "Z out of bound: %d, z bound=%d!\n", Z, gv->fluid_grid->z_dim); exit(1);
          }

          elastic_force_x = fibernode->elastic_force_x * tmp_dist;
          elastic_force_y = fibernode->elastic_force_y * tmp_dist;
          elastic_force_z = fibernode->elastic_force_z * tmp_dist;

#ifdef VERIFY
          fprintf(oFile, "elastic_force at (%2d,%2d):(%2d,%2d,%2d): %.24f,%.24f,%.24f\n", 
                      i, j, X, Y, Z, elastic_force_x, elastic_force_y, elastic_force_z);
#endif

          //find ifd Fluid toProc
          fl_tid = global2task(X, Y, Z, gv, &ifd_fld_proc);
          msg = gv->ifd_msg_tid[ifd_fld_proc*total_threads + fl_tid];

#ifdef PERF
          t2 = Timer::get_cur_time();
#endif

          pthread_mutex_lock(&gv->lock_ifd_fluid_task_msg[ifd_fld_proc*total_threads + fl_tid]);

          // t4 = Timer::get_cur_time();
          ret = ifdmap.insert(std::make_pair(global_index, gv->ifd_last_pos[ifd_fld_proc*total_threads + fl_tid]));
          // t5 = Timer::get_cur_time();
          // t_insert += t5 - t4;

          // printf("Tid%d: Fluid(%d, %d, %d) ifd_last_pos[%d]=%d, IFD_SIZE=%d\n", 
          //   tid, X, Y, Z, ifd_fld_proc, gv->ifd_last_pos[ifd_fld_proc], ifdmap.size());
          // fflush(stdout);

          if (ret.second == true){
            
            ifd_msg_pos = gv->ifd_last_pos[ifd_fld_proc*total_threads + fl_tid];

            *((int*)(msg + ifd_msg_pos))                  = X;
            *((int*)(msg + ifd_msg_pos + sizeof(int)))    = Y;
            *((int*)(msg + ifd_msg_pos + sizeof(int)* 2)) = Z;
            *((double*)(msg + ifd_msg_pos + sizeof(int) * 3))                      = elastic_force_x;
            *((double*)(msg + ifd_msg_pos + sizeof(int) * 3 + sizeof(double)))     = elastic_force_y;
            *((double*)(msg + ifd_msg_pos + sizeof(int) * 3 + sizeof(double) * 2)) = elastic_force_z;

            gv->gv->ifd_last_pos[ifd_fld_proc*total_threads+fl_tid] += sizeof(int) * 3 + sizeof(double) * 3;
          }
          else{ //already exist
            // t4 = Timer::get_cur_time();
            // auto res = ifdmap.find(global_index);
            // ifd_msg_pos = res->second;

            ifd_msg_pos = ifdmap[global_index];

            // printf("Tid%d: find existing point Fluid(%d, %d, %d) to Proc%d, ifd_msg_pos=%d, IFD_SIZE=%d\n", 
            // tid, X, Y, Z, ifd_fld_proc, ifd_msg_pos, ifdmap.size());
            // fflush(stdout);
            // t5 = Timer::get_cur_time();
            // t_find += t5 - t4;

            *((double*)(msg + ifd_msg_pos + sizeof(int) * 3))                      += elastic_force_x;
            *((double*)(msg + ifd_msg_pos + sizeof(int) * 3 + sizeof(double)))     += elastic_force_y;
            *((double*)(msg + ifd_msg_pos + sizeof(int) * 3 + sizeof(double) * 2)) += elastic_force_z;

          }
          
          pthread_mutex_unlock(&gv->lock_ifd_fluid_task_msg[ifd_fld_proc]);

#ifdef PERF          
          t3 = Timer::get_cur_time();
          t_lock += t3 - t2;
#endif
          fibernode->ifd_msg_pos[ii][jj][kk] = ifd_msg_pos;


        }//for innerk ends
      } //for fiberclmn ends
    }//if fiber2thread ends
  }//for fiber row ends

#ifdef VERIFY
  fclose(oFile);
#endif

  t1 = Timer::get_cur_time();
  // wait until all fiber threads complete computation
  pthread_barrier_wait(&(gv->barr));

#ifdef PERF
  printf("Fiber task%d tid%d: T_prep_ifd_msg=%f, t_lock=%f, t_insert=%f, t_find=%f\n", 
    gv->taskid, tid, t1-t0, t_lock, t_insert, t_find);
  fflush(stdout);
#endif

  //fiber task thread 0 send the message out
  if (tid == 0){
    t0 = Timer::get_cur_time();
    for (int toProc = 0; toProc < num_fluid_tasks; toProc++){
      if (gv->ifd_last_pos[toProc] > 0){

        // gv->bufferpool_msg_size[toProc] = gv->ifd_last_pos[toProc];

        gv->influenced_macs[gv->num_influenced_macs] = toProc;

#ifdef IFD_FIBER2FLUID        
        printf("-COUNT- Fiber task%d: send msg to Fluid task%d, ifd_last_pos[%d]=%d, ifd_max_bufsize=%d\n",
          my_rank, toProc, toProc, gv->ifd_last_pos[toProc], gv->ifd_max_bufsize);
        fflush(stdout);
#endif
        assert(gv->ifd_last_pos[toProc] <= gv->ifd_max_bufsize);

        MPI_Send(gv->ifd_bufpool[toProc], gv->ifd_last_pos[toProc], MPI_CHAR, toProc, 0, MPI_COMM_WORLD);

        gv->num_influenced_macs++;
        // set ifd_last_pos to default value 0
        // gv->ifd_last_pos[toProc] = 0;
      }
      else{
        char stop = 100;

#ifdef IFD_FIBER2FLUID
        printf("-COUNT- Fiber task%d: send STOP msg to %d ifd_max_bufsize=%d\n",
          my_rank, toProc, gv->ifd_max_bufsize);
        fflush(stdout);
#endif

        MPI_Send(&stop, 1, MPI_CHAR, toProc, 0, MPI_COMM_WORLD);
      }
    }
    t1 = Timer::get_cur_time();

#ifdef PERF
    printf("Fiber task%dtid%d: T_SpreadForce_send_msg=%f\n", my_rank, tid, t1-t0);
    fflush(stdout);
#endif    
  }

  //reset
  if (tid == 0){
    gv->num_influenced_macs = 0;
  }
  

#ifdef DEBUG_PRINT
  printf("----- Fiber task:%d fiber_SpreadForce Exit! -----\n", my_rank);
#endif //DEBUG_PRINT
}
