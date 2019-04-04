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

extern std::vector<std::vector<IFDMap> > Ifdmap_proc_thd;
extern MsgMap msg_pos;

#if 0
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
  printf("Error: Cannot find velocity\n");
  exit(0);
}
#endif

void fiber_get_SpreadVelocity(LV lv){ //Fiber recv spread velocity from Fluid

#ifdef DEBUG_PRINT
  // printf("****Inside fiber_get_SpreadVelocity******\n");
  fflush(stdout);
#endif //DEBUG_PRINT

  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

  int    i = 0, j = 0; //inner i corresponding to fluid index
  int    ii, jj, kk;
  int    istart, jstart, kstart;
  double dx = 1.0, dy = 1.0, dz = 1.0; //fluid distance
  double rx = 0.0, ry = 0.0, rz = 0.0; //local temporary variable

  // double dx1, dy1, dz1;
  // dx1 = 1.0/dx;
  // dy1 = 1.0/dy;
  // dz1 = 1.0/dz;

  int    total_fibers_row, total_fibers_clmn;
  Fiber  *fiberarray;
  Fibernode *fibernode;

  double tmp; //spreading_velocity coefficient

  Fluidnode *nodes;
  int BI, BJ, BK;
  int li, lj, lk;
  int total_sub_grids, dim_x, dim_y, dim_z;
  long cube_idx;
  int node_idx;
  int P, Q, R, total_threads;

  /* Todo: Move the following variables to GV to save time */
  double c_x = PI / (2.0 * dx);
  double c_y = PI / (2.0 * dy);
  double c_z = PI / (2.0 * dz);
  double c_64 = 1.0 / 64.0;

  /*MPI changes*/
  int num_fluid_tasks = gv->num_fluid_tasks;
  int my_rank = gv->taskid;
  double elastic_force_x, elastic_force_y, elastic_force_z;
  MPI_Status status;
  //ifd
  int ifd_fld_proc;
  int X, Y, Z; //influenced fluid point coordintate
  // int ifd_fld_thd;
  // long global_index;
  int ifd_msg_pos;
  // std::pair<IFDMap::iterator, bool> ret;
  std::array<int, 3> arr3;
  std::array<int, 2> arr2;
  int fl_tid;

  double t0 = 0, t1 = 0, t_search = 0, t2 = 0, t3 = 0;

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

  int recv_cnt; // SpreadVelocity_recv_count
  double s1 = 0, s2 = 0, s3 = 0;
  double vel_x = 0, vel_y = 0, vel_z = 0;

#ifdef IFD_VEL_PERF
  t2 = Timer::get_cur_time();
#endif

  //fiber task thread 0 send the message out
  if (tid == 0){
#ifdef IFD_VEL_PERF
    t0 = Timer::get_cur_time();
#endif
    for (int fromProc = 0; fromProc < num_fluid_tasks; fromProc++){
      MPI_Recv(gv->ifd_send_msg[fromProc], gv->ifd_max_bufsize, MPI_CHAR, fromProc, 0, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_CHAR, &recv_cnt);

      if (recv_cnt > 1){

      	assert(recv_cnt == gv->ifd_last_pos[fromProc]);

#ifdef IFD_FLUID2FIBER
        printf("-COUNT- Fiber%d: recv velocity from Fluid%d, recv_cnt=%d, ifd_max_bufsize=%d\n",
          my_rank, fromProc, recv_cnt, gv->ifd_max_bufsize);
        fflush(stdout);
#endif
      }
      else{

#ifdef DEBUG_PRINT
        printf("-COUNT- Fiber%d: get EMPTY msg from Fluid%d ifd_max_bufsize=%d\n",
          my_rank, fromProc, gv->ifd_max_bufsize);
        fflush(stdout);
#endif //DEBUG_PRINT

      }
    }

#ifdef IFD_VEL_PERF
    t1 = Timer::get_cur_time();
    printf("Fiber task%d: T_recv_vel_msg=%f\n", my_rank, t1 - t0);
    fflush(stdout);
#endif
  }

  // other fiber threads wait until tid=0 complete receive
  pthread_barrier_wait(&(gv->barr));

#ifdef IFD_FLUID2FIBER_DUMP
  char fName[80];
  sprintf(fName, "Fiber%d_vel.dat", tid);
  FILE* oFile = fopen(fName, "w");
#endif

  for (i = 0; i < total_fibers_row; ++i){
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
      for (j = 0; j < total_fibers_clmn; ++j){
        s1=0;
        s2=0;
        s3=0;

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

          /*Used for calculating eqn 21 ....PN*/
          /*distance between the fiber node and all the fluid nodes within the infulential doman.*/
          rx = dx*X - fibernode->x;
          ry = dy*Y - fibernode->y;
          rz = dz*Z - fibernode->z;

          /* You should not use tmp variable, since multiplication sequence influence result*/
          // tmp = c_64 * (1.0 + cos(c_x*rx)) * (1.0 + cos(c_y*ry)) * (1.0 + cos(c_z*rz));

#if 0
          BI = inneri / cube_size; //BI should be inneri/cube_size but BI_end seems correct
          BJ = innerj / cube_size;
          BK = innerk / cube_size;
          li = inneri % cube_size;
          lj = innerj % cube_size;
          lk = innerk % cube_size;

          cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          //node_idx = (inneri%cube_size) *cube_size*cube_size +(innerj%cube_size) *cube_size +innerk%cube_size;
          node_idx = (li) * cube_size * cube_size +(lj) * cube_size + lk;
#endif
          //find ifd Fluid toProc
          fl_tid = global2task_and_thread(X, Y, Z, gv, &ifd_fld_proc);
          arr2[0] = ifd_fld_proc;
          arr2[1] = fl_tid;
          ifd_msg_pos = fibernode->ifd_msg_pos[ii][jj][kk] + msg_pos[arr2];

#if 0
          printf("Tid%d: fibernode(%d, %d)read local_ifd_msg_pos[%d][%d][%d]=%d, (%d, %d, %d):(%d,%d) --> %d, ifd_msg_pos=%d\n",
              tid, i, j, ii, jj, kk, fibernode->ifd_msg_pos[ii][jj][kk],
              X, Y, Z, ifd_fld_proc, fl_tid, msg_pos[arr2],
              ifd_msg_pos);
          fflush(stdout);
#endif
          vel_x = 0;
          vel_y = 0;
          vel_z = 0;

          t0 = Timer::get_cur_time();

          // search_velocity(, gv->ifd_last_pos[ifd_fld_proc], inneri, innerj, innerk, &vel_x, &vel_y, &vel_z); //need to optimize

          vel_x = *((double*)(gv->ifd_send_msg[ifd_fld_proc] + ifd_msg_pos + sizeof(int) * 3));
          vel_y = *((double*)(gv->ifd_send_msg[ifd_fld_proc] + ifd_msg_pos + sizeof(int) * 3 + sizeof(double)));
          vel_z = *((double*)(gv->ifd_send_msg[ifd_fld_proc] + ifd_msg_pos + sizeof(int) * 3 + sizeof(double) * 2));

          t1 = Timer::get_cur_time();
          t_search += t1 - t0;

#ifdef IFD_FLUID2FIBER_DUMP
          fprintf(oFile, "velocity at (%2d,%2d):(%2d,%2d,%2d): %.24f,%.24f,%.24f\n",
                      i, j, X, Y, Z, vel_x, vel_y, vel_z);
#endif

          /*notice the difference for spreading and for interpolation here, cf is used for spreading*/
          /*factors in nondimensionalizing h^(-3) and dxdydz cancel each other*/
          s1 += vel_x * (1.0+cos(c_x*rx)) * (1.0+cos(c_y*ry)) * (1.0+cos(c_z*rz)) * c_64;
          s2 += vel_y * (1.0+cos(c_x*rx)) * (1.0+cos(c_y*ry)) * (1.0+cos(c_z*rz)) * c_64;
          s3 += vel_z * (1.0+cos(c_x*rx)) * (1.0+cos(c_y*ry)) * (1.0+cos(c_z*rz)) * c_64;

        }//for kk ends

        fibernode->x += gv->dt * s1;
        fibernode->y += gv->dt * s2;
        fibernode->z += gv->dt * s3;

      } //for fiberclmn ends
    }//if fiber2thread ends
  }//for fiber row ends

#ifdef IFD_FLUID2FIBER_DUMP
  fclose(oFile);
#endif

#if 0
  if(tid == 0){
  	printf("Fiber%dtid%d prepare to final update location, t_search=%f\n", my_rank, tid, t_search);
  }

  for (i = 0; i < total_fibers_row; ++i){
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
      for (j = 0; j < total_fibers_clmn; ++j){
        fiberarray[i].nodes[j].x += gv->dt * s1;
      	fiberarray[i].nodes[j].y += gv->dt * s2;
      	fiberarray[i].nodes[j].z += gv->dt * s3;
      } //for fiberclmn ends
    }//if fiber2thread ends
  }//for fiber row ends
#endif

  // wait until all fiber threads complete computation
  pthread_barrier_wait(&(gv->barr));

#ifdef IFD_VEL_PERF
  t3 = Timer::get_cur_time();
  printf("Fiber%dtid%d: T_get_Sprd_Vel=%f\n", gv->taskid, tid, t3 - t2);
  fflush(stdout);
#endif

  //reset
  if (tid == 0){
    gv->num_influenced_proc = 0;
    for (int i = 0; i < num_fluid_tasks; i++){
      // set ifd_last_pos to default value 0
      gv->ifd_last_pos[i] = 0;
      msg_pos.clear();

      for (j = 0; j < total_threads; j++){
        gv->ifd_last_pos_proc_thd[i][j] = 0;     // Initialize gv->ifd_fluid_thread_last_pos

        Ifdmap_proc_thd[i][j].clear(); //remove all elements from the map container
        // std::cout << Ifdmap_proc_thd[i][j].size() << "\n";

        // free gv->ifd_fluid_thread_msg[i][j]

      }
    }
  }

#ifdef DEBUG_PRINT
  printf("----- Fiber%dtid%d: fiber_get_SpreadVelocity Exit! -----\n", my_rank, tid);
#endif //DEBUG_PRINT
}
