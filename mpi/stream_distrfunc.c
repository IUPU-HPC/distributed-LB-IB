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
#include "lb.h"

// dir = aggr_stream_dir
void insert_msg (GV gv, int tid, int nextX, int nextY, int nextZ,
                        int dir, int toTid, int iPop, double df1_tosend){
  int my_rank = gv->taskid;

  int dim_x = gv->fluid_grid->x_dim;
  int dim_y = gv->fluid_grid->y_dim;
  int dim_z = gv->fluid_grid->z_dim;

  pthread_mutex_lock(&gv->lock_stream_thd_msg[dir][toTid]);

  *((int*)(gv->stream_thd_msg[dir][toTid] + gv->stream_thd_last_pos[dir][toTid]))                      = nextX;
  *((int*)(gv->stream_thd_msg[dir][toTid] + gv->stream_thd_last_pos[dir][toTid] + sizeof(int)))        = nextY;
  *((int*)(gv->stream_thd_msg[dir][toTid] + gv->stream_thd_last_pos[dir][toTid] + sizeof(int) * 2))    = nextZ;
  *((int*)(gv->stream_thd_msg[dir][toTid] + gv->stream_thd_last_pos[dir][toTid] + sizeof(int) * 3))    = iPop;
  *((double*)(gv->stream_thd_msg[dir][toTid] + gv->stream_thd_last_pos[dir][toTid] + sizeof(int) * 4)) = df1_tosend;

#if 0
  assertf(nextX < dim_x, "insert_msg_ERROR: nextX=%d, position=%d", nextX, gv->stream_thd_last_pos[dir][toTid]);
  assertf(nextY < dim_y, "insert_msg_ERROR: nextY=%d, position=%d", nextY, gv->stream_thd_last_pos[dir][toTid]);
  assertf(nextZ < dim_z, "insert_msg_ERROR: nextZ=%d, position=%d", nextZ, gv->stream_thd_last_pos[dir][toTid]);
  assertf(iPop < 19, "insert_msg_ERROR: iPop=%d, position=%d", iPop, gv->stream_thd_last_pos[dir][toTid]);
  assert(gv->stream_thd_last_pos[dir][toTid] < gv->stream_recv_max_bufsize);

  printf("Fluid%dtid%d: insert_msg on Sdir %2d, next(%d,%d,%d)[%2d]=%f, gv->stream_thd_last_pos[%d][%d]=%d\n",
              my_rank, tid, dir, nextX, nextY, nextZ, iPop, df1_tosend, dir, toTid, gv->stream_thd_last_pos[dir][toTid]);
#endif

  gv->stream_thd_last_pos[dir][toTid] += sizeof(int) * 4 + sizeof(double);

  pthread_mutex_unlock(&gv->lock_stream_thd_msg[dir][toTid]);
}

void check_msg(GV gv, int tid, char* stream_msg, int sendcnt, int dest, int dir){

  int my_rank = gv->taskid;
  int BI, BJ, BK, li, lj, lk, iPop;
  double df1_tosend;
  int toTid;
  int toProc;
  int position = 0;
  int cube_size = gv->cube_size;

  int dim_x = gv->fluid_grid->x_dim;
  int dim_y = gv->fluid_grid->y_dim;
  int dim_z = gv->fluid_grid->z_dim;

  while (position < sendcnt){

    int X = *((int*)(stream_msg + position));
    int Y = *((int*)(stream_msg + position + sizeof(int)));
    int Z = *((int*)(stream_msg + position + sizeof(int)* 2));
    iPop  = *((int*)(stream_msg + position + sizeof(int)* 3));
    df1_tosend = *((double*)(stream_msg + position + sizeof(int)* 4));

    assertf(X < dim_x, "check_msg_ERROR: nextX=%d, position=%d", X, gv->stream_last_pos[dir]);
    assertf(Y < dim_y, "check_msg_ERROR: nextY=%d, position=%d", Y, gv->stream_last_pos[dir]);
    assertf(Z < dim_z, "Fluid%d check_msg_ERROR: nextZ=%d, position=%d, dest=%d, dir=%d, dim_z=%d, sendcnt=%d",
      my_rank, Z, gv->stream_last_pos[dir], dest, dir, dim_z, sendcnt);
    assertf(iPop < 19, "check_msg_ERROR: iPop=%d, position=%d", iPop, gv->stream_last_pos[dir]);

    BI = X / cube_size;
    BJ = Y / cube_size;
    BK = Z / cube_size;

    li = X % cube_size;
    lj = Y % cube_size;
    lk = Z % cube_size;

    toTid = cube2thread_and_task(BI, BJ, BK, gv, &toProc);//toTid is thread id in the fluid machine

    // Need check my_rank == fluid_owner_mac?
    // assert(my_rank == toProc);
    assertf(dest == toProc, "Send_ERROR: Fluid%dtid%d_check_msg: dir %2d, iPop=%2d, (dest,me,toProc)=(%d,%d,%d), \
next(X,Y,Z)=(%d,%d,%d), (BI,BJ,BK)=(%d,%d,%d), (li,lj,lk)=(%d,%d,%d), df1_tosend=%f, pos=%d\n",
      my_rank, tid, dir, iPop, dest, my_rank, toProc, X, Y, Z, BI, BJ, BK, li, lj, lk, df1_tosend, position);

#if 0
    printf("Fluid%dtid%d: check_msg on Sdir %2d, next(%d,%d,%d)[%2d]=%f, position=%d \n",
              my_rank, tid, dir, X, Y, Z, iPop, df1_tosend, position);
#endif
    position += sizeof(int) * 4 + sizeof(double);
  }
}

void get_df2_from_stream_msg(GV gv, int tid, int stream_msg_count, int dir, int src){
  int my_rank = gv->taskid;

  int position = 0;
  int BI, BJ, BK, li, lj, lk, iPop;
  double df1_tosend;
  int toTid;
  int toProc;
  Fluidnode  *nodes_df2;
  int cube_df2_idx, node_df2_idx;

  int cube_size = gv->cube_size;
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;

  int X, Y, Z;

  while (position < stream_msg_count){

    X = *((int*)(gv->stream_recv_buf + position));
    Y = *((int*)(gv->stream_recv_buf + position + sizeof(int)));
    Z = *((int*)(gv->stream_recv_buf + position + sizeof(int)* 2));
    iPop  = *((int*)(gv->stream_recv_buf + position + sizeof(int)* 3));
    df1_tosend = *((double*)(gv->stream_recv_buf + position + sizeof(int)* 4));

    BI = X / cube_size;
    li = X % cube_size;

    BJ = Y / cube_size;
    lj = Y % cube_size;

    BK = Z / cube_size;
    lk = Z % cube_size;

    toTid = cube2thread_and_task(BI, BJ, BK, gv, &toProc);//toTid is thread id in the fluid machine

    // printf("Fluid%dtid%d_Get_df2: dir %2d, iPop=%2d, (src,me,toProc)=(%d,%d,%d), (X,Y,Z)=(%ld,%ld,%ld), (BI,BJ,BK)=(%ld,%ld,%ld), (li,lj,lk)=(%ld,%ld,%ld), df1_tosend=%f, pos=%d\n",
    //   my_rank, tid, dir, iPop, src, my_rank, toProc, X, Y, Z, BI, BJ, BK, li, lj, lk, df1_tosend, position);
    // fflush(stdout);

    // Need check my_rank == fluid_owner_mac?
    // assert(my_rank == toProc);
    assertf(my_rank == toProc, "RECV_ERROR: Fluid%dtid%d_Get_df2: dir %2d, iPop=%2d, (src,me,toProc)=(%d,%d,%d), \
(X,Y,Z)=(%ld,%ld,%ld), (BI,BJ,BK)=(%ld,%ld,%ld), (li,lj,lk)=(%ld,%ld,%ld), df1_tosend=%f, pos=%d\n",
      my_rank, tid, dir, iPop, src, my_rank, toProc, X, Y, Z, BI, BJ, BK, li, lj, lk, df1_tosend, position);

    if (tid == toTid){ // since stop message alonhg with data is sent to all fluid machines, so N-1 machine will recv wrong cube
#if 0
      printf("Fluid%dtid%d: update_DF2 on Sdir %2d, next(%d,%d,%d)[%2d]=%f\n",
              my_rank, tid, dir, X, Y, Z, iPop, df1_tosend);
#endif
      cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
      node_df2_idx  = li * cube_size * cube_size + lj * cube_size + lk;//only local k will be endg point of prev cube
      nodes_df2[node_df2_idx].df2[iPop]  = df1_tosend;
    }

    position += sizeof(int) * 4 + sizeof(double);
  }

}

void streaming_on_direction(GV gv, int tid, int dir, int dest, int src){

  int stream_msg_recv_cnt;
  int my_rank = gv->taskid;
  MPI_Status status;

  /* 1. Sendrecv msg*/
  if (tid == 0){

#ifdef CHECK_STREAM
    printf("#Pre_SndRecv Fluid%dtid%d: dir=%d, Recv_from_src=%d, Send_to_dest=%d, sendcnt=stream_last_pos[%d]=%d\n",
      my_rank, tid, dir, src, dest, dir, gv->stream_last_pos[dir]);
    fflush(stdout);
#endif //CHECK_STREAM

#if 0
    check_msg(gv, tid, gv->stream_msg[dir], gv->stream_last_pos[dir], dest, dir);
#endif

    // if (gv->stream_last_pos[dir] > 0) {
      MPI_Sendrecv(gv->stream_msg[dir], gv->stream_last_pos[dir], MPI_CHAR,
        dest, dir, gv->stream_recv_buf, gv->stream_recv_max_bufsize,
        MPI_CHAR, src, dir,
        gv->cartcomm, &status);
      MPI_Get_count(&status, MPI_CHAR, &gv->stream_msg_recv_cnt);

      // reset gv->stream_last_pos[dir]
      gv->stream_last_pos[dir] = 0;

#ifdef CHECK_STREAM
    printf("#Get_SndRecv Fluid%dtid%d: dir=%d, Recv_from_src=%d with gv->stream_msg_recv_cnt=%d\n",
      my_rank, tid, dir, src, gv->stream_msg_recv_cnt);
#endif
    // }
  }

  pthread_barrier_wait(&(gv->barr));

  /* 2. examine recv data*/
  //all threads in each fluid process start working on stream_recv_buf
#if 0
  printf("Fluid%dtid%d: Pass MPI_Sendrecv dir=%d, Recv_from_src=%d with gv->stream_msg_recv_cnt=%d\n",
    my_rank, tid, dir, src, gv->stream_msg_recv_cnt);
  fflush(stdout);
#endif

  // if (gv->stream_msg_recv_cnt > 0){
    get_df2_from_stream_msg(gv, tid, gv->stream_msg_recv_cnt, dir, src);

#ifdef CHECK_STREAM
  printf("#Pass_get_df2 Fluid%dtid%d: from_stream_msg src=%d, dir=%d, recv_cnt=%d\n",
    my_rank, tid, src, dir, gv->stream_msg_recv_cnt);
  fflush(stdout);
#endif //CHECK_STREAM
  // }

  pthread_barrier_wait(&(gv->barr));

  /* 3. reset gv->stream_msg_recv_cnt to 0*/
  if (tid == 0) {
    gv->stream_msg_recv_cnt = 0;
  }
}

int findDir(GV gv, int toProc){
  for(int aggr_stream_dir = 1; aggr_stream_dir < 19; aggr_stream_dir++){
    if (toProc == gv->streamDest[aggr_stream_dir])
      return aggr_stream_dir;
  }
}

void  stream_distrfunc(LV lv){//df2
  /*18 spaces for 18 fluid nodes(neighbours). Whn 18 neighbours update current node  they write in their own space. so no lock required */
  /*For MPi additional boundary conditions Refer Thesis figure*/
  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

#ifdef DEBUG_PRINT
  // printf("****************Inside stream_distfunc()*******at time %d\n", gv->time);
#endif //DEBUG_PRINT

  Fluidgrid  *fluidgrid;
  Fluidnode  *nodes_df1, *nodes_df2;

  long cube_df1_idx, cube_df2_idx;
  int node_df1_idx, node_df2_idx, cube_size;
  int starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;

  int BI, BJ, BK, li, lj, lk, X, Y, Z, iPop, aggr_stream_dir;
  int nextX, nextY, nextZ, nextBI, nextBJ, nextBK, nextli, nextlj, nextlk;

  // int P[3], cubes_per_task[3], T[3], cubes_per_thd[3], t[3], toTid, total_threads;
  int Px, Py, Pz, cubes_per_task_x, cubes_per_task_y, cubes_per_task_z,
      tx, ty, tz, cubes_per_thd_x, cubes_per_thd_y, cubes_per_thd_z,
      ti, tj, tk, toTid, total_threads;

  // P[0] = gv->num_fluid_task_x;
  // P[1] = gv->num_fluid_task_y;
  // P[2] = gv->num_fluid_task_z;

  // T[0] = gv->tx;
  // T[1] = gv->ty;
  // T[2] = gv->tz;

  Px = gv->num_fluid_task_x;
  Py = gv->num_fluid_task_y;
  Pz = gv->num_fluid_task_z;

  tx = gv->tx;
  ty = gv->ty;
  tz = gv->tz;

  total_threads = gv->threads_per_task;

  int rank_x, rank_y, rank_z;
  int Pyz = Py * Pz;

  fluidgrid = gv->fluid_grid;
  cube_size = gv->cube_size;
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;

  // cubes_per_task[0] = num_cubes_x / P[0]; // along x: how many cubes in each process
  // cubes_per_task[1] = num_cubes_y / P[1]; // along y: how many cubes in each process
  // cubes_per_task[2] = num_cubes_z / P[2]; // along y: how many cubes in each process

  // cubes_per_thd[0] = cubes_per_task[0] / T[0];
  // cubes_per_thd[1] = cubes_per_task[1] / T[1];
  // cubes_per_thd[2] = cubes_per_task[2] / T[2];

  cubes_per_task_x = num_cubes_x / Px; // along x: how many cubes in each process
  cubes_per_task_y = num_cubes_y / Py; // along y: how many cubes in each process
  cubes_per_task_z = num_cubes_z / Pz; // along y: how many cubes in each process

  cubes_per_thd_x = cubes_per_task_x / tx;
  cubes_per_thd_y = cubes_per_task_y / ty;
  cubes_per_thd_z = cubes_per_task_z / tz;
#if 0
  printf("cubes_per_task (%d, %d, %d) cubes_per_thd(%d, %d, %d)\n",
    cubes_per_task_x, cubes_per_task_y, cubes_per_task_z,
    cubes_per_thd_x, cubes_per_thd_y, cubes_per_thd_z);
#endif
  starting_x = starting_y = starting_z = 0;
  stopping_x = stopping_y = stopping_z = cube_size - 1;

  int fromProc, toProc, my_rank, send_tid;
  my_rank = gv->taskid;

  double t0 = 0, t1 = 0, t2 = 0, t3 = 0;
  double t_insert = 0, t_self = 0;


#ifdef STREAM_PERF
  t2 = Timer::get_cur_time();
#endif
  // step1: prepare message
  for (BI = 0; BI < num_cubes_x; ++BI)
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
  for (BK = 0; BK < num_cubes_z; ++BK){

    send_tid = cube2thread_and_task(BI, BJ, BK, gv, &fromProc);

    if (my_rank == fromProc && send_tid == tid){

      cube_df1_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes_df1 = gv->fluid_grid->sub_fluid_grid[cube_df1_idx].nodes;
      starting_x = starting_y = starting_z = 0;
      stopping_x = stopping_y = stopping_z = cube_size - 1;

      if (BI == 0) starting_x = 2;
      if (BI == num_cubes_x - 1) stopping_x = cube_size - 3;

      if (BJ == 0) starting_y = 2;
      if (BJ == num_cubes_y - 1) stopping_y = cube_size - 3;

      if (BK == 0) starting_z = 2;
      if (BK == num_cubes_z - 1) stopping_z = cube_size - 3;

      for (li = starting_x; li <= stopping_x; ++li)
      for (lj = starting_y; lj <= stopping_y; ++lj)
      for (lk = starting_z; lk <= stopping_z; ++lk){

        node_df1_idx = li * cube_size * cube_size + lj * cube_size + lk;

        X = BI * cube_size + li;
        Y = BJ * cube_size + lj;
        Z = BK * cube_size + lk;

        for (iPop = 0; iPop < 19; ++iPop){
          nextX = X + c[iPop][0];
          nextY = Y + c[iPop][1];
          nextZ = Z + c[iPop][2];

          nextBI = nextX / cube_size;
          nextli = nextX % cube_size;

          nextBJ = nextY / cube_size;
          nextlj = nextY % cube_size;

          nextBK = nextZ / cube_size;
          nextlk = nextZ % cube_size;

          // toProc = global2task(nextX, nextY, nextZ, gv);
          // toProc = my_rank; // Used in 1 fluid node test case to improve performance
          // toProc = (nextBI / num_cubes_per_proc_x) * Py * Pz
          //         + (nextBJ / num_cubes_per_proc_y) * Pz
          //         + nextBK / num_cubes_per_proc_z;
          // rank_x = nextBI / cubes_per_task[0];
          // rank_y = nextBJ / cubes_per_task[1];
          // rank_z = nextBK / cubes_per_task[2];
          rank_x = nextBI / cubes_per_task_x;
          rank_y = nextBJ / cubes_per_task_y;
          rank_z = nextBK / cubes_per_task_z;
          toProc = rank_x * Pyz + rank_y * Pz + rank_z;

          cube_df2_idx = nextBI * num_cubes_y * num_cubes_z + nextBJ * num_cubes_z + nextBK;
          nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = nextli * cube_size * cube_size + nextlj * cube_size + nextlk;

#if 0
          if(my_rank==0)
            printf("Fluid%dtid%d: iPop %2d, (X,Y,Z)=(%ld,%ld,%ld), (BI,BJ,BK)=(%ld,%ld,%ld), (li,lj,lk)=(%ld,%ld,%ld), df1_tosend=%f, \
next(X,Y,Z)=(%ld,%ld,%ld), next(BI,BJ,BK)=(%ld,%ld,%ld), next(li,lj,lk)=(%ld,%ld,%ld), node_df2_idx=%ld\n",
              my_rank, tid, iPop, X, Y, Z, BI, BJ, BK, li, lj, lk, nodes_df1[node_df1_idx].df1[iPop],
              nextX, nextY, nextZ, nextBI, nextBJ, nextBK, nextli, nextlj, nextlk, node_df2_idx);
#endif

          if(my_rank == toProc){
#if 0
            t0 = Timer::get_cur_time();
#endif
            nodes_df2[node_df2_idx].df2[iPop] = nodes_df1[node_df1_idx].df1[iPop];
#if 0
            t1 = Timer::get_cur_time();
            t_self += t1 - t0;
#endif
          }
          else{
            // t[0] = (nextBI % cubes_per_task[0]) / cubes_per_thd[0]; //i
            // t[1] = (nextBJ % cubes_per_task[1]) / cubes_per_thd[1]; //j
            // t[2] = (nextBK % cubes_per_task[2]) / cubes_per_thd[2]; //k

            ti = (nextBI % cubes_per_task_x) / cubes_per_thd_x; //i
            tj = (nextBJ % cubes_per_task_y) / cubes_per_thd_y; //j
            tk = (nextBK % cubes_per_task_z) / cubes_per_thd_z; //k

            //tid  = i * ty * tz + j * tz + k;
            // toTid = t[0] * T[1] * T[2] + t[1] * T[2] + t[2];
            toTid = ti * ty * tz + tj * tz + tk;
#if 0
            printf("Fluid%dTid%d: next(%d,%d,%d) nextB(%d,%d,%d) t(%d,%d,%d) toTid=%d\n",
              my_rank, tid,
              nextX, nextY, nextZ,
              nextBI, nextBJ, nextBK,
              ti, tj, tk,
              toTid);
#endif
            aggr_stream_dir = findDir(gv, toProc); //find the direction of streaming message
            assert(aggr_stream_dir > 0);
#if 0
            printf("Fluid%dtid%d: iPop%2d, Sdir %2d, (X,Y,Z)=(%ld,%ld,%ld), df1_tosend=%f, next(X,Y,Z)=(%ld,%ld,%ld)\n",
              my_rank, tid, iPop, dir, X, Y, Z, nodes_df1[node_df1_idx].df1[iPop],
              nextX, nextY, nextZ);
#endif

#if 0
            t0 = Timer::get_cur_time();
#endif
            insert_msg(gv, tid, nextX, nextY, nextZ,
                      aggr_stream_dir, toTid, iPop, nodes_df1[node_df1_idx].df1[iPop]);

#if 0
            t1 = Timer::get_cur_time();
            t_insert += t1 - t0;
#endif
          }

        }

        // //ksi==0 starts
        // //if (my_rank == fromProc){   //only the current machine needs to update its own node.
        //   nodes_df2[node_df2_idx].df2[0] = nodes_df1[node_df1_idx].df1[0];
        // //}
        // //ksi==0 ends

        // //ksi==1 starts
        // //surfaces[i+1].nodes[j*dim_z +k].df2[1]      = surfaces[i].nodes[j*dim_z +k].df1[1];
        // printf("Fluid task%dtid%d: Start Prepare streaming on Dir 1\n", my_rank, tid);
        // fflush(stdout);
        // if ((li + 1) <= cube_size - 1){//i.e the the computation domain is still in the current cube
        //     //same machine, same cube
        //     cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        //     nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
        //     node_df2_idx = (li + 1)*cube_size*cube_size + lj*cube_size + lk;
        //     nodes_df2[node_df2_idx].df2[1] = nodes_df1[node_df1_idx].df1[1];
        // }
        // else{
        //   cube2thread_and_task(BI + 1, BJ, BK, gv, &toProc);//finding the recv machine rank

        //   if(fromProc == toProc){ // streaming in different cube but in same machine

        //     cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        //     nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
        //     node_df2_idx = (0)*cube_size*cube_size + lj*cube_size + lk; // local i of next cube from 0
        //     nodes_df2[node_df2_idx].df2[1] = nodes_df1[node_df1_idx].df1[1];
        //   }
        //   else{ // streaming in different cube but in different machine

        //     df1_tosend = nodes_df1[node_df1_idx].df1[1];

        //     insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
        //   }
        // }
        // //ksi =1 ends

      }//lk
    }// if cube2thread ends
  }//BK

#ifdef STREAM_PERF
  t3 = Timer::get_cur_time();
#endif

#if 0
  printf("Fluid task%d: Pass Prepare streaming msg\n", my_rank);
  fflush(stdout);
#endif //CHECK_STREAM

  pthread_barrier_wait(&(gv->barr));

  // step2: memcpy message, iPop = aggr_stream_dir
  if(tid==0){
#ifdef STREAM_PERF
    t0 = Timer::get_cur_time();
#endif
    for (iPop = 1; iPop < 19; iPop++){
      for (toTid = 0; toTid < total_threads; ++toTid){
        char* src = gv->stream_thd_msg[iPop][toTid];
        int size = gv->stream_thd_last_pos[iPop][toTid];

        if (size > 0){
          int last_pos = gv->stream_last_pos[iPop];

          printf("Fluid%dTid%d: memcpy (%d, %d) --> (%d), size=%d\n",
            my_rank, tid, iPop, toTid, iPop, size);

          char* dest = (char*)(gv->stream_msg[iPop] + last_pos);
          memcpy(dest, src, size);

          gv->stream_last_pos[iPop] += size;
        }

        //reset
        gv->stream_thd_last_pos[iPop][toTid] = 0;
      }
    }
#ifdef STREAM_PERF
    t1 = Timer::get_cur_time();
    printf("Fluid%dTid%d: Pass stream accumulate copy, t_copy=%f\n", my_rank,tid, t1-t0);
#endif
  }

  pthread_barrier_wait(&(gv->barr));

  // step3: send message, iPop = aggr_stream_dir
  t0 = Timer::get_cur_time();

  for(iPop = 1; iPop < 19; iPop++){
#ifdef CHECK_STREAM
    if (tid == 0){
      printf("-COUNT- Fluid%dtid%d: start streaming msg iPop=%d, dest=%d, src=%d, sendcnt=%d\n",
        my_rank, tid, iPop, gv->streamDest[iPop], gv->streamSrc[iPop], gv->stream_last_pos[iPop]);
      fflush(stdout);
    }
#endif //CHECK_STREAM

    streaming_on_direction(gv, tid, iPop, gv->streamDest[iPop], gv->streamSrc[iPop]);

    // if(tid==0)
    //   MPI_Barrier(MPI_COMM_WORLD);
  }

  t1 = Timer::get_cur_time();

#ifdef STREAM_PERF
  printf("Fluid%dtid%d: T_prepare_stream=%f, T_self=%f, T_insert=%f, T_stream=%f\n",
    my_rank, tid, t3-t2, t_self, t_insert, t1-t0);
#endif
}
