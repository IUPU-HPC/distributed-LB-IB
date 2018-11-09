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
#include <assert.h>

extern int c[19][3];
extern double t[19];

void insert_msg(GV gv, long nextX, long nextY, long nextZ, long dir, double df1_tosend){

  pthread_mutex_lock(&gv->lock_stream_msg[dir]);

  *((long*)(gv->stream_msg[dir] + gv->stream_last_pos[dir]))                     = nextX;
  *((long*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(long)))      = nextY;
  *((long*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(long) * 2))  = nextZ;
  *((long*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(long) * 3))  = dir;
  *((double*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(long) * 4))  = df1_tosend;

  gv->stream_last_pos[dir] += sizeof(long) * 4 + sizeof(double);
  // printf("Fluid_Proc%d: gv->stream_last_pos[%ld]=%d\n", gv->taskid, dir, gv->stream_last_pos[dir]);
  // fflush(stdout);

  pthread_mutex_unlock(&gv->lock_stream_msg[dir]);
}

void get_df2_from_stream_msg(LV lv, int stream_msg_count){
  GV gv = lv->gv;
  int tid = lv->tid;
  int my_rank = gv->taskid;

  int position = 0;
  long BI, BJ, BK, li, lj, lk, dir;
  double df1_tosend;
  int toTid;
  int toProc;
  Fluidnode  *nodes_df2;
  long cube_df2_idx, node_df2_idx;

  int cube_size = gv->cube_size;
  long num_cubes_x = gv->fluid_grid->num_cubes_x;
  long num_cubes_y = gv->fluid_grid->num_cubes_y;
  long num_cubes_z = gv->fluid_grid->num_cubes_z;

  while (position < stream_msg_count){

	  long X = *((long*)(gv->ifd_recv_buf + position));
	  long Y = *((long*)(gv->ifd_recv_buf + position + sizeof(long)));
	  long Z = *((long*)(gv->ifd_recv_buf + position + sizeof(long)* 2));  
    dir = *((long*)(gv->stream_recv_buf + position + sizeof(long)* 3));
    df1_tosend = *((double*)(gv->stream_recv_buf + position + sizeof(long)* 4));

    position += sizeof(long) * 4 + sizeof(double);

    BI = X / cube_size;
    BJ = Y / cube_size;
    BK = Z / cube_size;

    li = X % cube_size;
    lj = Y % cube_size;
    lk = Z % cube_size;

    toTid = cube2thread_and_task(BI, BJ, BK, gv, &toProc);//toTid is thread id in the fluid machine

    printf("Fluid%dtid%d: (my_rank,toProc)=(%d, %d), (X,Y,Z)=(%ld, %ld, %ld), (BI,BJ,BK)=(%ld, %ld, %ld), (li, lj, lk)=(%ld, %ld, %ld)\n", 
      my_rank, tid, my_rank, toProc, X, Y, Z, BI, BJ, BK, li, lj, lk);
    fflush(stdout);
	
    // Need check my_rank == fluid_owner_mac?
	  assert(my_rank == toProc);
	
    if (tid == toTid){ // since stop message alonhg with data is sent to all fluid machines, so N-1 machine will recv wrong cube
      cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
      node_df2_idx  = li * cube_size * cube_size + lj * cube_size + lk;//only local k will be endg point of prev cube
      nodes_df2[node_df2_idx].df2[dir]  = df1_tosend;
    }
  }

  printf("Fluid%d: Pass get df2 from dir=%d stream_msg\n", gv->taskid, dir);
  fflush(stdout);
}

void streaming_on_direction(GV gv, LV lv, int dir, int dest, int src){
  int tid = lv->tid;
  int stream_msg_recv_cnt;
  int my_rank = gv->taskid;
  MPI_Status status;

  if (tid == 0){
    printf("Fluid%d: Prepare to MPI_Sendrecv to dest=%d gv->stream_last_pos[%d]=%d\n", my_rank, dest, dir, gv->stream_last_pos[dir]);
    fflush(stdout);

    MPI_Sendrecv(gv->stream_msg[dir], gv->stream_last_pos[dir], MPI_CHAR,
      dest, dir, gv->stream_recv_buf, gv->stream_recv_max_bufsize,
      MPI_CHAR, src, dir,
      gv->cartcomm, &status);
  }

  pthread_barrier_wait(&(gv->barr));

  //all threads in each fluid machine start working on stream_recv_buf
  MPI_Get_count(&status, MPI_CHAR, &stream_msg_recv_cnt);

  printf("Fluid%d: Pass MPI_Sendrecv dir=%d and start get df2 from stream_msg_recv_cnt%d\n", my_rank, dir, stream_msg_recv_cnt);
  fflush(stdout);

  // get_df2_from_stream_msg(lv, stream_msg_recv_cnt);

  pthread_barrier_wait(&(gv->barr));
}

int findDir(GV gv, int toProc){
  for(int iPop=1; iPop < 19; iPop++){
    if(toProc == gv->streamDest[iPop])
      return iPop;
  }
}

void  stream_distrfunc(GV gv, LV lv){//df2
  /*18 spaces for 18 fluid nodes(neighbours). Whn 18 neighbours update current node  they write in their own space. so no lock required */
  /*For MPi additional boundary conditions Refer Thesis figure*/
  int tid;
  // GV gv = lv->gv;
  tid = lv->tid;
#ifdef DEBUG_PRINT
  // printf("****************Inside stream_distfunc()*******at time %d\n", gv->time);
#endif //DEBUG_PRINT
  // int fiber_mac_rank = gv->num_macs - 1;

  Fluidgrid  *fluidgrid;
  Fluidnode  *nodes_df1, *nodes_df2;

  long cube_df1_idx, cube_df2_idx, node_df1_idx, node_df2_idx;
  int cube_size;
  int starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;

  long BI, BJ, BK, li, lj, lk, X, Y, Z;
  long nextBI, nextBJ, nextBK, nextli, nextlj, nextlk;

  int Px, Py, Pz;
  long num_cubes_per_proc_x, num_cubes_per_proc_y, num_cubes_per_proc_z;

  Px = gv->num_fluid_task_x;
  Py = gv->num_fluid_task_y;
  Pz = gv->num_fluid_task_z;

  fluidgrid = gv->fluid_grid;
  cube_size = gv->cube_size;
  long num_cubes_x = gv->fluid_grid->num_cubes_x;
  long num_cubes_y = gv->fluid_grid->num_cubes_y;
  long num_cubes_z = gv->fluid_grid->num_cubes_z;

  num_cubes_per_proc_x = num_cubes_x / Px; // along x: how many cubes in each process
  num_cubes_per_proc_y = num_cubes_y / Py; // along y: how many cubes in each process
  num_cubes_per_proc_z = num_cubes_z / Pz; // along y: how many cubes in each process

  starting_x = starting_y = starting_z = 0;
  stopping_x = stopping_y = stopping_z = cube_size - 1;

  int  fromProc, toProc, my_rank, send_tid;
  double df1_tosend;
  my_rank = gv->taskid;


  //prepare message
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

        node_df1_idx = li*cube_size*cube_size + lj*cube_size + lk;

        long X = BI * cube_size + li;
        long Y = BJ * cube_size + lj;
        long Z = BK * cube_size + lk;

        for (int iPop = 0; iPop < 19; ++iPop){
          long nextX = X + c[iPop][0];
          long nextY = Y + c[iPop][1];
          long nextZ = Z + c[iPop][2];

          int toProc = global2task(nextX, nextY, nextZ, gv);

          nextBI = nextX / cube_size;
          nextBJ = nextY / cube_size;
          nextBK = nextZ / cube_size;

          nextli = nextX % cube_size;
          nextlj = nextY % cube_size;
          nextlk = nextZ % cube_size;

          cube_df2_idx = nextBI * num_cubes_y * num_cubes_z + nextBJ * num_cubes_z + nextBK;
          nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = nextli * cube_size * cube_size + nextlj * cube_size + nextlk;

          if(my_rank == toProc){
            nodes_df2[node_df2_idx].df2[iPop] = nodes_df1[node_df1_idx].df1[iPop];
          }
          else{
            int dir = findDir(gv, toProc);
            assert(dir > 0);
            insert_msg(gv, nextX, nextY, nextZ, dir, df1_tosend);
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

  printf("Fluid task%d: Pass Prepare streaming msg\n", my_rank);
  fflush(stdout);

  pthread_barrier_wait(&(gv->barr));
  // if(tid==0)
  //   MPI_Barrier(MPI_COMM_WORLD);

  //send message
  for(int iPop =1; iPop < 19; iPop++){
    if (tid == 0){
      printf("Fluid task%d: start streaming msg iPop=%d, dest=%d, src=%d\n", 
        my_rank, iPop, gv->streamDest[iPop], gv->streamSrc[iPop]);
      fflush(stdout);
    }

    streaming_on_direction(gv, lv, iPop, gv->streamDest[iPop], gv->streamSrc[iPop]);
  }


#ifdef DEBUG_PRINT
  printf("****************stream_distfunc() Exit*******\n");
#endif //DEBUG_PRINT
}
