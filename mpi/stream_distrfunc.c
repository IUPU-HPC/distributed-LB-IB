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

extern int c[19][3];
extern double t[19];

void insert_msg(LV lv, int nextX, int nextY, int nextZ, int dir, int iPop, double df1_tosend){
  GV gv = lv->gv;
  int tid = lv->tid;

  int dim_x = gv->fluid_grid->x_dim;
  int dim_y = gv->fluid_grid->y_dim;
  int dim_z = gv->fluid_grid->z_dim;

  pthread_mutex_lock(&gv->lock_stream_msg[dir]);

  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir]))                     = nextX;
  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int)))       = nextY;
  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int) * 2))   = nextZ;
  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int) * 3))   = iPop;
  *((double*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int) * 4)) = df1_tosend;

  gv->stream_last_pos[dir] += sizeof(int) * 4 + sizeof(double);

#if 0
  assertf(nextX < dim_x, "insert_msg_ERROR: nextX=%d, position=%d", nextX, gv->stream_last_pos[dir]);
  assertf(nextY < dim_y, "insert_msg_ERROR: nextY=%d, position=%d", nextY, gv->stream_last_pos[dir]);
  assertf(nextZ < dim_z, "insert_msg_ERROR: nextZ=%d, position=%d", nextZ, gv->stream_last_pos[dir]);
  assertf(iPop < 19, "insert_msg_ERROR: iPop=%d, position=%d", iPop, gv->stream_last_pos[dir]);
  assert(gv->stream_last_pos[dir] < gv->stream_recv_max_bufsize);


  if(gv->stream_last_pos[dir] >= 413664 && gv->stream_last_pos[dir] <= 413760){
    printf("Fluid%dtid%d_insert_msg: stream_dir %d, iPop=%2d, next(x,y,z)=(%ld,%ld,%ld), df1=%f, gv->stream_last_pos[%ld]=%d\n", 
      gv->taskid, tid, dir, iPop, nextX, nextY, nextZ, df1_tosend, dir, gv->stream_last_pos[dir]);
    fflush(stdout);
  }
#endif  
  pthread_mutex_unlock(&gv->lock_stream_msg[dir]);
}

void check_msg(LV lv, char* stream_msg, int sendcnt, int dest, int dir){
  GV gv = lv->gv;
  int tid = lv->tid;
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

    position += sizeof(int) * 4 + sizeof(double);

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
  }  
}

void get_df2_from_stream_msg(LV lv, int stream_msg_count, int dir, int src){
  GV gv = lv->gv;
  int tid = lv->tid;
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

  while (position < stream_msg_count){

    int X = *((int*)(gv->stream_recv_buf + position));
    int Y = *((int*)(gv->stream_recv_buf + position + sizeof(int)));
    int Z = *((int*)(gv->stream_recv_buf + position + sizeof(int)* 2));  
    iPop  = *((int*)(gv->stream_recv_buf + position + sizeof(int)* 3));
    df1_tosend = *((double*)(gv->stream_recv_buf + position + sizeof(int)* 4));

    position += sizeof(int) * 4 + sizeof(double);

    BI = X / cube_size;
    BJ = Y / cube_size;
    BK = Z / cube_size;

    li = X % cube_size;
    lj = Y % cube_size;
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
      cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
      node_df2_idx  = li * cube_size * cube_size + lj * cube_size + lk;//only local k will be endg point of prev cube
      nodes_df2[node_df2_idx].df2[iPop]  = df1_tosend;
    }
  }

}

void streaming_on_direction(LV lv, int dir, int dest, int src){
  GV gv = lv->gv;
  int tid = lv->tid;
  int stream_msg_recv_cnt;
  int my_rank = gv->taskid;
  MPI_Status status;

  if (tid == 0){
    printf("Fluid%d: Prepare to MPI_Sendrecv on Dir=%d, Recv_from_src=%d, Send_to_dest=%d, sendcnt=gv->stream_last_pos[%d]=%d\n", 
      my_rank, dir, src, dest, dir, gv->stream_last_pos[dir]);
    fflush(stdout);

#if 0
    check_msg(lv, gv->stream_msg[dir], gv->stream_last_pos[dir], dest, dir);
#endif

    MPI_Sendrecv(gv->stream_msg[dir], gv->stream_last_pos[dir], MPI_CHAR,
      dest, dir, gv->stream_recv_buf, gv->stream_recv_max_bufsize,
      MPI_CHAR, src, dir,
      gv->cartcomm, &status);
    MPI_Get_count(&status, MPI_CHAR, &gv->stream_msg_recv_cnt);

    // reset gv->stream_last_pos[dir]
    gv->stream_last_pos[dir] = 0;
  }

  //all threads in each fluid machine start working on stream_recv_buf
  pthread_barrier_wait(&(gv->barr));

  printf("Fluid%dtid%d: Pass MPI_Sendrecv dir=%d, Recv_from_src=%d with gv->stream_msg_recv_cnt=%d\n", 
    my_rank, tid, dir, src, gv->stream_msg_recv_cnt);
  fflush(stdout);

  get_df2_from_stream_msg(lv, gv->stream_msg_recv_cnt, dir, src);

  printf("Fluid%dtid%d: Pass get_df2_from_stream_msg from dir=%d stream_msg\n", 
    my_rank, tid, dir);
  fflush(stdout);

  pthread_barrier_wait(&(gv->barr));
}

int findDir(GV gv, int toProc){
  for(int streamdir = 1; streamdir < 19; streamdir++){
    if (toProc == gv->streamDest[streamdir])
      return streamdir;
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

  int cube_df1_idx, cube_df2_idx, node_df1_idx, node_df2_idx;
  int cube_size;
  int starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;

  int BI, BJ, BK, li, lj, lk, X, Y, Z;
  int nextBI, nextBJ, nextBK, nextli, nextlj, nextlk;

  int Px, Py, Pz;
  int num_cubes_per_proc_x, num_cubes_per_proc_y, num_cubes_per_proc_z;

  Px = gv->num_fluid_task_x;
  Py = gv->num_fluid_task_y;
  Pz = gv->num_fluid_task_z;

  fluidgrid = gv->fluid_grid;
  cube_size = gv->cube_size;
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;

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

        node_df1_idx = li * cube_size * cube_size + lj * cube_size + lk;

        int X = BI * cube_size + li;
        int Y = BJ * cube_size + lj;
        int Z = BK * cube_size + lk;

        for (int iPop = 0; iPop < 19; ++iPop){
          int nextX = X + c[iPop][0];
          int nextY = Y + c[iPop][1];
          int nextZ = Z + c[iPop][2];

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

#if 0
          if(my_rank==0)
            printf("Fluid%dtid%d: iPop %2d, (X,Y,Z)=(%ld,%ld,%ld), (BI,BJ,BK)=(%ld,%ld,%ld), (li,lj,lk)=(%ld,%ld,%ld), df1_tosend=%f, \
next(X,Y,Z)=(%ld,%ld,%ld), next(BI,BJ,BK)=(%ld,%ld,%ld), next(li,lj,lk)=(%ld,%ld,%ld), node_df2_idx=%ld\n", 
              my_rank, tid, iPop, X, Y, Z, BI, BJ, BK, li, lj, lk, nodes_df1[node_df1_idx].df1[iPop],
              nextX, nextY, nextZ, nextBI, nextBJ, nextBK, nextli, nextlj, nextlk, node_df2_idx);
#endif

          if(my_rank == toProc){
            nodes_df2[node_df2_idx].df2[iPop] = nodes_df1[node_df1_idx].df1[iPop];
          }
          else{
            int dir = findDir(gv, toProc); //find the direction of streaming message 
            assert(dir > 0);
#if 0
            assertf(nextX != 128, "Error! Fluid%dtid%d: iPop %2d, Sdir %2d, (X,Y,Z)=(%ld,%ld,%ld), (BI,BJ,BK)=(%ld,%ld,%ld), (li,lj,lk)=(%ld,%ld,%ld), df1_tosend=%f, \
next(X,Y,Z)=(%ld,%ld,%ld), next(BI,BJ,BK)=(%ld,%ld,%ld), next(li,lj,lk)=(%ld,%ld,%ld), node_df2_idx=%ld\n", 
              my_rank, tid, iPop, dir, X, Y, Z, BI, BJ, BK, li, lj, lk, nodes_df1[node_df1_idx].df1[iPop],
              nextX, nextY, nextZ, nextBI, nextBJ, nextBK, nextli, nextlj, nextlk, node_df2_idx);
#endif
            insert_msg(lv, nextX, nextY, nextZ, dir, iPop, df1_tosend);
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
  for(int iPop = 1; iPop < 19; iPop++){
    if (tid == 0){
      printf("Fluid task%d: start streaming msg iPop=%d, dest=%d, src=%d, sendcnt=%d\n", 
        my_rank, iPop, gv->streamDest[iPop], gv->streamSrc[iPop], gv->stream_last_pos[iPop]);
      fflush(stdout);
    }

    streaming_on_direction(lv, iPop, gv->streamDest[iPop], gv->streamSrc[iPop]);

    // if(tid==0)
    //   MPI_Barrier(MPI_COMM_WORLD);
  }

#if 0
  printf("****************stream_distfunc() Exit*******\n");
#endif
}
