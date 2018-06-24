/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

int check_cube_on_proc_line_boundary(GV gv, int BI, int BJ, int BK, int Px, int Py, int Pz){
  int num_cubes_x, num_cubes_y, num_cubes_z;
  int num_cubes_per_proc_x, num_cubes_per_proc_y, num_cubes_per_proc_z;

  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  num_cubes_per_proc_x = num_cubes_x / Px; // along x: how many cubes in each process
  num_cubes_per_proc_y = num_cubes_y / Py; // along y: how many cubes in each process
  num_cubes_per_proc_z = num_cubes_z / Pz; // along y: how many cubes in each process

  if( (BI%num_cubes_per_proc_x==0 || BI%num_cubes_per_proc_x==(num_cubes_per_proc_x-1))
    && (BJ%num_cubes_per_proc_y==0 || BI%num_cubes_per_proc_y==(num_cubes_per_proc_y-1))
    )
    return 1;
  else if( (BI%num_cubes_per_proc_x==0 || BI%num_cubes_per_proc_x==(num_cubes_per_proc_x-1))
    && (BK%num_cubes_per_proc_z==0 || BI%num_cubes_per_proc_z==(num_cubes_per_proc_z-1))
    )
    return 1;
  else if( (BJ%num_cubes_per_proc_y==0 || BI%num_cubes_per_proc_y==(num_cubes_per_proc_y-1))
    && (BK%num_cubes_per_proc_z==0 || BI%num_cubes_per_proc_z==(num_cubes_per_proc_z-1))
    )
    return 1;
  else
    return 0;
}

void insert_msg(GV gv, int BI, int BJ, int BK, int li, int lj, int lk, int dir, double df1_tosend){

  pthread_mutex_lock(&gv->lock_stream_msg[dir]);

  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir]))                   = BI;
  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int)))     = BJ;
  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int)* 2))  = BK;
  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int)* 3))  = li;
  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int)* 4))  = lj;
  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int)* 5))  = lk;
  *((int*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int)* 6))  = dir;
  *((double*)(gv->stream_msg[dir] + gv->stream_last_pos[dir] + sizeof(int)* 7))  = df1_tosend;

  gv->stream_last_pos[dir] += sizeof(int)* 7 + sizeof(double);
  printf("Fluid_mac%d: gv->stream_last_pos[%d]=%d\n", gv->my_rank, dir, gv->stream_last_pos[dir]);
  fflush(stdout);

  pthread_mutex_unlock(&gv->lock_stream_msg[dir]);
}

void get_df2_from_stream_msg(LV lv, int stream_msg_count){
  GV gv = lv->gv;
  int tid = lv->tid;

  int position = 0;
  int BI, BJ, BK, li, lj, lk, dir;
  double df1_tosend;
  int owner_tid;
  int fluid_owner_mac;
  Fluidnode  *nodes_df2;
  int cube_df2_idx, node_df2_idx;

  int cube_size = gv->cube_size;
  int num_cubes_x = gv->num_cubes_x;
  int num_cubes_y = gv->num_cubes_y;
  int num_cubes_z = gv->num_cubes_z;

  while (position < stream_msg_count){

    BI = *((int*)(gv->stream_recv_buf + position));
    BJ = *((int*)(gv->stream_recv_buf + position + sizeof(int)));
    BK = *((int*)(gv->stream_recv_buf + position + sizeof(int)* 2));
    li = *((int*)(gv->stream_recv_buf + position + sizeof(int)* 3));
    lj = *((int*)(gv->stream_recv_buf + position + sizeof(int)* 4));
    lk = *((int*)(gv->stream_recv_buf + position + sizeof(int)* 5));
    dir = *((int*)(gv->stream_recv_buf + position + sizeof(int)* 6));
    df1_tosend = *((double*)(gv->stream_recv_buf + position + sizeof(int)* 7));

    position += sizeof(int)* 7 + sizeof(double);

    owner_tid = cube2thread_and_machine(BI, BJ, BK, gv, &fluid_owner_mac);//owner_tid is thread id in the fluid machine
    // Need check my_rank == fluid_owner_mac?
    if (tid == owner_tid){// since stop message alonhg with data is sent to all fluid machines, so N-1 machine will recv wrong cube

      cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
      node_df2_idx = li*cube_size*cube_size +lj*cube_size +lk;//only local k will be endg point of prev cube
      nodes_df2[node_df2_idx].df2[dir]  = df1_tosend;

    }
  }

  printf("Fluid%d: Pass get df2 from dir=%d stream_msg\n", gv->my_rank, dir);
  fflush(stdout);
}

void streaming_on_direction(GV gv, LV lv, int dir, int dest){
  int stream_msg_recv_cnt;
  int my_rank = gv->my_rank;
  MPI_Status status;

  printf("Fluid%d: Prepare to MPI_Sendrecv to dest=%d gv->stream_last_pos[%d]=%d\n", my_rank, dest, dir, gv->stream_last_pos[dir]);
  fflush(stdout);

  //MPI_PROC_NULL
  MPI_Sendrecv(gv->stream_msg[dir], gv->stream_last_pos[dir], MPI_CHAR,
    dest, dir, gv->stream_recv_buf, gv->stream_recv_max_bufsize,
    MPI_CHAR, my_rank, dir,
    MPI_COMM_WORLD, &status);

  pthread_barrier_wait(&(gv->barr));

  //all threads in each fluid machine start working on stream_recv_buf
  MPI_Get_count(&status, MPI_CHAR, &stream_msg_recv_cnt);

  printf("Fluid%d: Pass MPI_Sendrecv dir=%d and start get df2 from stream_msg_recv_cnt%d\n", my_rank, dir, stream_msg_recv_cnt);
  fflush(stdout);

  get_df2_from_stream_msg(lv, stream_msg_recv_cnt);

  pthread_barrier_wait(&(gv->barr));

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

  int        cube_df1_idx, cube_df2_idx, node_df1_idx, node_df2_idx, cube_size, num_cubes_x, num_cubes_y, num_cubes_z;
  int        starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;

  int        BI, BJ, BK, li, lj, lk;
  int        Px, Py, Pz;
  int        num_cubes_per_proc_x, num_cubes_per_proc_y, num_cubes_per_proc_z;


  Px = gv->Px;
  Py = gv->Py;
  Pz = gv->Pz;

  fluidgrid = gv->fluid_grid;
  cube_size = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  num_cubes_per_proc_x = num_cubes_x / Px; // along x: how many cubes in each process
  num_cubes_per_proc_y = num_cubes_y / Py; // along y: how many cubes in each process
  num_cubes_per_proc_z = num_cubes_z / Pz; // along y: how many cubes in each process


  starting_x = starting_y = starting_z = 0;
  stopping_x = stopping_y = stopping_z = cube_size - 1;

  int  temp_mac_rank_sender, temp_mac_rank_recv, my_rank, send_tid;
  double df1_tosend;
  my_rank = gv->my_rank;


  //prepare message
  for (BI = 0; BI < num_cubes_x; ++BI)
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
  for (BK = 0; BK < num_cubes_z; ++BK){

    send_tid = cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank_sender);

    if (my_rank == temp_mac_rank_sender && send_tid == tid){

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

        node_df1_idx = li*cube_size*cube_size + lj*cube_size + lk;//every machine knows node_df1_idx, cubes_df2, nodes_df2 and
        cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
        node_df2_idx = li*cube_size*cube_size + lj*cube_size + lk;
        //ksi==0 starts
        //if (my_rank == temp_mac_rank_sender){   //only the current machine needs to update its own node.
          nodes_df2[node_df2_idx].df2[0] = nodes_df1[node_df1_idx].df1[0];
        //}
        //ksi==0 ends

        //ksi==1 starts
        //surfaces[i+1].nodes[j*dim_z +k].df2[1]      = surfaces[i].nodes[j*dim_z +k].df1[1];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 1\n", my_rank, tid);
        fflush(stdout);
        if ((li + 1) <= cube_size - 1){//i.e the the computation domain is still in the current cube
            //same machine, same cube
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li + 1)*cube_size*cube_size + lj*cube_size + lk;
            nodes_df2[node_df2_idx].df2[1] = nodes_df1[node_df1_idx].df1[1];
        }
        else{
          cube2thread_and_machine(BI + 1, BJ, BK, gv, &temp_mac_rank_recv);//finding the recv machine rank

          if(temp_mac_rank_sender == temp_mac_rank_recv){ // streaming in different cube but in same machine

            cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (0)*cube_size*cube_size + lj*cube_size + lk; // local i of next cube from 0
            nodes_df2[node_df2_idx].df2[1] = nodes_df1[node_df1_idx].df1[1];
          }
          else{ // streaming in different cube but in different machine

            df1_tosend = nodes_df1[node_df1_idx].df1[1];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
          }
        }
        //ksi =1 ends


        //ksi==2 starts
        //surfaces[i-1].nodes[j*dim_z +k].df2[2]      = surfaces[i].nodes[j*dim_z +k].df1[2];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 2\n", my_rank, tid);
        fflush(stdout);
        if ((li - 1) >= 0){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li - 1)*cube_size*cube_size + lj*cube_size + lk;
            nodes_df2[node_df2_idx].df2[2] = nodes_df1[node_df1_idx].df1[2];

        }
        else{// computation domain is in the prev cube and proabaly prev machine only BI changes

          cube2thread_and_machine(BI - 1, BJ, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI - 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (cube_size - 1)*cube_size*cube_size + lj*cube_size + lk;//local i of prev cube will be ending point.
            nodes_df2[node_df2_idx].df2[2] = nodes_df1[node_df1_idx].df1[2];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[2];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 2, df1_tosend);
          }

        }
        //ksi==2 ends


        //ksi==3 starts
        //surfaces[i].nodes[j*dim_z +k+1].df2[3]      = surfaces[i].nodes[j*dim_z +k].df1[3];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 3\n", my_rank, tid);
        fflush(stdout);
        if ((lk + 1) <= cube_size - 1){//i.e the the computation domain is still in the current cube
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + lj*cube_size + lk + 1;
            nodes_df2[node_df2_idx].df2[3] = nodes_df1[node_df1_idx].df1[3];
        }
        else{//computation domain is in the next cube and probably next machine only BK changes
          cube2thread_and_machine(BI, BJ, BK + 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK + 1;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + lj*cube_size + 0;//local k of next cube will start from 0
            nodes_df2[node_df2_idx].df2[3] = nodes_df1[node_df1_idx].df1[3];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[3];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 3, df1_tosend);
          }

        }
        //ksi==3 ends


        //ksi==4 starts
        //surfaces[i].nodes[j*dim_z +k-1].df2[4]      = surfaces[i].nodes[j*dim_z +k].df1[4];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 4\n", my_rank, tid);
        fflush(stdout);
        if ((lk - 1) >= 0){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + lj*cube_size + lk - 1;
            nodes_df2[node_df2_idx].df2[4] = nodes_df1[node_df1_idx].df1[4];

        }
        else{// computation domain is in the prev cube and probably prev machine only BK changes
          cube2thread_and_machine(BI, BJ, BK - 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK - 1;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + lj*cube_size + cube_size - 1;//local k of prev cube will be the end point
            nodes_df2[node_df2_idx].df2[4] = nodes_df1[node_df1_idx].df1[4];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[4];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 4, df1_tosend);
          }

        }
        //ksi==4 ends

        //ksi==5 starts
        //surfaces[i].nodes[(j-1)*dim_z +k].df2[5]    = surfaces[i].nodes[j*dim_z +k].df1[5];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 5\n", my_rank, tid);
        fflush(stdout);
        if ((lj - 1) >= 0){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + (lj - 1)*cube_size + lk;
            nodes_df2[node_df2_idx].df2[5] = nodes_df1[node_df1_idx].df1[5];
        }
        else{// computation domain is in the prev cube and probably prev machine only BJ changes

          cube2thread_and_machine(BI, BJ - 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + (BJ - 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + (cube_size - 1)*cube_size + lk;//local j of prev cube will be the end point
            nodes_df2[node_df2_idx].df2[5] = nodes_df1[node_df1_idx].df1[5];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[5];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 5, df1_tosend);
          }

        }
        //ksi==5 ends

        //ksi==6 starts //surfaces[i].nodes[(j+1)*dim_z +k].df2[6]    = surfaces[i].nodes[j*dim_z +k].df1[6];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 6\n", my_rank, tid);
        fflush(stdout);
        if ((lj + 1) <= cube_size - 1){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + (lj + 1)*cube_size + lk;
            nodes_df2[node_df2_idx].df2[6] = nodes_df1[node_df1_idx].df1[6];

        }
        else{//computation domain is in the next cube and probably next machine only BJ changes

          cube2thread_and_machine(BI, BJ + 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + (BJ + 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + 0 * cube_size + lk;//local j of next cube will start from 0
            nodes_df2[node_df2_idx].df2[6] = nodes_df1[node_df1_idx].df1[6];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[6];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 6, df1_tosend);
          }

        }
        //ksi==6 ends


        //ksi==7 starts//surfaces[i+1].nodes[(j)*dim_z +k+1].df2[7]  = surfaces[i].nodes[j*dim_z +k].df1[7];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 7\n", my_rank, tid);
        fflush(stdout);
        if ((li + 1) <= cube_size - 1 && (lk + 1) <= cube_size - 1){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li + 1)*cube_size*cube_size + lj*cube_size + lk + 1;
            nodes_df2[node_df2_idx].df2[7] = nodes_df1[node_df1_idx].df1[7];

        }
        else if ((li + 1) >cube_size - 1 && (lk + 1) > cube_size - 1){//comp domain changes both BI and BK and probalby mahcines
          // up-right
          cube2thread_and_machine(BI + 1, BJ, BK + 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){ //different cube in same machine
            cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK + 1);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (0)*cube_size*cube_size + lj*cube_size + 0;//both local i and local k will be starting point of next cube
            nodes_df2[node_df2_idx].df2[7] = nodes_df1[node_df1_idx].df1[7];
          }
          else{ //different cube in different machines

            df1_tosend = nodes_df1[node_df1_idx].df1[7];

            if( (BI%num_cubes_per_proc_x==(num_cubes_per_proc_x-1)) && (BK%num_cubes_per_proc_z==(num_cubes_per_proc_z-1)) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 7, df1_tosend);
            else{
              if(BI%num_cubes_per_proc_x==(num_cubes_per_proc_x-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
              else if(BK%num_cubes_per_proc_z==(num_cubes_per_proc_z-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 3, df1_tosend);
              else{
                printf("Dir 7 Error!\n");
                fflush(stdout);
              }
            }

          }

        }
        else if ((li + 1) >cube_size - 1 && (lk + 1) <= cube_size - 1){//only BI changes and probably machine
          //right
          cube2thread_and_machine(BI + 1, BJ, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){ //different cubes in same machines

            cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (0)*cube_size*cube_size + lj*cube_size + lk + 1;//only local i will be starting point of next cube
            nodes_df2[node_df2_idx].df2[7] = nodes_df1[node_df1_idx].df1[7];

          }
          else{ //different cubes in different machines

            df1_tosend = nodes_df1[node_df1_idx].df1[7];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
          }

        }
        else{//only BK changes and probably machine
          // up
          cube2thread_and_machine(BI, BJ, BK + 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK + 1;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li + 1)*cube_size*cube_size + lj*cube_size + 0;//only local k will be starting point of next cube
            nodes_df2[node_df2_idx].df2[7] = nodes_df1[node_df1_idx].df1[7];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[7];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 3, df1_tosend);
          }

        }
        //ksi==7 ends


        //ksi==8 starts// surfaces[i-1].nodes[(j)*dim_z +k-1].df2[8]  = surfaces[i].nodes[j*dim_z +k].df1[8];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 8\n", my_rank, tid);
        fflush(stdout);
        if ((li - 1) >= 0 && (lk - 1) >= 0){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li - 1)*cube_size*cube_size + lj*cube_size + (lk - 1);
            nodes_df2[node_df2_idx].df2[8] = nodes_df1[node_df1_idx].df1[8];

        }
        else if ((li - 1) <0 && (lk - 1) <0){//comp domain changes both BI and BK and probably machines
          // down-left
          cube2thread_and_machine(BI - 1, BJ, BK - 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI - 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK - 1);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (cube_size - 1)*cube_size*cube_size + lj*cube_size + (cube_size - 1);//both local i and local k will be ending point of prev cube
            nodes_df2[node_df2_idx].df2[8] = nodes_df1[node_df1_idx].df1[8];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[8];

            if( (BI%num_cubes_per_proc_x==0) && (BK%num_cubes_per_proc_z==0) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 8, df1_tosend);
            else{
              if(BI%num_cubes_per_proc_x==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 2, df1_tosend);
              else if(BK%num_cubes_per_proc_z==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 4, df1_tosend);
              else{
                printf("Dir 8 Error!\n");
                fflush(stdout);
              }
            }
          }

        }
        else if ((li - 1)<0 && (lk - 1) >= 0){//only BI changes and probably machines
          //left
          cube2thread_and_machine(BI - 1, BJ, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){ // different cube in same machine

            cube_df2_idx = (BI - 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (cube_size - 1)*cube_size*cube_size + lj*cube_size + lk - 1;//only local i will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[8] = nodes_df1[node_df1_idx].df1[8];

          }
          else{ // different cube in different machine

            df1_tosend = nodes_df1[node_df1_idx].df1[8];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 2, df1_tosend);
          }

        }
        else{//only BK changes and prabably diff machine
          //down
          cube2thread_and_machine(BI, BJ, BK - 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK - 1);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li - 1)*cube_size*cube_size + lj*cube_size + cube_size - 1;//only local k will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[8] = nodes_df1[node_df1_idx].df1[8];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[8];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 4, df1_tosend);
          }

        }
        //ksi==8 ends


        //ksi==9 starts//surfaces[i+1].nodes[(j)*dim_z +k-1].df2[9]  = surfaces[i].nodes[j*dim_z +k].df1[9];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 9\n", my_rank, tid);
        fflush(stdout);
        if ((li + 1) <= cube_size - 1 && (lk - 1) >= 0){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li + 1)*cube_size*cube_size + lj*cube_size + (lk - 1);
            nodes_df2[node_df2_idx].df2[9] = nodes_df1[node_df1_idx].df1[9];

        }
        else if ((li + 1) >cube_size - 1 && (lk - 1)<0){//comp domain changes both BI and BK and probably machines
          //down-right
          cube2thread_and_machine(BI + 1, BJ, BK - 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK - 1;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (0)*cube_size*cube_size + lj*cube_size + cube_size - 1;// local i starting point of next cube and local k will be endng point of prev cube
            nodes_df2[node_df2_idx].df2[9] = nodes_df1[node_df1_idx].df1[9];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[9];

            if( (BI%num_cubes_per_proc_x==(num_cubes_x-1)) && (BK%num_cubes_per_proc_z==0) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 9, df1_tosend);
            else{
              if(BI%num_cubes_per_proc_x==(num_cubes_x-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
              else if(BK%num_cubes_per_proc_z==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 4, df1_tosend);
              else{
                printf("Dir 9 Error!\n");
                fflush(stdout);
              }
            }
          }
        }
        else if ((li + 1) >cube_size - 1 && (lk - 1) >= 0){//only BI changes and probably machine

          cube2thread_and_machine(BI + 1, BJ, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (0)*cube_size*cube_size + lj*cube_size + (lk - 1);//only local i will be starting point of next cube
            nodes_df2[node_df2_idx].df2[9] = nodes_df1[node_df1_idx].df1[9];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[9];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
          }

        }
        else{//only BK changes and probaly diff machine
          cube2thread_and_machine(BI, BJ, BK - 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK - 1;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li + 1)*cube_size*cube_size + lj*cube_size + cube_size - 1;//only local k will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[9] = nodes_df1[node_df1_idx].df1[9];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[9];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 4, df1_tosend);
          }

        }
        //ksi==9 ends


        //ksi==10 starts//surfaces[i-1].nodes[(j)*dim_z +k+1].df2[10] = surfaces[i].nodes[j*dim_z +k].df1[10];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 10\n", my_rank, tid);
        fflush(stdout);
        if ((li - 1) >= 0 && (lk + 1) <= cube_size - 1){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li - 1)*cube_size*cube_size + lj*cube_size + lk + 1;
            nodes_df2[node_df2_idx].df2[10] = nodes_df1[node_df1_idx].df1[10];

        }
        else if ((li - 1) <0 && (lk + 1) >cube_size - 1){//comp domain changes both BI and BK and probably machines
          cube2thread_and_machine(BI - 1, BJ, BK + 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = (BI - 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK + 1;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (cube_size - 1)*cube_size*cube_size + lj*cube_size + 0;//local i will be ending point of prev cube and local k strtng of next
            nodes_df2[node_df2_idx].df2[10] = nodes_df1[node_df1_idx].df1[10];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[10];

            if((BI%num_cubes_per_proc_x==0) && (BK%num_cubes_per_proc_z==(num_cubes_per_proc_z-1)) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 10, df1_tosend);
            else{
              if(BI%num_cubes_per_proc_x==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 2, df1_tosend);
              else if(BK%num_cubes_per_proc_z==(num_cubes_per_proc_z-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 3, df1_tosend);
              else{
                printf("Dir 10 Error!\n");
                fflush(stdout);
              }
            }

          }

        }
        else if ((li - 1) <0 && (lk + 1) <= cube_size - 1){ //only BI changes and probably machines
          cube2thread_and_machine(BI - 1, BJ, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI - 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (cube_size - 1)*cube_size*cube_size + lj*cube_size + lk + 1;//only local i will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[10] = nodes_df1[node_df1_idx].df1[10];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[10];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 2, df1_tosend);
          }


        }
        else{//only BK changes and probably machines, up

          cube2thread_and_machine(BI, BJ, BK + 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK + 1;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li - 1)*cube_size*cube_size + lj*cube_size + 0;//only local k will be strtn point of next cube
            nodes_df2[node_df2_idx].df2[10] = nodes_df1[node_df1_idx].df1[10];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[10];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 3, df1_tosend);
          }

        }
        //ksi==10 ends


        //ksi==11 starts//surfaces[i].nodes[(j-1)*dim_z +k+1].df2[11] = surfaces[i].nodes[j*dim_z +k].df1[11];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 11\n", my_rank, tid);
        fflush(stdout);
        if (lj - 1 >= 0 && lk + 1 <= cube_size - 1){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + (lj - 1)*cube_size + lk + 1;
            nodes_df2[node_df2_idx].df2[11] = nodes_df1[node_df1_idx].df1[11];

        }
        else if (lj - 1 <0 && lk + 1 >cube_size - 1){//comp domain changes both BJ and BK and probably machines, front-up

          cube2thread_and_machine(BI, BJ - 1, BK + 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + (BJ - 1) * num_cubes_z + BK + 1;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + (cube_size - 1)*cube_size + 0;//local j will be ending point of prev cube and local k strtng of next
            nodes_df2[node_df2_idx].df2[11] = nodes_df1[node_df1_idx].df1[11];

          }
          else{// front-up

            df1_tosend = nodes_df1[node_df1_idx].df1[11];

            if( (BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1)) && (BK%num_cubes_per_proc_z==(num_cubes_per_proc_z-1)) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 11, df1_tosend);
            else{
              if(BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 5, df1_tosend);
              else if(BK%num_cubes_per_proc_z==(num_cubes_per_proc_z-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 3, df1_tosend);
              else{
                printf("Dir 11 Error!\n");
                fflush(stdout);
              }
            }

          }

        }
        else if (lj - 1 <0 && lk + 1 <= cube_size - 1){//only BJ changes and probably machine, front

          cube2thread_and_machine(BI, BJ - 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + (BJ - 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + (cube_size - 1)*cube_size + lk + 1;//only local j will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[11] = nodes_df1[node_df1_idx].df1[11];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[11];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 5, df1_tosend);
          }

        }
        else{//only BK changes and probably machine, up

          cube2thread_and_machine(BI, BJ, BK + 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK + 1;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + (lj - 1)*cube_size + 0;//only local k will be strtn point of next cube
            nodes_df2[node_df2_idx].df2[11] = nodes_df1[node_df1_idx].df1[11];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[11];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 3, df1_tosend);
          }

        }
        //ksi==11 ends


        //ksi==12 starts  //surfaces[i].nodes[(j+1)*dim_z +k-1].df2[12] = surfaces[i].nodes[j*dim_z +k].df1[12];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 12\n", my_rank, tid);
        fflush(stdout);
        if ((lj + 1) <= cube_size - 1 && (lk - 1) >= 0){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + (lj + 1)*cube_size + lk - 1;
            nodes_df2[node_df2_idx].df2[12] = nodes_df1[node_df1_idx].df1[12];

        }
        else if ((lj + 1) >cube_size - 1 && (lk - 1)<0){//comp domain changes both BJ and BK and probably machine, back-down

          cube2thread_and_machine(BI, BJ + 1, BK - 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = (BI)* num_cubes_y * num_cubes_z + (BJ + 1) * num_cubes_z + (BK - 1);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + 0 * cube_size + (cube_size - 1);// local j starting point of next cube and local k will be endng point of prev cube
            nodes_df2[node_df2_idx].df2[12] = nodes_df1[node_df1_idx].df1[12];

          }
          else{// back-down

            df1_tosend = nodes_df1[node_df1_idx].df1[12];

            if( (BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1)) && (BK%num_cubes_per_proc_z==0) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 12, df1_tosend);
            else{
              if(BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 6, df1_tosend);
              else if(BK%num_cubes_per_proc_z==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 4, df1_tosend);
              else{
                printf("Dir 12 Error!\n");
                fflush(stdout);
              }
            }

          }

        }
        else if ((lj + 1) >cube_size - 1 && (lk - 1) >= 0){//only BJ changes and probably machine
          cube2thread_and_machine(BI, BJ + 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI)* num_cubes_y * num_cubes_z + (BJ + 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + 0 * cube_size + lk - 1;//only local j will be starting point of next cube
            nodes_df2[node_df2_idx].df2[12] = nodes_df1[node_df1_idx].df1[12];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[12];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 6, df1_tosend);
          }

        }
        else{//only BK changes and probably maxhine
          cube2thread_and_machine(BI, BJ, BK - 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK - 1);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + (lj + 1)*cube_size + cube_size - 1;//only local k will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[12] = nodes_df1[node_df1_idx].df1[12];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[12];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 4, df1_tosend);
          }

        }
        //ksi==12 ends


        //ksi==13 starts //surfaces[i].nodes[(j+1)*dim_z +k+1].df2[13] = surfaces[i].nodes[j*dim_z +k].df1[13];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 13\n", my_rank, tid);
        fflush(stdout);
        if ((lj + 1) <= cube_size - 1 && (lk + 1) <= cube_size - 1){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + (lj + 1)*cube_size + lk + 1;
            nodes_df2[node_df2_idx].df2[13] = nodes_df1[node_df1_idx].df1[13];

        }
        else if ((lj + 1) >cube_size - 1 && (lk + 1) >cube_size - 1){//comp domain changes both BI and BK and probably maxhine, back-up

          cube2thread_and_machine(BI, BJ + 1, BK + 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = (BI)* num_cubes_y * num_cubes_z + (BJ + 1) * num_cubes_z + (BK + 1);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + 0 * cube_size + 0;//both local j and local k will be starting point of next cube
            nodes_df2[node_df2_idx].df2[13] = nodes_df1[node_df1_idx].df1[13];

          }
          else{// back-up

            df1_tosend = nodes_df1[node_df1_idx].df1[13];

            if( (BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1)) && (BK%num_cubes_per_proc_z==(num_cubes_z-1)) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 13, df1_tosend);
            else{
              if(BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 6, df1_tosend);
              else if(BK%num_cubes_per_proc_z==(num_cubes_z-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 3, df1_tosend);
              else{
                printf("Dir 13 Error!\n");
                fflush(stdout);
              }
            }

          }
        }
        else if ((lj + 1) >cube_size - 1 && (lk + 1) <= cube_size - 1){//only BJ changes and probably maxhine, back

          cube2thread_and_machine(BI, BJ + 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = (BI)* num_cubes_y * num_cubes_z + (BJ + 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + 0 * cube_size + lk + 1;//only local j will be starting point of next cube
            nodes_df2[node_df2_idx].df2[13] = nodes_df1[node_df1_idx].df1[13];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[13];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 6, df1_tosend);
          }

        }
        else{//only BK changes and probably maxhine, up
          cube2thread_and_machine(BI, BJ, BK + 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK + 1);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + (lj + 1)*cube_size + 0;//only local k will be starting point of next cube
            nodes_df2[node_df2_idx].df2[13] = nodes_df1[node_df1_idx].df1[13];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[13];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 3, df1_tosend);
          }

        }
        //ksi==13 ends


        //ksi==14 starts//surfaces[i].nodes[(j-1)*dim_z +k-1].df2[14] = surfaces[i].nodes[j*dim_z +k].df1[14];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 14\n", my_rank, tid);
        fflush(stdout);
        if ((lj - 1) >= 0 && (lk - 1) >= 0){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + (lj - 1)*cube_size + (lk - 1);
            nodes_df2[node_df2_idx].df2[14] = nodes_df1[node_df1_idx].df1[14];

        }
        else if ((lj - 1) <0 && (lk - 1) <0){//comp domain changes both BJ and BK and robably maxhine, front-down

          cube2thread_and_machine(BI, BJ - 1, BK - 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI)* num_cubes_y * num_cubes_z + (BJ - 1) * num_cubes_z + (BK - 1);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + (cube_size - 1)*cube_size + (cube_size - 1);//both local j and local k will be ending point of prev cube
            nodes_df2[node_df2_idx].df2[14] = nodes_df1[node_df1_idx].df1[14];
          }
          else{// front-down

            df1_tosend = nodes_df1[node_df1_idx].df1[14];

            if( (BJ%num_cubes_per_proc_y==0) && (BK%num_cubes_per_proc_z==0) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 14, df1_tosend);
            else{
              if(BJ%num_cubes_per_proc_y==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 5, df1_tosend);
              else if(BK%num_cubes_per_proc_z==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 4, df1_tosend);
              else{
                printf("Dir 14 Error!\n");
                fflush(stdout);
              }
            }

          }

        }
        else if ((lj - 1)<0 && (lk - 1) >= 0){//only BJ changes and probably maxhine, front

          cube2thread_and_machine(BI, BJ - 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI)* num_cubes_y * num_cubes_z + (BJ - 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li)*cube_size*cube_size + (cube_size - 1)*cube_size + lk - 1;//only local j will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[14] = nodes_df1[node_df1_idx].df1[14];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[14];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 5, df1_tosend);
          }

        }
        else{//only BK changes and probably machine, down
          cube2thread_and_machine(BI, BJ, BK - 1, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK - 1);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = li*cube_size*cube_size + (lj - 1)*cube_size + cube_size - 1;//only local k will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[14] = nodes_df1[node_df1_idx].df1[14];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[14];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 4, df1_tosend);
          }

        }
        //ksi==14 ends

        //ksi==15 starts//surfaces[i+1].nodes[(j-1)*dim_z +k].df2[15] = surfaces[i].nodes[j*dim_z +k].df1[15];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 15\n", my_rank, tid);
        fflush(stdout);
        if ((li + 1) <= cube_size - 1 && (lj - 1) >= 0){//i.e the the computation domain is still in the current cube, fornt-right
          // if (my_rank == temp_mac_rank_sender){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li + 1)*cube_size*cube_size + (lj - 1)*cube_size + (lk);
            nodes_df2[node_df2_idx].df2[15] = nodes_df1[node_df1_idx].df1[15];
          // }
        }
        else if ((li + 1) >cube_size - 1 && (lj - 1)<0){//comp domain changes both BI and BJ and probably maxhine, front-right
          cube2thread_and_machine(BI + 1, BJ - 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + (BJ - 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (0)*cube_size*cube_size + (cube_size - 1)*cube_size + lk;// local i starting point of next cube and local j will be endng point of prev cube
            nodes_df2[node_df2_idx].df2[15] = nodes_df1[node_df1_idx].df1[15];
          }
          else{// front-right

            df1_tosend = nodes_df1[node_df1_idx].df1[15];

            if( (BI%num_cubes_per_proc_x==(num_cubes_per_proc_x-1)) && (BJ%num_cubes_per_proc_y==0) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 15, df1_tosend);
            else{
              if(BI%num_cubes_per_proc_x==(num_cubes_per_proc_x-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
              else if(BJ%num_cubes_per_proc_y==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 5, df1_tosend);
              else{
                printf("Dir 15 Error!\n");
                fflush(stdout);
              }
            }

          }

        }
        else if ((li + 1) >cube_size - 1 && (lj - 1) >= 0){//only BI changes and probably maxhine, right

          cube2thread_and_machine(BI + 1, BJ, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (0)*cube_size*cube_size + (lj - 1)*cube_size + lk;//only local i will be starting point of next cube
            nodes_df2[node_df2_idx].df2[15] = nodes_df1[node_df1_idx].df1[15];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[15];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
          }

        }
        else{//only BJ changes and probably maxhine
          cube2thread_and_machine(BI, BJ - 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + (BJ - 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li + 1)*cube_size*cube_size + (cube_size - 1)*cube_size + lk;//only local k will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[15] = nodes_df1[node_df1_idx].df1[15];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[15];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 5, df1_tosend);
          }

        }
        //ksi==15 ends

        //ksi==16 starts //surfaces[i-1].nodes[(j+1)*dim_z +k].df2[16] = surfaces[i].nodes[j*dim_z +k].df1[16];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 16\n", my_rank, tid);
        fflush(stdout);
        if ((li - 1) >= 0 && (lj + 1) <= cube_size - 1){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li - 1)*cube_size*cube_size + (lj + 1)*cube_size + lk;
            nodes_df2[node_df2_idx].df2[16] = nodes_df1[node_df1_idx].df1[16];

        }
        else if ((li - 1) <0 && (lj + 1) >cube_size - 1){//comp domain changes both BI and BJ and probably maxhine, back-left

          cube2thread_and_machine(BI - 1, BJ + 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = (BI - 1) * num_cubes_y * num_cubes_z + (BJ + 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (cube_size - 1)*cube_size*cube_size + 0 * cube_size + lk;//local i will be ending point of prev cube and local j strtng of next
            nodes_df2[node_df2_idx].df2[16] = nodes_df1[node_df1_idx].df1[16];

          }
          else{// back-left

            df1_tosend = nodes_df1[node_df1_idx].df1[16];

            if( (BI%num_cubes_per_proc_x==0) && (BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1)) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 16, df1_tosend);
            else{
              if(BI%num_cubes_per_proc_x==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 2, df1_tosend);
              else if(BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 6, df1_tosend);
              else{
                printf("Dir 16 Error!\n");
                fflush(stdout);
              }
            }

          }

        }
        else if ((li - 1) <0 && (lj + 1) <= cube_size - 1){//only BI changes and probably maxhine, left

          cube2thread_and_machine(BI - 1, BJ, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = (BI - 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (cube_size - 1)*cube_size*cube_size + (lj + 1)*cube_size + lk;//only local i will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[16] = nodes_df1[node_df1_idx].df1[16];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[16];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 2, df1_tosend);
          }

        }
        else{//only BJ changes and probably machine, back
          cube2thread_and_machine(BI, BJ + 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + (BJ + 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li - 1)*cube_size*cube_size + 0 * cube_size + lk;//only local k will be strtn point of next cube
            nodes_df2[node_df2_idx].df2[16] = nodes_df1[node_df1_idx].df1[16];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[16];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 6, df1_tosend);
          }

        }
        //ksi==16 ends

        //ksi==17 starts //surfaces[i+1].nodes[(j+1)*dim_z +k].df2[17] = surfaces[i].nodes[j*dim_z +k].df1[17];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 17\n", my_rank, tid);
        fflush(stdout);
        if ((li + 1) <= cube_size - 1 && (lj + 1) <= cube_size - 1){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li + 1)*cube_size*cube_size + (lj + 1)*cube_size + lk;
            nodes_df2[node_df2_idx].df2[17] = nodes_df1[node_df1_idx].df1[17];

        }
        else if ((li + 1) >cube_size - 1 && (lj + 1) >cube_size - 1){//comp domain changes both BI and BJ and probably maxchines, back-right
          cube2thread_and_machine(BI + 1, BJ + 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + (BJ + 1) * num_cubes_z + (BK);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (0)*cube_size*cube_size + (0)*cube_size + lk;//both local i and local j will be starting point of next cube
            nodes_df2[node_df2_idx].df2[17] = nodes_df1[node_df1_idx].df1[17];
          }
          else{// front-right

            df1_tosend = nodes_df1[node_df1_idx].df1[17];

            if( (BI%num_cubes_per_proc_x==(num_cubes_per_proc_x-1)) && (BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1)) )
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 17, df1_tosend);
            else{
              if(BI%num_cubes_per_proc_x==(num_cubes_per_proc_x-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
              else if(BJ%num_cubes_per_proc_y==(num_cubes_per_proc_y-1))
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 6, df1_tosend);
              else{
                printf("Dir 17 Error!\n");
                fflush(stdout);
              }
            }

          }

        }
        else if ((li + 1) >cube_size - 1 && (lj + 1) <= cube_size - 1){//only BI changes and probably machiine
          cube2thread_and_machine(BI + 1, BJ, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI + 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (0)*cube_size*cube_size + (lj + 1)*cube_size + lk;//only local i will be starting point of next cube
            nodes_df2[node_df2_idx].df2[17] = nodes_df1[node_df1_idx].df1[17];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[17];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 1, df1_tosend);
          }

        }
        else{//only BJ changes and probably machine, back
          cube2thread_and_machine(BI, BJ + 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + (BJ + 1) * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li + 1)*cube_size*cube_size + 0 * cube_size + lk;//only local k will be starting point of next cube
            nodes_df2[node_df2_idx].df2[17] = nodes_df1[node_df1_idx].df1[17];

          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[17];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 6, df1_tosend);
          }

        }
        //ksi==17 ends


        //ksi==18 starts  //surfaces[i-1].nodes[(j-1)*dim_z +k].df2[18] = surfaces[i].nodes[j*dim_z +k].df1[18];
        printf("Fluid mac%dtid%d: Start Prepare streaming on Dir 18\n", my_rank, tid);
        fflush(stdout);
        if ((li - 1) >= 0 && (lj - 1) >= 0){//i.e the the computation domain is still in the current cube

            cube_df2_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li - 1)*cube_size*cube_size + (lj - 1)*cube_size + (lk);
            nodes_df2[node_df2_idx].df2[18] = nodes_df1[node_df1_idx].df1[18];

        }
        else if ((li - 1) <0 && (lj - 1) <0){//comp domain changes both BI and BJ and proabaly machine, front-left
          cube2thread_and_machine(BI - 1, BJ - 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI - 1) * num_cubes_y * num_cubes_z + (BJ - 1) * num_cubes_z + (BK);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (cube_size - 1)*cube_size*cube_size + (cube_size - 1)*cube_size + lk;//both local i and local j will be ending point of prev cube
            nodes_df2[node_df2_idx].df2[18] = nodes_df1[node_df1_idx].df1[18];
          }
          else{// front-left

            df1_tosend = nodes_df1[node_df1_idx].df1[18];

            if( (BI%num_cubes_per_proc_x==0) && (BJ%num_cubes_per_proc_y==0))
              insert_msg(gv, BI, BJ, BK, li, lj, lk, 18, df1_tosend);
            else{
              if(BI%num_cubes_per_proc_x==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 2, df1_tosend);
              else if(BJ%num_cubes_per_proc_y==0)
                insert_msg(gv, BI, BJ, BK, li, lj, lk, 5, df1_tosend);
              else{
                printf("Dir 18 Error!\n");
                fflush(stdout);
              }
            }

          }

        }
        else if ((li - 1)<0 && (lj - 1) >= 0){//only BI changes and prabably machine, left

          cube2thread_and_machine(BI - 1, BJ, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = (BI - 1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (cube_size - 1)*cube_size*cube_size + (lj - 1)*cube_size + lk;//only local i will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[18] = nodes_df1[node_df1_idx].df1[18];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[18];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 2, df1_tosend);
          }

        }
        else{//only BJ changes and probably machine, front
          cube2thread_and_machine(BI, BJ - 1, BK, gv, &temp_mac_rank_recv);

          if(temp_mac_rank_sender == temp_mac_rank_recv){
            cube_df2_idx = BI * num_cubes_y * num_cubes_z + (BJ - 1) * num_cubes_z + (BK);
            nodes_df2 = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
            node_df2_idx = (li - 1)*cube_size*cube_size + (cube_size - 1)*cube_size + lk;//only local k will be endg point of prev cube
            nodes_df2[node_df2_idx].df2[18] = nodes_df1[node_df1_idx].df1[18];
          }
          else{

            df1_tosend = nodes_df1[node_df1_idx].df1[18];

            insert_msg(gv, BI, BJ, BK, li, lj, lk, 5, df1_tosend);
          }

        }
        //ksi==18 ends

      }//lk
    }// if cube2thread ends
  }//BK

  printf("Fluid mac%d: Pass Prepare streaming msg\n", my_rank);
  fflush(stdout);

  pthread_barrier_wait(&(gv->barr));
  // if(tid==0)
  //   MPI_Barrier(MPI_COMM_WORLD);

  //send message
  if (tid == 0){
    //add condition 1d,2d,3d
    int direction, i, dest; //int flag;
    if(Px==1 && Py==1 && Pz==1){
      printf("Fluid mac%d: don't need to stream!\n", my_rank);
      fflush(stdout);
    }
    else if (check_1d(gv->Px, gv->Py, gv->Pz, &direction)){//1d

      printf("Fluid mac%d: Enter 1D streaming msg!\n", my_rank);
      fflush(stdout);

      if(direction == X_transfer_1D){//1,2

        printf("Fluid mac%d: Enter X_transfer_1D streaming msg!\n", my_rank);
        fflush(stdout);

        //dir 1
        streaming_on_direction(gv, lv, 1, dest);

        //dir 2
        streaming_on_direction(gv, lv, 2, dest);
      }
      else if(direction == Z_transfer_1D){
        //dir 3
        streaming_on_direction(gv, lv, 3, dest);

        //dir 4
        streaming_on_direction(gv, lv, 4, dest);

      }
      else if(direction == Y_transfer_2D){
        //dir 5
        streaming_on_direction(gv, lv, 5, dest);

        //dir 6
        streaming_on_direction(gv, lv, 6, dest);
      }
      else{
        printf("1d streaming Error!\n");
        fflush(stdout);
      }
    }
    else if(check_2d(gv->Px, gv->Py, gv->Pz, &direction)) {//2d

      if(direction==X_transfer_2D){// 1, 2, 5, 6, 15~18

        //dir 1, 2
        for(i=1; i<=2; i++)
          streaming_on_direction(gv, lv, i, gv->dest_mac[i]);

        //dir 5, 6
        for(i=5; i<=6; i++)
          streaming_on_direction(gv, lv, i, gv->dest_mac[i]);

        //dir 15~18
        for(i=15; i<=18; i++)
          streaming_on_direction(gv, lv, i, gv->dest_mac[i]);

      }
      else if(direction==Y_transfer_2D){// 3, 4, 5, 6, 11~14

        //dir 3~6
        for(i=3; i<=6; i++)
          streaming_on_direction(gv, lv, i, gv->dest_mac[i]);

        //dir 11~14
        for(i=11; i<=14; i++)
          streaming_on_direction(gv, lv, i, gv->dest_mac[i]);

      }
      else if(direction==Z_transfer_2D){// 1, 2, 3, 4, 7~10
        //dir 1~4
        for(i=1; i<=4; i++)
          streaming_on_direction(gv, lv, i, gv->dest_mac[i]);

        //dir 7~10
        for(i=7; i<=10; i++)
          streaming_on_direction(gv, lv, i, gv->dest_mac[i]);

      }
      else{
        printf("2D Error!\n");
      }

    }
    else{//3d
      //dir 1~18
      for(i=7; i<=10; i++)
        streaming_on_direction(gv, lv, i, gv->dest_mac[i]);
    }
  }




#ifdef DEBUG_PRINT
  printf("****************stream_distfunc() Exit*******\n");
#endif //DEBUG_PRINT
}
