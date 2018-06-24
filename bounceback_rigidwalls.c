/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/
#include "do_thread.h"

void bounceback_rigidwalls(LV lv){
   //printf("******************Inside bounceback_rigidwalls************ at time  %d\n");
  int tid;
  GV gv = lv->gv;
  tid   = lv->tid;

  Fluidgrid     *fluidgrid;
  Fluidnode     *nodes;


  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube
  int           dim_z;

  fluidgrid      = gv->fluid_grid;
  dim_z          = fluidgrid->z_dim;

  cube_size   = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;

   int my_rank, temp_mac_rank;
   my_rank = gv->my_rank;

   /* 1.1 half-way bounce back on the bottom */
   //k = gv->kb;
   lk = gv->kb; BK = 0;

   for(BI =0; BI < num_cubes_x; ++BI)
    for(BJ =0; BJ < num_cubes_y; ++BJ){
      if(cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank) ==tid){
        if(my_rank == temp_mac_rank){
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_x = starting_y = 0;
      stopping_x = stopping_y = cube_size-1;
      if(BI == 0) starting_x = gv->ib+1;//3

      if(BI == num_cubes_x-1) stopping_x = cube_size-4;//gv->ie-1

      if(BJ == 0) starting_y = gv->jb;//2

      if(BJ == num_cubes_y-1) stopping_y = cube_size-3;//je

     for(li=starting_x; li <= stopping_x; ++li)
      for(lj=starting_y; lj <= stopping_y; ++lj){
            node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. //lk =gv->kb
         nodes[node_idx].df2[3] =nodes[node_idx].df1[4];
         nodes[node_idx].df2[7] =nodes[node_idx].df1[8];
         nodes[node_idx].df2[10]=nodes[node_idx].df1[9];
         nodes[node_idx].df2[11]=nodes[node_idx].df1[12];
         nodes[node_idx].df2[13]=nodes[node_idx].df1[14];


       }
      }// if mahcine chk ends
      }//if cube2thread ends
     }

   /* 1.2 bounce back on the top*/
 //k=gv->ke;
  lk = cube_size-3;
  BK = num_cubes_z -1;//lk=gv->ke;//BI=gv->ib+1;BI<=gv->ie-1;BJ=gv->jb;BJ<=gv->je
   for(BI =0; BI < num_cubes_x; ++BI)
    for(BJ =0; BJ < num_cubes_y; ++BJ){
       if(cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank) ==tid){
        if(my_rank == temp_mac_rank){
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_x = starting_y = 0;
      stopping_x = stopping_y = cube_size-1;
      if(BI == 0) starting_x = gv->ib+1;//3

      if(BI == num_cubes_x-1) stopping_x = cube_size-4;//gv->ie-1

      if(BJ == 0) starting_y = gv->jb;//2

      if(BJ == num_cubes_y-1) stopping_y = cube_size-3;//je

     for(li=starting_x; li <= stopping_x; ++li)
      for(lj=starting_y; lj <= stopping_y; ++lj){
            node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. //lk =gv->ke
       nodes[node_idx].df2[4] =nodes[node_idx].df1[3];
       nodes[node_idx].df2[8] =nodes[node_idx].df1[7];
       nodes[node_idx].df2[9] =nodes[node_idx].df1[10];
       nodes[node_idx].df2[12]=nodes[node_idx].df1[11];
       nodes[node_idx].df2[14]=nodes[node_idx].df1[13];
         }
        }// if mahcine chk ends
      }//if cube2 thread ends
     }

   /* 1.3 bounce back on the front*/
  //j=gv->jb;
       //BJ = gv->jb;BI=gv->ib+1;BI<=gv->ie-1;BK=gv->kb;BK<=gv->ke
   lj = gv->jb; BJ =0;
   for(BI =0; BI < num_cubes_x; ++BI)
    for(BK =0; BK < num_cubes_z; ++BK){
      if(cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank) ==tid){
        if(my_rank == temp_mac_rank){
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_x = starting_z = 0;
      stopping_x  = stopping_z = cube_size-1;
      if(BI == 0) starting_x = gv->ib+1;//3

      if(BI == num_cubes_x-1) stopping_x = cube_size-4;//gv->ie-1

      if(BK == 0) starting_z = gv->kb;//2

      if(BK == num_cubes_x-1) stopping_z = cube_size-3;//gv->ke

     for(li=starting_x; li <= stopping_x; ++li)
      for(lk=starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. //lj =gv->jb =2
           nodes[node_idx].df2[6] =nodes[node_idx].df1[5];
           nodes[node_idx].df2[12]=nodes[node_idx].df1[11];
           nodes[node_idx].df2[13]=nodes[node_idx].df1[14];
           nodes[node_idx].df2[16]=nodes[node_idx].df1[15];
           nodes[node_idx].df2[17]=nodes[node_idx].df1[18];
          }
      }//if machine chek ends
     }//if cube2 thread ends
   }

   /* 1.4 bounce back on the rear*/
  //j=gv->je;
       //BJ = gv->je;BI=gv->ib+1;BI<=gv->ie-1;BK=gv->kb;BK<=gv->ke
   lj = cube_size-3; BJ = num_cubes_y-1;//gv->je;
   for(BI =0; BI < num_cubes_x; ++BI)
    for(BK =0; BK < num_cubes_z; ++BK){
      if(cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank) ==tid){
        if(my_rank == temp_mac_rank){
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_x = starting_z = 0;
      stopping_x = stopping_z = cube_size-1;
      if(BI == 0) starting_x = gv->ib+1;//3

      if(BI == num_cubes_x-1) stopping_x = cube_size-4;//gv->ie-1

      if(BK == 0) starting_z = gv->kb;//2

      if(BK == num_cubes_x-1) stopping_z = cube_size-3;//gv->ke

     for(li=starting_x; li <= stopping_x; ++li)
      for(lk=starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. //lj =gv->je
             nodes[node_idx].df2[5] =nodes[node_idx].df1[6];
             nodes[node_idx].df2[11]=nodes[node_idx].df1[12];
             nodes[node_idx].df2[14]=nodes[node_idx].df1[13];
             nodes[node_idx].df2[15]=nodes[node_idx].df1[16];
             nodes[node_idx].df2[18]=nodes[node_idx].df1[17];
          }
         }//if machine chk ends
       }// if cube2thread ends
     }


   //printf("******************bounceback_rigidwalls Exit************\n ");

 }
