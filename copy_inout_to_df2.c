/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

void copy_inout_to_df2(LV lv){
  //printf("***********************Inside copy_inout_to_df2***********\n");
  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

  Fluidgrid     *fluidgrid;
  Fluidnode     *nodes;

  int           ksi;
  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_y, starting_z, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube
  int           dim_z;

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;
  int my_rank, temp_mac_rank;
  my_rank = gv->my_rank;

  fluidgrid = gv->fluid_grid;
  dim_z = fluidgrid->z_dim;

  cube_size = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;


  //i = gv->ib;BJ=gv->jb-1; BJ<=gv->je+1;BK=gv->kb-1; BK<=gv->ke+1;
  li = gv->ib; BI = 0;
  for (BJ = 0; BJ < num_cubes_y; BJ++)
    for (BK = 0; BK < num_cubes_x; BK++){
      if (cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank) == tid){
        if (my_rank == temp_mac_rank){
          cube_idx = BI *num_cubes_y*num_cubes_z + BJ *num_cubes_z + BK;
          nodes = fluidgrid->sub_fluid_grid[cube_idx].nodes;
          starting_y = starting_z = 0;
          stopping_y = stopping_z = cube_size - 1;
          if (BJ == 0) starting_y = gv->jb - 1;

          if (BJ == num_cubes_y - 1) stopping_y = cube_size - 2;//gv->je+1

          if (BK == 0) starting_z = gv->kb - 1;

          if (BK == num_cubes_z - 1) stopping_z = cube_size - 2;//ke +1

          for (lj = starting_y; lj <= stopping_y; ++lj)
            for (lk = starting_z; lk <= stopping_z; ++lk){
              node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. li =gv->ib
              for (ksi = 0; ksi <= 18; ksi++)
                nodes[node_idx].df2[ksi] = gv->fluid_grid->inlet->nodes[((BJ*cube_size + lj)*dim_z + BK*cube_size + lk)].df_inout[0][ksi];//is this correct?
            }//lk
        }//if machine chek
      }//if cube2thread ends
    }


  //i=gv->ib-1;

  li = gv->ib - 1; BI = 0;
  for (BJ = 0; BJ < num_cubes_y; BJ++)
    for (BK = 0; BK < num_cubes_x; BK++){
      if (cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank) == tid){
        if (my_rank == temp_mac_rank){
          cube_idx = BI *num_cubes_y*num_cubes_z + BJ *num_cubes_z + BK;
          nodes = fluidgrid->sub_fluid_grid[cube_idx].nodes;
          starting_y = starting_z = 0;
          stopping_y = stopping_z = cube_size - 1;
          if (BJ == 0) starting_y = gv->jb - 1;

        if (BJ == num_cubes_y - 1) stopping_y = cube_size - 2;//gv->je+1

        if (BK == 0) starting_z = gv->kb - 1;

        if (BK == num_cubes_z - 1) stopping_z = cube_size - 2;//ke+1

        for (lj = starting_y; lj <= stopping_y; ++lj)
          for (lk = starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. li =gv->ib-1
          for (ksi = 0; ksi <= 18; ksi++)
            nodes[node_idx].df2[ksi] = gv->fluid_grid->inlet->nodes[((BJ*cube_size + lj)*dim_z + BK*cube_size + lk)].df_inout[0][ksi];
        }//lk
      }//if machine chek
    }//if cube2thread ends
  }



  /* outlet */
  //i=gv->ie;
  li = cube_size - 3; BI = num_cubes_x - 1;
  for (BJ = 0; BJ < num_cubes_y; BJ++)
    for (BK = 0; BK < num_cubes_x; BK++){
      if (cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank) == tid){
        if (my_rank == temp_mac_rank){
          cube_idx = BI *num_cubes_y*num_cubes_z + BJ *num_cubes_z + BK;
          nodes = fluidgrid->sub_fluid_grid[cube_idx].nodes;
          starting_y = starting_z = 0;
          stopping_y = stopping_z = cube_size - 1;
          if (BJ == 0) starting_y = gv->jb - 1;

        if (BJ == num_cubes_y - 1) stopping_y = cube_size - 2;//gv->je+1

        if (BK == 0) starting_z = gv->kb - 1;

        if (BK == num_cubes_z - 1) stopping_z = cube_size - 2;//ke +1

        for (lj = starting_y; lj <= stopping_y; ++lj)
          for (lk = starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk;//local node index inside a cube. li =gv->ie
          for (ksi = 0; ksi <= 18; ksi++)
            nodes[node_idx].df2[ksi] = gv->fluid_grid->outlet->nodes[((BJ*cube_size + lj)*dim_z + BK*cube_size + lk)].df_inout[1][ksi];
        }//lk
      }//if machine chek
    }//if cube2thread ends
  }

  //i=gv->ie+1;
  li = cube_size - 2; BI = num_cubes_x - 1;
  for (BJ = 0; BJ < num_cubes_y; BJ++)
    for (BK = 0; BK < num_cubes_x; BK++){
      if (cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank) == tid){
        if (my_rank == temp_mac_rank){
          cube_idx = BI *num_cubes_y*num_cubes_z + BJ *num_cubes_z + BK;
          nodes = fluidgrid->sub_fluid_grid[cube_idx].nodes;
          starting_y = starting_z = 0;
          stopping_y = stopping_z = cube_size - 1;
          if (BJ == 0) starting_y = gv->jb - 1;

        if (BJ == num_cubes_y - 1) stopping_y = cube_size - 2;//gv->je+1

        if (BK == 0) starting_z = gv->kb - 1;

        if (BK == num_cubes_z - 1) stopping_z = cube_size - 2;//ke +1

        for (lj = starting_y; lj <= stopping_y; ++lj)
          for (lk = starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk;//local node index inside a cube. li =gv->ie+1
          for (ksi = 0; ksi <= 18; ksi++)
            nodes[node_idx].df2[ksi] = gv->fluid_grid->outlet->nodes[((BJ*cube_size + lj)*dim_z + BK*cube_size + lk)].df_inout[1][ksi];
        }//lk
      }//if machine chek
    }//if cube2thread ends
  }

  //printf("*********************** copy_inout_to_df2 Exit ***********\n");
}
