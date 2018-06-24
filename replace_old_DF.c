/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

void replace_old_DF(LV lv){
  //printf("Inside replace_oldDF\n");
  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

  Fluidgrid     *fluidgrid;
  Fluidnode     *nodes;

  int           ksi;
  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;
  int my_rank = gv->my_rank;
  int temp_mac_rank;

  fluidgrid = gv->fluid_grid;

  cube_size = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  /* replacing the old d.f. values by the newly computed ones */
  for (BI = 0; BI < num_cubes_x; ++BI)
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
  for (BK = 0; BK < num_cubes_z; ++BK){
    if (cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank) == tid){
      if (my_rank == temp_mac_rank){
        cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        nodes = fluidgrid->sub_fluid_grid[cube_idx].nodes;
        starting_x = starting_y = starting_z = 0;
        stopping_x = stopping_y = stopping_z = cube_size - 1;
        if (BI == 0) starting_x = 1;//ib-1

        if (BI == num_cubes_x - 1) stopping_x = cube_size - 2;//ie+1

        if (BJ == 0) starting_y = 1;//jb-1

        if (BJ == num_cubes_y - 1) stopping_y = cube_size - 2;//je+1

        if (BK == 0) starting_z = 1;//kb-1

        if (BK == num_cubes_z - 1) stopping_z = cube_size - 2;//ke+1

        for (li = starting_x; li <= stopping_x; ++li)
        for (lj = starting_y; lj <= stopping_y; ++lj)
        for (lk = starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk;
          for (ksi = 0; ksi <= 18; ksi++)
            nodes[node_idx].df1[ksi] = nodes[node_idx].df2[ksi];

        }//lk
      }//if machine check
    }//if cube2thread
  }//BK

  //printf("replace_oldDF EXIT\n");
}
