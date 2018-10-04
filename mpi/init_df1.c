/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

void init_df1(GV gv) {
  Fluidnode     *nodes;
  int           ksi;
  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube
  int           temp_mac_rank;

  cube_size = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;
  /*PTHREAD_Change*/
  for (BI = 0; BI < num_cubes_x; ++BI)
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
  for (BK = 0; BK < num_cubes_z; ++BK){
    cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank);
    if (gv->my_rank == temp_mac_rank){ //MPI changes
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
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
        node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
        for (ksi = 0; ksi <= 18; ksi++)
          nodes[node_idx].df1[ksi] = nodes[node_idx].dfeq[ksi];
      }//local k
    }//if machine check ends
  }//for BK
  /*PTHREAD_Change*/
  //printf("****************init_df1 EXIT ***************\n");
}
