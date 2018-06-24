/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

/* Called only once before simulation starts */
void init_df_inout(GV gv){
  //printf("****************Inside copy_df_to_inout*************\n");
  Fluidgrid     *fluidgrid;
  Fluidnode     *nodes;

  int           ksi, dim_z;
  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_y, starting_z, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube
  int           temp_mac_rank;

  fluidgrid = gv->fluid_grid;
  dim_z = fluidgrid->z_dim;

  cube_size = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;
  //printf("*************ENTRY OF init_df_inout\n");
  // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);

  /* inlet, copy once */
  li = gv->ib; BI = 0;
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
  for (BK = 0; BK < num_cubes_z; ++BK){
    cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank);
    if (gv->my_rank == temp_mac_rank){ //MPI changes
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_y = starting_z = 0;
      stopping_y = stopping_z = cube_size - 1;
      for (lj = starting_y; lj <= stopping_y; ++lj)
      for (lk = starting_z; lk <= stopping_z; ++lk){
        node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
        for (ksi = 0; ksi <= 18; ksi++){
          fluidgrid->inlet->nodes[(BJ*cube_size + lj)*dim_z + BK*cube_size + lk].df_inout[0][ksi] =
            nodes[node_idx].dfeq[ksi];
        }
      }
    }//if machine check ends
  }

  /* outlet, copy once */
  li = cube_size - 1; BI = num_cubes_x - 1;//ie
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
  for (BK = 0; BK < num_cubes_z; ++BK){
    cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank);
    if (gv->my_rank == temp_mac_rank){ //MPI changes
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_y = starting_z = 0;
      stopping_y = stopping_z = cube_size - 1;
      for (lj = starting_y; lj <= stopping_y; ++lj)
      for (lk = starting_z; lk <= stopping_z; ++lk){
        node_idx = li * cube_size * cube_size + lj * cube_size + lk;
        for (ksi = 0; ksi <= 18; ksi++){
          fluidgrid->outlet->nodes[(BJ*cube_size + lj)*dim_z + BK*cube_size + lk].df_inout[1][ksi] =
            nodes[node_idx].dfeq[ksi];
        }
      }
    }//if machine check ends
  }
  // printf("****************init_df_inout: EXIT*************\n");
}

