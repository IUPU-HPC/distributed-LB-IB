/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

void init_eqlbrmdistrfuncDF0(GV gv){/*stored in dfeq*/
  Fluidgrid     *fluidgrid;
  Fluidnode     *nodes;

  int           ksi;
  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube

  fluidgrid = gv->fluid_grid;

  cube_size = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;
  int           temp_mac_rank;

  /*PTHREAD_Change*/

  //For each cube BI, BJ ,BK, li, lj, lk 0 to cube size-1
  for (BI = 0; BI < num_cubes_x; ++BI){
    for (BJ = 0; BJ < num_cubes_y; ++BJ){
      for (BK = 0; BK < num_cubes_z; ++BK){
        cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank);
        if (gv->my_rank == temp_mac_rank){ //MPI changes
          cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes = fluidgrid->sub_fluid_grid[cube_idx].nodes;
          starting_x = starting_y = starting_z = 0;
          stopping_x = stopping_y = stopping_z = cube_size - 1;
          for (li = starting_x; li <= stopping_x; ++li){
            for (lj = starting_y; lj <= stopping_y; ++lj){
              for (lk = starting_z; lk <= stopping_z; ++lk){
                node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
                for (ksi = 0; ksi<19; ++ksi){
                  if (ksi == 0){
                    nodes[node_idx].dfeq[ksi] =
                      1.0 / 3.0 * nodes[node_idx].rho
                      * (1.0 - 1.5 *
                      (nodes[node_idx].vel_x * nodes[node_idx].vel_x
                      + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                      + nodes[node_idx].vel_z * nodes[node_idx].vel_z));
                  }
                  else if (ksi >0 && ksi<7) {
                    nodes[node_idx].dfeq[ksi] =
                      1.0 / 18.0 * nodes[node_idx].rho
                      *(1.0 + 3.0 *
                      (gv->c[ksi][0] * nodes[node_idx].vel_x
                      + gv->c[ksi][1] * nodes[node_idx].vel_y
                      + gv->c[ksi][2] * nodes[node_idx].vel_z)
                      + 4.5 *(pow(gv->c[ksi][0] * nodes[node_idx].vel_x
                      + gv->c[ksi][1] * nodes[node_idx].vel_y
                      + gv->c[ksi][2] * nodes[node_idx].vel_z, 2))
                      - 1.5 *(nodes[node_idx].vel_x * nodes[node_idx].vel_x
                      + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                      + nodes[node_idx].vel_z * nodes[node_idx].vel_z));
                  }
                  else {
                    nodes[node_idx].dfeq[ksi] =
                      1.0 / 36.0 * nodes[node_idx].rho
                      *(1.0 + 3.0*
                      (gv->c[ksi][0] * nodes[node_idx].vel_x
                      + gv->c[ksi][1] * nodes[node_idx].vel_y
                      + gv->c[ksi][2] * nodes[node_idx].vel_z)
                      + 4.5 *(pow(gv->c[ksi][0] * nodes[node_idx].vel_x
                      + gv->c[ksi][1] * nodes[node_idx].vel_y
                      + gv->c[ksi][2] * nodes[node_idx].vel_z, 2))
                      - 1.5 *(nodes[node_idx].vel_x *nodes[node_idx].vel_x
                      + nodes[node_idx].vel_y *nodes[node_idx].vel_y
                      + nodes[node_idx].vel_z *nodes[node_idx].vel_z));
                  }
                } //for loop ksi ends
              }//for loop local z
            }//for loop local y
          } // for loop local x
        }//if machine check ends
      }//for BK
    }//for BJ
  }//FOR BI
  /*PTHREAD_Change*/
  // printf("****************init_eqlbrmdistrfuncDF0 EXIT ***************\n");
}
