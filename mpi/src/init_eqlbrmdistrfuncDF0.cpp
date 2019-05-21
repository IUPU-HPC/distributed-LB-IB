/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"
#include "lb.h"

void init_eqlbrmdistrfuncDF0(GV gv){/*stored in dfeq*/
  Fluidgrid     *fluidgrid;

  int starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone

  fluidgrid = gv->fluid_grid;

  int cube_size = gv->cube_size;
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;
  int tmp_task;

  /*PTHREAD_Change*/

  //For each cube BI, BJ ,BK, li, lj, lk 0 to cube size-1
  // for (int BI = 0; BI < num_cubes_x; ++BI){
  //   for (int BJ = 0; BJ < num_cubes_y; ++BJ){
  //     for (int BK = 0; BK < num_cubes_z; ++BK){
  for (int BI = gv->start_B[0]; BI < gv->stop_B[0]; ++BI)
    for (int BJ = gv->start_B[1]; BJ < gv->stop_B[1]; ++BJ)
      for (int BK = gv->start_B[2]; BK < gv->stop_B[2]; ++BK){
        tmp_task = cube2task(BI, BJ, BK, gv);
        if (gv->taskid == tmp_task){ //MPI changes
          long cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          Fluidnode *nodes = fluidgrid->sub_fluid_grid[cube_idx].nodes;
          starting_x = starting_y = starting_z = 0;
          stopping_x = stopping_y = stopping_z = cube_size - 1;
          for (int li = starting_x; li <= stopping_x; ++li){
            for (int lj = starting_y; lj <= stopping_y; ++lj){
              for (int lk = starting_z; lk <= stopping_z; ++lk){
                int node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
                computeEquilibrium(nodes+node_idx);

                // for (int ksi = 0; ksi<19; ++ksi){
                //   if (ksi == 0){
                //     nodes[node_idx].dfeq[ksi] =
                //       1.0 / 3.0 * nodes[node_idx].rho
                //       * (1.0 - 1.5 *
                //       (nodes[node_idx].vel_x * nodes[node_idx].vel_x
                //       + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                //       + nodes[node_idx].vel_z * nodes[node_idx].vel_z));
                //   }
                //   else if (ksi >0 && ksi<7) {
                //     nodes[node_idx].dfeq[ksi] =
                //       1.0 / 18.0 * nodes[node_idx].rho
                //       *(1.0 + 3.0 *
                //       (gv->c[ksi][0] * nodes[node_idx].vel_x
                //       + gv->c[ksi][1] * nodes[node_idx].vel_y
                //       + gv->c[ksi][2] * nodes[node_idx].vel_z)
                //       + 4.5 *(pow(gv->c[ksi][0] * nodes[node_idx].vel_x
                //       + gv->c[ksi][1] * nodes[node_idx].vel_y
                //       + gv->c[ksi][2] * nodes[node_idx].vel_z, 2))
                //       - 1.5 *(nodes[node_idx].vel_x * nodes[node_idx].vel_x
                //       + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                //       + nodes[node_idx].vel_z * nodes[node_idx].vel_z));
                //   }
                //   else {
                //     nodes[node_idx].dfeq[ksi] =
                //       1.0 / 36.0 * nodes[node_idx].rho
                //       *(1.0 + 3.0*
                //       (gv->c[ksi][0] * nodes[node_idx].vel_x
                //       + gv->c[ksi][1] * nodes[node_idx].vel_y
                //       + gv->c[ksi][2] * nodes[node_idx].vel_z)
                //       + 4.5 *(pow(gv->c[ksi][0] * nodes[node_idx].vel_x
                //       + gv->c[ksi][1] * nodes[node_idx].vel_y
                //       + gv->c[ksi][2] * nodes[node_idx].vel_z, 2))
                //       - 1.5 *(nodes[node_idx].vel_x *nodes[node_idx].vel_x
                //       + nodes[node_idx].vel_y *nodes[node_idx].vel_y
                //       + nodes[node_idx].vel_z *nodes[node_idx].vel_z));
                //   }
                // } //for loop ksi ends

              }//for loop local z
            }//for loop local y
          } // for loop local x
        }//if machine check ends
      }//for BK

  /*PTHREAD_Change*/
  // printf("****************init_eqlbrmdistrfuncDF0 EXIT ***************\n");
}
