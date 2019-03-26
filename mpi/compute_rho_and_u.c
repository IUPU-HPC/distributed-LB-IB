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
 * @file: compute_rho_and_u.c
 *
 * @date:
 */

#include "do_thread.h"

extern int c[19][3];

void compute_rho_and_u(LV lv){
  //printf("**********Inside compute_rho**********\n");
  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

  Fluidgrid *fluidgrid;
  Fluidnode *nodes;

  int ksi;
  double s1, s2, s3, s4; //sum of df2, df2[x,y,z] * c[ksi][x,y,z]
  int BI, BJ, BK; //to identify the Sub grids
  long cube_idx;
  int starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z; //To identify buffer zone
  int li, lj, lk, node_idx; //local access point inside cube

  fluidgrid = gv->fluid_grid;

  int cube_size = gv->cube_size;
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;
  int my_rank, toProc, ownertid;
  my_rank = gv->taskid;
  s1 = s2 = s3 = s4 = 0;

  for (BI = 0; BI < num_cubes_x; ++BI)
    for (BJ = 0; BJ < num_cubes_y; ++BJ)//for computing womega near bdy
      for (BK = 0; BK < num_cubes_z; ++BK){
        ownertid = cube2thread_and_task(BI, BJ, BK, gv, &toProc);
        if (my_rank == toProc && ownertid == tid){

            cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
            nodes = fluidgrid->sub_fluid_grid[cube_idx].nodes;

            starting_x = starting_y = starting_z = 0;
            stopping_x = stopping_y = stopping_z = cube_size - 1;

            if (BI == 0) starting_x = 3;//ib+1
            if (BI == num_cubes_x - 1) stopping_x = cube_size - 4;//ie-1

            if (BJ == 0) starting_y = 2;//jb
            if (BJ == num_cubes_y - 1) stopping_y = cube_size - 3;//je

            if (BK == 0) starting_z = 2;//kb
            if (BK == num_cubes_z - 1) stopping_z = cube_size - 3;//ke

            s1 = s2 = s3 = s4 = 0;

            for (li = starting_x; li <= stopping_x; ++li)
              for (lj = starting_y; lj <= stopping_y; ++lj)
                for (lk = starting_z; lk <= stopping_z; ++lk){
                  node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.

                  s1 = nodes[node_idx].df2[0];
                  s2 = c[0][0] * nodes[node_idx].df2[0];
                  s3 = c[0][1] * nodes[node_idx].df2[0];
                  s4 = c[0][2] * nodes[node_idx].df2[0];

                  for (ksi = 1; ksi <= 18; ksi++){
                    s1 += nodes[node_idx].df2[ksi];
                    s2 += c[ksi][0] * nodes[node_idx].df2[ksi];
                    s3 += c[ksi][1] * nodes[node_idx].df2[ksi];
                    s4 += c[ksi][2] * nodes[node_idx].df2[ksi];
                  }

                  nodes[node_idx].rho = s1; /* Eqn 11 from paper...PN*/
                  nodes[node_idx].vel_x = (s2 + 0.5 * gv->dt * nodes[node_idx].elastic_force_x) / s1;/*Eqn 12*/
                  nodes[node_idx].vel_y = (s3 + 0.5 * gv->dt * nodes[node_idx].elastic_force_y) / s1;/*Eqn 12*/
                  nodes[node_idx].vel_z = (s4 + 0.5 * gv->dt * (nodes[node_idx].elastic_force_z + nodes[node_idx].rho * gv->g_l)) / s1;

                }
        } //if machine check ends
    }//For BK

  //printf("**********compute_rho Exit**********\n");
}
