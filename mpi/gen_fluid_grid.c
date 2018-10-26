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

void gen_fluid_grid(Fluidgrid *fluid_grid, int cube_size, int taskid, GV gv){
  long total_sub_grids;
  long dim_x = fluid_grid->x_dim;
  long dim_y = fluid_grid->y_dim;
  long dim_z = fluid_grid->z_dim;

  /*Calculate total no of cubes*/
  total_sub_grids = (dim_x * dim_y * dim_z) / pow(cube_size, 3);
  long num_cubes_x = gv->fluid_grid->num_cubes_x = dim_x / cube_size;
  long num_cubes_y = gv->fluid_grid->num_cubes_y = dim_y / cube_size;
  long num_cubes_z = gv->fluid_grid->num_cubes_z = dim_z / cube_size;

  /*Allocate Memory for sub_fluid grid*/
  fluid_grid->sub_fluid_grid = (Sub_Fluidgrid*) malloc(sizeof(Sub_Fluidgrid) * total_sub_grids);

  /* allocate memory for the surface structure */
  fluid_grid->inlet  = (Fluidsurface*) malloc(sizeof(Fluidsurface));
  fluid_grid->outlet = (Fluidsurface*) malloc(sizeof(Fluidsurface));

  // MPI change: only fluid machine will start malloc space
  /* allocate memory for the surface's fluid nodes */
  fluid_grid->inlet->nodes  = (Fluidnode*)calloc(dim_y * dim_z, sizeof(Fluidnode));
  fluid_grid->outlet->nodes = (Fluidnode*)calloc(dim_y * dim_z, sizeof(Fluidnode));

  for (long BI = 0; BI < num_cubes_x; ++BI)
    for (long BJ = 0; BJ < num_cubes_y; ++BJ)
      for (long BK = 0; BK < num_cubes_z; ++BK){
        int tmp_taskid = cube2task(BI, BJ, BK, gv); //MPI changes
        if (taskid == tmp_taskid){
          long cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          Fluidnode *nodes = fluid_grid->sub_fluid_grid[cube_idx].nodes =
            (Fluidnode*)calloc(cube_size * cube_size * cube_size, sizeof(Fluidnode));
          for (long li = 0; li < cube_size; ++li)
            for (long lj = 0; lj < cube_size; ++lj)
              for (long lk = 0; lk < cube_size; ++lk){
                long node_idx = li * cube_size * cube_size + lj * cube_size + lk;
                nodes[node_idx].rho = gv->rho_l;//same as rho_l
                nodes[node_idx].vel_x = gv->u_l;
                nodes[node_idx].vel_y = 0.0;
                nodes[node_idx].vel_z = 0.0;
              }
        }
  }
}