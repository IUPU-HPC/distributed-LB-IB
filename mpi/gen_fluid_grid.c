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
  long BI, BJ, BK;      //Thread Block index in whole fluid domain
  long cube_idx;
  int temp_mac_rank;

  /*Calculate total no of cubes*/
  total_sub_grids = (fluid_grid->dim_x * fluid_grid->dim_y * fluid_grid->dim_z) / pow(cube_size, 3);
  long num_cubes_x = fluid_grid->dim_x / cube_size;
  long num_cubes_y = fluid_grid->dim_y / cube_size;
  long num_cubes_z = fluid_grid->dim_z / cube_size;

  /*Allocate Memory for sub_fluid grid*/
  fluid_grid->sub_fluid_grid = (Sub_Fluidgrid*) malloc(sizeof(Sub_Fluidgrid) * total_sub_grids);

  /* allocate memory for the surface structure */
  fluid_grid->inlet  = (Fluidsurface*) malloc(sizeof(Fluidsurface));
  fluid_grid->outlet = (Fluidsurface*) malloc(sizeof(Fluidsurface));

  // MPI change: only fluid machine will start malloc space
  /* allocate memory for the surface's fluid nodes */
  fluid_grid->inlet->nodes  = (Fluidnode*)calloc(fluid_grid->dim_z * fluid_grid->dim_y, sizeof(Fluidnode));
  fluid_grid->outlet->nodes = (Fluidnode*)calloc(fluid_grid->dim_z * fluid_grid->dim_y, sizeof(Fluidnode));

  for (BI = 0; BI < num_cubes_x; ++BI)
    for (BJ = 0; BJ < num_cubes_y; ++BJ)
      for (BK = 0; BK < num_cubes_z; ++BK){
        cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank); //MPI changes
        if (taskid == temp_mac_rank){
          cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          fluid_grid->sub_fluid_grid[cube_idx].nodes =
            (Fluidnode*)calloc(cube_size * cube_size * cube_size, sizeof(Fluidnode));
        }
  }

}