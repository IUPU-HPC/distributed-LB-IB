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

void print_fiber_sub_grid(GV gv, int start_y, int start_z,
                          int end_y, int end_z) {
  /*Assuming one fiber sheet!!*/
  Fiber     *fiber_array;
  int        i, j, k;
  Fibernode *node;

  for (k = 0; k < gv->fiber_shape->num_sheets; ++k){
    fiber_array = gv->fiber_shape->sheets[k].fibers;
    printf("Fiber_sheets[%d] (i,j): {cord_x, cord_y, cord_z} || SF_x, SF_y, SF_z || BF_X, BF_y, BF_z || EF_X, EF_y, EF_z\n", k);
    for (i = start_y; i <= end_y; ++i) {
      for (j = start_z; j <= end_z; ++j) {
        node = fiber_array[i].nodes + j;
        printf("(%2d,%2d):{%f,%f,%f} || %.6f,%.24f,%.24f || %.6f,%.24f,%.24f || %.6f,%.24f,%.24f\n",
          i, j, node->x, node->y, node->z,
          node->stretch_force_x, node->stretch_force_y, node->stretch_force_z,
          node->bend_force_x, node->bend_force_y, node->bend_force_z,
          node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
      }
      printf("\n");
    }
  }
}

// save one lattice population to disk
void save_fiber_sub_grid(GV gv, int start_y, int start_z,
                         int end_y, int end_z, char fName[]) {
  FILE* oFile = fopen(fName, "w");

  /*Assuming one fiber sheet!!*/
  Fiber     *fiber_array;
  int        i, j, k;
  Fibernode *node;
  for (k = 0; k < gv->fiber_shape->num_sheets; ++k){
    fiber_array = gv->fiber_shape->sheets[k].fibers;
    fprintf(oFile, "Fiber_sheets[%d] (i,j): {cord_x, cord_y, cord_z} || SF_x, SF_y, SF_z || BF_x, BF_y, BF_z || EF_x, EF_y, EF_z\n", k);
    for (i = start_y; i <= end_y; ++i) {
      for (j = start_z; j <= end_z; ++j) {
        node = fiber_array[i].nodes + j;
        fprintf(oFile, "(%2d,%2d):{%f,%f,%f} || %.6f,%.24f,%.24f || %.6f,%.24f,%.24f || %.6f,%.24f,%.24f\n",
          i, j, node->x, node->y, node->z,
          node->stretch_force_x, node->stretch_force_y, node->stretch_force_z,
          node->bend_force_x, node->bend_force_y, node->bend_force_z,
          node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
      }
    }
  }
  fclose(oFile);
}

// Description: print all the fluid nodes info within all the cubes pointed by fluid coordinate
// Input: fluid coordinate
void print_fluid_sub_grid(GV gv, int start_x, int start_y, int start_z,
                          int end_x, int end_y, int end_z) {// start, stop are closed: inclusive

  Fluidgrid *grid;
  Sub_Fluidgrid *sub_grid;
  int li, lj, lk, BI, BJ, BK, BI_start, BJ_start, BK_start, BI_end, BJ_end, BK_end;
  int ksi, cube_idx, num_cubes_y, num_cubes_z;
  int temp_taskid;
  Fluidnode *node;
  int cube_size = gv->cube_size;

  grid = gv->fluid_grid;
  sub_grid = grid->sub_fluid_grid;
  num_cubes_y = gv->fluid_grid->num_cubes_y;
  num_cubes_z = gv->fluid_grid->num_cubes_z;

  BI_start = start_x / cube_size;
  BJ_start = start_y / cube_size;
  BK_start = start_z / cube_size;

  BI_end = end_x / cube_size;
  BJ_end = end_y / cube_size;
  BK_end = end_z / cube_size;

  /*li     = start_x%cube_size;
  lj     = start_y%cube_size;
  lk     = start_z%cube_size;
  li_end = end_x%cube_size;
  lj_end = end_y%cube_size;
  lk_end = end_z%cube_size;*/

  // printf("Rank %d Enter print_fluid_sub_grid\n", gv->taskid);

  for (BI = BI_start; BI <= BI_end; ++BI)
    for (BJ = BJ_start; BJ <= BJ_end; ++BJ)
      for (BK = BK_start; BK <= BK_end; ++BK) {
        temp_taskid = cube2task(BI, BJ, BK, gv);
        cube2task(BI, BJ, BK, gv);
        if (gv->taskid == temp_taskid){ //MPI changes

          // printf("temp_taskid = %d\n", temp_taskid);
          printf("(BI,BJ,BK): {vel_x, vel_y, vel_z} || {G0, DF1, DF2}|| rho || {ElasticF_x, y, z}\n");

          cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          for (li = 0; li < cube_size; li++)
            for (lj = 0; lj < cube_size; lj++)
              for (lk = 0; lk < cube_size; lk++){
                node = &sub_grid[cube_idx].nodes[li*cube_size*cube_size + lj*cube_size + lk];
                //if( li == start_x%cube_size && lj ==start_y%cube_size && lk ==start_z%cube_size ){
                printf("For cube <%d,%d,%d> Rank %d\n", BI, BJ, BK, gv->taskid);
                for (ksi = 0; ksi < 19; ksi++){
                  printf("Rank-%d- (%d,%d,%d, %d):{%.12f,%.12f,%.12f} || {%.12f,%.12f,%.12f} || %.12f || {%.12f,%.12f,%.12f}\n",
                    gv->taskid, li, lj, lk, ksi, node->vel_x, node->vel_y, node->vel_z,
                    node->dfeq[ksi], node->df1[ksi], node->df2[ksi],
                    node->rho,
                    node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
                }
          }
        }
  }

  printf("\n");

}

// Input: cube block coordinate
void print_fluid_cube(GV gv, int BI_start, int BJ_start, int BK_start,
  int BI_end, int BJ_end, int BK_end) {// start, stop are closed: inclusive

  Fluidgrid *grid;
  Sub_Fluidgrid *sub_grid;
  int li, lj, lk, BI, BJ, BK;
  int ksi, cube_idx, num_cubes_y, num_cubes_z;
  int temp_taskid;
  Fluidnode *node;
  int cube_size = gv->cube_size;

  grid = gv->fluid_grid;
  sub_grid = grid->sub_fluid_grid;
  num_cubes_y = gv->fluid_grid->num_cubes_y;
  num_cubes_z = gv->fluid_grid->num_cubes_z;

  // printf("Rank %d Enter print_fluid_cube\n", gv->taskid);

  for (BI = BI_start; BI <= BI_end; ++BI)
    for (BJ = BJ_start; BJ <= BJ_end; ++BJ)
      for (BK = BK_start; BK <= BK_end; ++BK) {
      temp_taskid = cube2task(BI, BJ, BK, gv);
      if (gv->taskid == temp_taskid){ //MPI changes
        printf("(BI,BJ,BK): {vel_x, vel_y, vel_z} || {G0, DF1, DF2}|| rho || {ElasticF_x, y, z}\n");

        cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        for (li = 0; li < cube_size; li++)
          for (lj = 0; lj < cube_size; lj++)
            for (lk = 0; lk < cube_size; lk++){
              node = &sub_grid[cube_idx].nodes[li*cube_size*cube_size + lj*cube_size + lk];
              printf("For cube <%d,%d,%d>, Rank %d\n", BI, BJ, BK, gv->taskid);
              for (ksi = 0; ksi < 19; ksi++){
                printf("Rank-%d- (%d,%d,%d, %d):{%.12f,%.12f,%.12f} || {%.12f,%.12f,%.12f} || %.12f || {%.12f,%.12f,%.12f}\n",
                  gv->taskid, li, lj, lk, ksi, node->vel_x, node->vel_y, node->vel_z,
                  node->dfeq[ksi], node->df1[ksi], node->df2[ksi],
                  node->rho,
                  node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
              }
            }
      }
    }
    printf("\n");
}

// save one subgrid of fluid to disk
void save_fluid_sub_grid(GV gv, int start_x, int start_y, int start_z,
                         int end_x, int end_y, int end_z, char fName[]) {
  FILE* oFile = fopen(fName, "w");

  Fluidgrid *grid;
  Sub_Fluidgrid *sub_grid;
  int li, lj, lk, BI, BJ, BK, BI_start, BJ_start, BK_start, BI_end, BJ_end, BK_end;
  int ksi, cube_idx, num_cubes_y, num_cubes_z;
  int temp_taskid;
  Fluidnode *node;
  int cube_size = gv->cube_size;

  grid = gv->fluid_grid;
  sub_grid = grid->sub_fluid_grid;
  num_cubes_y = gv->fluid_grid->num_cubes_y;
  num_cubes_z = gv->fluid_grid->num_cubes_z;

  BI_start = start_x / cube_size;
  BJ_start = start_y / cube_size;
  BK_start = start_z / cube_size;

  BI_end = end_x / cube_size;
  BJ_end = end_y / cube_size;
  BK_end = end_z / cube_size;

  /*li     = start_x%cube_size;
  lj     = start_y%cube_size;
  lk     = start_z%cube_size;
  li_end = end_x%cube_size;
  lj_end = end_y%cube_size;
  lk_end = end_z%cube_size;*/

  // printf("Rank %d Enter print_fluid_sub_grid\n", gv->taskid);

  for (BI = BI_start; BI <= BI_end; ++BI)
    for (BJ = BJ_start; BJ <= BJ_end; ++BJ)
      for (BK = BK_start; BK <= BK_end; ++BK) {
        temp_taskid = cube2task(BI, BJ, BK, gv);
        cube2task(BI, BJ, BK, gv);
        if (gv->taskid == temp_taskid){ //MPI changes

          // printf("temp_taskid = %d\n", temp_taskid);
          fprintf(oFile, "(BI,BJ,BK): {vel_x, vel_y, vel_z} || {G0, DF1, DF2}|| rho || {ElasticF_x, y, z}\n");

          cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          for (li = 0; li < cube_size; li++)
            for (lj = 0; lj < cube_size; lj++)
              for (lk = 0; lk < cube_size; lk++){
                node = &sub_grid[cube_idx].nodes[li*cube_size*cube_size + lj*cube_size + lk];
                //if( li == start_x%cube_size && lj ==start_y%cube_size && lk ==start_z%cube_size ){
                fprintf(oFile, "For cube <%d,%d,%d> Rank %d\n", BI, BJ, BK, gv->taskid);
                int X = BI * cube_size + li;
                int Y = BJ * cube_size + lj;
                int Z = BK * cube_size + lk;
                for (ksi = 0; ksi < 19; ksi++){
                  fprintf(oFile, "(%d,%d,%d, %2d):{%.6f,%.6f,%.6f} || {%.24f,%.24f,%.6f} || %.6f || {%.6f,%.6f,%.6f}\n",
                    X, Y, Z, ksi, node->vel_x, node->vel_y, node->vel_z,
                    node->dfeq[ksi], node->df1[ksi], node->df2[ksi],
                    node->rho,
                    node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
                }
              }
        }
  }

  fclose(oFile);
}