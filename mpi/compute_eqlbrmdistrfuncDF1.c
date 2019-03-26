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

extern int c[19][3];
extern double t[19];

// initialize a node to its local equilibrium term
void computeEquilibrium(Fluidnode* node) {
  int iPop;
  double rho = node->rho;
  double ux = node->vel_x;
  double uy = node->vel_y;
  double uz = node->vel_z;

  double uSqr = ux*ux + uy*uy + uz*uz;
  for (iPop = 0; iPop < 19; ++iPop) {
      node->dfeq[iPop] =
          computeEquilibrium(iPop, rho, ux, uy, uz, uSqr);
  }
}

void computeDF1(GV gv, Fluidnode* node){
  int iPop;
  double rho = node->rho;
  double ux = node->vel_x;
  double uy = node->vel_y;
  double uz = node->vel_z;
  double fx = node->elastic_force_x;
  double fy = node->elastic_force_y;
  double fz = node->elastic_force_z;

  // double uSqr = ux*ux + uy*uy + uz*uz;
  double c_u;

  for (iPop = 0; iPop < 19; ++iPop) {
    c_u = c[iPop][0]*ux + c[iPop][1]*uy + c[iPop][2]*uz;

    node->df1[iPop] =
          node->df1[iPop] * (1.0 - 1.0/ gv->tau)
        + 1.0 / gv->tau * node->dfeq[iPop]
        + gv->dt * (1.0 - 1.0 / (2.0 * gv->tau)) * t[iPop]
        * (
            (   (c[iPop][0] - ux) * fx
              + (c[iPop][1] - uy) * fy
              + (c[iPop][2] - uz) * (fz + rho * gv->g_l)
            ) / (gv->cs_l*gv->cs_l)
            +
            c_u / pow(gv->cs_l, 4) 
            * (c[iPop][0]*fx + c[iPop][1]*fy + c[iPop][2]*(fz + rho * gv->g_l))
          );
  }
}

void compute_eqlbrmdistrfuncDF1(LV lv){
#ifdef DEBUG_PRINT
  // printf("****************Inside compute_eqlbrmdistrfuncDF1  *************\n");
#endif //DEBUG_PRINT
  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

  Fluidgrid* fluidgrid = gv->fluid_grid;
  Fluidnode* nodes;

  int ksi;
  int BI, BJ, BK; //to identify the Sub grids
  long cube_idx;
  int starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  int li, lj, lk, node_idx;//local access point inside cube

  int  total_sub_grids, dim_x, dim_y, dim_z;
  dim_x = gv->fluid_grid->x_dim;
  dim_y = gv->fluid_grid->y_dim;
  dim_z = gv->fluid_grid->z_dim;
  int cube_size = gv->cube_size;
  int num_cubes_x = gv->fluid_grid->num_cubes_x;
  int num_cubes_y = gv->fluid_grid->num_cubes_y;
  int num_cubes_z = gv->fluid_grid->num_cubes_z;
  total_sub_grids = (dim_x*dim_y*dim_z) / pow(cube_size, 3);

  int toProc, my_rank, owner_tid;
  my_rank = gv->taskid;

  /*Pthread Changes*/ /*OPTIMISE using Loop unrolling*/
  for (BI = 0; BI < num_cubes_x; ++BI)
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
  for (BK = 0; BK < num_cubes_z; ++BK){
    owner_tid = cube2thread_and_task(BI, BJ, BK, gv, &toProc);

    if (tid == owner_tid && my_rank == toProc){

    cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
    nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
    starting_x = starting_y = starting_z = 0;
    stopping_x = stopping_y = stopping_z = cube_size - 1;

    if (BI == 0) starting_x = 2;//ib
    if (BI == num_cubes_x - 1) stopping_x = cube_size - 3;//ie
    if (BJ == 0) starting_y = 2;//jb
    if (BJ == num_cubes_y - 1) stopping_y = cube_size - 3;//je
    if (BK == 0) starting_z = 2;//kb
    if (BK == num_cubes_z - 1) stopping_z = cube_size - 3;//ke

    for (li = starting_x; li <= stopping_x; ++li)
    for (lj = starting_y; lj <= stopping_y; ++lj)
    for (lk = starting_z; lk <= stopping_z; ++lk){
      node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
      computeEquilibrium(nodes+node_idx);
      computeDF1(gv, nodes+node_idx);

      // for (ksi = 0; ksi <= 18; ksi++){
      //   if (ksi == 0){
      //     nodes[node_idx].dfeq[ksi] =
      //       1.0 / 3.0 * nodes[node_idx].rho
      //       * (1.0 - 1.5 *
      //       (nodes[node_idx].vel_x * nodes[node_idx].vel_x
      //       + nodes[node_idx].vel_y * nodes[node_idx].vel_y
      //       + nodes[node_idx].vel_z * nodes[node_idx].vel_z));

      //     nodes[node_idx].df1[ksi] =
      //       nodes[node_idx].df1[ksi] * (1.0 - 1.0 / gv->tau)
      //       + 1.0 / gv->tau
      //       *nodes[node_idx].dfeq[ksi]
      //       + gv->dt*(1.0 - 1.0 / (2.0*gv->tau)) / 3.0
      //       *(
      //       ((gv->c[ksi][0] - nodes[node_idx].vel_x)*nodes[node_idx].elastic_force_x
      //       + (gv->c[ksi][1] - nodes[node_idx].vel_y)*nodes[node_idx].elastic_force_y
      //       + (gv->c[ksi][2] - nodes[node_idx].vel_z)*(nodes[node_idx].elastic_force_z
      //       + nodes[node_idx].rho*gv->g_l)) / (gv->cs_l*gv->cs_l)
      //       + (gv->c[ksi][0] * nodes[node_idx].vel_x
      //       + gv->c[ksi][1] * nodes[node_idx].vel_y
      //       + gv->c[ksi][2] * nodes[node_idx].vel_z) / pow(gv->cs_l, 4)
      //       *(gv->c[ksi][0] * nodes[node_idx].elastic_force_x
      //       + gv->c[ksi][1] * nodes[node_idx].elastic_force_y
      //       + gv->c[ksi][2] * (nodes[node_idx].elastic_force_z
      //       + nodes[node_idx].rho*gv->g_l)));
      //   } //ksi==0
      //   else if (ksi >= 1 && ksi <= 6){
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
      //       + nodes[node_idx].vel_y *nodes[node_idx].vel_y
      //       + nodes[node_idx].vel_z * nodes[node_idx].vel_z));

      //     nodes[node_idx].df1[ksi]
      //       = nodes[node_idx].df1[ksi] * (1.0 - 1.0 / gv->tau) + 1.0 / gv->tau
      //       *nodes[node_idx].dfeq[ksi] + gv->dt*(1.0 - 1.0 / (2.0*gv->tau)) / 18.0
      //       *(
      //       ((gv->c[ksi][0] - nodes[node_idx].vel_x)*nodes[node_idx].elastic_force_x
      //       + (gv->c[ksi][1] - nodes[node_idx].vel_y)*nodes[node_idx].elastic_force_y
      //       + (gv->c[ksi][2] - nodes[node_idx].vel_z)*(nodes[node_idx].elastic_force_z
      //       + nodes[node_idx].rho*gv->g_l)) / (gv->cs_l*gv->cs_l)
      //       + (gv->c[ksi][0] * nodes[node_idx].vel_x
      //       + gv->c[ksi][1] * nodes[node_idx].vel_y
      //       + gv->c[ksi][2] * nodes[node_idx].vel_z) / pow(gv->cs_l, 4)
      //       *(gv->c[ksi][0] * nodes[node_idx].elastic_force_x
      //       + gv->c[ksi][1] * nodes[node_idx].elastic_force_y
      //       + gv->c[ksi][2] * (nodes[node_idx].elastic_force_z
      //       + nodes[node_idx].rho*gv->g_l)));
      //   } //ksi between 1 and 7 ends
      //   else{
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

      //     nodes[node_idx].df1[ksi] =
      //       nodes[node_idx].df1[ksi] * (1.0 - 1.0 / gv->tau)
      //       + 1.0 / gv->tau * nodes[node_idx].dfeq[ksi]
      //       + gv->dt*(1.0 - 1.0 / (2.0*gv->tau)) / 36.0
      //       * (
      //       ((gv->c[ksi][0] - nodes[node_idx].vel_x) * nodes[node_idx].elastic_force_x
      //       + (gv->c[ksi][1] - nodes[node_idx].vel_y) * nodes[node_idx].elastic_force_y
      //       + (gv->c[ksi][2] - nodes[node_idx].vel_z) * (nodes[node_idx].elastic_force_z
      //       + nodes[node_idx].rho*gv->g_l)) / (gv->cs_l*gv->cs_l)
      //       + (gv->c[ksi][0] * nodes[node_idx].vel_x
      //       + gv->c[ksi][1] * nodes[node_idx].vel_y
      //       + gv->c[ksi][2] * nodes[node_idx].vel_z) / pow(gv->cs_l, 4)
      //       * (gv->c[ksi][0] * nodes[node_idx].elastic_force_x
      //       + gv->c[ksi][1] * nodes[node_idx].elastic_force_y
      //       + gv->c[ksi][2] * (nodes[node_idx].elastic_force_z
      //       + nodes[node_idx].rho*gv->g_l)));
      //   }//ksi between 7 and 18 ends
      // }//ksi loop

        }//lk loop
      }// if machine check
  }//For BK

#if 0
  printf("****************compute_eqlbrmdistrfuncDF1 EXIT my_rankis :%d ***************\n", gv->my_rank);
#endif //DEBUG_PRINT

}
