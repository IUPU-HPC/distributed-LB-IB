/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

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
  long BI, BJ, BK; //to identify the Sub grids
  long cube_idx;
  long starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  long li, lj, lk, node_idx;//local access point inside cube

  long  total_sub_grids, dim_x, dim_y, dim_z;
  dim_x = gv->fluid_grid->x_dim;
  dim_y = gv->fluid_grid->y_dim;
  dim_z = gv->fluid_grid->z_dim;
  int cube_size = gv->cube_size;
  long num_cubes_x = gv->fluid_grid->num_cubes_x;
  long num_cubes_y = gv->fluid_grid->num_cubes_y;
  long num_cubes_z = gv->fluid_grid->num_cubes_z;
  total_sub_grids = (dim_x*dim_y*dim_z) / pow(cube_size, 3);

  int temp_mac_rank, my_rank;
  my_rank = gv->taskid;

  /*Pthread Changes*/ /*OPTIMISE using Loop unrolling*/
  for (BI = 0; BI < num_cubes_x; ++BI)
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
  for (BK = 0; BK < num_cubes_z; ++BK){
    if (cube2thread_and_task(BI, BJ, BK, gv, &temp_mac_rank) == tid){
      if (temp_mac_rank == my_rank){
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
          for (ksi = 0; ksi <= 18; ksi++){
            if (ksi == 0){
              nodes[node_idx].dfeq[ksi] =
                1.0 / 3.0 * nodes[node_idx].rho
                * (1.0 - 1.5 *
                (nodes[node_idx].vel_x * nodes[node_idx].vel_x
                + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                + nodes[node_idx].vel_z * nodes[node_idx].vel_z));

              nodes[node_idx].df1[ksi] =
                nodes[node_idx].df1[ksi] * (1.0 - 1.0 / gv->tau)
                + 1.0 / gv->tau
                *nodes[node_idx].dfeq[ksi]
                + gv->dt*(1.0 - 1.0 / (2.0*gv->tau)) / 3.0
                *(
                ((gv->c[ksi][0] - nodes[node_idx].vel_x)*nodes[node_idx].elastic_force_x
                + (gv->c[ksi][1] - nodes[node_idx].vel_y)*nodes[node_idx].elastic_force_y
                + (gv->c[ksi][2] - nodes[node_idx].vel_z)*(nodes[node_idx].elastic_force_z
                + nodes[node_idx].rho*gv->g_l)) / (gv->cs_l*gv->cs_l)
                + (gv->c[ksi][0] * nodes[node_idx].vel_x
                + gv->c[ksi][1] * nodes[node_idx].vel_y
                + gv->c[ksi][2] * nodes[node_idx].vel_z) / pow(gv->cs_l, 4)
                *(gv->c[ksi][0] * nodes[node_idx].elastic_force_x
                + gv->c[ksi][1] * nodes[node_idx].elastic_force_y
                + gv->c[ksi][2] * (nodes[node_idx].elastic_force_z
                + nodes[node_idx].rho*gv->g_l)));
            } //ksi==0
            else if (ksi >= 1 && ksi <= 6){
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
                + nodes[node_idx].vel_y *nodes[node_idx].vel_y
                + nodes[node_idx].vel_z * nodes[node_idx].vel_z));

              nodes[node_idx].df1[ksi]
                = nodes[node_idx].df1[ksi] * (1.0 - 1.0 / gv->tau) + 1.0 / gv->tau
                *nodes[node_idx].dfeq[ksi] + gv->dt*(1.0 - 1.0 / (2.0*gv->tau)) / 18.0
                *(
                ((gv->c[ksi][0] - nodes[node_idx].vel_x)*nodes[node_idx].elastic_force_x
                + (gv->c[ksi][1] - nodes[node_idx].vel_y)*nodes[node_idx].elastic_force_y
                + (gv->c[ksi][2] - nodes[node_idx].vel_z)*(nodes[node_idx].elastic_force_z
                + nodes[node_idx].rho*gv->g_l)) / (gv->cs_l*gv->cs_l)
                + (gv->c[ksi][0] * nodes[node_idx].vel_x
                + gv->c[ksi][1] * nodes[node_idx].vel_y
                + gv->c[ksi][2] * nodes[node_idx].vel_z) / pow(gv->cs_l, 4)
                *(gv->c[ksi][0] * nodes[node_idx].elastic_force_x
                + gv->c[ksi][1] * nodes[node_idx].elastic_force_y
                + gv->c[ksi][2] * (nodes[node_idx].elastic_force_z
                + nodes[node_idx].rho*gv->g_l)));
            } //ksi between 1 and 7 ends
            else{
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

              nodes[node_idx].df1[ksi] =
                nodes[node_idx].df1[ksi] * (1.0 - 1.0 / gv->tau)
                + 1.0 / gv->tau * nodes[node_idx].dfeq[ksi]
                + gv->dt*(1.0 - 1.0 / (2.0*gv->tau)) / 36.0
                * (
                ((gv->c[ksi][0] - nodes[node_idx].vel_x) * nodes[node_idx].elastic_force_x
                + (gv->c[ksi][1] - nodes[node_idx].vel_y) * nodes[node_idx].elastic_force_y
                + (gv->c[ksi][2] - nodes[node_idx].vel_z) * (nodes[node_idx].elastic_force_z
                + nodes[node_idx].rho*gv->g_l)) / (gv->cs_l*gv->cs_l)
                + (gv->c[ksi][0] * nodes[node_idx].vel_x
                + gv->c[ksi][1] * nodes[node_idx].vel_y
                + gv->c[ksi][2] * nodes[node_idx].vel_z) / pow(gv->cs_l, 4)
                * (gv->c[ksi][0] * nodes[node_idx].elastic_force_x
                + gv->c[ksi][1] * nodes[node_idx].elastic_force_y
                + gv->c[ksi][2] * (nodes[node_idx].elastic_force_z
                + nodes[node_idx].rho*gv->g_l)));
            }//ksi between 7 and 18 ends
          }//ksi loop
        }//lk loop
      }// if machine check
    }//if cube2thread ends
  }//For BK
#ifdef DEBUG_PRINT
  // printf("****************compute_eqlbrmdistrfuncDF1 EXIT my_rankis :%d ***************\n", gv->my_rank);
#endif //DEBUG_PRINT

}
