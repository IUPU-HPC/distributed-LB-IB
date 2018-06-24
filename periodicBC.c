/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

void periodicBC(LV lv){
  /* use periodic boundary condition on y and z directions,as a result of this, we are actually
  *            doing simulation on a array of micro-channels */
  /*along y-direction,j=1 & n2-1 */
  //printf("***************Inside periodicBC********\n");
  int tid;
  GV gv = lv->gv;
  tid = lv->tid;

  int           ksi;
  Fluidgrid     *fluidgrid;
  Fluidnode     *nodes_first, *nodes_last;
  int           BI, BJ, BJ_first, BJ_last, BK, BK_first, BK_last;
  int           num_cubes_x, num_cubes_y, num_cubes_z, cube_size;
  int           li, lj, lk;
  int           cube_idx_first, cube_idx_last, node_idx_first, node_idx_last;
  int           starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;
  int my_rank = gv->my_rank;
  int temp_mac_rank;

  fluidgrid = gv->fluid_grid;

  cube_size = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  BJ_first = 0; BJ_last = num_cubes_y - 1;
  for (BI = 0; BI < num_cubes_x; ++BI)
  for (BK = 0; BK < num_cubes_z; ++BK){
    //if(cube2thread(BI, BJ_first, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R)==tid){
    if (cube2thread_and_machine(BI, BJ_first, BK, gv, &temp_mac_rank) == tid){
      if (my_rank == temp_mac_rank){
        cube_idx_first = BI * num_cubes_y * num_cubes_z + BJ_first * num_cubes_z + BK;
        nodes_first = fluidgrid->sub_fluid_grid[cube_idx_first].nodes;


        cube_idx_last = BI * num_cubes_y * num_cubes_z + BJ_last * num_cubes_z + BK;
        nodes_last = fluidgrid->sub_fluid_grid[cube_idx_last].nodes;


        starting_x = starting_z = 0;
        stopping_x = stopping_z = cube_size - 1;
        if (BI == 0) starting_x = 3;//ib+1
        if (BI == num_cubes_x - 1) stopping_x = cube_size - 4;//ie-1
        if (BK == 0) starting_z = gv->kb;//kb-1
        if (BK == num_cubes_z - 1) stopping_z = cube_size - 3;//ke

        for (li = starting_x; li <= stopping_x; ++li)
        for (lk = starting_z; lk <= stopping_z; ++lk){
          node_idx_first = li * cube_size * cube_size + (gv->jb - 1) * cube_size + lk;
          node_idx_last = li * cube_size * cube_size + (cube_size - 2) * cube_size + lk;//je+1
          for (ksi = 0; ksi <= 18; ksi++){
            nodes_first[node_idx_first].df2[ksi] = nodes_last[li*cube_size*cube_size + (cube_size - 3)*cube_size + lk].df2[ksi];//gv->je*dim_z +k
            nodes_last[node_idx_last].df2[ksi] = nodes_first[li*cube_size*cube_size + (gv->jb)*cube_size + lk].df2[ksi];//gv->jb*dim_z +k
          }//ksi
        }//lk
      }//if machine chek
    }//if cube 2thread
  }//BK

  /* along z-direction, z=1 & n3-1 */

  BK_first = 0; BK_last = num_cubes_z - 1;
  for (BI = 0; BI < num_cubes_x; ++BI)
  for (BJ = 0; BJ < num_cubes_y; ++BJ){
    //if(cube2thread(BI, BJ, BK_first, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) == tid){
    if (cube2thread_and_machine(BI, BJ, BK_first, gv, &temp_mac_rank) == tid){
      if (my_rank == temp_mac_rank){
        cube_idx_first = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK_first;
        cube_idx_last = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK_last;
        nodes_first = fluidgrid->sub_fluid_grid[cube_idx_first].nodes;
        nodes_last = fluidgrid->sub_fluid_grid[cube_idx_last].nodes;
        starting_x = starting_y = 0;
        stopping_x = stopping_y = cube_size - 1;

        if (BI == 0) starting_x = 3;//ib+1
        if (BI == num_cubes_x - 1) stopping_x = cube_size - 4;//ie-1
        if (BJ == 0) starting_y = gv->jb;//jb
        if (BJ == num_cubes_z - 1) stopping_y = cube_size - 3;//ke

        for (li = starting_x; li <= stopping_x; ++li)
        for (lj = starting_y; lj <= stopping_y; ++lj){
          node_idx_first = li * cube_size * cube_size + lj * cube_size + gv->kb - 1;
          node_idx_last = li * cube_size * cube_size + lj * cube_size + cube_size - 2;//ke+1
          for (ksi = 0; ksi <= 18; ksi++){
            nodes_first[node_idx_first].df2[ksi] = nodes_last[li*cube_size*cube_size + lj*cube_size + cube_size - 3].df2[ksi];//j*dim_z +gv->ke

            nodes_last[node_idx_last].df2[ksi] = nodes_first[li*cube_size*cube_size + lj*cube_size + gv->kb].df2[ksi];//j*dim_z +gv->kb

          }//ksi
        }//lj
      }//if machine chek
    }//if cube 2thread
  }//BJ

  // printf("***************periodicBC Exit********\n");
}
