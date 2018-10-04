/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

void print_fiber_sub_grid(GV gv, int start_y, int start_z,
  int end_y, int end_z) {
  /*Assuming one fiber sheet!!*/
  Fiber     *fiber_array;
  int        i, j;
  Fibernode *node;

  fiber_array = gv->fiber_shape->sheets[0].fibers;

  printf("(i,j): {cord_x, cord_y, cord_z} || SF_x,SF_y, SF_z || BF_X, BF_y, BF_z\n");
  for (i = start_y; i <= end_y; ++i) {
    for (j = start_z; j <= end_z; ++j) {
      node = fiber_array[i].nodes + j;
      printf("(%d, %d):{%f,%f,%f} || %.12f,%.12f,%.12f || %.12f,%.12f,%.12f\n",
        i, j, node->x, node->y, node->z,
        node->stretch_force_x, node->stretch_force_y, node->stretch_force_z,
        node->bend_force_x, node->bend_force_y, node->bend_force_z);
    }
    printf("\n");
  }

}

// Description: print all the fluid nodes info within all the cubes pointed by fluid coordinate
// Input: fluid coordinate
void print_fluid_sub_grid(GV gv, int start_x, int start_y, int start_z,
  int end_x, int end_y, int end_z,
  int cube_size) {// start, stop are closed: inclusive

  Fluidgrid *grid;
  Sub_Fluidgrid *sub_grid;
  int        li, lj, lk, BI, BJ, BK, BI_start, BJ_start, BK_start, BI_end, BJ_end, BK_end;
  int        ksi, cube_idx, num_cubes_y, num_cubes_z;
  int        temp_mac_rank;
  Fluidnode *node;

  grid = gv->fluid_grid;
  sub_grid = grid->sub_fluid_grid;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

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

  // printf("Rank %d Enter print_fluid_sub_grid\n", gv->my_rank);

  for (BI = BI_start; BI <= BI_end; ++BI)
  for (BJ = BJ_start; BJ <= BJ_end; ++BJ)
  for (BK = BK_start; BK <= BK_end; ++BK) {
    cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank);
    if (gv->my_rank == temp_mac_rank){ //MPI changes

      // printf("temp_mac_rank = %d\n", temp_mac_rank);
      printf("(BI,BJ,BK): {vel_x, vel_y, vel_z} ||  {G0, DF1, DF2}|| rho || {ElasticF_x, y, z}\n");

      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      for (li = 0; li < cube_size; li++)
      for (lj = 0; lj < cube_size; lj++)
      for (lk = 0; lk < cube_size; lk++){
        node = &sub_grid[cube_idx].nodes[li*cube_size*cube_size + lj*cube_size + lk];
        //if( li == start_x%cube_size && lj ==start_y%cube_size && lk ==start_z%cube_size ){
        printf("For cube <%d,%d,%d> Rank %d\n", BI, BJ, BK, gv->my_rank);
        for (ksi = 0; ksi < 19; ksi++){
          printf("Rank-%d- (%d,%d,%d, %d):{%.12f,%.12f,%.12f} || {%.12f,  %.12f,%.12f} || %.12f || {%.12f,%.12f,%.12f} \n",
            gv->my_rank, li, lj, lk, ksi, node->vel_x, node->vel_y, node->vel_z,
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
  int BI_end, int BJ_end, int BK_end,
  int cube_size) {// start, stop are closed: inclusive

  Fluidgrid *grid;
  Sub_Fluidgrid *sub_grid;
  int        li, lj, lk, BI, BJ, BK;
  int        ksi, cube_idx, num_cubes_y, num_cubes_z;
  int        temp_mac_rank;
  Fluidnode *node;

  grid = gv->fluid_grid;
  sub_grid = grid->sub_fluid_grid;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  // printf("Rank %d Enter print_fluid_cube\n", gv->my_rank);


  for (BI = BI_start; BI <= BI_end; ++BI)
  for (BJ = BJ_start; BJ <= BJ_end; ++BJ)
  for (BK = BK_start; BK <= BK_end; ++BK) {
    cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank);
    if (gv->my_rank == temp_mac_rank){ //MPI changes
      printf("(BI,BJ,BK): {vel_x, vel_y, vel_z} ||  {G0, DF1, DF2}|| rho || {ElasticF_x, y, z}\n");

      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      for (li = 0; li < cube_size; li++)
      for (lj = 0; lj < cube_size; lj++)
      for (lk = 0; lk < cube_size; lk++){
        node = &sub_grid[cube_idx].nodes[li*cube_size*cube_size + lj*cube_size + lk];
        printf("For cube <%d,%d,%d>, Rank %d\n", BI, BJ, BK, gv->my_rank);
        for (ksi = 0; ksi < 19; ksi++){
          printf("Rank-%d- (%d,%d,%d, %d):{%.12f,%.12f,%.12f} || {%.12f,  %.12f,%.12f} || %.12f || {%.12f,%.12f,%.12f} \n",
            gv->my_rank, li, lj, lk, ksi, node->vel_x, node->vel_y, node->vel_z,
            node->dfeq[ksi], node->df1[ksi], node->df2[ksi],
            node->rho,
            node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
        }
      }
    }

  }
  printf("\n");

}


double get_cur_time() {
  struct timeval   tv;
  struct timezone  tz;
  double cur_time;

  gettimeofday(&tv, &tz);
  cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

  return cur_time;
}
