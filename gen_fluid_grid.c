/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

Fluidgrid*  gen_fluid_grid(int dim_x, int dim_y, int dim_z, int cube_size, GV gv){
  Fluidgrid     *fluid_grid;
  int            total_sub_grids, num_cubes_x, num_cubes_y, num_cubes_z;
  int            BI, BJ, BK;      //Thread Block index in whole fluid domain
  int            cube_idx;
  int            temp_mac_rank;


  /*Calculate total no of cubes*/
  total_sub_grids = (dim_x*dim_y*dim_z) / pow(cube_size, 3);
  num_cubes_x = gv->num_cubes_x;  //should be sub_grid_dimx for later
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  /* allocate memory for the entire fluid grid.*/
  fluid_grid = (Fluidgrid*)malloc(sizeof(Fluidgrid));
  fluid_grid->x_dim = dim_x;
  fluid_grid->y_dim = dim_y;
  fluid_grid->z_dim = dim_z;

  /*Allocate Memory for sub_fluid grid*/
  fluid_grid->sub_fluid_grid = (Sub_Fluidgrid*)malloc(sizeof(Sub_Fluidgrid)*total_sub_grids);

  /* allocate memory for the surface structure */
  fluid_grid->inlet = (Fluidsurface*)malloc(sizeof(Fluidsurface));
  fluid_grid->outlet = (Fluidsurface*)malloc(sizeof(Fluidsurface));

  // MPI change: only fluid machine will start malloc space
  if (gv->my_rank != gv->num_macs - 1){
    /* allocate memory for the surface's fluid nodes */
    fluid_grid->inlet->nodes = (Fluidnode*)calloc(dim_z * dim_y, sizeof(Fluidnode));
    fluid_grid->outlet->nodes = (Fluidnode*)calloc(dim_z * dim_y, sizeof(Fluidnode));

    /**/
    for (BI = 0; BI < num_cubes_x; ++BI)
    for (BJ = 0; BJ < num_cubes_y; ++BJ)
    for (BK = 0; BK < num_cubes_z; ++BK){
      cube2thread_and_machine(BI, BJ, BK, gv, &temp_mac_rank);//MPI changes
      if (gv->my_rank == temp_mac_rank){
        cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        fluid_grid->sub_fluid_grid[cube_idx].nodes =
          (Fluidnode*)calloc(cube_size*cube_size*cube_size, sizeof(Fluidnode));
      }
    }

  }

  return fluid_grid;
}
