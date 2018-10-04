/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

/* Assming ONLY one sheet at this moment */
Fibershape*  gen_fiber_shape(double w, double h, int num_cols, int num_rows,
  double x_orig, double y_orig, double z_orig, GV gv){
  int         i, row, col;
  Fibershape  *shape;
  Fibersheet  *sheet;
  Fiber       *fiber;
  int         fiber_mac_rank = gv->num_macs - 1;

  double      zdist = w / (num_cols - 1);
  double      ydist = h / (num_rows - 1);
  /*Fiber sheet is kept on yz plane with exactly midway between fluid grid , therfore added +30 */

  /* allocate memory for fiber_shape.*/
  shape = (Fibershape*)malloc(sizeof(Fibershape));

  // allocate memory for fiber sheets. assumption: Only One Sheet For Now!!!!!
  shape->num_sheets = 1;
  shape->sheets = (Fibersheet*)malloc(sizeof(Fibersheet)* shape->num_sheets);


  /* allocate memory for fibers in each sheet */
  for (i = 0; i < shape->num_sheets; ++i) {
    sheet = shape->sheets + i;  //The i-th fiber sheet
    sheet->num_cols = num_cols;
    sheet->num_rows = num_rows;
    sheet->width = w;
    sheet->height = h;
    sheet->x_orig = x_orig;
    sheet->y_orig = y_orig;
    sheet->z_orig = z_orig;

    // fiber_mac malloc nodes
    if (gv->my_rank == fiber_mac_rank){
      sheet->fibers = (Fiber*)malloc(num_rows * sizeof(*sheet->fibers));
      for (row = 0; row < num_rows; ++row) {
        fiber = sheet->fibers + row;       // The row-th fiber
        fiber->nodes = (Fibernode*)malloc(num_cols * sizeof(*fiber->nodes));
        /* Fill in the contents of the points in each fiber */
        for (col = 0; col < num_cols; ++col){
          /*z-coordinate of sheet is from (center of z-direction - sheet height) to (center of z-direction);
          y-coordinate of sheet is  from (center of y-direction - half of sheet width) to (center of y-direction+half of sheet width) */

          // fiber->nodes[col].x = x_orig;               //perpendicular to x-axis.
          fiber->nodes[col].x = x_orig + (row - 10)*0.4;
          fiber->nodes[col].y = ydist * row + y_orig; //ydist * i + 30;
          fiber->nodes[col].z = zdist * col + z_orig; //zdist * j + 30;// middle of the fluid grid
          //If fiber nodes' coordinates are given by a user, remove the above 3 lines.
        }
      }
      /*For Tethered Points*/
      for (row = num_rows - 1; row > num_rows - 3; --row) {
        for (col = 0; col < num_cols; ++col){
          //TODO: Tether points
          fiber->nodes[col].x = fiber->nodes[col].x;
          fiber->nodes[col].y = fiber->nodes[col].y;
          fiber->nodes[col].z = fiber->nodes[col].z;
        }
      }
    } // endif fiber_mac malloc
  }
  return shape;
}
