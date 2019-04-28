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

void gen_fiber_sheet(Fibersheet* sheet){
  
  double width = sheet->width;
  double height = sheet->height; 
  int num_cols = sheet->num_cols;
  int num_rows = sheet->num_rows;;
  double x_orig = sheet->x_orig;
  double y_orig = sheet->y_orig;
  double z_orig = sheet->z_orig;

  int row, col;
  Fiber *fiber;

  double zdist = width  / (num_cols - 1);
  double ydist = height / (num_rows - 1);
  /* Fiber sheet is kept on yz plane with exactly midway between fluid grid */

  /* allocate memory for fibers in each sheet */

  // fiber_mac malloc nodes
  sheet->fibers = (Fiber*) calloc(num_rows, sizeof(Fiber));
  for (row = 0; row < num_rows; ++row) {
    fiber = sheet->fibers + row;       // The row-th fiber
    fiber->nodes = (Fibernode*) calloc(num_cols, sizeof(Fibernode));
    /* Fill in the contents of the points in each fiber */
    for (col = 0; col < num_cols; ++col){
      //If fiber nodes' coordinates are given by a user, remove 3 lines below.
      fiber->nodes[col].x = x_orig; // perpendicular to x-axis
      // fiber->nodes[col].x = x_orig + (row - 10)*0.4; // not perpendicular to x-axis, with some angle
      fiber->nodes[col].y = ydist * row + y_orig;
      fiber->nodes[col].z = zdist * col + z_orig;
#if 0      
      printf("(%d, %d)  (%f, %f, %f)\n", 
        row, col, fiber->nodes[col].x, fiber->nodes[col].y, fiber->nodes[col].z);
#endif        
    }
  }
  /*For Tethered Points*/
  for (row = num_rows - 1; row > num_rows - 3; --row) {
    for (col = 0; col < num_cols; ++col){
      //TODO: Tether points have user given position
      fiber->nodes[col].x = fiber->nodes[col].x;
      fiber->nodes[col].y = fiber->nodes[col].y;
      fiber->nodes[col].z = fiber->nodes[col].z;
    }
  }
}