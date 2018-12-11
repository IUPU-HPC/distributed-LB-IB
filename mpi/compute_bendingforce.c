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

void  compute_bendingforce(LV lv) {
#ifdef DEBUG_PRINT
  // printf("****Inside compute_bendingforce*******\n");
#endif //DEBUG_PRINT
  int         i, j;
  double      ds1, ds2;
  double      bending_const;
  Fibershape *fiber_shape; //aka IB shape
  Fiber*      fibers;      //an array of fibers
  double      fiber_sheet_w, fiber_sheet_h;
  int         total_fibers_clmn;
  int         total_fibers_row;
  Fibernode  *nodes;
  int         tid;
  int         total_threads;

  GV gv = lv->gv;
  tid = lv->tid;

  total_threads = gv->threads_per_task;
  fiber_shape = gv->fiber_shape;
  /* Assuming one sheet for now */
  fibers = fiber_shape->sheets[0].fibers;
  total_fibers_row = fiber_shape->sheets[0].num_rows;
  total_fibers_clmn = fiber_shape->sheets[0].num_cols;
  fiber_sheet_w = fiber_shape->sheets[0].width;  //z-direction
  fiber_sheet_h = fiber_shape->sheets[0].height; //y-direction

  ds1 = fiber_sheet_w / (total_fibers_clmn - 1); //gap between fiber nodes along the z-dim.
  ds2 = fiber_sheet_h / (total_fibers_row - 1); //gap between two neighboring fibers along y-dim.
  bending_const = gv->Kb_l / pow(ds1, 3);  //KB and alpha1 pow 3

  // printf("fiber_sheet_w = %f, fiber_sheet_h = %f, ds1 = %f, ds2 = %f, bending_const = %f\n", 
  //   fiber_sheet_w, fiber_sheet_h, ds1, ds2, bending_const);

  /* compute force along the fibers (row)*/
  for (i = 0; i < total_fibers_row; ++i) {
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){

      nodes = fibers[i].nodes; //an array of fiber nodes in the i-th fiber ???
      //For the 1st point
      nodes[0].bend_force_x = bending_const * (-nodes[2].x + 2 * nodes[1].x - nodes[0].x);
      nodes[0].bend_force_y = bending_const * (-nodes[2].y + 2 * nodes[1].y - nodes[0].y);
      nodes[0].bend_force_z = bending_const * (-nodes[2].z + 2 * nodes[1].z - nodes[0].z);

      //For the 2nd Point
      nodes[1].bend_force_x = bending_const * (-nodes[3].x + 4 * nodes[2].x - 5 * nodes[1].x +
        2 * nodes[0].x);
      nodes[1].bend_force_y = bending_const * (-nodes[3].y + 4 * nodes[2].y - 5 * nodes[1].y +
        2 * nodes[0].y);
      nodes[1].bend_force_z = bending_const * (-nodes[3].z + 4 * nodes[2].z - 5 * nodes[1].z +
        2 * nodes[0].z);

      //For the middle points
      for (j = 2; j < total_fibers_clmn - 2; ++j){
        nodes[j].bend_force_x = bending_const * (-nodes[j+2].x + 4.0 * nodes[j+1].x
          - 6.0 * nodes[j].x + 4.0 * nodes[j-1].x
          - nodes[j-2].x);
        nodes[j].bend_force_y = bending_const * (-nodes[j+2].y + 4.0 * nodes[j+1].y
          - 6.0 * nodes[j].y + 4.0 * nodes[j-1].y
          - nodes[j-2].y);
        nodes[j].bend_force_z = bending_const * (-nodes[j+2].z + 4.0 * nodes[j+1].z
          - 6.0 * nodes[j].z + 4.0 * nodes[j-1].z
          - nodes[j-2].z);
      }

      //For the last but one point
      j = total_fibers_clmn - 2;
      nodes[j].bend_force_x = bending_const * (2 * nodes[j+1].x
        - 5 * nodes[j].x
        + 4 * nodes[j-1].x
        - nodes[j].x);
      nodes[j].bend_force_y = bending_const * (2 * nodes[j+1].y
        - 5 * nodes[j].y
        + 4 * nodes[j-1].y
        - nodes[j - 2].y);
      nodes[j].bend_force_z = bending_const * (2 * nodes[j+1].z
        - 5 * nodes[j].z
        + 4 * nodes[j-1].z
        - nodes[j-2].z);
      //For last ponit
      j = total_fibers_clmn - 1;
      nodes[j].bend_force_x = bending_const * (-nodes[j].x
        + 2 * nodes[j-1].x
        - nodes[j-2].x);
      nodes[j].bend_force_y = bending_const * (-nodes[j].y
        + 2 * nodes[j-1].y
        - nodes[j-2].y);
      nodes[j].bend_force_z = bending_const * (-nodes[j].z
        + 2 * nodes[j-1].z
        - nodes[j-2].z);
    }//if fiber2thread ends
  }

  pthread_barrier_wait(&(gv->barr));

#if 0
      if (tid == 0){
        printf("Inside BF -- Printing for Corner Points(z,y) : 0,0 \n");
        print_fiber_sub_grid(gv, 0, 0, 0, 0);
        printf("Inside BF -- Printing for Corner Points(z,y) : 51,0 \n");
        print_fiber_sub_grid(gv, 0, 51, 0, 51);
        printf("Inside BF -- Printing for Corner Points(z,y) : 0,51 \n");
        print_fiber_sub_grid(gv, 51, 0, 51, 0);
        printf("Inside BF -- Printing for Corner Points (z,y): 51,51 \n");
        print_fiber_sub_grid(gv, 51, 51, 51, 51);
      }
      pthread_barrier_wait(&(gv->barr));
#endif //DEBUG_PRINT

// #ifdef DEBUG_PRINT
//   printf("After barrier in BF\n");
// #endif //DEBUG_PRINT
  /* computing bending force along direction normal to fibers (column), bending_const may be different */
  Fibersheet* sheet = fiber_shape->sheets + 0;
  bending_const = gv->Kb_l / pow(ds2, 3);                 //KB and alpha1 pow 3
  for (j = 0; j < total_fibers_clmn; ++j){
    if (fiber2thread(0, total_fibers_row, total_threads) == tid){
      //For first fiber
      sheet->fibers[0].nodes[j].bend_force_x += bending_const*(-sheet->fibers[2].nodes[j].x
        + 2 * (sheet->fibers[1].nodes[j].x)
        - (sheet->fibers[0].nodes[j].x));
      sheet->fibers[0].nodes[j].bend_force_y += bending_const*(-sheet->fibers[2].nodes[j].y
        + 2 * (sheet->fibers[1].nodes[j].y)
        - (sheet->fibers[0].nodes[j].y));
      sheet->fibers[0].nodes[j].bend_force_z += bending_const*(-sheet->fibers[2].nodes[j].z
        + 2 * (sheet->fibers[1].nodes[j].z)
        - (sheet->fibers[0].nodes[j].z));
    } //if fiber2thread ends
    //For 2nd fiber
    if (fiber2thread(1, total_fibers_row, total_threads) == tid){
      sheet->fibers[1].nodes[j].bend_force_x += bending_const*(-sheet->fibers[3].nodes[j].x
        + 4 * (sheet->fibers[2].nodes[j].x)
        - 5 * (sheet->fibers[1].nodes[j].x)
        + 2 * (sheet->fibers[0].nodes[j].x));
      sheet->fibers[1].nodes[j].bend_force_y += bending_const*(-sheet->fibers[3].nodes[j].y
        + 4 * (sheet->fibers[2].nodes[j].y)
        - 5 * (sheet->fibers[1].nodes[j].y)
        + 2 * (sheet->fibers[0].nodes[j].y));
      sheet->fibers[1].nodes[j].bend_force_z += bending_const*(-sheet->fibers[3].nodes[j].z
        + 4 * (sheet->fibers[2].nodes[j].z)
        - 5 * (sheet->fibers[1].nodes[j].z)
        + 2 * (sheet->fibers[0].nodes[j].z));
    }//if fiber2thread ends

    // For the middle fibers
    for (i = 2; i < total_fibers_row - 2; ++i){
      if (fiber2thread(i, total_fibers_row, total_threads) == tid){
        sheet->fibers[i].nodes[j].bend_force_x += bending_const*(-sheet->fibers[i+2].nodes[j].x
          + 4.0 * (sheet->fibers[i+1].nodes[j].x)
          - 6.0 * (sheet->fibers[i].nodes[j].x) +
          4.0 * (sheet->fibers[i-1].nodes[j].x)
          - (sheet->fibers[i-2].nodes[j].x));
        sheet->fibers[i].nodes[j].bend_force_y += bending_const*(-sheet->fibers[i+2].nodes[j].y
          + 4.0*(sheet->fibers[i+1].nodes[j].y)
          - 6.0*(sheet->fibers[i].nodes[j].y)
          + 4.0*(sheet->fibers[i-1].nodes[j].y)
          - (sheet->fibers[i-2].nodes[j].y));
        sheet->fibers[i].nodes[j].bend_force_z += bending_const*(-sheet->fibers[i+2].nodes[j].z
          + 4.0*(sheet->fibers[i+1].nodes[j].z)
          - 6.0*(sheet->fibers[i].nodes[j].z)
          + 4.0*(sheet->fibers[i-1].nodes[j].z)
          - (sheet->fibers[i-2].nodes[j].z));
      }//if fiber2thread ends
    }

    //For last but one fiber
    i = total_fibers_row - 2;
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
      sheet->fibers[i].nodes[j].bend_force_x += bending_const * (2 * (sheet->fibers[i+1].nodes[j].x)
        - 5 * (sheet->fibers[i].nodes[j].x)
        + 4 * (sheet->fibers[i-1].nodes[j].x)
        - (sheet->fibers[i-2].nodes[j].x));
      sheet->fibers[i].nodes[j].bend_force_y += bending_const * (2 * (sheet->fibers[i+1].nodes[j].y)
        - 5 * (sheet->fibers[i].nodes[j].y)
        + 4 * (sheet->fibers[i-1].nodes[j].y)
        - (sheet->fibers[i-2].nodes[j].y));
      sheet->fibers[i].nodes[j].bend_force_z += bending_const  * (2 * (sheet->fibers[i+1].nodes[j].z)
        - 5 * (sheet->fibers[i].nodes[j].z)
        + 4 * (sheet->fibers[i-1].nodes[j].z)
        - (sheet->fibers[i - 2].nodes[j].z));
    }
    //For the last fiber
    i = total_fibers_row - 1;
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
      sheet->fibers[i].nodes[j].bend_force_x += bending_const * (-(sheet->fibers[i].nodes[j].x)
        + 2 * (sheet->fibers[i-1].nodes[j].x)
        - (sheet->fibers[i-2].nodes[j].x));
      sheet->fibers[i].nodes[j].bend_force_y += bending_const * (-(sheet->fibers[i].nodes[j].y)
        + 2 * (sheet->fibers[i-1].nodes[j].y)
        - (sheet->fibers[i-2].nodes[j].y));
      sheet->fibers[i].nodes[j].bend_force_z += bending_const * (-(sheet->fibers[i].nodes[j].z)
        + 2 * (sheet->fibers[i-1].nodes[j].z)
        - (sheet->fibers[i-2].nodes[j].z));
    }
  }

#ifdef DEBUG_PRINT
  // printf("**** compute_bendingforce Exit*******\n");
#endif //DEBUG_PRINT
}
