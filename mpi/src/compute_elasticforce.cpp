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

void compute_elasticforce(LV lv) { //emabarassingly parallel
#ifdef DEBUG_PRINT
   printf("****Inside compute_elasticforce*******\n");
#endif //DEBUG_PRINT

  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  int total_threads = gv->threads_per_task;

  int i, j;
  // TODO: Modify to multi sheets
  int total_fibers_row  = gv->fiber_shape->sheets[0].num_rows;
  int total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  Fiber *fiberarray     = gv->fiber_shape->sheets[0].fibers;

  // Adding up Streching Force and Bending Force
  for (i = 0; i < total_fibers_row; ++i){ 
    if (fiber2thread(i, total_fibers_row, total_threads) == tid ){ 
      for (j = 0; j < total_fibers_clmn; ++j){
        fiberarray[i].nodes[j].elastic_force_x = fiberarray[i].nodes[j].stretch_force_x 
                                               + fiberarray[i].nodes[j].bend_force_x;
        fiberarray[i].nodes[j].elastic_force_y = fiberarray[i].nodes[j].stretch_force_y 
                                               + fiberarray[i].nodes[j].bend_force_y;
        fiberarray[i].nodes[j].elastic_force_z = fiberarray[i].nodes[j].stretch_force_z 
                                               + fiberarray[i].nodes[j].bend_force_z;
#if 0                                               
        printf("tid%d: (%2d, %2d) || (%.6f,%.24f,%.24f)\n", 
          tid, i, j, fiberarray[i].nodes[j].elastic_force_x, fiberarray[i].nodes[j].elastic_force_y, fiberarray[i].nodes[j].elastic_force_z);
#endif          
      } // end of total_fluidgridpoints i.e fibers along y axis 
    } // if fiber2thread ends
   } // end of total_fibers i.e fibers along x axis 
   
#ifdef DEBUG_PRINT
     printf("**** compute_elasticforce Exit*******\n");
#endif //DEBUG_PRINT
}
