/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

void   compute_elasticforce(LV lv) { //emabarassingly parallel
#ifdef DEBUG_PRINT
  // printf("****Inside compute_elasticforce*******\n");
#endif //DEBUG_PRINT
  int tid;
  GV gv = lv->gv;
  tid = lv->tid;
  int total_threads = gv->total_threads;

  int i = 0, j = 0;
  int total_fibers_row = gv->fiber_shape->sheets[0].num_rows;
  int total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  Fiber *fiberarray = gv->fiber_shape->sheets[0].fibers;


  // Adding up Streching Force and Bending Force
  for (i = 0; i<total_fibers_row; ++i){
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
      for (j = 0; j<total_fibers_clmn; ++j){
        fiberarray[i].nodes[j].elastic_force_x = fiberarray[i].nodes[j].stretch_force_x
          + fiberarray[i].nodes[j].bend_force_x;
        fiberarray[i].nodes[j].elastic_force_y = fiberarray[i].nodes[j].stretch_force_y
          + fiberarray[i].nodes[j].bend_force_y;
        fiberarray[i].nodes[j].elastic_force_z = fiberarray[i].nodes[j].stretch_force_z
          + fiberarray[i].nodes[j].bend_force_z;
      }//end of total_fluidgridpoints i.e fibers along y axis
    }//if fiber2thread ends
  }//end of total_fibers i.e fibers along x axis

  pthread_barrier_wait(&(gv->barr));

#ifdef DEBUG_PRINT
  // printf("**** compute_elasticforce Exit*******\n");
#endif //DEBUG_PRINT
}
