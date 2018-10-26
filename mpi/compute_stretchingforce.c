/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

void compute_stretchingforce(LV lv) {
  #ifdef DEBUG_PRINT
  // printf("****Inside compute_stretchingforce*******\n");
  #endif //DEBUG_PRINT
  int tid;
  GV gv = lv->gv;
  tid   = lv->tid;
  int total_threads = gv->threads_per_task;


  int i =0;
  int j = 0;

  double ds1 = gv->fiber_shape->sheets[0].width / (gv->fiber_shape->sheets[0].num_cols - 1);
  double ds2 = gv->fiber_shape->sheets[0].height / (gv->fiber_shape->sheets[0].num_rows - 1);  // corresponds to ydist on the fibersheet.

  double stretching_const = gv->Ks_l/(ds1 * ds1);

  /* for convenience to save space. */
  int    total_fibers_row  = gv->fiber_shape->sheets[0].num_rows;
  int    total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  Fiber* fiberarray        = gv->fiber_shape->sheets[0].fibers;


  double tangx_right, tangx_left, tangy_right, tangy_left, tangz_right, tangz_left, dist_right, dist_left;
  double tangx_top, tangx_bottom, tangy_top, tangy_bottom, tangz_top, tangz_bottom, dist_top, dist_bottom;


  //computing stretching force along fibers (row)
  for(i=0; i<total_fibers_row; ++i){  /* for all fibers  */
   if(fiber2thread(i, total_fibers_row, total_threads) == tid ){
    /* general cases for the points in the middle */
    for(j=1; j < total_fibers_clmn-1; ++j){ /* for all points for a given fiber*/
      tangx_right = fiberarray[i].nodes[j+1].x - fiberarray[i].nodes[j].x;
      tangy_right = fiberarray[i].nodes[j+1].y - fiberarray[i].nodes[j].y;
      tangz_right = fiberarray[i].nodes[j+1].z - fiberarray[i].nodes[j].z;
      tangx_left = fiberarray[i].nodes[j-1].x -fiberarray[i].nodes[j].x;
      tangy_left = fiberarray[i].nodes[j-1].y -fiberarray[i].nodes[j].y;
      tangz_left = fiberarray[i].nodes[j-1].z -fiberarray[i].nodes[j].z;
      dist_right = sqrt(tangx_right * tangx_right + tangy_right * tangy_right + tangz_right * tangz_right);
      dist_left = sqrt(tangx_left * tangx_left + tangy_left * tangy_left + tangz_left * tangz_left);
      fiberarray[i].nodes[j].stretch_force_x = stretching_const*((dist_right-ds1)* (tangx_right/dist_right)+ (dist_left-ds1)* (tangx_left/dist_left));
      fiberarray[i].nodes[j].stretch_force_y = stretching_const*((dist_right-ds1)* (tangy_right/dist_right)+ (dist_left-ds1)* (tangy_left/dist_left));
      fiberarray[i].nodes[j].stretch_force_z = stretching_const*((dist_right-ds1)* (tangz_right/dist_right)+ (dist_left-ds1)* (tangz_left/dist_left));
    }


    /* for the first point */
    tangx_right = fiberarray[i].nodes[1].x - fiberarray[i].nodes[0].x;
    tangy_right = fiberarray[i].nodes[1].y - fiberarray[i].nodes[0].y;
    tangz_right = fiberarray[i].nodes[1].z - fiberarray[i].nodes[0].z;
    dist_right = sqrt(tangx_right * tangx_right + tangy_right * tangy_right + tangz_right * tangz_right);
    fiberarray[i].nodes[0].stretch_force_x = stretching_const*((dist_right-ds1)* (tangx_right/dist_right));
    fiberarray[i].nodes[0].stretch_force_y = stretching_const*((dist_right-ds1)* (tangy_right/dist_right));
    fiberarray[i].nodes[0].stretch_force_z = stretching_const*( (dist_right-ds1)* (tangz_right/dist_right));


    /* for the last point */
    tangx_left = fiberarray[i].nodes[total_fibers_clmn-1].x -fiberarray[i].nodes[total_fibers_clmn-2].x;
    tangy_left = fiberarray[i].nodes[total_fibers_clmn-1].y -fiberarray[i].nodes[total_fibers_clmn-2].y;
    tangz_left = fiberarray[i].nodes[total_fibers_clmn-1].z -fiberarray[i].nodes[total_fibers_clmn-2].z;
    dist_left = sqrt(tangx_left * tangx_left + tangy_left * tangy_left + tangz_left * tangz_left);
    fiberarray[i].nodes[total_fibers_clmn-1].stretch_force_x =stretching_const*( (dist_left-ds1)* (tangx_left/dist_left));
    fiberarray[i].nodes[total_fibers_clmn-1].stretch_force_y =stretching_const*( (dist_left-ds1)* (tangy_left/dist_left));
    fiberarray[i].nodes[total_fibers_clmn-1].stretch_force_z =stretching_const*((dist_left-ds1)* (tangz_left/dist_left));
   }//if fiber2thread ends
  } //end of total_fibers along row i.e fibers along z axis


     pthread_barrier_wait(&(gv->barr));

  // COmputing Streching Force along the direction normal to fibers
  stretching_const = gv->Ks_l/(ds2 * ds2);
  for(j=0; j<total_fibers_clmn; ++j){  //for the j-th point
    for(i=1; i<total_fibers_row-1; ++i){ //for the i-th fiber, middle fibers
      if(fiber2thread(i,total_fibers_row,total_threads) == tid ){
      tangx_top = fiberarray[i+1].nodes[j].x -fiberarray[i].nodes[j].x;
      tangy_top = fiberarray[i+1].nodes[j].y -fiberarray[i].nodes[j].y;
      tangz_top = fiberarray[i+1].nodes[j].z -fiberarray[i].nodes[j].z;
      tangx_bottom = fiberarray[i-1].nodes[j].x -fiberarray[i].nodes[j].x;
      tangy_bottom = fiberarray[i-1].nodes[j].y -fiberarray[i].nodes[j].y;
      tangz_bottom = fiberarray[i-1].nodes[j].z -fiberarray[i].nodes[j].z;
      dist_top = sqrt(tangx_top * tangx_top + tangy_top * tangy_top + tangz_top * tangz_top);
      dist_bottom = sqrt(tangx_bottom * tangx_bottom + tangy_bottom * tangy_bottom + tangz_bottom * tangz_bottom);
      fiberarray[i].nodes[j].stretch_force_x += stretching_const*((dist_top-ds2)* (tangx_top/dist_top)+ (dist_bottom-ds2)* (tangx_bottom/dist_bottom));
      fiberarray[i].nodes[j].stretch_force_y += stretching_const*((dist_top-ds2)* (tangy_top/dist_top)+ (dist_bottom-ds2)* (tangy_bottom/dist_bottom));
      fiberarray[i].nodes[j].stretch_force_z += stretching_const*((dist_top-ds2)* (tangz_top/dist_top)+ (dist_bottom-ds2)* (tangz_bottom/dist_bottom));
    }//if fiber2thread ends
   }

    /* For the first fiber on the sheet boundary */
    if(fiber2thread(0,total_fibers_row,total_threads) == tid ){ //i =0 :first fiber node
    tangx_top = fiberarray[1].nodes[j].x -fiberarray[0].nodes[j].x;
    tangy_top = fiberarray[1].nodes[j].y -fiberarray[0].nodes[j].y;
    tangz_top = fiberarray[1].nodes[j].z -fiberarray[0].nodes[j].z;
    dist_top = sqrt(tangx_top * tangx_top + tangy_top * tangy_top + tangz_top * tangz_top);
    fiberarray[0].nodes[j].stretch_force_x +=stretching_const*((dist_top-ds2)* (tangx_top/dist_top));
    fiberarray[0].nodes[j].stretch_force_y +=stretching_const*((dist_top-ds2)* (tangy_top/dist_top));
    fiberarray[0].nodes[j].stretch_force_z += stretching_const*((dist_top-ds2)* (tangz_top/dist_top));
  }//if fiber2thread ends

    /* For the last fiber on the sheet boundary */
  if(fiber2thread(total_fibers_row-1,total_fibers_row,total_threads) == tid ){ //i = total_fibers_row-1 :last fiber node
    tangx_bottom = fiberarray[total_fibers_row-1].nodes[j].x -fiberarray[total_fibers_row-2].nodes[j].x;
    tangy_bottom = fiberarray[total_fibers_row-1].nodes[j].y -fiberarray[total_fibers_row-2].nodes[j].y;
    tangz_bottom = fiberarray[total_fibers_row-1].nodes[j].z -fiberarray[total_fibers_row-2].nodes[j].z;
    dist_bottom = sqrt(tangx_bottom * tangx_bottom + tangy_bottom * tangy_bottom + tangz_bottom * tangz_bottom);
    fiberarray[total_fibers_row-1].nodes[j].stretch_force_x += stretching_const*((dist_bottom-ds2)* (tangx_bottom/dist_bottom));
    fiberarray[total_fibers_row-1].nodes[j].stretch_force_y += stretching_const*((dist_bottom-ds2)* (tangy_bottom/dist_bottom));
    fiberarray[total_fibers_row-1].nodes[j].stretch_force_z += stretching_const*((dist_bottom-ds2)* (tangz_bottom/dist_bottom));
   }
  }

  pthread_barrier_wait(&(gv->barr));

  #ifdef DEBUG_PRINT
// printf("**** compute_stretchingforce Exit*******\n");
 #endif //DEBUG_PRINT
}
