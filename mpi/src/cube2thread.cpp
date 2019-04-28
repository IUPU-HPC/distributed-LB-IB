/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

int cube2thread(int BI, int BJ, int BK, int num_cubes_x, int num_cubes_y, int num_cubes_z, int P, int Q, int R){
 int tid, i, j, k;
 int Block_size_x, Block_size_y, Block_size_z;

 Block_size_x  =  num_cubes_x/P;
 Block_size_y  =  num_cubes_y/Q;
 Block_size_z  =  num_cubes_z/R;


 i = BI / Block_size_x;
 j = BJ / Block_size_y;
 k = BK / Block_size_z;

 tid  = i*Q*R + j*R + k;

 return tid;
 /*return 0;//for 1 thread*/
}
