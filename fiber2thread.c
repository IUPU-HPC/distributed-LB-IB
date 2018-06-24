/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

int fiber2thread(int fiber_row_i, int num_fibers, int num_threads){
  int tid;
  int Block_size;

  Block_size = num_fibers/num_threads;

  if( num_fibers%num_threads !=0  && Block_size * num_threads < num_fibers)
  {
    Block_size = Block_size +1;
  }
  tid        = fiber_row_i/Block_size;//ith fiber row.
 return tid;
  //return 0;
}
