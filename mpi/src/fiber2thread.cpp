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

int fiber2thread(int fiber_row_i, int num_fibers, int num_threads){
  int tid;
  int Block_size;

/* Original code by Prateek. Aim to get the ceiling.
  Block_size = num_fibers / num_threads;

  if ( num_fibers % num_threads != 0  && Block_size * num_threads < num_fibers){
    Block_size = Block_size + 1;
  }
*/

  Block_size = (num_fibers + num_fibers - 1) / num_threads;

  tid  = fiber_row_i / Block_size; //ith fiber row.

  return tid;
}
