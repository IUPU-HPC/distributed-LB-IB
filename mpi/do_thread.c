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
#include "timer.h"

static inline double get_cur_time()
{
  struct timeval tv;
  double t;

  gettimeofday(&tv,NULL);
  t = tv.tv_sec + tv.tv_usec / 1e6;
  return t;
}

void* do_thread(void* v){
  LV lv = (LV)v;
  GV gv = lv->gv;
  int barrcheck;
  int tid = lv->tid;
  int num_fluid_tasks = gv->num_fluid_tasks;
  int my_rank = gv->taskid;
  double t0, t1;
  double t2=0, t3=0, t4=0, t5=0, t6=0, tail=0;
  double t00, t11, t4_1=0, t4_2=0, t4_3=0, t4_4=0;
  double startTime=0, endTime=0;
  char filename[80];

  startTime = get_cur_time();

  while (gv->time <= gv->timesteps) {

    t0 = get_cur_time();

    /******************* Part-1: Fiber compute force locally *******************/
    if (my_rank >= num_fluid_tasks){
#ifdef DEBUG_PRINT
      if (tid==0){
        printf("\n\nTask%d Tid%d: Start time step %ld ...\n", my_rank, tid, gv->time);
      }
#endif //DEBUG_PRINT

      compute_bendingforce(lv);
      pthread_barrier_wait(&(gv->barr));
#if 0 //Verify results
      if (tid == 0){
        sprintf(filename, "Fiber%d_bending_force%d.dat", my_rank, gv->time);
        save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif

      compute_stretchingforce(lv);
      pthread_barrier_wait(&(gv->barr));
#if 0 //Verify results
      if (tid == 0){
        sprintf(filename, "Fiber%d_stretching_force%d.dat", my_rank, gv->time);
        save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif

      compute_elasticforce(lv);
      pthread_barrier_wait(&(gv->barr));
#if 0 //Verify results
      if (tid == 0){
        my_rank = 1;
        sprintf(filename, "Fiber%d_elasitc_force%d.dat", my_rank, gv->time);
        save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif      

#ifdef DEBUG_PRINT
    printf("Fiber%dtid%d finish calculate local force\n", my_rank, tid);
#endif //DEBUG_PRINT
    }
    else{ //fluid tasks wait here
      // pthread_barrier_wait(&(gv->barr));
    }

    t1 = get_cur_time();
    t2 += t1 - t0;

    if(tid == 0){
      MPI_Barrier(MPI_COMM_WORLD);
    }
    /**************************************Part-1: end **************************************/

    /******************* Part-2: Fiber spread force to influential fluid *******************/
    t0 = get_cur_time();

    if (my_rank >= num_fluid_tasks){

#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fiber%dtid%d: Start fiber_SpreadForce\n", my_rank, tid);
        fflush(stdout);
      }
#endif //DEBUG_PRINT

      fiber_SpreadForce(lv);

    }
    else{
#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fluid%dtid%d: Start fluid_get_SpreadForce\n", my_rank, tid);
        fflush(stdout);
      }
#endif //DEBUG_PRINT

      fluid_get_SpreadForce(lv);

    }

    t1 = get_cur_time();
    t3 += t1 - t0;

    pthread_barrier_wait(&(gv->barr));    
    if(tid==0){
      MPI_Barrier(MPI_COMM_WORLD);
    }

#ifdef DEBUG_PRINT
    if(tid==0){
      printf("Task%d: After fluid_get_SpreadForce\n", my_rank);
      fflush(stdout);
    }
#endif //DEBUG_PRINT

#ifdef VERIFY //Verify results
    if (my_rank < num_fluid_tasks){
      if (tid == 0){
        sprintf(filename, "Fluid%d_get_SpreadForce_step%d.dat", my_rank, gv->time);
        save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        // save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
      }
      pthread_barrier_wait(&(gv->barr));
    }
    else{
      if (tid == 0){
        sprintf(filename, "Fiber%d_SpreadForce_step%d.dat", my_rank, gv->time);
        save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
      }
      pthread_barrier_wait(&(gv->barr));
    }
#endif
    /**************************************Part-2: end **************************************/   

    /******************************** Part-3: Fluid streaming *******************************/
    // Fluid tasks
    t0 = get_cur_time();

    if(my_rank < num_fluid_tasks){

      t00 = get_cur_time();
#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fluid%d: Start compute_DF1\n", my_rank);
        fflush(stdout);
      }
#endif //DEBUG_PRINT

      compute_eqlbrmdistrfuncDF1(lv);
      pthread_barrier_wait(&(gv->barr));
#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fluid%d: After compute DF1\n", my_rank);
        fflush(stdout);
      }
#endif //DEBUG_PRINT
#ifdef VERIFY //Verify results
      if (tid == 0){
        sprintf(filename, "Fluid%d_compute_DF1_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif
      t11 = get_cur_time();
      t4_1 += t11 - t00;

      t00 = get_cur_time();
#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fluid%d: Start streaming\n", my_rank);
        fflush(stdout);
      }
#endif //DEBUG_PRINT
      stream_distrfunc(lv);
      pthread_barrier_wait(&(gv->barr));
      // if(tid==0)
      //   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      // if(tid==0){
        printf("Fluid%dtid%d: After streaming\n", my_rank, tid);
        fflush(stdout);
      // }
#endif //DEBUG_PRINT
#ifdef VERIFY //Verify results
      if (tid == 0){
        sprintf(filename, "Fluid%d_streaming_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif
      t11 = get_cur_time();
      t4_2 += t11 - t00;

      t00 = get_cur_time();
      bounceback_rigidwalls(lv);
      pthread_barrier_wait(&(gv->barr));
      // if (tid == 0)
      //   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      printf("Fluid%d: After bounceback_rigidwalls\n", my_rank);
#endif //DEBUG_PRINT
#ifdef VERIFY //Verify results
      if (tid == 0){
        sprintf(filename, "Fluid%d_bounceback_rigidwalls_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif
      t11 = get_cur_time();
      t4_3 += t11 - t00;

      t00 = get_cur_time();
      compute_rho_and_u(lv);
      pthread_barrier_wait(&(gv->barr));
      // if (tid == 0)
      //   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      printf("Fluid%d: After compute_rho_and_u\n", my_rank);
#endif //DEBUG_PRINT
#ifdef VERIFY //Verify results
      if (tid == 0){
        sprintf(filename, "Fluid%d_rho_and_u_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif
      t11 = get_cur_time();
      t4_4 += t11 - t00;

    }

    t1 = get_cur_time();
    t4 += t1 - t0;

    // pthread_barrier_wait(&(gv->barr));
    if (tid == 0)
      MPI_Barrier(MPI_COMM_WORLD);
    pthread_barrier_wait(&(gv->barr));

    /**************************************Part-3: end **************************************/ 

    /******************************** Part-4: Fluid to Fiber velocity *******************************/
    // move_fiber
    if (my_rank < num_fluid_tasks){
      t0 = get_cur_time();
      fluid_SpreadVelocity(lv);
      t1 = get_cur_time();
      t5 += t1 - t0;
    }
    else{
      t0 = get_cur_time();
      fiber_get_SpreadVelocity(lv);
      t1 = get_cur_time();
      t5 += t1 - t0;
    }

    pthread_barrier_wait(&(gv->barr));
    // if (tid == 0)
    //   MPI_Barrier(MPI_COMM_WORLD);
    // pthread_barrier_wait(&(gv->barr));

#ifdef DEBUG_PRINT
    printf("Task%dtid%d: After moving fibersheet \n", my_rank, tid);
#endif //DEBUG_PRINT

#ifdef VERIFY //Verify results
    if (my_rank < num_fluid_tasks){
      if (tid == 0){
        sprintf(filename, "Fluid%d_moveFiber_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        // save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
      }
      pthread_barrier_wait(&(gv->barr));
    }
    else{
      if (tid == 0){
        sprintf(filename, "Fiber%d_moveFiber_step%d.dat", my_rank, gv->time);
        save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
      }
      pthread_barrier_wait(&(gv->barr));
    }
#endif

    /**************************************Part-4: end **************************************/

    /******************************** Part-5: Fluid wrap up *******************************/
    // Fluid tasks
    t0 = get_cur_time();

    if(my_rank < num_fluid_tasks){
      copy_inout_to_df2(lv);
      pthread_barrier_wait(&(gv->barr));
      // if (tid == 0) 
          // MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      printf("Fluid%dtid%d: After copy_inout_to_df2 \n", my_rank, tid);
#endif //DEBUG_PRINT

      replace_old_DF(lv);
      pthread_barrier_wait(&(gv->barr));
      // if (tid == 0) 
      //   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      printf("Fluid%dtid%d: After replace_old_DF \n", my_rank, tid);
#endif //DEBUG_PRINT

      // periodicBC(lv);
      pthread_barrier_wait(&(gv->barr));
      // if (tid == 0) 
      //   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      printf("Fluid%dtid%d: After PeriodicBC\n", my_rank, tid);
#endif //DEBUG_PRINT
    }

    t1 = get_cur_time();
    t6 += t1 - t0;
    /**************************************Part-5: end **************************************/

#ifdef SAVE
    if (tid == 0 && (gv->time % gv->dump == 0)){
      if (my_rank < num_fluid_tasks){
#if 0      
      sprintf(filename, "Fluid%d_dump_step%d.dat", my_rank, gv->time);
      // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
      save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
#endif
      }  
      else{
        // sprintf(filename, "Fiber%d_dump_step%d.dat", my_rank, gv->time);
        // save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
        sprintf(filename, "Fiber%d_dump_step%d.bin", my_rank, gv->time);
        save_fiber_sub_grid_binary(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
      }
    }
    pthread_barrier_wait(&(gv->barr));
#endif

#if 1
    if (my_rank == 0 && tid == 0){
      printf("End of time step %d\n", gv->time);
      fflush(stdout);
    }
#endif


    if (tid == 0)//does it requires gv->my_rank==0 i.e only one machine updating the counter value
      gv->time += gv->dt;

#if 0
    if ( (my_rank >= num_fluid_tasks) && (tid == 0)){
      printf("Timesteps:%ld complete!\n", gv->time - 1);
      printf("Printing for Corner Points(z,y) : 0,0 \n");
      print_fiber_sub_grid(gv, 0, 0, 0, 0);
      printf("Printing for Corner Points(z,y) : 51,0 \n");
      print_fiber_sub_grid(gv, 0, 51, 0, 51);
      printf("Printing for Corner Points(z,y) : 0,51 \n");
      print_fiber_sub_grid(gv, 51, 0, 51, 0);
      printf("Printing for Corner Points (z,y): 51,51 \n");
      print_fiber_sub_grid(gv, 51, 51, 51, 51);
    }
#endif //DEBUG_PRINT

    t0 = get_cur_time();
    if (tid == 0){
      MPI_Barrier(MPI_COMM_WORLD);
    }
    pthread_barrier_wait(&(gv->barr));
    t1 = get_cur_time();
    tail += t1 - t0;

  }

  endTime = get_cur_time();

  // if(tid == 0){
    if(my_rank >= num_fluid_tasks){
      printf("Fiber%dtid%d: T_comp_F=%.3f, T_Sprd_F=%.3f, T_get_Sprd_Vel=%.3f, T_barr=%.3f, total=%.3f\n",
        my_rank, tid, t2, t3, t5, tail, endTime-startTime);
      fflush(stdout);
    }
    else{
      printf("Fluid%dtid%d: T_wait=%.3f, T_get_Sprd_F=%.3f, \
T_stream=%.3f, T_4_1=%.3f, T_4_2=%.3f, T4_3=%.3f, T4_4=%.3f, \
T_Sprd_Vel=%.3f, T_tail=%.3f, T_barr=%.3f, total=%.3f\n",
        my_rank, tid, t2, t3, 
        t4, t4_1, t4_2, t4_3, t4_4,
        t5, t6, tail, endTime-startTime);
      fflush(stdout);
    }
  // }
  return NULL;
}