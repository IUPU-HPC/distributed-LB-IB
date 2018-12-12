/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
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
  double t2=0, t2_1=0, t3=0, t3_1=0, t4=0, t5=0, t6=0, tail=0;
  double startTime=0, endTime=0;
  char filename[80];

  Timer::time_start();

  while (gv->time <= gv->timesteps) {

    t0 = get_cur_time();

    // Fiber tasks
    if (my_rank >= num_fluid_tasks){
#ifdef DEBUG_PRINT
      if(tid==0){
        printf("\n\nTask%d Tid%d: Start time step %ld ...\n", my_rank, tid, gv->time);
      }
#endif //DEBUG_PRINT
      t0 = get_cur_time();

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

      t1 = get_cur_time();
      t2_1 += t1 - t0;

#ifdef DEBUG_PRINT
    printf("Fiber%dtid%d finish calculate local force\n", my_rank, tid);
#endif //DEBUG_PRINT
    }
    else{ //fluid tasks wait here
      // pthread_barrier_wait(&(gv->barr));
    }

    //check whether barrier long is becuase of this
    if(tid==0)
      MPI_Barrier(MPI_COMM_WORLD);

    t1 = get_cur_time();
    t2 += t1 - t0;

    if (my_rank >= num_fluid_tasks){

#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fiber%dtid%d: Start fiber_SpreadForce\n", my_rank, tid);
        fflush(stdout);
      }
#endif //DEBUG_PRINT

      t0 = get_cur_time();
      fiber_SpreadForce(lv);
      t1 = get_cur_time();
      t3_1 += t1- t0;
    }
    else{
#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fluid%dtid%d: Start fluid_get_SpreadForce\n", my_rank, tid);
        fflush(stdout);
      }
#endif //DEBUG_PRINT

      t0 = get_cur_time();
      fluid_get_SpreadForce(lv);
      t1 = get_cur_time();
      t4 += t1- t0;
    }

    pthread_barrier_wait(&(gv->barr));
    if(tid==0)
      MPI_Barrier(MPI_COMM_WORLD);
    t1 = get_cur_time();
    t3 += t1- t0;
    t4 += t1- t0;
#ifdef DEBUG_PRINT
    if(tid==0){
      printf("Task%d: After fluid_get_SpreadForce\n", my_rank);
      fflush(stdout);
    }
#endif //DEBUG_PRINT

#if 1 //Verify results
    if (my_rank < num_fluid_tasks){
      if (tid == 0){
        sprintf(filename, "Fluid%d_get_SpreadForce_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
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

    // Fluid tasks
    if(my_rank < num_fluid_tasks){

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
#if 1 //Verify results
      if (tid == 0){
        sprintf(filename, "Fluid%d_compute_DF1_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif

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
#if 0 //Verify results
      if (tid == 0){
        sprintf(filename, "Fluid%d_streaming_step%d.dat", my_rank, gv->time);
        save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif

      bounceback_rigidwalls(lv);
      pthread_barrier_wait(&(gv->barr));
      // if (tid == 0)
      //   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      printf("Fluid%d: After bounceback_rigidwalls\n", my_rank);
#endif //DEBUG_PRINT
#if 1 //Verify results
      if (tid == 0){
        sprintf(filename, "Fluid%d_bounceback_rigidwalls_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif

      compute_rho_and_u(lv);
      pthread_barrier_wait(&(gv->barr));
      // if (tid == 0)
      //   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      printf("Fluid%d: After compute_rho_and_u\n", my_rank);
#endif //DEBUG_PRINT
#if 1 //Verify results
      if (tid == 0){
        sprintf(filename, "Fluid%d_rho_and_u_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
      }
      pthread_barrier_wait(&(gv->barr));
#endif
    }

// VERIFY
    // pthread_barrier_wait(&(gv->barr));
    if (tid == 0)
      MPI_Barrier(MPI_COMM_WORLD);
    pthread_barrier_wait(&(gv->barr));

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
      t6 += t1 - t0;
    }

    pthread_barrier_wait(&(gv->barr));
    // if (tid == 0)
    //   MPI_Barrier(MPI_COMM_WORLD);
    // pthread_barrier_wait(&(gv->barr));

#ifdef DEBUG_PRINT
    printf("Task%dtid%d: After moving fibersheet \n", my_rank, tid);
#endif //DEBUG_PRINT

#ifdef SAVE //Verify results
    if (my_rank < num_fluid_tasks){
      if (tid == 0){
        sprintf(filename, "Fluid%d_moveFiber_step%d.dat", my_rank, gv->time);
        // save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
        save_fluid_sub_grid(gv, 19, 20, 10, 22, 25, 33, filename);
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

    // Fluid tasks
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

// VERIFY

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

  double time_elapsed = Timer::time_end();

  if(my_rank >= num_fluid_tasks){
    printf("Fiber%dtid%d: compute_force=%.3f || %.3f, fiber_SpreadForce=%.3f || %.3f, fiber_get_SpreadVelocity=%.3f, tail=%.3f, total=%.3f\n",
      my_rank, tid, t2, t2_1, t3, t3_1, t6, tail, endTime-startTime);
    fflush(stdout);
  }
  else{
    printf("Fluid%dtid%d: t2=%.3f, fluid_get_SpreadForce=%.3f, fluid_SpreadVelocity=%.3f, tail=%.3f, total=%.3f\n",
      my_rank, tid, t2, t4, t5, tail, endTime-startTime);
    fflush(stdout);
  }
  
  return NULL;
}