/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

int check_1d(int x, int y, int z, int* dir){

  if ( (y == 1) && (z == 1) ){
    if(dir!=NULL)
      *dir = X_transfer_1D;
    return 1;
  }
  else if( (x == 1) && (z == 1) ){
    if(dir!=NULL)
      *dir = Y_transfer_1D;
    return 1;
  }
  else if ( (x == 1) && (y == 1) ){
    if(dir!=NULL)
      *dir = Z_transfer_1D;
    return 1;
  }
  else
    return 0;
}

int check_2d(int x, int y, int z, int* dir){
  if(x == 1){
    if(dir!=NULL)
      *dir = Y_transfer_2D;
    return 1;
  }
  else if (y == 1) {
    if(dir!=NULL)
      *dir = Z_transfer_2D;
    return 1;
  }
  else if (z == 1) {
    if(dir!=NULL)
      *dir = X_transfer_2D;
    return 1;
  }
  else
    return 0;
}

void* do_thread(void* v){
  LV lv = (LV)v;
  GV gv = lv->gv;
  int barrcheck;
  int tid = lv->tid;
  int fiber_mac_rank = gv->num_macs - 1;
  int my_rank = gv->my_rank;
  double t0, t1;
  double t2=0, t2_1=0, t3=0, t3_1=0, t4=0, t5=0, t6=0, tail=0;
  double startTime=0, endTime=0;

  startTime = get_cur_time();

  while (gv->time <= gv->TIME_STOP) {

    t0 = get_cur_time();
    if (my_rank == fiber_mac_rank){
// #ifdef DEBUG_PRINT
      if(tid==0)
        printf("\n\n\n\nStart time step %d ...\n", gv->time);
// #endif //DEBUG_PRINT
      // t0 = get_cur_time();
      compute_bendingforce(lv);

      compute_stretchingforce(lv);

      compute_elasticforce(lv);

      t1 = get_cur_time();
      t2_1 += t1 - t0;

#ifdef DEBUG_PRINT
    // printf("Force calculation ::IB completes by  myrank :%d thread id:%d\n", my_rank, tid);
#endif //DEBUG_PRINT
    }
    else{ //fluid wait here
      pthread_barrier_wait(&(gv->barr));
    }

    //check whether barrier long is becuase of this
    if(tid==0)
      MPI_Barrier(MPI_COMM_WORLD);
    t1 = get_cur_time();
    t2 += t1 - t0;


    if (my_rank == fiber_mac_rank){
      t0 = get_cur_time();
      fiber_SpreadForce(lv);
      t1 = get_cur_time();
      t3_1 += t1- t0;
    }
    else{
#ifdef DEBUG_PRINT
      if(tid==0){
        printf("rank%d: Start get_influence_domain\n", my_rank);
        fflush(stdout);
      }
#endif //DEBUG_PRINT
      t0 = get_cur_time();
      fluid_get_SpreadForce(lv);

      // t1 = get_cur_time();
      // t4 += t1- t0;
    }

    pthread_barrier_wait(&(gv->barr));
    if(tid==0)
      MPI_Barrier(MPI_COMM_WORLD);
    t1 = get_cur_time();
    t3 += t1- t0;
    t4 += t1- t0;
#ifdef DEBUG_PRINT
    if(tid==0){
      printf("mac%d: After get_influence_domain\n", my_rank);
      fflush(stdout);
    }
#endif //DEBUG_PRINT


    if(my_rank != fiber_mac_rank){ //:TODO why this if doesnot work? it hangs all the process

#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fluid mac%d: Start compute DF1\n", my_rank);
        fflush(stdout);
      }
#endif //DEBUG_PRINT
      compute_eqlbrmdistrfuncDF1(lv);
      pthread_barrier_wait(&(gv->barr));
      // if(tid==0)
      //   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fluid mac%d: After compute DF1\n", my_rank);
        fflush(stdout);
      }
#endif //DEBUG_PRINT


#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fluid mac%d: Start streaming\n", my_rank);
        fflush(stdout);
      }
#endif //DEBUG_PRINT
      stream_distrfunc(gv, lv);
      // pthread_barrier_wait(&(gv->barr));
      // if(tid==0)
      //   MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINT
      if(tid==0){
        printf("Fluid mac%d: After streaming\n", my_rank);
        fflush(stdout);
      }

#endif //DEBUG_PRINT

//       //Comment Here
//       bounceback_rigidwalls(lv);
//       pthread_barrier_wait(&(gv->barr));
//       MPI_Barrier(MPI_COMM_WORLD);
// #ifdef DEBUG_PRINT
//       printf("After bounceback\n");
// #endif //DEBUG_PRINT

//       compute_rho_and_u(lv);
//       pthread_barrier_wait(&(gv->barr));
//       MPI_Barrier(MPI_COMM_WORLD);
// #ifdef DEBUG_PRINT
//       printf("After compute rho and u \n");
// #endif //DEBUG_PRINT
    }
    pthread_barrier_wait(&(gv->barr));
    if(tid==0)
        MPI_Barrier(MPI_COMM_WORLD);

    // move_fiber
    if (my_rank != fiber_mac_rank){
      t0 = get_cur_time();
      fluid_SpreadVelocity(lv);
      // t1 = get_cur_time();
      // t5 += t1 - t0;
    }
    else{
      t0 = get_cur_time();
      fiber_get_SpreadVelocity(lv);
      // t1 = get_cur_time();
      // t6 += t1 - t0;
    }
    pthread_barrier_wait(&(gv->barr));
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = get_cur_time();
    t5 += t1 - t0;
    t6 += t1 - t0;
#ifdef DEBUG_PRINT
    printf("rank%d: After moving fibersheet \n", my_rank);
#endif //DEBUG_PRINT

//     copy_inout_to_df2(lv);
//     pthread_barrier_wait(&(gv->barr));
//     if (tid == 0) MPI_Barrier(MPI_COMM_WORLD);
// #ifdef DEBUG_PRINT
//     printf("After  copy_inout_to_df2 \n");
// #endif //DEBUG_PRINT

//     replace_old_DF(lv);
//     pthread_barrier_wait(&(gv->barr));
//     if (tid == 0) MPI_Barrier(MPI_COMM_WORLD);
// #ifdef DEBUG_PRINT
//     printf("After  replace_old_DF \n");
// #endif //DEBUG_PRINT

//     periodicBC(lv);
//     pthread_barrier_wait(&(gv->barr));
//     if (tid == 0) MPI_Barrier(MPI_COMM_WORLD);
// #ifdef DEBUG_PRINT
//     printf("After PeriodicBC\n");
// #endif //DEBUG_PRINT


    if (tid == 0)//does it requires gv->my_rank==0 i.e only one machine updating the counter value
      gv->time += gv->dt;

#ifdef DEBUG_PRINT
    if ( (gv->my_rank == fiber_mac_rank) && (tid == 0)){
      printf("after Running for Timesteps:%d:\n", gv->time - 1);
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
    pthread_barrier_wait(&(gv->barr));
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = get_cur_time();
    tail += t1 - t0;
  }

  endTime = get_cur_time();

  if(my_rank == fiber_mac_rank){
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
