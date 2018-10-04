/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

void moveFiberSheet(LV lv){

  int tid;
  GV gv = lv->gv;

#ifdef DEBUG_PRINT
  printf("************** Inside moveFibersheet at time %d**********\n", gv->time);
#endif DEBUG_PRINT
  // int lookup_fluid_start_x = 20, lookup_fluid_start_y = gv->jb, lookup_fluid_start_z = 10;
  // int lookup_fluid_end_x = 20, lookup_fluid_end_y = gv->je, lookup_fluid_end_z = 20;

  tid = lv->tid;
  int total_threads = gv->total_threads;

  double       s1, s2, s3, dx, dy, dz, dx1, dy1, dz1;
  int          gi, gj, gk, i2, j2, k2, istart, jstart, kstart, jf, kf;
  int          total_fibers_row, total_fibers_clmn;
  Fiber        *fiberarray;


  int i;
  dx = 1.0;
  dy = 1.0;
  dz = 1.0;

  dx1 = 1.0 / dx;
  dy1 = 1.0 / dy;
  dz1 = 1.0 / dz;

  double  PI = 3.14159265358979;
  double  c_x = PI / (2.0*dx);
  double  c_y = PI / (2.0*dy);
  double  c_z = PI / (2.0*dz);
  double  c_64 = 1.0 / (6.4e1);
  double   rx = 0.0;
  double   ry = 0.0;
  double   rz = 0.0;

  Fluidnode *nodes;
  int        BI, BJ, BK, cube_size, cube_idx, node_idx;
  int        num_cubes_x, num_cubes_y, num_cubes_z;
  int        li, lj, lk;

  istart = 0;
  jstart = 0;
  kstart = 0;
  s1 = 0;
  s2 = 0;
  s3 = 0;


  total_fibers_row = gv->fiber_shape->sheets[0].num_rows;
  total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  fiberarray = gv->fiber_shape->sheets[0].fibers;

  cube_size = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  //double s1_buffer[total_fibers_row*64 +1], s2_buffer[total_fibers_row*64 +1], s3_buffer[total_fibers_row*64 +1];
  double *s1_buffer = malloc(((total_fibers_clmn* total_fibers_row) * 64 + 1)* sizeof(double));
  double *s2_buffer = malloc(((total_fibers_clmn* total_fibers_row) * 64 + 1)* sizeof(double));
  double *s3_buffer = malloc(((total_fibers_clmn* total_fibers_row) * 64 + 1)* sizeof(double));

  int position, pos_from_fluid;
  int bufferSize = (gv->num_macs) * 1000;
  char *buffer = malloc(bufferSize*sizeof(buffer));
  // char buffer[1000], buffer_for_fluid[3*(total_fibers_row*64 +1)];
  int itr_stop, itr_stop_from_fluid, owner_mac_rank;

  int my_rank = gv->my_rank;
  int fiber_mac_rank = gv->num_macs - 1;
  MPI_Status status;
  s1_buffer[0] = -1; s2_buffer[0] = -2; s3_buffer[0] = -3;
  int fluid_mac_rank;
  char recv_history[gv->num_macs - 1];
  for (i = 0; i <= gv->num_macs - 1; i++){
    recv_history[i] = '0';
  }

  if (my_rank == fiber_mac_rank){
    for (jf = 0; jf < total_fibers_row; jf++) {
      if (fiber2thread(jf, total_fibers_row, total_threads) == tid){
        for (kf = 0; kf < total_fibers_clmn; kf++){

          //finding influential domain again!!!!!
          istart = floor(fiberarray[jf].nodes[kf].x*dx1 - 2) + 1;
          jstart = floor(fiberarray[jf].nodes[kf].y*dy1 - 2) + 1;
          kstart = floor(fiberarray[jf].nodes[kf].z*dz1 - 2) + 1;

          for (i2 = 0; i2 <= 3; i2++){
            gi = istart + i2;
            rx = dx * gi - fiberarray[jf].nodes[kf].x;/*dimensionless*/
            for (j2 = 0; j2 <= 3; j2++){
              gj = jstart + j2;
              ry = dy * gj - fiberarray[jf].nodes[kf].y;
              for (k2 = 0; k2 <= 3; k2++){
                gk = kstart + k2;
                rz = dz * gk - fiberarray[jf].nodes[kf].z;

                BI = gi / cube_size; BJ = gj / cube_size; BK = gk / cube_size;
                li = gi%cube_size; lj = gj%cube_size; lk = gk%cube_size;
                cube2thread_and_machine(BI, BJ, BK, gv, &owner_mac_rank);

                if (k2 >= 3 && j2 >= 3 && i2 >= 3 && jf == total_fibers_row - 1 && kf == total_fibers_clmn - 1){
                  itr_stop = 1;
                  for (fluid_mac_rank = 0; fluid_mac_rank <= gv->num_macs - 1; fluid_mac_rank++){
                    position = 0;
                    MPI_Pack(&itr_stop, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                    MPI_Pack(&jf, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);//sending fiber node rowindex
                    MPI_Pack(&kf, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);//sending fiber node clmnindex
                    MPI_Pack(&BI, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                    MPI_Pack(&BJ, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                    MPI_Pack(&BK, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                    MPI_Pack(&li, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                    MPI_Pack(&lj, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                    MPI_Pack(&lk, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                    MPI_Pack(&rx, 1, MPI_DOUBLE, buffer, bufferSize, &position, MPI_COMM_WORLD);
                    MPI_Pack(&ry, 1, MPI_DOUBLE, buffer, bufferSize, &position, MPI_COMM_WORLD);
                    MPI_Pack(&rz, 1, MPI_DOUBLE, buffer, bufferSize, &position, MPI_COMM_WORLD);

                    MPI_Send(buffer, bufferSize, MPI_PACKED, fluid_mac_rank, 0, MPI_COMM_WORLD);
                  }
                }
                else{
                  itr_stop = 0;
                  position = 0;
                  MPI_Pack(&itr_stop, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                  MPI_Pack(&jf, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);//sending fiber node rowindex
                  MPI_Pack(&kf, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);//sending fiber node clmnindex
                  MPI_Pack(&BI, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                  MPI_Pack(&BJ, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                  MPI_Pack(&BK, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                  MPI_Pack(&li, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                  MPI_Pack(&lj, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                  MPI_Pack(&lk, 1, MPI_INT, buffer, bufferSize, &position, MPI_COMM_WORLD);
                  MPI_Pack(&rx, 1, MPI_DOUBLE, buffer, bufferSize, &position, MPI_COMM_WORLD);
                  MPI_Pack(&ry, 1, MPI_DOUBLE, buffer, bufferSize, &position, MPI_COMM_WORLD);
                  MPI_Pack(&rz, 1, MPI_DOUBLE, buffer, bufferSize, &position, MPI_COMM_WORLD);

                  MPI_Send(buffer, bufferSize, MPI_PACKED, owner_mac_rank, 0, MPI_COMM_WORLD);
                }

              }
            }
          }


        }// for fiber clmn ends
      }//if fiber2thread ends
    }//for fiber row ends
    //recvng partial sums from all fluid machine


  }
  else {

    //Srring up s buffers for all fluid machines to be 0. another way was inside while commented right now.....
    for (jf = 0; jf < total_fibers_row; jf++)
      for (kf = 0; kf < total_fibers_clmn; kf++){
        s1_buffer[jf*total_fibers_clmn + kf + 1] = 0.0;
        s2_buffer[jf*total_fibers_clmn + kf + 1] = 0.0;
        s3_buffer[jf*total_fibers_clmn + kf + 1] = 0.0;
      }
    //s1_buffer[total_fibers_row*64 +1]=0.0;   s2_buffer[total_fibers_row*64 +1]=0.0;     s3_buffer[total_fibers_row*64 +1]=0.0;
    //if(my_rank== gv->owner_fluid_mac){
    while (true){
      //s1_buffer[70]=0;   s2_buffer[70]=0;     s3_buffer[70]=0;
      MPI_Recv(buffer, bufferSize, MPI_PACKED, fiber_mac_rank, 0, MPI_COMM_WORLD, &status);
      position = 0;
      MPI_Unpack(buffer, bufferSize, &position, &itr_stop, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &jf, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &kf, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &BI, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &BJ, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &BK, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &li, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &lj, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &lk, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &rx, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &ry, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buffer, bufferSize, &position, &rz, 1, MPI_DOUBLE, MPI_COMM_WORLD);

      cube2thread_and_machine(BI, BJ, BK, gv, &fluid_mac_rank);
      recv_history[fluid_mac_rank] = '1';
      if (my_rank == fluid_mac_rank){
        cube_idx = BI *num_cubes_y *num_cubes_z + BJ *num_cubes_z + BK;
        nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
        node_idx = (li)*cube_size*cube_size + (lj)*cube_size + lk;
        s1_buffer[0] = -1; s2_buffer[0] = -2; s3_buffer[0] = -3;

        s1_buffer[jf*total_fibers_clmn + kf + 1] +=
          nodes[node_idx].vel_x*(1.0 + cos(c_x*rx))*(1.0 + cos(c_y*ry))*(1.0 + cos(c_z*rz))*c_64;

        s2_buffer[jf*total_fibers_clmn + kf + 1] +=
          nodes[node_idx].vel_y*(1.0 + cos(c_x*rx))*(1.0 + cos(c_y*ry))*(1.0 + cos(c_z*rz))*c_64;

        s3_buffer[jf*total_fibers_clmn + kf + 1] +=
          nodes[node_idx].vel_z*(1.0 + cos(c_x*rx))*(1.0 + cos(c_y*ry))*(1.0 + cos(c_z*rz))*c_64;

      }
      if (itr_stop == 1)
        break;
      //Settinging up s buffers for all fluid machines to be 0 which are not involved in the above computations.
      /* if(my_rank!=fluid_mac_rank && recv_history[my_rank]!='0'){
      for(jf=0;jf<total_fibers_row;jf++)
      for(kf=0; kf<total_fibers_clmn;kf++){
      s1_buffer[jf*total_fibers_clmn +kf +1]=0.0;
      s2_buffer[jf*total_fibers_clmn +kf +1]=0.0;
      s3_buffer[jf*total_fibers_clmn +kf +1]=0.0;
      }
      }*/
    }//while ends
    //}//to match equal number of send and recvs.
    MPI_Barrier(MPI_COMM_WORLD);
    int bufferlen = (total_fibers_clmn* total_fibers_row) * 64 + 1;

    //Sending buffers
    MPI_Send(s1_buffer, bufferlen, MPI_DOUBLE, fiber_mac_rank, 0, MPI_COMM_WORLD);
    MPI_Send(s2_buffer, bufferlen, MPI_DOUBLE, fiber_mac_rank, 0, MPI_COMM_WORLD);
    MPI_Send(s3_buffer, bufferlen, MPI_DOUBLE, fiber_mac_rank, 0, MPI_COMM_WORLD);


  }
  if (my_rank == fiber_mac_rank) MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == fiber_mac_rank){
    int bufferlen = (total_fibers_clmn* total_fibers_row) * 64 + 1;
    for (i = 0; i < (gv->num_macs - 1); ++i){//reciving from all fluid machines
      MPI_Recv(s1_buffer, bufferlen, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(s2_buffer, bufferlen, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(s3_buffer, bufferlen, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      //MPI_Barrier(MPI_COMM_WORLD);
      for (jf = 0; jf < total_fibers_row; jf++) {
        if (fiber2thread(jf, total_fibers_row, total_threads) == tid){
          for (kf = 0; kf < total_fibers_clmn; kf++){

            if (s1_buffer[0] == -1)
              fiberarray[jf].nodes[kf].x += gv->dt * s1_buffer[jf*total_fibers_row + kf + 1];// +1  since at 0 there is a flag

            if (s2_buffer[0] == -2)
              fiberarray[jf].nodes[kf].y += gv->dt * s2_buffer[jf*total_fibers_row + kf + 1];// +1  since at 0 there is a flag

            if (s3_buffer[0] == -3)
              fiberarray[jf].nodes[kf].z += gv->dt *  s3_buffer[jf*total_fibers_row + kf + 1];// +1  since at 0 there is a flag


            // printf("recvng");
          }
        }
      }
    }
  }

#ifdef DEBUG_PRINT
  printf("**************  moveFibersheet   Exit**********\n");
#endif DEBUG_PRINT
}
