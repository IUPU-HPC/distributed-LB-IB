/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

extern int c[19][3];
extern double t[19];
// static const int c[19][3] = { //{x, y, z}
//             { 0, 0, 0},

//             { 1, 0, 0}, {-1, 0, 0}, { 0, 0, 1}, //1, 2, 3
//             { 0, 0,-1}, { 0,-1, 0}, { 0, 1, 0}, //4, 5, 6
//             { 1, 0, 1}, {-1, 0,-1}, { 1, 0,-1}, //7, 8, 9

//             {-1, 0, 1}, { 0,-1, 1}, { 0, 1,-1}, //10, 11, 12
//             { 0, 1, 1}, { 0,-1,-1}, { 1,-1, 0}, //13, 14, 15
//             {-1, 1, 0}, { 1, 1, 0}, {-1,-1, 0}  //16, 17, 18
// }; // stores 19 different velocity directions for ksi

// static const double t[19] = {
//             1./3.,

//             1./18., 1./18., 1./18., //1~6
//             1./18., 1./18., 1./18.,

//             1./36., 1./36., 1./36., //7~18
//             1./36., 1./36., 1./36.,
//             1./36., 1./36., 1./36.,
//             1./36., 1./36., 1./36.
// };

// compute local equilibrium from rho and u
double computeEquilibrium(int iPop, double rho,
                          double ux, double uy, double uz, double uSqr)
{
    double c_u = c[iPop][0]*ux + c[iPop][1]*uy + c[iPop][2]*uz;
    return rho * t[iPop] * (
               1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr
           );
}

  // initialize a node to its local equilibrium term
void iniEquilibrium(Fluidnode* node) {
    int iPop;
    double rho = node->rho;
    double ux = node->vel_x;
    double uy = node->vel_y;
    double uz = node->vel_z;

    double uSqr = ux*ux + uy*uy + uz*uz;
    for (iPop = 0; iPop < 19; ++iPop) {
        node->dfeq[iPop] =
            computeEquilibrium(iPop, rho, ux, uy, uz, uSqr);
    }
}

void init_eqlbrmdistrfuncDF0(GV gv){/*stored in dfeq*/
  Fluidgrid     *fluidgrid;

  long starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone


  fluidgrid = gv->fluid_grid;

  int cube_size = gv->cube_size;
  long num_cubes_x = gv->fluid_grid->num_cubes_x;
  long num_cubes_y = gv->fluid_grid->num_cubes_y;
  long num_cubes_z = gv->fluid_grid->num_cubes_z;
  int  tmp_task;

  /*PTHREAD_Change*/

  //For each cube BI, BJ ,BK, li, lj, lk 0 to cube size-1
  for (long BI = 0; BI < num_cubes_x; ++BI){
    for (long BJ = 0; BJ < num_cubes_y; ++BJ){
      for (long BK = 0; BK < num_cubes_z; ++BK){
        tmp_task = cube2task(BI, BJ, BK, gv);
        if (gv->taskid == tmp_task){ //MPI changes
          long cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          Fluidnode *nodes = fluidgrid->sub_fluid_grid[cube_idx].nodes;
          starting_x = starting_y = starting_z = 0;
          stopping_x = stopping_y = stopping_z = cube_size - 1;
          for (long li = starting_x; li <= stopping_x; ++li){
            for (long lj = starting_y; lj <= stopping_y; ++lj){
              for (long lk = starting_z; lk <= stopping_z; ++lk){
                long node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
                iniEquilibrium(nodes+node_idx);

                // for (int ksi = 0; ksi<19; ++ksi){
                //   if (ksi == 0){
                //     nodes[node_idx].dfeq[ksi] =
                //       1.0 / 3.0 * nodes[node_idx].rho
                //       * (1.0 - 1.5 *
                //       (nodes[node_idx].vel_x * nodes[node_idx].vel_x
                //       + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                //       + nodes[node_idx].vel_z * nodes[node_idx].vel_z));
                //   }
                //   else if (ksi >0 && ksi<7) {
                //     nodes[node_idx].dfeq[ksi] =
                //       1.0 / 18.0 * nodes[node_idx].rho
                //       *(1.0 + 3.0 *
                //       (gv->c[ksi][0] * nodes[node_idx].vel_x
                //       + gv->c[ksi][1] * nodes[node_idx].vel_y
                //       + gv->c[ksi][2] * nodes[node_idx].vel_z)
                //       + 4.5 *(pow(gv->c[ksi][0] * nodes[node_idx].vel_x
                //       + gv->c[ksi][1] * nodes[node_idx].vel_y
                //       + gv->c[ksi][2] * nodes[node_idx].vel_z, 2))
                //       - 1.5 *(nodes[node_idx].vel_x * nodes[node_idx].vel_x
                //       + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                //       + nodes[node_idx].vel_z * nodes[node_idx].vel_z));
                //   }
                //   else {
                //     nodes[node_idx].dfeq[ksi] =
                //       1.0 / 36.0 * nodes[node_idx].rho
                //       *(1.0 + 3.0*
                //       (gv->c[ksi][0] * nodes[node_idx].vel_x
                //       + gv->c[ksi][1] * nodes[node_idx].vel_y
                //       + gv->c[ksi][2] * nodes[node_idx].vel_z)
                //       + 4.5 *(pow(gv->c[ksi][0] * nodes[node_idx].vel_x
                //       + gv->c[ksi][1] * nodes[node_idx].vel_y
                //       + gv->c[ksi][2] * nodes[node_idx].vel_z, 2))
                //       - 1.5 *(nodes[node_idx].vel_x *nodes[node_idx].vel_x
                //       + nodes[node_idx].vel_y *nodes[node_idx].vel_y
                //       + nodes[node_idx].vel_z *nodes[node_idx].vel_z));
                //   }
                // } //for loop ksi ends

              }//for loop local z
            }//for loop local y
          } // for loop local x
        }//if machine check ends
      }//for BK
    }//for BJ
  }//FOR BI
  /*PTHREAD_Change*/
  // printf("****************init_eqlbrmdistrfuncDF0 EXIT ***************\n");
}
