/* Copyright 2018 Indiana University Purdue University Indianapolis 
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
 * Author: Yuankun Fu (Purdue University, fu121@purdue.edu)
 *         Prateek Nagar (Purdue University, pnagar@iupui.edu)
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>
#include <stdbool.h>

#define DEBUG_PRINT
#define SAVE
#define PI 3.14159265358979
/*
 * ASSUMPTIONS: Cube Size >=4
 * There are no irregular cubes
 * The boundary cubes will be always be first and last cubes for eg BI = BJ = BK = 0 
 * The number of threads n = P*Q*R where n = 2^i, P,Q,R = 2^j, fluid_elem_x = fluid_elem_y, fluid_elem_z = 2^k :FOR SIMPLE Distribution
 * Using Block Distribution on Threads for both 1d (FIber to thread:fiber2thread ) and 3D (Fluid_Sub_cube to thread :cube2thread)
 */


//#define DEBUG_PRINT
int lookup_fluid_start_x = 20, lookup_fluid_start_y = 11, lookup_fluid_start_z = 20;
int lookup_fluid_end_x = 20, lookup_fluid_end_y = 11, lookup_fluid_end_z = 20;
//int lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z ;//for Debug purpose
//pthread_mutex_t mutex_EF_forFluid;
//pthread_mutex_t mutex_EF_forVel;

/* Defining immersed flexible structures consisting of fibers */


/* A node on a fiber */
typedef struct fiber_node_t {
  double x, y, z;
  double bend_force_x, bend_force_y, bend_force_z;
  double stretch_force_x, stretch_force_y, stretch_force_z;
  double elastic_force_x, elastic_force_y, elastic_force_z; // Bendingforce + Stretchingforce
} Fibernode;

/* a single fiber consisting of a number of fiber nodes */ 
typedef struct fiber_t {
  Fibernode* nodes; /* 8 bytes, pointing to an array of fiber nodes*/
  int num_nodes; /* 4 bytes, how many grid points on a fiber*/
} Fiber;

/* one fiber sheet */
typedef struct fiber_sheet_t {
  Fiber* fibers; // i.e., an array of fibers
  double width, height; // Todo: NOT NEEDED. floating point values, width is the dimension in z, height is the dimension in y.
  int num_cols, num_rows;  
  /* bottom left corner: located at <min_y, min_z> */
  double x_orig; // starting point for fibersheet x0 = 20
  double y_orig; // starting point for fibersheet y0 = 21.5
  double z_orig; // starting point for fibersheer z0 = 11.5
} Fibersheet; /* initial configuration of the sheet is given by users!!!*/


/* a set of fiber sheets compose a general IB structure*/
typedef struct fiber_shape_t {
  Fibersheet *sheets;
  int num_sheets; // right now, only one rectangle sheet. Will be more general later.
} Fibershape; // IBShape 


/* Defining a fluid grid */
 
/* A single fluid node */
typedef struct fluid_node_t {
  double vel_x, vel_y, vel_z; // u,v and w given initally as 1 and later computed using rho from eqn 12
  double rho, pressure; // add pressure, a scalar given initially as 1.0 and for next time step computed from df2[] using eqn 11
  double df1[19]; // old value, corresponds to 19 different values of distribution function along 19 different discrete velocities ξj given on FIg1 computed from  df0[] using eq 10
  double df2[19]; // new value, after streaming
  double dfeq[19]; /*???*/            // added since getting errroneous value for timestep>1  That Is, df0 //Todo: no needed
  double df_inout[2][19]; /*???*/    // df_inout[0] for inlet and df_inout[1] for outlet    //TODO: Too Expensive!!
  double elastic_force_x, elastic_force_y, elastic_force_z;  // bending + stretching force, which is the little f in equation 19.
} Fluidnode;

/* A 2-D plane (or surface) of fluid nodes */
typedef struct fluid_surface_t {
  Fluidnode* nodes; //Pointing to a 2-D array of fluid nodes; nodes or 2DarrayOfNodes???
                    //We assume 1st surface is inlet, last one is outlet. Surfaces between is the actual fluid grid.
} Fluidsurface;

/*PTHREAD_Change*/
/* A sub fluid grid adding up together to form the entire fluid grid
The bigger fluid grid composed of smaller cubes of sub fluid grid */
typedef struct sub_fluid_grid_t {
  int sub_coordx_start, sub_coordx_end; // To identify the cube across X : accessed by FluidGrid->SubFluidGrid->sub_x_dim
  int sub_coordy_start, sub_coordy_end; // To identify the cube along Y direction (or #columns) accessed by FluidGrid->SubFluidGrid->sub_y_dim
  int sub_coordz_start, sub_coordz_end; // To identify the cube along Z direction (or #rows)    accessed by FluidGrid->SubFluidGrid->sub_z_dim
  Fluidnode* nodes; // element
  bool isboundary; // To check if this subgrid is part of inlet, outlet or boundary to handle LBM boundary cases
  // =k input from the user, shluld be dim_x, y, z :Some cube might not be k*K*K
  int grid_dim_x;
  int grid_dim_y;
  int grid_dim_z;
//int subgridID;//To identify individual subgrid identified by identified by globalI/cube_size, globalJ/cube_size, globalK/cube_size
} Sub_Fluidgrid;
/*PTHREAD_Change*/

/* A fluid grid is simply a stack of fluid surfaces, where surface itself is a 2D matrix */       

typedef struct fluid_grid_t {
  int x_dim;  // change name to surfaces??  //How many surfaces a long the X direction.
  int y_dim;  // Surface dimension along Y direction (or #columns)
  int z_dim;  // Surface dimension along Z direction (or #rows)
  Fluidsurface*  inlet; // with constant values
  Fluidsurface*  outlet; // with constant values
  /*Added following for Pthread version*/ /*PTHREAD_Change*/
  Sub_Fluidgrid* sub_fluid_grid;
} Fluidgrid;

/* Global info */
typedef struct gv_t {
  /* The immersed structure */
  Fibershape* fiber_shape;
  /* The fluid grid */  
  Fluidgrid* fluid_grid;
  /* Constant parameters used in LBM-IB */
  double tau, nu_l, u_l, rho_l, L_l, M_l, g_l, Ks_l, Kb_l, Kst_l;
  double Re, cs_l, Ma, Kn, Kshat, Ksthat, Kbhat, Mhat, Fr;
  int dt, time, time1, TIME_STOP, TIME_WR, TIME_WR1, N_WR; 
  double c[19][3]; // stores 19 different velocities for ksi
  int ib, ie, je, jb, ke, kb; // For Fluid Grid's actual computation part

/*PTHREAD_Change*/
  pthread_t *threads;
  //pthread_mutex_t lock_Fluid;
  /*Lock for every thread for optimisation*/
  pthread_mutex_t *lock_Fluid;

  int total_threads;
  int cube_size; // Tuning factor k having dimension of the subgrid
  int num_cubes_x, num_cubes_y, num_cubes_z;

  int P,Q,R; // number of threads = P*Q*R :: Used for thread distribtuion

  pthread_barrier_t barr; // to put a bariier after every routine and also in some parts of Fiber's force compuatation
/*PTHREAD_Change*/
}* GV; 


typedef struct lv_t{
  int tid;
  GV  gv;	
}* LV; 


/* Debug help function */
//print range of cubes ip: global i, j,k , iend, jend ,kend
void print_fluid_sub_grid(GV gv, int start_x, int start_y, int start_z, 
                          int end_x, int end_y, int end_z,
                          int cube_size);

//print one cube ip BI, BJ ,BK
void print_fluid_cube(GV gv, int BI, int BJ, int BK, 
                      int BI_end, int BJ_end, int BK_end,
                      int cube_size);

/* Debug help function */
void print_fiber_sub_grid(GV gv, int start_y, int start_z,
         			      int end_y, int end_z);
/* Debug help function */
double get_cur_time();

/* Requires following 2 more Debug Help Functions
 * Cube Distribution Function: Map a cube to Threadid 0 to N-1
 * Distribtion is Arbitary Surface or Block Cyclic 
 * Distribuion is Block based where every thread is working on a block  
 * # of threads = P*Q*R--
 */
int cube2thread(int BI, int BJ, int BK, 
                int num_cubes_x, int num_cubes_y, int num_cubes_z,
                int P, int Q, int R);
 
/*Fiber Distribution Function: Map a fiber row to each thread
  fiber_row -ith row*/
int fiber2thread(int fiber_row, int num_fibers, int num_threads);

/* Initiallize the global GV */ /*PTHREAD_Change*/
void init_gv(GV, Fibershape*, Fluidgrid*, int cube_size); //TODO: pass all necessary information from users to initialize GV.
 
/* Create a fiber shape (user defined arbitray shape)*/
Fibershape* gen_fiber_shape(GV gv, double width, double height, 
                            int num_cols, int num_rows, 
                            double x0, double y0, double z0);

/* Create a 3D fluid grid */ /*PTHREAD_Change*/
Fluidgrid* gen_fluid_grid(GV gv, int elem_x, int elem_y, int elem_z, int sub_fluid_grid_dim_k);

/* Methods for IB computations */
void compute_bendingforce(LV);    // eqn 18
void compute_stretchingforce(LV); // eqn 16
void compute_elasticforce(LV);    // F for eqn 19
void get_influentialdomain_fluid_grid_and_SpreadForce(LV); // eqn 19

/* Methods for LBM computations */
void init_eqlbrmdistrfuncDF0(GV);    // eqn13 for df0 or dfeq
void init_df_inout(GV);              // used by inlet, outlet, before the simuation starts.
void compute_eqlbrmdistrfuncDF1(LV); // eqn10 for df1
void stream_distrfunc(LV);           // eqn10 for df2
void bounceback_rigidwalls(LV);      // after streaming, special handling of boundary fluid nodes  
void compute_rho_and_u(LV);          // eqn 11 and eqn 12 for rho and velocity
void moveFiberSheet(LV);             // eqn 22 to move fiber sheet, also  where eqn 20 is implemented.
void copy_inout_to_df2(LV);          // computed in every time step
void periodicBC(LV);                 // boundary condition again?? It is used but does not affect result.?
void replace_old_DF(LV);

/*TODO: move everything above to a separate header file */

void* do_thread(void* v); //IB-LBM simulation moved to do_thread called via ptheread_create

void print_fluid_sub_grid(GV gv, int start_x, int start_y, int start_z, 
                          int end_x, int end_y, int end_z,
                          int cube_size) { // start, stop are closed: inclusive
  
  Fluidgrid *grid;
  Sub_Fluidgrid *sub_grid;
  int li, lj, lk, BI, BJ, BK, BI_start, BJ_start, BK_start, BI_end, BJ_end, BK_end;
  int ksi, cube_idx, num_cubes_y, num_cubes_z;

  Fluidnode *node;
  
  grid        = gv->fluid_grid;
  sub_grid    = grid->sub_fluid_grid;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  BI_start = start_x/cube_size;
  BJ_start = start_y/cube_size;
  BK_start = start_z/cube_size;

  BI_end = end_x/cube_size;
  BJ_end = end_y/cube_size;
  BK_end = end_z/cube_size;
  
  /*li     = start_x%cube_size;
  lj     = start_y%cube_size;
  lk     = start_z%cube_size;
  li_end = end_x%cube_size;
  lj_end = end_y%cube_size;
  lk_end = end_z%cube_size;*/

  printf("(BI,BJ,BK): {vel_x, vel_y, vel_z} ||  {G0, DF1, DF2}|| rho || {ElasticF_x, y, z}\n");
 
  for (BI = BI_start; BI <= BI_end; ++BI) 
    for (BJ = BJ_start; BJ <= BJ_end; ++BJ) 
      for (BK = BK_start; BK <= BK_end; ++BK) {
        cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        for (li = 0; li < cube_size ;li++)
          for (lj = 0; lj < cube_size; lj++)
            for (lk = 0; lk < cube_size; lk++){
              node = &sub_grid[cube_idx].nodes[li*cube_size*cube_size+ lj*cube_size+lk];
              //if( li == start_x%cube_size && lj ==start_y%cube_size && lk ==start_z%cube_size ){
              printf("For cube <%d, %d, %d>\n", BI, BJ, BK);
              for (ksi = 0 ; ksi < 19; ksi++){
                printf("(%d, %d, %d, %d):{%.12f, %.12f, %.12f} || {%.12f, %.12f, %.12f} || %.12f || {%.12f, %.12f, %.12f}\n", 
                        li, lj, lk, ksi, node->vel_x, node->vel_y, node->vel_z,
                        node->dfeq[ksi], node->df1[ksi],node->df2[ksi] ,
                        node->rho,
                        node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
              }
              //} 
            }
       }   

   printf("\n");
}


void print_fluid_cube(GV gv, int BI_start, int BJ_start, int BK_start,
                      int BI_end, int BJ_end, int BK_end,
                      int cube_size) { // start, stop are closed: inclusive
  
  Fluidgrid *grid;
  Sub_Fluidgrid *sub_grid;
  int li, lj, lk, BI, BJ, BK;
  int ksi, cube_idx, num_cubes_y, num_cubes_z;

  Fluidnode *node;
  
  grid        = gv->fluid_grid;
  sub_grid    = grid->sub_fluid_grid;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;
 
  printf("(BI,BJ,BK): {vel_x, vel_y, vel_z} ||  {G0, DF1, DF2}|| rho || {ElasticF_x, y, z}\n");
  for (BI = BI_start; BI <= BI_end; ++BI) 
    for (BJ = BJ_start; BJ <= BJ_end; ++BJ) 
      for(BK = BK_start; BK <= BK_end; ++BK) {
        cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        for (li = 0; li < cube_size ;li++)
          for (lj = 0; lj < cube_size; lj++)
            for (lk = 0; lk < cube_size; lk++){
              node = &sub_grid[cube_idx].nodes[li*cube_size*cube_size+ lj*cube_size+lk];
              printf("FOr cube <%d,%d,%d>\n",BI, BJ, BK);
              for (ksi = 0 ; ksi < 19; ksi++){
                printf("(%d, %d, %d, %d):{%.12f, %.12f, %.12f} || {%.12f, %.12f, %.12f} || %.12f || {%.12f, %.12f, %.12f}\n", 
                       li, lj, lk, ksi, node->vel_x, node->vel_y, node->vel_z,
                       node->dfeq[ksi], node->df1[ksi],node->df2[ksi],
                       node->rho,
                       node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
              }
            } 
      }

  printf("\n");
}

// save one subgrid of fluid to disk
void save_fluid_sub_grid(GV gv, int start_x, int start_y, int start_z,
                         int end_x, int end_y, int end_z, char fName[]) {
  FILE* oFile = fopen(fName, "w");

  Fluidgrid *grid;
  Sub_Fluidgrid *sub_grid;
  int li, lj, lk, BI, BJ, BK, BI_start, BJ_start, BK_start, BI_end, BJ_end, BK_end;
  int ksi, cube_idx, num_cubes_y, num_cubes_z;
  int temp_taskid;
  Fluidnode *node;
  int cube_size = gv->cube_size;

  grid        = gv->fluid_grid;
  sub_grid    = grid->sub_fluid_grid;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  BI_start = start_x / cube_size;
  BJ_start = start_y / cube_size;
  BK_start = start_z / cube_size;

  BI_end = end_x / cube_size;
  BJ_end = end_y / cube_size;
  BK_end = end_z / cube_size;

  /*li     = start_x%cube_size;
  lj     = start_y%cube_size;
  lk     = start_z%cube_size;
  li_end = end_x%cube_size;
  lj_end = end_y%cube_size;
  lk_end = end_z%cube_size;*/

  int my_rank = 0;

  for (BI = BI_start; BI <= BI_end; ++BI)
    for (BJ = BJ_start; BJ <= BJ_end; ++BJ)
      for (BK = BK_start; BK <= BK_end; ++BK) {
        fprintf(oFile, "(BI,BJ,BK): {vel_x, vel_y, vel_z} || {G0, DF1, DF2}|| rho || {ElasticF_x, y, z}\n");

        cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        for (li = 0; li < cube_size; li++)
          for (lj = 0; lj < cube_size; lj++)
            for (lk = 0; lk < cube_size; lk++){
              node = &sub_grid[cube_idx].nodes[li*cube_size*cube_size + lj*cube_size + lk];

              fprintf(oFile, "For cube <%d,%d,%d> Rank %d\n", BI, BJ, BK, my_rank);
              int X = BI * cube_size + li;
              int Y = BJ * cube_size + lj;
              int Z = BK * cube_size + lk;
              for (ksi = 0; ksi < 19; ksi++){
                fprintf(oFile, "(%d,%d,%d, %2d):{%.6f,%.6f,%.6f} || {%.12f,%.12f,%.12f} || %.6f || {%.12f,%.12f,%.24f}\n",
                  X, Y, Z, ksi, node->vel_x, node->vel_y, node->vel_z,
                  node->dfeq[ksi], node->df1[ksi], node->df2[ksi],
                  node->rho,
                  node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
              }
            }
  }

  fclose(oFile);
}

void print_fiber_sub_grid(GV gv, int start_y, int start_z,
			                    int end_y, int end_z) {
  /*Assuming one fiber sheet!!*/
  Fiber     *fiber_array;
  int        i, j;
  Fibernode *node;
 
  // TODO: 1 fiber now
  fiber_array = gv->fiber_shape->sheets[0].fibers;

  printf("(i,j): {cord_x, cord_y, cord_z} || SF_x,SF_y, SF_z || BF_X, BF_y, BF_z\n");

  for (i = start_y; i <= end_y; ++i) {
    for (j = start_z; j <= end_z; ++j) {
      node = fiber_array[i].nodes + j;
      printf("(%d, %d):{%f, %f, %f} || %.12f, %.12f, %.12f || %.12f, %.12f, %.12f\n", 
              i, j,  node->x, node->y, node->z,
	            node->stretch_force_x, node->stretch_force_y, node->stretch_force_z, 
              node->bend_force_x, node->bend_force_y, node->bend_force_z);
    }
    printf("\n");
  }
}

// save one lattice population to disk
void save_fiber_sub_grid(GV gv, int start_y, int start_z,
                         int end_y, int end_z, char fName[]) {
  FILE* oFile = fopen(fName, "w");

    /*Assuming one fiber sheet!!*/
  Fiber     *fiber_array;
  int        i, j, k;
  Fibernode *node;
  
  k = 0;
  fiber_array = gv->fiber_shape->sheets[k].fibers;

  fprintf(oFile, "Fiber_sheets[%d] (i,j): {cord_x, cord_y, cord_z} || SF_x, SF_y, SF_z || BF_x, BF_y, BF_z || EF_x, EF_y, EF_z\n", k);
  for (i = start_y; i <= end_y; ++i) {
    for (j = start_z; j <= end_z; ++j) {
      node = fiber_array[i].nodes + j;
      fprintf(oFile, "(%2d,%2d):{%f,%f,%f} || %.6f,%.24f,%.24f || %.6f,%.24f,%.24f || %.6f,%.24f,%.24f\n",
        i, j, node->x, node->y, node->z,
        node->stretch_force_x, node->stretch_force_y, node->stretch_force_z,
        node->bend_force_x, node->bend_force_y, node->bend_force_z,
        node->elastic_force_x, node->elastic_force_y, node->elastic_force_z);
    }
  }

  fclose(oFile);
}

double get_cur_time() {
  struct timeval   tv;
  struct timezone  tz;
  double cur_time;

  gettimeofday(&tv, &tz);
  cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

  return cur_time;
}

void init_gv(GV gv, Fibershape* fiber_shape, Fluidgrid* fluid_grid, int cube_size) {
  int i, j, k;
  int node_idx, subgrid_idx;
  int dim_x, dim_y, dim_z;
  int total_sub_grids;
  Sub_Fluidgrid *sub_fluid_grid;
  // int num_cubes_x, num_cubes_y, num_cubes_z;

  if(fiber_shape == NULL) {
    fprintf(stderr, "Error(%s): Null fiber shape!\n", __func__); exit(1);
  }
  if(fluid_grid == NULL) {
    fprintf(stderr, "Error(%s): Null fluid grid!\n", __func__); exit(1);
  }
  
  gv->fiber_shape = fiber_shape;
  gv->fluid_grid  = fluid_grid;
  sub_fluid_grid  = fluid_grid->sub_fluid_grid;
  gv->cube_size   = cube_size;
  //nodes = gv->fluid_grid->sub_fluid_grid->nodes;

  dim_x = gv->fluid_grid->x_dim;
  dim_y = gv->fluid_grid->y_dim;
  dim_z = gv->fluid_grid->z_dim;

  total_sub_grids = (dim_x*dim_y*dim_z) / pow(cube_size,3);

  /*gv->num_cubes_y     = total_sub_grids/cube_size;
  gv->num_cubes_z     = total_sub_grids/cube_size;
  gv->num_cubes_x     = total_sub_grids/cube_size;*/

  gv->num_cubes_x = dim_x/cube_size;
  gv->num_cubes_y = dim_y/cube_size;
  gv->num_cubes_z = dim_z/cube_size;
  // num_cubes_x     = gv->num_cubes_x;
  // num_cubes_y     = gv->num_cubes_y;
  // num_cubes_z     = gv->num_cubes_z;

  gv->ib = 2;       // 0 and 1 are used for buffer zone
  gv->ie = dim_x-3; // dim_x-1 (i.e., the last node), dim_x-2 are used for buffer zone 
  gv->jb = 2;
  gv->je = dim_y-3;
  gv->kb = 2;
  gv->ke = dim_z-3;

  gv->Re    = 1.5e2;
  gv->rho_l = 1.0e0;
  // gv->u_l   = 0.001; /* choice of u_l should make Ma <0.1 */
  gv->u_l = 0.08;
  
  //Todo: move it to fiber sheet!!
  gv->L_l   = 2.0e1;       /*the dimensionless char LENGTH, should shortest fiber or width of the sheet should be fibersheet_w*/
  
  gv->nu_l  = gv->u_l * gv->L_l / gv->Re; //viscosity
  gv->tau   = 0.5 * (6.0 * gv->nu_l + 1.0);
  gv->Fr    = 1.0e10;        /* huge number of Fr means gra is not significant */
  //NO IB if Kbhat and Kshat = 0.0
  gv->Kbhat = 0.005; /* Dimensionless flexure Modulus Kb(.0001-.05) Prateek */
  gv->Kshat = 2.0e1; /*Streching Compression Coefficient of fiber Kst(20) Prateek */
  gv->Ksthat= 0.0e0; //TODO:add Ksthat value    /*relative stiffness of the virtual spring tethered to the fixed point*/
  gv->cs_l  = 1.0 / (sqrt(3)); // non-dimensional cs is the speed of the sound of the model and c is taken to be 1....Why other values of j(from paper not considered) Prateek

  /* calculate the values of M, Kb, Ks, and g in LB world */
  gv->g_l = gv->u_l * gv->u_l / (gv->Fr * gv->L_l);  /*How...? PN*/
  gv->g_l = 0.0; // gravity is equal to zero.

  gv->Kb_l  = gv->Kbhat * gv->rho_l * gv->u_l * gv->u_l * pow(gv->L_l,4); // For bending Force
  //Not used any more.gv->Kst_l = gv->Ksthat * gv->rho_l * gv->u_l * gv->u_l * gv->L_l;       // For stretching Force may be used fo tethered points only
  gv->Ks_l = gv->Kshat * gv->rho_l * gv->u_l * gv->u_l * gv->L_l;
  //printf("GV->Ks_l:::::%f\n",gv->Ks_l);
  
  /* discrete particle velocity, 19 velocities: see Figure 1 */
  gv->c[0][0]=0.0;
  gv->c[0][1]=0.0;
  gv->c[0][2]=0;   //direction 0
  gv->c[1][0]=1.0;
  gv->c[1][1]=0.0;
  gv->c[1][2]=0.0; //direction 1
  gv->c[2][0]=-1.0;
  gv->c[2][1]=0.0;
  gv->c[2][2]=0.0;
  gv->c[3][0]=0.0;
  gv->c[3][1]=0.0;
  gv->c[3][2]=1;
  gv->c[4][0]=0.0;
  gv->c[4][1]=0.0;
  gv->c[4][2]=-1;
  gv->c[5][0]=0.0;
  gv->c[5][1]=-1.0;
  gv->c[5][2]=0.0;
  gv->c[6][0]=0.0;
  gv->c[6][1]=1.0;
  gv->c[6][2]=0.0;
  gv->c[7][0]=1.0;
  gv->c[7][1]=0.0;
  gv->c[7][2]=1;
  gv->c[8][0]=-1.0;
  gv->c[8][1]=0.0;
  gv->c[8][2]=-1;
  gv->c[9][0]=1.0;
  gv->c[9][1]=0.0;
  gv->c[9][2]=-1.0;
  gv->c[10][0]=-1.0;
  gv->c[10][1]=0.0;
  gv->c[10][2]=1.0;
  gv->c[11][0]=0.0;
  gv->c[11][1]=-1.0;
  gv->c[11][2]=1.0;
  gv->c[12][0]=0.0;
  gv->c[12][1]=1.0;
  gv->c[12][2]=-1.0;
  gv->c[13][0]=0.0;
  gv->c[13][1]=1.0;
  gv->c[13][2]=1.0;
  gv->c[14][0]=0.0;
  gv->c[14][1]=-1.0;
  gv->c[14][2]=-1.0;
  gv->c[15][0]=1.0;
  gv->c[15][1]=-1.0;
  gv->c[15][2]=0.0;
  gv->c[16][0]=-1.0;
  gv->c[16][1]=1.0;
  gv->c[16][2]=0.0;
  gv->c[17][0]=1.0;
  gv->c[17][1]=1.0;
  gv->c[17][2]=0.0;
  gv->c[18][0]=-1.0;
  gv->c[18][1]=-1.0;
  gv->c[18][2]=0.0; //velocity 18

  // initiallization of rho and vel for fluid grid points
  /*PTHREAD_Change*/
  for (subgrid_idx = 0; subgrid_idx < total_sub_grids; subgrid_idx++){
    for (i = 0; i < cube_size; ++i){ // dim_x = NX and n1 = NX -1 therfore, changed from dim_x to dim_x-1
      for (j = 0; j < cube_size; ++j){    
        for (k = 0; k < cube_size; ++k){ 
          node_idx =  i * cube_size * cube_size + j * cube_size + k;
          sub_fluid_grid[subgrid_idx].nodes[node_idx].rho = gv->rho_l;//same as rho_l
          sub_fluid_grid[subgrid_idx].nodes[node_idx].vel_x = gv->u_l;
          sub_fluid_grid[subgrid_idx].nodes[node_idx].vel_y = 0.0;
          sub_fluid_grid[subgrid_idx].nodes[node_idx].vel_z = 0.0;
        }
      }
    }
  }
 /*PTHREAD_Change*/
  
  //gv->N_WR = gv->TIME_STOP / gv->TIME_WR + 1;//param for Cd need to add later....

  gv->dt   = 1;      // Time step size, always equal 1!    
  gv->time = gv->dt; // Time step#, starting from 1.

  /******* Shift fiber sheet for one time step ************/
  //TODO: for all the fiber sheets, shift. 
  for (i = 0; i < gv->fiber_shape->sheets[0].num_rows; i++) {
    for (j = 0; j < gv->fiber_shape->sheets[0].num_cols; j++){
      gv->fiber_shape->sheets[0].fibers[i].nodes[j].x += gv->u_l*gv->dt;
      //TODO: if v_l, w_l !=0, then ADD .y, .z = v_l, w_l times gv->dt.
     }
  }
}

/* Assming ONLY one sheet at this moment */
Fibershape*  gen_fiber_shape(GV gv, double w, double h, int num_cols, int num_rows,
                             double x_orig, double y_orig, double z_orig){
  int         i, row, col;
  Fibershape  *shape = gv->fiber_shape;
  Fibersheet  *sheet;
  Fiber       *fiber;
  
  double      zdist =  w / (num_cols - 1);
  double      ydist =  h / (num_rows - 1);
  /*Fiber sheet is kept on yz plane with exactly midway between fluid grid , therfore added +30 */

  /* allocate memory for fiber_shape.*/
  // shape = (Fibershape*) malloc(sizeof(*shape));

  // TODO: Modify to multi sheets
  // allocate memory for fiber sheets. assumption: Only One Sheet For Now!!!!!
  // shape->sheets = (Fibersheet*) malloc(sizeof(*shape->sheets) * shape->num_sheets); 

  /* allocate memory for fibers in each sheet */
  for (i = 0; i < shape->num_sheets; ++i) {
    sheet = shape->sheets + i;  //The i-th fiber sheet
    sheet->num_cols = num_cols;
    sheet->num_rows = num_rows;
    sheet->width    = w;
    sheet->height   = h;
    sheet->x_orig   = x_orig;
    sheet->y_orig   = y_orig;
    sheet->z_orig   = z_orig;
    sheet->fibers   = (Fiber*) calloc(num_rows, sizeof(*sheet->fibers));
    for (row = 0; row < num_rows; ++row) {	
      fiber = sheet->fibers + row;       // The row-th fiber
      fiber->nodes = (Fibernode*) calloc(num_cols, sizeof(*fiber->nodes));
      /* Fill in the contents of the points in each fiber */
      for (col = 0; col < num_cols; ++col){
        /* z-coordinate of sheet is from (center of z-direction - sheet height) to (center of z-direction);
           y-coordinate of sheet is from (center of y-direction - half of sheet width) to (center of y-direction+half of sheet width) */
        fiber->nodes[col].x = x_orig;               // perpendicular to x-axis.
        fiber->nodes[col].y = ydist * row + y_orig; // ydist * i + 30;
        fiber->nodes[col].z = zdist * col + z_orig; // zdist * j + 30;// middle of the fluid grid
        // If fiber nodes' coordinates are given by a user, remove the above 3 lines.
      }	
    }
    /* For Tethered Points */
    for (row = num_rows-1; row > num_rows-3; --row) {
      fiber = sheet->fibers + row;       // The row-th fiber
      for (col = 0; col < num_cols; ++col){
        //TODO: Tether points
        fiber->nodes[col].x = fiber->nodes[col].x;
        fiber->nodes[col].y = fiber->nodes[col].y;
        fiber->nodes[col].z = fiber->nodes[col].z;
      }
    }
  }

  return shape;
}		


Fluidgrid* gen_fluid_grid(GV gv, int dim_x, int dim_y, int dim_z, int cube_size){
  
  Fluidgrid *fluid_grid = gv->fluid_grid;
  int idx, total_sub_grids; //total_sub_grids_x, total_sub_grids_y, total_sub_grids_z; 

  /* calculate total # of cubes */
  total_sub_grids   = (dim_x * dim_y * dim_z) / (cube_size * cube_size * cube_size);
  // total_sub_grids_x = dim_x / cube_size; // shoould be sub_grid_dimx for later
  // total_sub_grids_y = dim_y / cube_size;
  // total_sub_grids_z = dim_z / cube_size;

  /* allocate memory for the entire fluid grid */
  fluid_grid = (Fluidgrid *) malloc(sizeof(*fluid_grid));
  fluid_grid->x_dim = dim_x;
  fluid_grid->y_dim = dim_y;
  fluid_grid->z_dim = dim_z;

  /* allocate Memory for sub_fluid grid */
  fluid_grid->sub_fluid_grid = (Sub_Fluidgrid*) malloc(sizeof(*fluid_grid->sub_fluid_grid)*total_sub_grids);
  
  /* allocate memory for the surface structure */
  fluid_grid->inlet  = (Fluidsurface*) malloc(sizeof(*fluid_grid->inlet));
  fluid_grid->outlet = (Fluidsurface*) malloc(sizeof(*fluid_grid->outlet));

   /* allocate memory for the surface's fluid nodes */
  fluid_grid->inlet->nodes  =  (Fluidnode*) calloc(dim_z * dim_y, sizeof(*fluid_grid->inlet->nodes));
  fluid_grid->outlet->nodes =  (Fluidnode*) calloc(dim_z * dim_y, sizeof(*fluid_grid->outlet->nodes));

  for (idx = 0; idx < total_sub_grids; ++idx){
    fluid_grid->sub_fluid_grid[idx].nodes = 
      (Fluidnode*) calloc(cube_size * cube_size * cube_size, sizeof(*fluid_grid->sub_fluid_grid->nodes));
  }  
  
  return fluid_grid;
}
		
	
void compute_bendingforce(LV lv) { 
  #ifdef DEBUG_PRINT
    printf("Tid%d: ****Inside compute_bendingforce*******\n", lv->tid);
  #endif //DEBUG_PRINT

  int         i, j;
  double      ds1, ds2;
  double      bending_const;
  Fibershape  *fiber_shape; // aka IB shape
  Fiber*      fibers;       // an array of fibers
  double      fiber_sheet_w, fiber_sheet_h;
  int         total_fibers_clmn;
  int         total_fibers_row;
  Fibernode   *nodes;
  int         tid;
  int         total_threads;
  
  GV gv = lv->gv;   
  tid   = lv->tid;
  // pthread_barrier_t barr = gv->barr;
  // int barrcheck;
  
  total_threads       = gv->total_threads;
  fiber_shape         = gv->fiber_shape;

  // TODO: Modify to multi sheets
  /* Assuming one sheet for now */
  fibers              = fiber_shape->sheets[0].fibers;
  total_fibers_row    = fiber_shape->sheets[0].num_rows;
  total_fibers_clmn   = fiber_shape->sheets[0].num_cols;
  fiber_sheet_w       = fiber_shape->sheets[0].width;  // z-direction
  fiber_sheet_h       = fiber_shape->sheets[0].height; // y-direction 

  ds1                 = fiber_sheet_w / (total_fibers_clmn -1); // gap between fiber nodes along the z-dim.
  ds2                 = fiber_sheet_h / (total_fibers_row-1);   // gap between two neighboring fibers along y-dim.
  bending_const       = gv->Kb_l / pow(ds2, 3);           // KB and alpha1 pow 3

  /* compute force along the fibers (row) */
  for (i = 0; i < total_fibers_row; ++i) {
    if (fiber2thread(i, total_fibers_row, total_threads) == tid){
    //printf("Tid:%d\n",tid);
    nodes = fibers[i].nodes; // an array of fiber nodes in the i-th fiber ???

    // For the 1st point
    nodes[0].bend_force_x = bending_const * (-nodes[2].x + 2 * nodes[1].x - nodes[0].x);
    nodes[0].bend_force_y = bending_const * (-nodes[2].y + 2 * nodes[1].y - nodes[0].y);
    nodes[0].bend_force_z = bending_const * (-nodes[2].z + 2 * nodes[1].z - nodes[0].z);

    // For the 2nd Point
    nodes[1].bend_force_x = bending_const * (-nodes[3].x + 4 * nodes[2].x - 5 * nodes[1].x 
	    				    + 2 * nodes[0].x);
    nodes[1].bend_force_y = bending_const * (-nodes[3].y + 4 * nodes[2].y - 5 * nodes[1].y 
					        + 2 * nodes[0].y);
    nodes[1].bend_force_z = bending_const * (-nodes[3].z + 4 * nodes[2].z - 5 * nodes[1].z 
					        + 2 * nodes[0].z);

    // For the middle points
    for(j = 2; j < total_fibers_clmn - 2; ++j){
      nodes[j].bend_force_x = bending_const * (-nodes[j+2].x + 4.0 * nodes[j+1].x 
					          - 6.0 * nodes[j].x + 4.0 * nodes[j-1].x
					          - nodes[j-2].x);
      nodes[j].bend_force_y = bending_const * (-nodes[j+2].y + 4.0 * nodes[j+1].y 
					          - 6.0 * nodes[j].y + 4.0 * nodes[j-1].y 
					          - nodes[j-2].y);
      nodes[j].bend_force_z = bending_const * (-nodes[j+2].z + 4.0 * nodes[j+1].z 
					          - 6.0 * nodes[j].z + 4.0 * nodes[j-1].z 
					          - nodes[j-2].z);
    }

    // For the last but one point
    j = total_fibers_clmn-2;
    nodes[j].bend_force_x = bending_const * (2 * nodes[j+1].x 
					        - 5 * nodes[j].x 
					        + 4 * nodes[j-1].x 
					        - nodes[j].x);
    nodes[j].bend_force_y = bending_const * (2 * nodes[j+1].y 
					        - 5 * nodes[j].y 
					        + 4 * nodes[j-1].y 
					        - nodes[j-2].y);
    nodes[j].bend_force_z = bending_const * (2 * nodes[j+1].z 
					        - 5 * nodes[j].z 
					        + 4 * nodes[j-1].z 
					        - nodes[j-2].z);
    // For last ponit
    j = total_fibers_clmn-1;
    nodes[j].bend_force_x = bending_const * (-nodes[j].x 
					        + 2 * nodes[j-1].x
					        - nodes[j-2].x);
    nodes[j].bend_force_y = bending_const * (-nodes[j].y 
					        + 2 * nodes[j-1].y
					        - nodes[j-2].y);
    nodes[j].bend_force_z = bending_const * (-nodes[j].z 
					        + 2 * nodes[j-1].z
					        - nodes[j-2].z);
    } // if fiber2thread ends
  } 
  
  pthread_barrier_wait(&(gv->barr));
     
  #ifdef DEBUG_PRINT
    printf("Tid%d: After barrier in BF\n", tid);
  #endif //DEBUG_PRINT
  /* computing bending force along direction normal to fibers (column), bending_const may be different */									
  Fibersheet* sheet = fiber_shape->sheets + 0;
  bending_const     = gv->Kb_l / pow(ds2, 3); // KB and alpha1 pow 3

  for (j = 0; j < total_fibers_clmn; ++j){	
    if (fiber2thread(0, total_fibers_row, total_threads) == tid){	
      // For first fiber
      sheet->fibers[0].nodes[j].bend_force_x += bending_const*(-sheet->fibers[2].nodes[j].x 
                                              + 2 * (sheet->fibers[1].nodes[j].x) 
                                              - (sheet->fibers[0].nodes[j].x));
      sheet->fibers[0].nodes[j].bend_force_y += bending_const*(-sheet->fibers[2].nodes[j].y 
                                              + 2 * (sheet->fibers[1].nodes[j].y) 
							                  - (sheet->fibers[0].nodes[j].y));
      sheet->fibers[0].nodes[j].bend_force_z += bending_const*(-sheet->fibers[2].nodes[j].z 
							                  + 2 * (sheet->fibers[1].nodes[j].z) 
							                  - (sheet->fibers[0].nodes[j].z));
    } // if fiber2thread ends

    // For 2nd fiber
    if (fiber2thread(1, total_fibers_row, total_threads) == tid ){
      sheet->fibers[1].nodes[j].bend_force_x += bending_const*(-sheet->fibers[3].nodes[j].x 
							                  + 4 * (sheet->fibers[2].nodes[j].x) 
							                  - 5 * (sheet->fibers[1].nodes[j].x) 
							                  + 2 * (sheet->fibers[0].nodes[j].x));
      sheet->fibers[1].nodes[j].bend_force_y += bending_const*(-sheet->fibers[3].nodes[j].y 
							                  + 4 * (sheet->fibers[2].nodes[j].y) 
							                  - 5 * (sheet->fibers[1].nodes[j].y) 
							                  + 2 * (sheet->fibers[0].nodes[j].y));
      sheet->fibers[1].nodes[j].bend_force_z += bending_const*(-sheet->fibers[3].nodes[j].z 
							                  + 4 * (sheet->fibers[2].nodes[j].z) 
							                  - 5 * (sheet->fibers[1].nodes[j].z) 
							                  + 2 * (sheet->fibers[0].nodes[j].z));
    } // if fiber2thread ends
    
    // For the middle fibers
    for (i = 2; i < total_fibers_row-2; ++i){
      if (fiber2thread(i, total_fibers_row, total_threads) == tid){
        sheet->fibers[i].nodes[j].bend_force_x += bending_const*(-sheet->fibers[i+2].nodes[j].x 
             							        + 4.0 * (sheet->fibers[i+1].nodes[j].x) 
							                    - 6.0 * (sheet->fibers[i].nodes[j].x)
							                    + 4.0 * (sheet->fibers[i-1].nodes[j].x)
							                    - (sheet->fibers[i-2].nodes[j].x));
        sheet->fibers[i].nodes[j].bend_force_y += bending_const*(-sheet->fibers[i+2].nodes[j].y 
                                                + 4.0 * (sheet->fibers[i+1].nodes[j].y) 
                                                - 6.0 * (sheet->fibers[i].nodes[j].y) 
                                                + 4.0 * (sheet->fibers[i-1].nodes[j].y)
                                                - (sheet->fibers[i-2].nodes[j].y));
        sheet->fibers[i].nodes[j].bend_force_z += bending_const*(-sheet->fibers[i+2].nodes[j].z 
                                                + 4.0 * (sheet->fibers[i+1].nodes[j].z) 
                                                - 6.0 * (sheet->fibers[i].nodes[j].z) 
                                                + 4.0 * (sheet->fibers[i-1].nodes[j].z)
                                                - (sheet->fibers[i-2].nodes[j].z));
     } // if fiber2thread ends
    }

    // For last but one fiber
    i = total_fibers_row-2;
    if (fiber2thread(i, total_fibers_row, total_threads) == tid ){
      sheet->fibers[i].nodes[j].bend_force_x += bending_const * (2 * sheet->fibers[i+1].nodes[j].x 
            							      - 5 * (sheet->fibers[i].nodes[j].x) 
			            				      + 4 * (sheet->fibers[i-1].nodes[j].x)
						            	      - (sheet->fibers[i-2].nodes[j].x));
      sheet->fibers[i].nodes[j].bend_force_y += bending_const * (2*(sheet->fibers[i+1].nodes[j].y) 
							                  - 5 * (sheet->fibers[i].nodes[j].y) 
							                  + 4 * (sheet->fibers[i-1].nodes[j].y)
							                  - (sheet->fibers[i-2].nodes[j].y));
      sheet->fibers[i].nodes[j].bend_force_z += bending_const  * (2*(sheet->fibers[i+1].nodes[j].z) 
								              - 5 * (sheet->fibers[i].nodes[j].z) 
								              + 4 * (sheet->fibers[i-1].nodes[j].z)
								              - (sheet->fibers[i-2].nodes[j].z));
   }
    // For the last fiber
    i = total_fibers_row-1;
    if (fiber2thread(i, total_fibers_row, total_threads) == tid ){
      sheet->fibers[i].nodes[j].bend_force_x += bending_const * (- (sheet->fibers[i].nodes[j].x) 
						            	      + 2 * (sheet->fibers[i-1].nodes[j].x)
							                  - (sheet->fibers[i-2].nodes[j].x));
      sheet->fibers[i].nodes[j].bend_force_y += bending_const * (- (sheet->fibers[i].nodes[j].y) 
							                  + 2 * (sheet->fibers[i-1].nodes[j].y)
							                  - (sheet->fibers[i-2].nodes[j].y));
      sheet->fibers[i].nodes[j].bend_force_z += bending_const * (- (sheet->fibers[i].nodes[j].z) 
							                  + 2 * (sheet->fibers[i-1].nodes[j].z)
							                  - (sheet->fibers[i-2].nodes[j].z));
   }
  }

  #ifdef DEBUG_PRINT
    printf("**** Tid%d: compute_bendingforce Exit*******\n", tid);
  #endif //DEBUG_PRINT
}
     
void compute_stretchingforce(LV lv) { 
#ifdef DEBUG_PRINT
  printf("****Inside compute_stretchingforce*******\n");
#endif //DEBUG_PRINT

  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  int total_threads = gv->total_threads;
  
  //pthread_barrier_t barr = gv->barr;

  int i, j; 
  
  double ds1 = gv->fiber_shape->sheets[0].width / (gv->fiber_shape->sheets[0].num_cols - 1); 
  double ds2 = gv->fiber_shape->sheets[0].height / (gv->fiber_shape->sheets[0].num_rows - 1);  // corresponds to ydist on the fibersheet.  

  double stretching_const = gv->Ks_l/(ds1 * ds1);	

  // TODO: Modify to multi sheets
  /* for convenience to save space. */
  int    total_fibers_row  = gv->fiber_shape->sheets[0].num_rows;
  int    total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  Fiber* fiberarray        = gv->fiber_shape->sheets[0].fibers;
  //double fibersheet_w      = gv->fiber_shape->sheets[0].width;
  //double fibersheet_h      = gv->fiber_shape->sheets[0].height;	
 
  double tangx_right, tangx_left, tangy_right, tangy_left, tangz_right, tangz_left, dist_right, dist_left;
  double tangx_top, tangx_bottom, tangy_top, tangy_bottom, tangz_top, tangz_bottom, dist_top, dist_bottom;

  //computing stretching force along fibers (row)
  for (i = 0; i < total_fibers_row; ++i){	/* for all fibers  */	
    if (fiber2thread(i, total_fibers_row, total_threads) == tid ){
    /* general cases for the points in the middle */
    for (j = 1; j < total_fibers_clmn-1; ++j){ /* for all points for a given fiber*/
      tangx_right = fiberarray[i].nodes[j+1].x - fiberarray[i].nodes[j].x;
      tangy_right = fiberarray[i].nodes[j+1].y - fiberarray[i].nodes[j].y;
      tangz_right = fiberarray[i].nodes[j+1].z - fiberarray[i].nodes[j].z;
      tangx_left = fiberarray[i].nodes[j-1].x -fiberarray[i].nodes[j].x;
      tangy_left = fiberarray[i].nodes[j-1].y -fiberarray[i].nodes[j].y;
      tangz_left = fiberarray[i].nodes[j-1].z -fiberarray[i].nodes[j].z;
      dist_right = sqrt(tangx_right * tangx_right + tangy_right * tangy_right + tangz_right * tangz_right);
      dist_left = sqrt(tangx_left * tangx_left + tangy_left * tangy_left + tangz_left * tangz_left);
      fiberarray[i].nodes[j].stretch_force_x = stretching_const * ((dist_right-ds1) * (tangx_right/dist_right)
                                               + (dist_left - ds1) * (tangx_left/dist_left));
      fiberarray[i].nodes[j].stretch_force_y = stretching_const * ((dist_right-ds1) * (tangy_right/dist_right)
                                               + (dist_left - ds1) * (tangy_left/dist_left));
      fiberarray[i].nodes[j].stretch_force_z = stretching_const * ((dist_right-ds1) * (tangz_right/dist_right)
                                               + (dist_left - ds1) * (tangz_left/dist_left));
    }

    /* for the first point */
    tangx_right = fiberarray[i].nodes[1].x - fiberarray[i].nodes[0].x;
    tangy_right = fiberarray[i].nodes[1].y - fiberarray[i].nodes[0].y;
    tangz_right = fiberarray[i].nodes[1].z - fiberarray[i].nodes[0].z;
    dist_right = sqrt(tangx_right * tangx_right + tangy_right * tangy_right + tangz_right * tangz_right);
    fiberarray[i].nodes[0].stretch_force_x = stretching_const * ((dist_right - ds1) * (tangx_right/dist_right));
    fiberarray[i].nodes[0].stretch_force_y = stretching_const * ((dist_right - ds1) * (tangy_right/dist_right));
    fiberarray[i].nodes[0].stretch_force_z = stretching_const * ((dist_right - ds1) * (tangz_right/dist_right));    

    /* for the last point */
    tangx_left = fiberarray[i].nodes[total_fibers_clmn-1].x - fiberarray[i].nodes[total_fibers_clmn-2].x;
    tangy_left = fiberarray[i].nodes[total_fibers_clmn-1].y - fiberarray[i].nodes[total_fibers_clmn-2].y;
    tangz_left = fiberarray[i].nodes[total_fibers_clmn-1].z - fiberarray[i].nodes[total_fibers_clmn-2].z;
    dist_left = sqrt(tangx_left * tangx_left + tangy_left * tangy_left + tangz_left * tangz_left);
    fiberarray[i].nodes[total_fibers_clmn-1].stretch_force_x = stretching_const * ((dist_left - ds1) * (tangx_left/dist_left));
    fiberarray[i].nodes[total_fibers_clmn-1].stretch_force_y = stretching_const * ((dist_left - ds1) * (tangy_left/dist_left));
    fiberarray[i].nodes[total_fibers_clmn-1].stretch_force_z = stretching_const * ((dist_left - ds1) * (tangz_left/dist_left));
    } //if fiber2thread ends
  } //end of total_fibers along row i.e fibers along z axis

  //Need a barrier
  //Error Chek#5
  pthread_barrier_wait(&(gv->barr));
  /* barrcheck = pthread_barrier_wait(&barr);
  if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD){
      fprintf(stderr,"Could not wait on barrier\n");
      exit(1);
  }*/

  // Computing Streching Force along the direction normal to fibers
  stretching_const = gv->Ks_l/(ds2 * ds2);	
  for (j = 0; j < total_fibers_clmn; ++j){  //for the j-th point	
    for (i = 1; i < total_fibers_row-1; ++i){ //for the i-th fiber, middle fibers
      if (fiber2thread(i, total_fibers_row, total_threads) == tid ){
        tangx_top = fiberarray[i+1].nodes[j].x - fiberarray[i].nodes[j].x;
        tangy_top = fiberarray[i+1].nodes[j].y - fiberarray[i].nodes[j].y;
        tangz_top = fiberarray[i+1].nodes[j].z - fiberarray[i].nodes[j].z;
        tangx_bottom = fiberarray[i-1].nodes[j].x - fiberarray[i].nodes[j].x;
        tangy_bottom = fiberarray[i-1].nodes[j].y - fiberarray[i].nodes[j].y;
        tangz_bottom = fiberarray[i-1].nodes[j].z - fiberarray[i].nodes[j].z;
        dist_top = sqrt(tangx_top * tangx_top + tangy_top * tangy_top + tangz_top * tangz_top);
        dist_bottom = sqrt(tangx_bottom * tangx_bottom + tangy_bottom * tangy_bottom + tangz_bottom * tangz_bottom);
        fiberarray[i].nodes[j].stretch_force_x += stretching_const * ((dist_top - ds2) * (tangx_top / dist_top)
                                                  + (dist_bottom - ds2) * (tangx_bottom / dist_bottom));
        fiberarray[i].nodes[j].stretch_force_y += stretching_const*((dist_top - ds2) * (tangy_top / dist_top)
                                                  + (dist_bottom - ds2) * (tangy_bottom/dist_bottom));
        fiberarray[i].nodes[j].stretch_force_z += stretching_const*((dist_top - ds2) * (tangz_top / dist_top)
                                                  + (dist_bottom - ds2) * (tangz_bottom / dist_bottom));
    } // if fiber2thread ends
   }
    
    /* For the first fiber on the sheet boundary */
    if (fiber2thread(0, total_fibers_row, total_threads) == tid ){ // i =0 :first fiber node 
      tangx_top = fiberarray[1].nodes[j].x - fiberarray[0].nodes[j].x;
      tangy_top = fiberarray[1].nodes[j].y - fiberarray[0].nodes[j].y;
      tangz_top = fiberarray[1].nodes[j].z - fiberarray[0].nodes[j].z;
      dist_top = sqrt(tangx_top * tangx_top + tangy_top * tangy_top + tangz_top * tangz_top);
      fiberarray[0].nodes[j].stretch_force_x += stretching_const * ((dist_top - ds2) * (tangx_top/dist_top));
      fiberarray[0].nodes[j].stretch_force_y += stretching_const * ((dist_top - ds2) * (tangy_top/dist_top));
      fiberarray[0].nodes[j].stretch_force_z += stretching_const * ((dist_top - ds2) * (tangz_top/dist_top));
    } // if fiber2thread ends

    /* For the last fiber on the sheet boundary */
    if (fiber2thread(total_fibers_row-1, total_fibers_row, total_threads) == tid ){ // i = total_fibers_row-1 :last fiber node 
      tangx_bottom = fiberarray[total_fibers_row-1].nodes[j].x - fiberarray[total_fibers_row-2].nodes[j].x;
      tangy_bottom = fiberarray[total_fibers_row-1].nodes[j].y - fiberarray[total_fibers_row-2].nodes[j].y;
      tangz_bottom = fiberarray[total_fibers_row-1].nodes[j].z - fiberarray[total_fibers_row-2].nodes[j].z;
      dist_bottom = sqrt(tangx_bottom * tangx_bottom + tangy_bottom * tangy_bottom + tangz_bottom * tangz_bottom);
      fiberarray[total_fibers_row-1].nodes[j].stretch_force_x += stretching_const * ((dist_bottom - ds2) * (tangx_bottom / dist_bottom));
      fiberarray[total_fibers_row-1].nodes[j].stretch_force_y += stretching_const * ((dist_bottom - ds2) * (tangy_bottom / dist_bottom));
      fiberarray[total_fibers_row-1].nodes[j].stretch_force_z += stretching_const * ((dist_bottom - ds2) * (tangz_bottom / dist_bottom));
    }
  }

#ifdef DEBUG_PRINT
  printf("**** compute_stretchingforce Exit*******\n");
#endif //DEBUG_PRINT
}


void compute_elasticforce(LV lv) { //emabarassingly parallel
#ifdef DEBUG_PRINT
   printf("****Inside compute_elasticforce*******\n");
#endif //DEBUG_PRINT

  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  int total_threads = gv->total_threads;

  int i, j;
  // TODO: Modify to multi sheets
  int total_fibers_row  = gv->fiber_shape->sheets[0].num_rows;
  int total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  Fiber *fiberarray     = gv->fiber_shape->sheets[0].fibers;

  // Adding up Streching Force and Bending Force
  for (i = 0; i < total_fibers_row; ++i){	
    if (fiber2thread(i, total_fibers_row, total_threads) == tid ){ 
      for (j = 0; j < total_fibers_clmn; ++j){
        fiberarray[i].nodes[j].elastic_force_x = fiberarray[i].nodes[j].stretch_force_x 
	                                             + fiberarray[i].nodes[j].bend_force_x;
        fiberarray[i].nodes[j].elastic_force_y = fiberarray[i].nodes[j].stretch_force_y 
	                                             + fiberarray[i].nodes[j].bend_force_y;
        fiberarray[i].nodes[j].elastic_force_z = fiberarray[i].nodes[j].stretch_force_z 
	                                             + fiberarray[i].nodes[j].bend_force_z;
      } // end of total_fluidgridpoints i.e fibers along y axis	
    } // if fiber2thread ends
   } // end of total_fibers i.e fibers along x axis	
   
#ifdef DEBUG_PRINT
     printf("**** compute_elasticforce Exit*******\n");
#endif //DEBUG_PRINT
}


/* Step1: Influential domain for force-spreading and velocity interpolation. */
/* Step2: Do actual force spreading. */
void get_influentialdomain_fluid_grid_and_SpreadForce(LV lv){//Fiber influences fluid
  #ifdef DEBUG_PRINT
    printf("****Inside get_influentaildomain_fluid_grid******\n");
  #endif //DEBUG_PRINT

  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  //pthread_mutex_t lock_Fluid = gv->lock_Fluid;

  int     i=0, j=0, inneri=0, innerj=0, innerk=0; // inner i corresponding to fluid index
  int     istart,istop, jstart, jstop, kstart, kstop;
  double  dx= 1.0, dy = 1.0, dz= 1.0;
  double  rx = 1.0, ry = 1.0, rz =1.0;

  int     total_fibers_row, total_fibers_clmn;
  Fiber   *fiberarray;

  double  tmp_dist;

  Fluidnode     *nodes;
  int     BI, BJ, BK;
  int     li,lj,lk;
  int     total_sub_grids, cube_size, dim_x, dim_y, dim_z;
  int     num_cubes_x, num_cubes_y, num_cubes_z, cube_idx, node_idx;
  int     P, Q, R , total_threads;
  //printf("*************ENTRY OF getDomain and spreadforce\n");
  // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, 
                            //lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
  
  /* TODO: Move the following variables to GV to save time */
  double  c_x= PI / (2.0 * dx);
  double  c_y= PI / (2.0 * dy);
  double  c_z= PI / (2.0 * dz);
  double  cf = 1.0 / (64.0 * dx * dy * dz); //cf[jf]=1.0/(64.0*dx*dy*dz);

  int owner_id;//only owner will lock 

  // TODO: Modify to multi sheets
  total_fibers_row  = gv->fiber_shape->sheets[0].num_rows;
  total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  fiberarray        = gv->fiber_shape->sheets[0].fibers;
   /*Pthread chages*/

  dim_x           = gv->fluid_grid->x_dim;
  dim_y           = gv->fluid_grid->y_dim;
  dim_z           = gv->fluid_grid->z_dim;
  cube_size       = gv->cube_size;
  total_sub_grids = (dim_x*dim_y*dim_z) /pow(cube_size,3);
  num_cubes_x     = gv->num_cubes_x;
  num_cubes_y     = gv->num_cubes_y;
  num_cubes_z     = gv->num_cubes_z;

  P = gv->P;
  Q = gv->Q;
  R = gv->R;
  total_threads = P*Q*R;//gv->total_threads
  
  //pthread_mutex_t lock_Fluid[total_threads] = gv->lock_Fluid[total_threads];

  // Annuling Forces on Fluid grid ::SHOULD be HERE!! 
  for (BI = 0 ; BI < num_cubes_x; BI++)
    for (BJ = 0; BJ < num_cubes_y; BJ++)
      for (BK = 0; BK < num_cubes_z; BK++){
        if (cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q, R) ==tid){
          cube_idx = BI*num_cubes_y*num_cubes_z +BJ*num_cubes_z +BK;
          for (li = 0; li < cube_size; li++)
            for (lj = 0; lj < cube_size; lj++)
              for (lk = 0; lk < cube_size; lk++){
                node_idx = li* cube_size * cube_size + lj * cube_size + lk;
                gv->fluid_grid->sub_fluid_grid[cube_idx].nodes[node_idx].elastic_force_x  = 0.0e0;
                gv->fluid_grid->sub_fluid_grid[cube_idx].nodes[node_idx].elastic_force_y  = 0.0e0;
                gv->fluid_grid->sub_fluid_grid[cube_idx].nodes[node_idx].elastic_force_z  = 0.0e0;
              }
        } // if cube2thread ends 
      }  

      /*for(subgrid_idx=0; subgrid_idx< total_sub_grids; subgrid_idx++)
       for(li=0; li< cube_size; li++)
        for(lj= 0; lj< cube_size; lj++)
         for(lk=0; lk< cube_size; lk++ ){
            node_idx = li* cube_size * cube_size + lj * cube_size + lk;
            gv->fluid_grid->sub_fluid_grid[subgrid_idx].nodes[node_idx].elastic_force_x  = 0.0e0;
            gv->fluid_grid->sub_fluid_grid[subgrid_idx].nodes[node_idx].elastic_force_y  = 0.0e0;
            gv->fluid_grid->sub_fluid_grid[subgrid_idx].nodes[node_idx].elastic_force_z  = 0.0e0;
         }*/
    /*Pthread chages*/
    //printf("*************Afer Anuuling\n");
      // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);

  for (i = 0; i < total_fibers_row; ++i){
    if (fiber2thread(i, total_fibers_row, total_threads)==tid){
      for (j = 0; j < total_fibers_clmn; ++j){
        /* To find influential domain */
        istart = floor(fiberarray[i].nodes[j].x/dx - 2) + 1; // x dimension
        istop  = istart + 3; 
        jstart = floor(fiberarray[i].nodes[j].y/dy - 2) + 1; // y dimension
        jstop  = jstart + 3;
        kstart = floor(fiberarray[i].nodes[j].z/dz - 2) + 1; // z dimension
        kstop  = kstart + 3;

        /* Spreading force computation */
     
        for (inneri = istart; inneri <= istop; inneri++) // x direction
          for (innerj = jstart; innerj <= jstop; innerj++) // y direction
            for (innerk = kstart; innerk <= kstop; innerk++){ // z direction
              /* Used for calculating eqn 21 ....PN */
              // distance between the fiber node and all the fluid nodes within the infulential doman.
              rx = dx * inneri - fiberarray[i].nodes[j].x;
              ry = dy * innerj - fiberarray[i].nodes[j].y;
              rz = dz * innerk - fiberarray[i].nodes[j].z;
              /*annul ffx,ffy,ffz before spread for each time step Done Above */
              BI = inneri / cube_size; //BI should be inneri/cube_size but BI_end seems correct
              BJ = innerj / cube_size; 
              BK = innerk / cube_size; 
              li = inneri % cube_size;
              lj = innerj % cube_size;
              lk = innerk % cube_size;
                
              cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
              //node_idx = (inneri%cube_size) *cube_size*cube_size +(innerj%cube_size) *cube_size +innerk%cube_size;
              node_idx = (li) * cube_size * cube_size +(lj) * cube_size + lk;
              nodes    = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;     	    
              /* Eqn 19 for Step2 ....PN*/
               
              tmp_dist = cf * (1.0e0 + cos(c_x * rx)) * (1.0e0 + cos(c_y * ry)) * (1.0e0 + cos(c_z * rz));
              // TODO: remove check
              if (inneri < 0 || inneri >= gv->fluid_grid->x_dim) {
                fprintf(stderr, "inneri out of bound: %d, x bound=%d!\n", inneri, gv->fluid_grid->x_dim); exit(1);
              }
              if (innerj < 0 || innerj >= gv->fluid_grid->y_dim) {
                fprintf(stderr, "innerj out of bound: %d, y bound=%d!\n", innerj, gv->fluid_grid->y_dim); exit(1);
              }
              if (innerk < 0 || innerk >= gv->fluid_grid->z_dim) {
                fprintf(stderr, "innerk out of bound: %d, z bound=%d!\n", innerk, gv->fluid_grid->z_dim); exit(1);
              }
              /* if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
                && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size ){
                    // printf("*************Before calculating\n");
                    // print_fluid_sub_grid(gv, lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, 
                                            lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
                 }*/
      
              // Need more locks--12/22
              // if node_idx belongs to thread k then lock thread k's lock gv->locksFLuid[tid_k];
              owner_id = cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R);
          
              pthread_mutex_lock(&(gv->lock_Fluid[owner_id]));
              nodes[node_idx].elastic_force_x += fiberarray[i].nodes[j].elastic_force_x * tmp_dist;
              nodes[node_idx].elastic_force_y += fiberarray[i].nodes[j].elastic_force_y * tmp_dist;
              nodes[node_idx].elastic_force_z += fiberarray[i].nodes[j].elastic_force_z * tmp_dist;
              pthread_mutex_unlock(&(gv->lock_Fluid[owner_id]));
              //unlock

              /*if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
              pthread_mutex_lock(&(gv->lock_Fluid[tid]));
              nodes[node_idx].elastic_force_x += fiberarray[i].nodes[j].elastic_force_x * tmp_dist;
              nodes[node_idx].elastic_force_y += fiberarray[i].nodes[j].elastic_force_y * tmp_dist;
              nodes[node_idx].elastic_force_z += fiberarray[i].nodes[j].elastic_force_z * tmp_dist;
              pthread_mutex_unlock(&(gv->lock_Fluid[tid]));
              //unlock
              }
              else{
              nodes[node_idx].elastic_force_x += fiberarray[i].nodes[j].elastic_force_x * tmp_dist;
              nodes[node_idx].elastic_force_y += fiberarray[i].nodes[j].elastic_force_y * tmp_dist;
              nodes[node_idx].elastic_force_z += fiberarray[i].nodes[j].elastic_force_z * tmp_dist; 
              }*/
            }
      } // end of total_fluidgridpoints i.e fibers along y axis
    } // if fiber2thread ends
  } // end of total_fibers i.e fibers along x axis
  
  //printf("*************Exiting gtinf and spred\n");
  /* print_fluid_sub_grid(gv, lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, 
                          lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
  */  
  // printf("****   get_influentaildomain_fluid_grid Exit******\n");
}	
 
 
void init_eqlbrmdistrfuncDF0(GV gv){/*stored in dfreq*/
  Fluidgrid *fluid_grid;
  Fluidnode *nodes;
 
  int       ksi;
  int       BI, BJ, BK; // to identify the Sub grids
  int       cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int       starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z; // To identify buffer zone
  int       li, lj, lk, node_idx; // local access point inside cube

  fluid_grid      = gv->fluid_grid;
     
  cube_size   = gv->cube_size; 
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;
  

/*PTHREAD_Change*/

  //For each cube BI, BJ ,BK, li, lj, lk 0 to cube size-1
  for (BI = 0; BI < num_cubes_x; ++BI){
    for (BJ = 0; BJ < num_cubes_y; ++BJ){
      for (BK = 0; BK < num_cubes_z; ++BK){
        cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        nodes = fluid_grid->sub_fluid_grid[cube_idx].nodes;
        starting_x = starting_y = starting_z = 0;
        stopping_x = stopping_y = stopping_z = cube_size-1;
        for (li = starting_x; li <= stopping_x; ++li){
          for (lj = starting_y; lj <= stopping_y; ++lj){
            for (lk = starting_z; lk <= stopping_z; ++lk){
              node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
              for (ksi = 0; ksi < 19; ++ksi){
                if (ksi == 0){
                  nodes[node_idx].dfeq[ksi] = 
                        1.0/3.0 * nodes[node_idx].rho
                        * (1.0 - 1.5 *              
                            (nodes[node_idx].vel_x * nodes[node_idx].vel_x 
                            + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                            + nodes[node_idx].vel_z * nodes[node_idx].vel_z));
                } else if (ksi >0 && ksi<7) {
                  nodes[node_idx].dfeq[ksi]=
                        1.0/18.0 * nodes[node_idx].rho 
                        *(1.0 + 3.0 *
                            (gv->c[ksi][0]*nodes[node_idx].vel_x
                                + gv->c[ksi][1]*nodes[node_idx].vel_y
                                + gv->c[ksi][2]*nodes[node_idx].vel_z)
                            + 4.5 * (pow(gv->c[ksi][0]*nodes[node_idx].vel_x
                                + gv->c[ksi][1]*nodes[node_idx].vel_y
                                + gv->c[ksi][2]*nodes[node_idx].vel_z, 2))
                            - 1.5 * (nodes[node_idx].vel_x * nodes[node_idx].vel_x
                                + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                                + nodes[node_idx].vel_z * nodes[node_idx].vel_z));
               } else {
                 nodes[node_idx].dfeq[ksi]=
                        1.0/36.0 * nodes[node_idx].rho 
                        *(1.0 + 3.0*
                            (gv->c[ksi][0]*nodes[node_idx].vel_x
                            +gv->c[ksi][1]*nodes[node_idx].vel_y
                            +gv->c[ksi][2]*nodes[node_idx].vel_z)
                         + 4.5 *(pow(gv->c[ksi][0]*nodes[node_idx].vel_x
                            + gv->c[ksi][1]*nodes[node_idx].vel_y
                            + gv->c[ksi][2]*nodes[node_idx].vel_z,2))
                        - 1.5 * (nodes[node_idx].vel_x *nodes[node_idx].vel_x
                            +nodes[node_idx].vel_y *nodes[node_idx].vel_y
                            +nodes[node_idx].vel_z *nodes[node_idx].vel_z));
               }
             } // for loop ksi ends
           }// for loop local z
         }// for loop local y
       } // for loop local x
     }//for BK
    }//for BJ
  }//FOR BI

/*PTHREAD_Change*/

  //printf("****************init_eqlbrmdistrfuncDF0 EXIT ***************\n");
}


/* Called only once before simulation starts */
void init_df_inout(GV gv){
  //printf("****************Inside copy_df_to_inout*************\n");
  Fluidgrid *fluid_grid;
  Fluidnode *nodes;
 
  int ksi, dim_z;
  int BI, BJ, BK; // to identify the Sub grids
  int cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int starting_y, starting_z, stopping_y, stopping_z; // To identify buffer zone
  int li, lj, lk, node_idx; // local access point inside cube

  fluid_grid = gv->fluid_grid;
  dim_z     = fluid_grid->z_dim;
     
  cube_size   = gv->cube_size; 
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;
  //printf("*************ENTRY OF init_df_inout\n");
  // print_fluid_sub_grid(gv, lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, 
                          //lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);

  /* inlet, copy once */ //TODO: dfeq could be removed
  li = gv->ib;
  BI=0;
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
    for (BK = 0; BK < num_cubes_z; ++BK){
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_y = starting_z = 0;
      stopping_y = stopping_z = cube_size-1;
      for (lj = starting_y; lj <= stopping_y; ++lj)
        for (lk=starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
          for (ksi = 0; ksi <= 18; ksi++){ 
            fluid_grid->inlet->nodes[(BJ * cube_size + lj) * dim_z + BK * cube_size + lk].df_inout[0][ksi] =
                nodes[node_idx].dfeq[ksi];
          }
        }    
    }

  /* outlet, copy once */ //TODO: dfeq could be removed 
  li = cube_size-1;
  BI = num_cubes_x-1;//ie
  for (BJ = 0; BJ < num_cubes_y; ++BJ)
    for (BK = 0; BK < num_cubes_z; ++BK){
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_y = starting_z = 0;
      stopping_y = stopping_z = cube_size-1;
      for (lj = starting_y; lj <= stopping_y; ++lj)
        for (lk = starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk;
          for (ksi = 0; ksi <= 18; ksi++){
            fluid_grid->outlet->nodes[(BJ * cube_size + lj) * dim_z + BK * cube_size + lk].df_inout[1][ksi] =
                nodes[node_idx].dfeq[ksi];
          }
        }
     }
  //printf("****************copy_df_to_inout: EXIT*************\n");  
}


void compute_eqlbrmdistrfuncDF1(LV lv){
  #ifdef DEBUG_PRINT
printf("****************Inside compute_eqlbrmdistrfuncDF1*************\n");
#endif //DEBUG_PRINT
  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;

  Fluidgrid *fluid_grid;
  Fluidnode *nodes;
 
  int ksi;
  int BI, BJ, BK; // to identify the Sub grids
  int cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z; // To identify buffer zone
  int li, lj, lk, node_idx; // local access point inside cube

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;

  fluid_grid = gv->fluid_grid;
     
  cube_size   = gv->cube_size; 
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  //printf("*************ENTRY OF DF1\n");
  //print_fluid_sub_grid(gv, lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, 
                        // lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
  /*Pthread Changes*/ /*OPTIMISE using Loop unrolling*/
  for (BI = 0; BI < num_cubes_x; ++BI)
    for (BJ = 0; BJ < num_cubes_y; ++BJ)
      for (BK = 0; BK < num_cubes_z; ++BK){
        if (cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
          cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
          starting_x = starting_y = starting_z = 0;
          stopping_x = stopping_y = stopping_z = cube_size-1;
          
          if (BI == 0) starting_x = 2;//ib
          
          if(BI == num_cubes_x-1) stopping_x = cube_size-3;//ie
           
          if(BJ == 0) starting_y = 2;//jb
          
          if(BJ == num_cubes_y-1) stopping_y = cube_size-3;//je
          
          if(BK == 0) starting_z = 2;//kb
          
          if(BK == num_cubes_z-1) stopping_z = cube_size-3;//ke
          
          for (li = starting_x; li <= stopping_x; ++li)
            for (lj = starting_y; lj <= stopping_y; ++lj)
              for (lk = starting_z; lk <= stopping_z; ++lk){
                node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
                for (ksi = 0; ksi <= 18; ksi++){
                  if (ksi == 0){
                    nodes[node_idx].dfeq[ksi] = 
                        1.0 / 3.0 * nodes[node_idx].rho
                        * (1.0 - 1.5 *              
                            (nodes[node_idx].vel_x * nodes[node_idx].vel_x 
                            + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                            + nodes[node_idx].vel_z * nodes[node_idx].vel_z));

                    nodes[node_idx].df1[ksi] = 
                        nodes[node_idx].df1[ksi] * (1.0 - 1.0/ gv->tau)
                        + 1.0 / gv->tau * nodes[node_idx].dfeq[ksi]
                        + gv->dt * (1.0 - 1.0 / (2.0 * gv->tau)) / 3.0
                          * (
                            ( (gv->c[ksi][0] - nodes[node_idx].vel_x) * nodes[node_idx].elastic_force_x
                            + (gv->c[ksi][1] - nodes[node_idx].vel_y) * nodes[node_idx].elastic_force_y
                            + (gv->c[ksi][2] - nodes[node_idx].vel_z) * (nodes[node_idx].elastic_force_z         
                            + nodes[node_idx].rho * gv->g_l) ) / (gv->cs_l*gv->cs_l)
                        + ( gv->c[ksi][0] * nodes[node_idx].vel_x 
                          + gv->c[ksi][1] * nodes[node_idx].vel_y
                          + gv->c[ksi][2] * nodes[node_idx].vel_z) / pow(gv->cs_l,4) 
                            * ( gv->c[ksi][0] * nodes[node_idx].elastic_force_x 
                              + gv->c[ksi][1] * nodes[node_idx].elastic_force_y
                              + gv->c[ksi][2] * ( nodes[node_idx].elastic_force_z
                                                + nodes[node_idx].rho * gv->g_l)));
                  } //ksi==0
                  else if (ksi >= 1 && ksi <= 6){
                    nodes[node_idx].dfeq[ksi] =
                        1.0 / 18.0 * nodes[node_idx].rho 
                        *(1.0 + 3.0 *
                            ( gv->c[ksi][0] * nodes[node_idx].vel_x
                            + gv->c[ksi][1] * nodes[node_idx].vel_y
                            + gv->c[ksi][2] * nodes[node_idx].vel_z)
                         + 4.5 * ( pow(gv->c[ksi][0] * nodes[node_idx].vel_x
                                 + gv->c[ksi][1] * nodes[node_idx].vel_y
                                 + gv->c[ksi][2] * nodes[node_idx].vel_z, 2))
                         - 1.5 * ( nodes[node_idx].vel_x * nodes[node_idx].vel_x
                                 + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                                 + nodes[node_idx].vel_z * nodes[node_idx].vel_z));

                    nodes[node_idx].df1[ksi] =
                        nodes[node_idx].df1[ksi]*(1.0 - 1.0 / gv->tau) 
                        + 1.0 / gv->tau * nodes[node_idx].dfeq[ksi]
                        + gv->dt * (1.0 - 1.0 / (2.0 * gv->tau)) / 18.0
                        * (
                          ( (gv->c[ksi][0] - nodes[node_idx].vel_x) * nodes[node_idx].elastic_force_x
                          + (gv->c[ksi][1] - nodes[node_idx].vel_y) * nodes[node_idx].elastic_force_y
                          + (gv->c[ksi][2] - nodes[node_idx].vel_z) * 
                                    (nodes[node_idx].elastic_force_z + nodes[node_idx].rho*gv->g_l) ) / (gv->cs_l*gv->cs_l)
                          + ( gv->c[ksi][0] * nodes[node_idx].vel_x 
                            + gv->c[ksi][1] * nodes[node_idx].vel_y
                            + gv->c[ksi][2] * nodes[node_idx].vel_z) / pow(gv->cs_l,4)
                              * ( gv->c[ksi][0] * nodes[node_idx].elastic_force_x 
                                + gv->c[ksi][1] * nodes[node_idx].elastic_force_y
                                + gv->c[ksi][2] * ( nodes[node_idx].elastic_force_z
                                                  + nodes[node_idx].rho * gv->g_l) ) );
                  } //ksi between 1 and 7 ends
                  else {
                     nodes[node_idx].dfeq[ksi] =
                        1.0 / 36.0 * nodes[node_idx].rho 
                        * (1.0 + 3.0 *
                             ( gv->c[ksi][0] * nodes[node_idx].vel_x
                             + gv->c[ksi][1] * nodes[node_idx].vel_y
                             + gv->c[ksi][2] * nodes[node_idx].vel_z)
                          + 4.5 * ( pow( gv->c[ksi][0] * nodes[node_idx].vel_x
                                       + gv->c[ksi][1] * nodes[node_idx].vel_y
                                       + gv->c[ksi][2] * nodes[node_idx].vel_z, 2))
                          - 1.5 * ( nodes[node_idx].vel_x * nodes[node_idx].vel_x
                                  + nodes[node_idx].vel_y * nodes[node_idx].vel_y
                                  + nodes[node_idx].vel_z * nodes[node_idx].vel_z));

                     nodes[node_idx].df1[ksi] =
                        nodes[node_idx].df1[ksi] * (1.0 - 1.0 / gv->tau)
                        + 1.0 / gv->tau * nodes[node_idx].dfeq[ksi]
                        + gv->dt * (1.0 - 1.0 / (2.0 * gv->tau)) / 36.0
                        * (
                          ( (gv->c[ksi][0] - nodes[node_idx].vel_x) * nodes[node_idx].elastic_force_x
                          + (gv->c[ksi][1] - nodes[node_idx].vel_y) * nodes[node_idx].elastic_force_y
                          + (gv->c[ksi][2] - nodes[node_idx].vel_z) 
                                             * ( nodes[node_idx].elastic_force_z
                                               + nodes[node_idx].rho * gv->g_l) ) / (gv->cs_l*gv->cs_l)
                          + ( gv->c[ksi][0] * nodes[node_idx].vel_x 
                            + gv->c[ksi][1] * nodes[node_idx].vel_y
                            + gv->c[ksi][2] * nodes[node_idx].vel_z) / pow(gv->cs_l,4)
                              * ( gv->c[ksi][0] * nodes[node_idx].elastic_force_x 
                                + gv->c[ksi][1] * nodes[node_idx].elastic_force_y
                                + gv->c[ksi][2] * ( nodes[node_idx].elastic_force_z
                                                  + nodes[node_idx].rho*gv->g_l) ) );
                  }//ksi between 7 and 18 ends
                }//ksi loop
             }//lk loop
        }//if cube2thread ends
      }//For BK

  // printf("****************compute_eqlbrmdistrfuncDF1 EXIT ***************\n");  
}

void  stream_distrfunc(LV lv){//df2
  /*18 spaces for 18 fluid nodes(neighbours). Whn 18 neighbours update current node  they write in their own space. so no lock required */
  //printf("****************Inside stream_distfunc()*******at time %d\n",gv->time);
  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  
  Fluidgrid *fluid_grid;
  Fluidnode *nodes_df1, *nodes_df2;
  
  int cube_df1_idx, cube_df2_idx, node_df1_idx, node_df2_idx, cube_size, num_cubes_x, num_cubes_y, num_cubes_z;
  int starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;
  
  int BI, BJ, BK, li, lj, lk;
  
  fluid_grid   = gv->fluid_grid;
  cube_size   = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;
  
  starting_x = starting_y = starting_z = 0;
  stopping_x = stopping_y = stopping_z = cube_size-1;
  // printf("*************ENTRY OF Streaming\n");
  // print_fluid_sub_grid(gv, lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, 
                         //lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
  
  for (BI = 0; BI < num_cubes_x; ++BI)
    for (BJ = 0; BJ < num_cubes_y; ++BJ)
      for (BK = 0; BK < num_cubes_z; ++BK){
        if (cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q, R) ==tid){
          cube_df1_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df1    = gv->fluid_grid->sub_fluid_grid[cube_df1_idx].nodes;
          starting_x = starting_y = starting_z = 0;
          stopping_x = stopping_y = stopping_z = cube_size-1;

          if(BI == 0)   starting_x = 2;
             
          if(BI == num_cubes_x-1) stopping_x = cube_size-3;
               
          if(BJ == 0) starting_y = 2;
             
          if(BJ == num_cubes_y-1) stopping_y = cube_size-3;
              
          if(BK == 0) starting_z = 2;
              
          if(BK == num_cubes_z-1) stopping_z = cube_size-3;
          
     for(li=starting_x; li <= stopping_x; ++li)
      for(lj=starting_y; lj <= stopping_y; ++lj)
       for(lk=starting_z; lk <= stopping_z; ++lk){
       // surfaces[i].nodes[j*dim_z +k].df2[0]        = surfaces[i].nodes[j*dim_z +k].df1[0];
        node_df1_idx = li*cube_size*cube_size +lj*cube_size +lk;
//ksi==0 starts
        cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
        nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
        node_df2_idx = li*cube_size*cube_size +lj*cube_size +lk;
        nodes_df2[node_df2_idx].df2[0]  = nodes_df1[node_df1_idx].df1[0];
//ksi==0 ends
//ksi==1 starts
//surfaces[i+1].nodes[j*dim_z +k].df2[1]      = surfaces[i].nodes[j*dim_z +k].df1[1]; 
        if((li+1) <= cube_size-1){//i.e the the computation domain is still in the current cube
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = (li+1)*cube_size*cube_size +lj*cube_size +lk;
         nodes_df2[node_df2_idx].df2[1]  = nodes_df1[node_df1_idx].df1[1];
        }else{// computation domain is in the next cube only BI changes
         cube_df2_idx  = (BI+1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = (0)*cube_size*cube_size +lj*cube_size +lk; // local i of next cube from 0
         nodes_df2[node_df2_idx].df2[1]  = nodes_df1[node_df1_idx].df1[1];
        }
//ksi==1 ends 
//ksi==2 starts
//surfaces[i-1].nodes[j*dim_z +k].df2[2]      = surfaces[i].nodes[j*dim_z +k].df1[2];
        if((li-1) >=0){//i.e the the computation domain is still in the current cube
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = (li-1)*cube_size*cube_size +lj*cube_size +lk;
         nodes_df2[node_df2_idx].df2[2]  = nodes_df1[node_df1_idx].df1[2];
        }else{// computation domain is in the prev cube only BI changes
         cube_df2_idx  = (BI-1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = (cube_size-1)*cube_size*cube_size +lj*cube_size +lk;//local i of prev cube will be ending point.
         nodes_df2[node_df2_idx].df2[2]  = nodes_df1[node_df1_idx].df1[2];
        }
//ksi==2 ends      
//ksi==3 starts
//surfaces[i].nodes[j*dim_z +k+1].df2[3]      = surfaces[i].nodes[j*dim_z +k].df1[3];
        if((lk+1) <=cube_size-1){//i.e the the computation domain is still in the current cube
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = li*cube_size*cube_size +lj*cube_size +lk+1;
         nodes_df2[node_df2_idx].df2[3]  = nodes_df1[node_df1_idx].df1[3];
        }else{//computation domain is in the next cube only BK changes
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK+1;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = li*cube_size*cube_size +lj*cube_size +0;//local k of next cube will start from 0
         nodes_df2[node_df2_idx].df2[3]  = nodes_df1[node_df1_idx].df1[3];
        }	
//ksi==3 ends
//ksi==4 starts
//surfaces[i].nodes[j*dim_z +k-1].df2[4]      = surfaces[i].nodes[j*dim_z +k].df1[4];
        if((lk-1) >= 0){//i.e the the computation domain is still in the current cube
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = li*cube_size*cube_size +lj*cube_size +lk-1;
         nodes_df2[node_df2_idx].df2[4]  = nodes_df1[node_df1_idx].df1[4];
        }else{// computation domain is in the prev cube only BK changes
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK-1;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = li*cube_size*cube_size +lj*cube_size +cube_size-1;//local k of prev cube will be the end point
         nodes_df2[node_df2_idx].df2[4]  = nodes_df1[node_df1_idx].df1[4];
        }
//ksi==4 ends	
//ksi==5 starts 
//surfaces[i].nodes[(j-1)*dim_z +k].df2[5]    = surfaces[i].nodes[j*dim_z +k].df1[5];
       if((lj-1) >= 0){//i.e the the computation domain is still in the current cube
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = li*cube_size*cube_size +(lj-1)*cube_size +lk;
         nodes_df2[node_df2_idx].df2[5]  = nodes_df1[node_df1_idx].df1[5];
        }else{// computation domain is in the prev cube only BJ changes
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + (BJ-1) * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = li*cube_size*cube_size +(cube_size-1)*cube_size +lk;//local j of prev cube will be the end point
         nodes_df2[node_df2_idx].df2[5]  = nodes_df1[node_df1_idx].df1[5];
        }
//ksi==5 ends	
//ksi==6 starts //surfaces[i].nodes[(j+1)*dim_z +k].df2[6]    = surfaces[i].nodes[j*dim_z +k].df1[6];
        if((lj+1) <=cube_size-1){//i.e the the computation domain is still in the current cube
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = li*cube_size*cube_size +(lj+1)*cube_size +lk;
         nodes_df2[node_df2_idx].df2[6]  = nodes_df1[node_df1_idx].df1[6];
        }else{//computation domain is in the next cube only BJ changes
         cube_df2_idx  = BI * num_cubes_y * num_cubes_z + (BJ+1) * num_cubes_z + BK;
         nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
         node_df2_idx = li*cube_size*cube_size +0*cube_size +lk;//local j of next cube will start from 0
         nodes_df2[node_df2_idx].df2[6]  = nodes_df1[node_df1_idx].df1[6];
        }
//ksi==6 ends
//ksi==7 starts//surfaces[i+1].nodes[(j)*dim_z +k+1].df2[7]  = surfaces[i].nodes[j*dim_z +k].df1[7];
        if((li+1) <=cube_size-1 && (lk+1) <=cube_size-1){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li+1)*cube_size*cube_size +lj*cube_size +lk+1;
          nodes_df2[node_df2_idx].df2[7]  = nodes_df1[node_df1_idx].df1[7];
        }else if((li+1) >cube_size-1 && (lk+1) >cube_size-1){//comp domain changes both BI and BK
          cube_df2_idx  = (BI+1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK+1);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (0)*cube_size*cube_size +lj*cube_size +0;//both local i and local k will be starting point of next cube
          nodes_df2[node_df2_idx].df2[7]  = nodes_df1[node_df1_idx].df1[7];
        }else if((li+1) >cube_size-1 && (lk+1) <=cube_size-1){//only BI changes
          cube_df2_idx  = (BI+1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (0)*cube_size*cube_size +lj*cube_size +lk+1;//only local i will be starting point of next cube
          nodes_df2[node_df2_idx].df2[7]  = nodes_df1[node_df1_idx].df1[7];
        }else{//only BK changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK+1;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li+1)*cube_size*cube_size +lj*cube_size +0;//only local k will be starting point of next cube
          nodes_df2[node_df2_idx].df2[7]  = nodes_df1[node_df1_idx].df1[7];
        }
//ksi==7 ends
//ksi==8 starts// surfaces[i-1].nodes[(j)*dim_z +k-1].df2[8]  = surfaces[i].nodes[j*dim_z +k].df1[8];
        if((li-1) >=0 && (lk-1) >=0){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li-1)*cube_size*cube_size +lj*cube_size +(lk-1);
          nodes_df2[node_df2_idx].df2[8]  = nodes_df1[node_df1_idx].df1[8];
        }else if((li-1) <0 && (lk-1) <0){//comp domain changes both BI and BK
          cube_df2_idx  = (BI-1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK-1);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (cube_size-1)*cube_size*cube_size +lj*cube_size + (cube_size-1);//both local i and local k will be ending point of prev cube
          nodes_df2[node_df2_idx].df2[8]  = nodes_df1[node_df1_idx].df1[8];
        }else if((li-1)<0 && (lk-1) >=0){//only BI changes
          cube_df2_idx  = (BI-1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (cube_size-1)*cube_size*cube_size +lj*cube_size +lk-1;//only local i will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[8]  = nodes_df1[node_df1_idx].df1[8];
        }else{//only BK changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK-1);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li-1)*cube_size*cube_size +lj*cube_size +cube_size-1;//only local k will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[8]  = nodes_df1[node_df1_idx].df1[8];
        }
//ksi==8 ends
//ksi==9 starts//surfaces[i+1].nodes[(j)*dim_z +k-1].df2[9]  = surfaces[i].nodes[j*dim_z +k].df1[9]; 
       if((li+1) <=cube_size-1 && (lk-1) >=0){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li+1)*cube_size*cube_size +lj*cube_size +(lk-1);
          nodes_df2[node_df2_idx].df2[9]  = nodes_df1[node_df1_idx].df1[9];
        }else if((li+1) >cube_size-1 && (lk-1)<0){//comp domain changes both BI and BK
          cube_df2_idx  = (BI+1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK-1;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (0)*cube_size*cube_size +lj*cube_size +cube_size-1;// local i starting point of next cube and local k will be endng point of prev cube
          nodes_df2[node_df2_idx].df2[9]  = nodes_df1[node_df1_idx].df1[9];
        }else if((li+1) >cube_size-1 && (lk-1) >=0){//only BI changes
          cube_df2_idx  = (BI+1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (0)*cube_size*cube_size +lj*cube_size +(lk-1);//only local i will be starting point of next cube
          nodes_df2[node_df2_idx].df2[9]  = nodes_df1[node_df1_idx].df1[9];
        }else{//only BK changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK-1;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li+1)*cube_size*cube_size +lj*cube_size +cube_size-1;//only local k will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[9]  = nodes_df1[node_df1_idx].df1[9];
        }
//ksi==9 ends
//ksi==10 starts//surfaces[i-1].nodes[(j)*dim_z +k+1].df2[10] = surfaces[i].nodes[j*dim_z +k].df1[10];
        if((li-1) >=0 && (lk+1) <=cube_size-1){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li-1)*cube_size*cube_size +lj*cube_size +lk+1;
          nodes_df2[node_df2_idx].df2[10]  = nodes_df1[node_df1_idx].df1[10];
        }else if((li-1) <0 && (lk+1) >cube_size-1){//comp domain changes both BI and BK
          cube_df2_idx  = (BI-1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK+1;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (cube_size-1)*cube_size*cube_size +lj*cube_size + 0;//local i will be ending point of prev cube and local k strtng of next 
          nodes_df2[node_df2_idx].df2[10]  = nodes_df1[node_df1_idx].df1[10];
        }else if((li-1) <0 && (lk+1) <=cube_size-1){//only BI changes
          cube_df2_idx  = (BI-1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (cube_size-1)*cube_size*cube_size +lj*cube_size +lk+1;//only local i will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[10]  = nodes_df1[node_df1_idx].df1[10];
        }else{//only BK changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK+1;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li-1)*cube_size*cube_size +lj*cube_size +0;//only local k will be strtn point of next cube
          nodes_df2[node_df2_idx].df2[10]  = nodes_df1[node_df1_idx].df1[10];
        }
//ksi==10 ends
//ksi==11 starts//surfaces[i].nodes[(j-1)*dim_z +k+1].df2[11] = surfaces[i].nodes[j*dim_z +k].df1[11];
        if(lj-1 >=0 && lk+1 <=cube_size-1){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +(lj-1)*cube_size +lk+1;
          nodes_df2[node_df2_idx].df2[11]  = nodes_df1[node_df1_idx].df1[11];
        }else if(lj-1 <0 && lk+1 >cube_size-1){//comp domain changes both BJ and BK
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + (BJ-1) * num_cubes_z + BK+1;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +(cube_size-1)*cube_size + 0;//local j will be ending point of prev cube and local k strtng of next 
          nodes_df2[node_df2_idx].df2[11]  = nodes_df1[node_df1_idx].df1[11];
        }else if(lj-1 <0 && lk+1 <=cube_size-1){//only BJ changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + (BJ-1) * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +(cube_size-1)*cube_size +lk+1;//only local j will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[11]  = nodes_df1[node_df1_idx].df1[11];
        }else{//only BK changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK+1;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = li*cube_size*cube_size +(lj-1)*cube_size +0;//only local k will be strtn point of next cube
          nodes_df2[node_df2_idx].df2[11]  = nodes_df1[node_df1_idx].df1[11];
        }
//ksi==11 ends
//ksi==12 starts  //surfaces[i].nodes[(j+1)*dim_z +k-1].df2[12] = surfaces[i].nodes[j*dim_z +k].df1[12];
        if((lj+1) <=cube_size-1 && (lk-1) >=0){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +(lj+1)*cube_size +lk-1;
          nodes_df2[node_df2_idx].df2[12]  = nodes_df1[node_df1_idx].df1[12];
        }else if((lj+1) >cube_size-1 && (lk-1)<0){//comp domain changes both BJ and BK
          cube_df2_idx  = (BI) * num_cubes_y * num_cubes_z + (BJ+1) * num_cubes_z + (BK-1);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +0*cube_size +(cube_size-1);// local j starting point of next cube and local k will be endng point of prev cube
          nodes_df2[node_df2_idx].df2[12]  = nodes_df1[node_df1_idx].df1[12];
        }else if((lj+1) >cube_size-1 && (lk-1) >=0){//only BJ changes
          cube_df2_idx  = (BI) * num_cubes_y * num_cubes_z + (BJ+1) * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +0*cube_size +lk-1;//only local j will be starting point of next cube
          nodes_df2[node_df2_idx].df2[12]  = nodes_df1[node_df1_idx].df1[12];
        }else{//only BK changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK-1);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = li*cube_size*cube_size +(lj+1)*cube_size +cube_size-1;//only local k will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[12]  = nodes_df1[node_df1_idx].df1[12];
        }
//ksi==12 ends
//ksi==13 starts //surfaces[i].nodes[(j+1)*dim_z +k+1].df2[13] = surfaces[i].nodes[j*dim_z +k].df1[13];
        if((lj+1) <=cube_size-1 && (lk+1) <=cube_size-1){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +(lj+1)*cube_size +lk+1;
          nodes_df2[node_df2_idx].df2[13]  = nodes_df1[node_df1_idx].df1[13];
        }else if((lj+1) >cube_size-1 && (lk+1) >cube_size-1){//comp domain changes both BI and BK
          cube_df2_idx  = (BI) * num_cubes_y * num_cubes_z + (BJ+1) * num_cubes_z + (BK+1);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +0*cube_size +0;//both local j and local k will be starting point of next cube
          nodes_df2[node_df2_idx].df2[13]  = nodes_df1[node_df1_idx].df1[13];
        }else if((lj+1) >cube_size-1 && (lk+1) <=cube_size-1){//only BJ changes
          cube_df2_idx  = (BI) * num_cubes_y * num_cubes_z + (BJ+1) * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +0*cube_size +lk+1;//only local j will be starting point of next cube
          nodes_df2[node_df2_idx].df2[13]  = nodes_df1[node_df1_idx].df1[13];
        }else{//only BK changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK+1);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = li*cube_size*cube_size +(lj+1)*cube_size +0;//only local k will be starting point of next cube
          nodes_df2[node_df2_idx].df2[13]  = nodes_df1[node_df1_idx].df1[13];
        }
//ksi==13 ends
//ksi==14 starts//surfaces[i].nodes[(j-1)*dim_z +k-1].df2[14] = surfaces[i].nodes[j*dim_z +k].df1[14];
        if((lj-1) >=0 && (lk-1) >=0){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +(lj-1)*cube_size +(lk-1);
          nodes_df2[node_df2_idx].df2[14]  = nodes_df1[node_df1_idx].df1[14];
        }else if((lj-1) <0 && (lk-1) <0){//comp domain changes both BJ and BK
          cube_df2_idx  = (BI) * num_cubes_y * num_cubes_z + (BJ-1) * num_cubes_z + (BK-1);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +(cube_size-1)*cube_size + (cube_size-1);//both local j and local k will be ending point of prev cube
          nodes_df2[node_df2_idx].df2[14]  = nodes_df1[node_df1_idx].df1[14];
        }else if((lj-1)<0 && (lk-1) >=0){//only BJ changes
          cube_df2_idx  = (BI) * num_cubes_y * num_cubes_z + (BJ-1) * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li)*cube_size*cube_size +(cube_size-1)*cube_size +lk-1;//only local j will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[14]  = nodes_df1[node_df1_idx].df1[14];
        }else{//only BK changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + (BK-1);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = li*cube_size*cube_size +(lj-1)*cube_size +cube_size-1;//only local k will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[14]  = nodes_df1[node_df1_idx].df1[14];
        }
//ksi==14 ends
//ksi==15 starts//surfaces[i+1].nodes[(j-1)*dim_z +k].df2[15] = surfaces[i].nodes[j*dim_z +k].df1[15];
        if((li+1) <=cube_size-1 && (lj-1) >=0){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li+1)*cube_size*cube_size +(lj-1)*cube_size +(lk);
          nodes_df2[node_df2_idx].df2[15]  = nodes_df1[node_df1_idx].df1[15];
        }else if((li+1) >cube_size-1 && (lj-1)<0){//comp domain changes both BI and BJ
          cube_df2_idx  = (BI+1) * num_cubes_y * num_cubes_z + (BJ-1) * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (0)*cube_size*cube_size +(cube_size-1)*cube_size +lk;// local i starting point of next cube and local j will be endng point of prev cube
          nodes_df2[node_df2_idx].df2[15]  = nodes_df1[node_df1_idx].df1[15];
        }else if((li+1) >cube_size-1 && (lj-1) >=0){//only BI changes
          cube_df2_idx  = (BI+1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (0)*cube_size*cube_size +(lj-1)*cube_size +lk;//only local i will be starting point of next cube
          nodes_df2[node_df2_idx].df2[15]  = nodes_df1[node_df1_idx].df1[15];
        }else{//only BJ changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + (BJ-1) * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li+1)*cube_size*cube_size +(cube_size-1)*cube_size +lk;//only local k will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[15]  = nodes_df1[node_df1_idx].df1[15];
        }
//ksi==15 ends
//ksi==16 starts //surfaces[i-1].nodes[(j+1)*dim_z +k].df2[16] = surfaces[i].nodes[j*dim_z +k].df1[16];
        if((li-1) >=0 && (lj+1) <=cube_size-1){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li-1)*cube_size*cube_size +(lj+1)*cube_size +lk;
          nodes_df2[node_df2_idx].df2[16]  = nodes_df1[node_df1_idx].df1[16];
        }else if((li-1) <0 && (lj+1) >cube_size-1){//comp domain changes both BI and BJ
          cube_df2_idx  = (BI-1) * num_cubes_y * num_cubes_z + (BJ+1) * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (cube_size-1)*cube_size*cube_size +0*cube_size + lk;//local i will be ending point of prev cube and local j strtng of next 
          nodes_df2[node_df2_idx].df2[16]  = nodes_df1[node_df1_idx].df1[16];
        }else if((li-1) <0 && (lj+1) <=cube_size-1){//only BI changes
          cube_df2_idx  = (BI-1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (cube_size-1)*cube_size*cube_size +(lj+1)*cube_size +lk;//only local i will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[16]  = nodes_df1[node_df1_idx].df1[16];
        }else{//only BJ changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + (BJ+1) * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li-1)*cube_size*cube_size +0*cube_size +lk;//only local k will be strtn point of next cube
          nodes_df2[node_df2_idx].df2[16]  = nodes_df1[node_df1_idx].df1[16];
        }
//ksi==16 ends
//ksi==17 starts //surfaces[i+1].nodes[(j+1)*dim_z +k].df2[17] = surfaces[i].nodes[j*dim_z +k].df1[17];
        if((li+1) <=cube_size-1 && (lj+1) <=cube_size-1){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li+1)*cube_size*cube_size +(lj+1)*cube_size +lk;
          nodes_df2[node_df2_idx].df2[17]  = nodes_df1[node_df1_idx].df1[17];
        }else if((li+1) >cube_size-1 && (lj+1) >cube_size-1){//comp domain changes both BI and BJ
          cube_df2_idx  = (BI+1) * num_cubes_y * num_cubes_z + (BJ+1) * num_cubes_z + (BK);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (0)*cube_size*cube_size +(0)*cube_size +lk;//both local i and local j will be starting point of next cube
          nodes_df2[node_df2_idx].df2[17]  = nodes_df1[node_df1_idx].df1[17];
        }else if((li+1) >cube_size-1 && (lj+1) <=cube_size-1){//only BI changes
          cube_df2_idx  = (BI+1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (0)*cube_size*cube_size +(lj+1)*cube_size +lk;//only local i will be starting point of next cube
          nodes_df2[node_df2_idx].df2[17]  = nodes_df1[node_df1_idx].df1[17];
        }else{//only BJ changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + (BJ+1) * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li+1)*cube_size*cube_size +0*cube_size +lk;//only local k will be starting point of next cube
          nodes_df2[node_df2_idx].df2[17]  = nodes_df1[node_df1_idx].df1[17];
        }
//ksi==17 ends
//ksi==18 starts  //surfaces[i-1].nodes[(j-1)*dim_z +k].df2[18] = surfaces[i].nodes[j*dim_z +k].df1[18];
        if((li-1) >=0 && (lj-1) >=0){//i.e the the computation domain is still in the current cube
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li-1)*cube_size*cube_size +(lj-1)*cube_size +(lk);
          nodes_df2[node_df2_idx].df2[18]  = nodes_df1[node_df1_idx].df1[18];
        }else if((li-1) <0 && (lj-1) <0){//comp domain changes both BI and BJ
          cube_df2_idx  = (BI-1) * num_cubes_y * num_cubes_z + (BJ-1) * num_cubes_z + (BK);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (cube_size-1)*cube_size*cube_size +(cube_size-1)*cube_size + lk;//both local i and local j will be ending point of prev cube
          nodes_df2[node_df2_idx].df2[18]  = nodes_df1[node_df1_idx].df1[18];
        }else if((li-1)<0 && (lj-1) >=0){//only BI changes
          cube_df2_idx  = (BI-1) * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (cube_size-1)*cube_size*cube_size +(lj-1)*cube_size +lk;//only local i will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[18]  = nodes_df1[node_df1_idx].df1[18];
        }else{//only BJ changes
          cube_df2_idx  = BI * num_cubes_y * num_cubes_z + (BJ-1) * num_cubes_z + (BK);
          nodes_df2     = gv->fluid_grid->sub_fluid_grid[cube_df2_idx].nodes;
          node_df2_idx = (li-1)*cube_size*cube_size +(cube_size-1)*cube_size +lk;//only local k will be endg point of prev cube
          nodes_df2[node_df2_idx].df2[18]  = nodes_df1[node_df1_idx].df1[18];
        }
//ksi==18 ends	
     }//lk
    }// if cube2thread ends
   }//BK
     
      //printf("****************stream_distfunc() Exit*******");
 }


 void bounceback_rigidwalls(LV lv){
  //printf("******************Inside bounceback_rigidwalls************ at time  %d\n");
  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  
  Fluidgrid *fluid_grid;
  Fluidnode *nodes;
  
  int BI, BJ, BK; //to identify the Sub grids
  int cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  int li, lj, lk, node_idx;//local access point inside cube
  int dim_z;

  fluid_grid      = gv->fluid_grid;
  dim_z          = fluid_grid->z_dim;
     
  cube_size   = gv->cube_size; 
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;

//printf("*************ENTRY OF Bounceback\n");
      //print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
   /* 1.1 half-way bounce back on the bottom */
   //k = gv->kb;
   lk = gv->kb; BK = 0;
    
   for(BI =0; BI < num_cubes_x; ++BI)
    for(BJ =0; BJ < num_cubes_y; ++BJ){
      if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_x = starting_y = 0;
      stopping_x = stopping_y = cube_size-1;
      if(BI == 0) starting_x = gv->ib+1;//3
      
      if(BI == num_cubes_x-1) stopping_x = cube_size-4;//gv->ie-1
      
      if(BJ == 0) starting_y = gv->jb;//2
       
      if(BJ == num_cubes_y-1) stopping_y = cube_size-3;//je
       
     for(li=starting_x; li <= stopping_x; ++li)
      for(lj=starting_y; lj <= stopping_y; ++lj){
            node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. //lk =gv->kb
         nodes[node_idx].df2[3] =nodes[node_idx].df1[4];
         nodes[node_idx].df2[7] =nodes[node_idx].df1[8];
         nodes[node_idx].df2[10]=nodes[node_idx].df1[9];
         nodes[node_idx].df2[11]=nodes[node_idx].df1[12];
         nodes[node_idx].df2[13]=nodes[node_idx].df1[14];
        // printf("1.1 half-way bounce back on the bottom\n");
         if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
      && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size ){
         //   print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
         }
       }
      }//if cube2thread ends
     }
   
   /* 1.2 bounce back on the top*/
 //k=gv->ke;
  lk = cube_size-3; 
  BK = num_cubes_z -1;//lk=gv->ke;//BI=gv->ib+1;BI<=gv->ie-1;BJ=gv->jb;BJ<=gv->je
   for(BI =0; BI < num_cubes_x; ++BI)
    for(BJ =0; BJ < num_cubes_y; ++BJ){
       if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_x = starting_y = 0;
      stopping_x = stopping_y = cube_size-1;
      if(BI == 0) starting_x = gv->ib+1;//3
       
      if(BI == num_cubes_x-1) stopping_x = cube_size-4;//gv->ie-1
       
      if(BJ == 0) starting_y = gv->jb;//2
       
      if(BJ == num_cubes_y-1) stopping_y = cube_size-3;//je
       
     for(li=starting_x; li <= stopping_x; ++li)
      for(lj=starting_y; lj <= stopping_y; ++lj){
            node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. //lk =gv->ke
       nodes[node_idx].df2[4] =nodes[node_idx].df1[3];
       nodes[node_idx].df2[8] =nodes[node_idx].df1[7];
       nodes[node_idx].df2[9] =nodes[node_idx].df1[10];
       nodes[node_idx].df2[12]=nodes[node_idx].df1[11];
       nodes[node_idx].df2[14]=nodes[node_idx].df1[13];
      //  printf("1.2 bounce back on the top\n");
         if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
      && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size ){
          //  print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
         }

       }
      }//if cube2 thread ends
     } 
   
   /* 1.3 bounce back on the front*/
  //j=gv->jb;
       //BJ = gv->jb;BI=gv->ib+1;BI<=gv->ie-1;BK=gv->kb;BK<=gv->ke
   lj = gv->jb; BJ =0;
   for(BI =0; BI < num_cubes_x; ++BI)
    for(BK =0; BK < num_cubes_z; ++BK){
      if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_x = starting_z = 0;
      stopping_x  = stopping_z = cube_size-1;
      if(BI == 0) starting_x = gv->ib+1;//3
       
      if(BI == num_cubes_x-1) stopping_x = cube_size-4;//gv->ie-1
       
      if(BK == 0) starting_z = gv->kb;//2
       
      if(BK == num_cubes_x-1) stopping_z = cube_size-3;//gv->ke
       
     for(li=starting_x; li <= stopping_x; ++li)
      for(lk=starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. //lj =gv->jb =2
           nodes[node_idx].df2[6] =nodes[node_idx].df1[5];
           nodes[node_idx].df2[12]=nodes[node_idx].df1[11];
           nodes[node_idx].df2[13]=nodes[node_idx].df1[14];
           nodes[node_idx].df2[16]=nodes[node_idx].df1[15];
           nodes[node_idx].df2[17]=nodes[node_idx].df1[18];
         //  printf("1.2 bounce back on the front\n");
         if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
      && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size ){
          //  print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
          }
        }
       }//if cube2 thread ends
   }

   /* 1.4 bounce back on the rear*/
  //j=gv->je;
       //BJ = gv->je;BI=gv->ib+1;BI<=gv->ie-1;BK=gv->kb;BK<=gv->ke
   lj = cube_size-3; BJ = num_cubes_y-1;//gv->je;
   for(BI =0; BI < num_cubes_x; ++BI)
    for(BK =0; BK < num_cubes_z; ++BK){
      if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
      cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
      nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_x = starting_z = 0;
      stopping_x = stopping_z = cube_size-1;
      if(BI == 0) starting_x = gv->ib+1;//3
      
      if(BI == num_cubes_x-1) stopping_x = cube_size-4;//gv->ie-1
       
      if(BK == 0) starting_z = gv->kb;//2
       
      if(BK == num_cubes_x-1) stopping_z = cube_size-3;//gv->ke
       
     for(li=starting_x; li <= stopping_x; ++li)
      for(lk=starting_z; lk <= stopping_z; ++lk){
          node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. //lj =gv->je 
             nodes[node_idx].df2[5] =nodes[node_idx].df1[6];
             nodes[node_idx].df2[11]=nodes[node_idx].df1[12];
             nodes[node_idx].df2[14]=nodes[node_idx].df1[13];
             nodes[node_idx].df2[15]=nodes[node_idx].df1[16];
             nodes[node_idx].df2[18]=nodes[node_idx].df1[17];
             //printf("1.2 bounce back on the rear\n");
         if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
      && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size ){
            //print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
         }
       }
      }// if cube2thread ends
     }
// printf("EXIT bounceback\n");
   //       print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
         
   //printf("******************bounceback_rigidwalls Exit************\n ");
 
 }


 void compute_rho_and_u(LV lv){
   //printf("**********Inside compute_rho**********\n");
  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  
  Fluidgrid     *fluid_grid;
  Fluidnode     *nodes;
 
  int           ksi;
  double        s1, s2, s3, s4;
  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube
  
  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;

  fluid_grid      = gv->fluid_grid;
     
  cube_size   = gv->cube_size; 
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;
   
s1 =s2= s3 = s4 =0;
//printf("*************ENTRY OF RHO AND U\n");
//print_fluid_sub_grid(gv, 19,10,19, 19,10,21,gv->cube_size);
  // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
   for(BI =0; BI < num_cubes_x; ++BI)
   for(BJ =0; BJ < num_cubes_y; ++BJ)//for computing womega near bdy
    for(BK =0; BK < num_cubes_z; ++BK){
      if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
     cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
     nodes = fluid_grid->sub_fluid_grid[cube_idx].nodes;
     starting_x = starting_y = starting_z = 0;
     stopping_x = stopping_y = stopping_z = cube_size-1;
     if(BI == 0) starting_x = 3;//ib+1
      
     if(BI == num_cubes_x-1) stopping_x = cube_size-4;//ie-1
      
     if(BJ == 0) starting_y = 2;//jb
    
     if(BJ == num_cubes_y-1) stopping_y = cube_size-3;//je
     
     if(BK == 0) starting_z = 2;//kb
     
     if(BK == num_cubes_z-1) stopping_z = cube_size-3;//ke
     s1 =s2= s3 = s4 =0; 
    for(li=starting_x; li <= stopping_x; ++li)
      for(lj=starting_y; lj <= stopping_y; ++lj)
       for(lk=starting_z; lk <= stopping_z; ++lk){
        node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.

      s1=nodes[node_idx].df2[0];
      s2=gv->c[0][0]*nodes[node_idx].df2[0];
      s3=gv->c[0][1]*nodes[node_idx].df2[0];
      s4=gv->c[0][2]*nodes[node_idx].df2[0];
     if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
      && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size ){
     // printf("*************Start of KSI\n");
     // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
     }
     for (ksi=1;ksi<=18;ksi++)
       {
         s1 += nodes[node_idx].df2[ksi];
         s2 += gv->c[ksi][0]*nodes[node_idx].df2[ksi];
         s3 += gv->c[ksi][1]*nodes[node_idx].df2[ksi];
         s4 += gv->c[ksi][2]*nodes[node_idx].df2[ksi];
       }
        if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
      && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size ){
     // printf("*************END of KSI\n");
     // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
     }
    // if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
     // && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size )
      //printf(" B4....s2 , s3, s4  NODES :%.12f, %.12f, %.12f",s2,s3,s4);
     nodes[node_idx].rho=s1;/* Eqn 11 from paper...PN*/
     nodes[node_idx].vel_x=(s2+0.5*gv->dt * nodes[node_idx].elastic_force_x) / s1;/*Eqn 12*/
     nodes[node_idx].vel_y=(s3+0.5*gv->dt * nodes[node_idx].elastic_force_y) /s1;/*Eqn 12*/
     nodes[node_idx].vel_z=(s4+0.5*gv->dt * (nodes[node_idx].elastic_force_z +nodes[node_idx].rho*gv->g_l))/s1;
     /*if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
      && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size ){
      printf("*************EXIT OF RHO AND U\n");
      print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
      printf("s2 , s3, s4  NODES :%.12f, %.12f, %.12f",s2,s3,s4);
     }*/
     }
    } //if cube2thread ends
  }//For BK  
  //printf("*************EXIT OF RHO AND U\n");
  //print_fluid_sub_grid(gv, 19,10,19, 19,10,21,gv->cube_size);
   //printf("**********compute_rho Exit**********\n");
}

void moveFiberSheet(LV lv){
  #ifdef DEBUG_PRINT
  printf("************** Inside moveFibersheet**********\n");
  #endif //DEBUG_PRINT

  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  int total_threads = gv->total_threads;

  double       s1, s2, s3, dx, dy, dz, dx1, dy1, dz1;
  int          gi, gj, gk, i2, j2, k2, istart, jstart, kstart, jf, kf;
  int          total_fibers_row, total_fibers_clmn;
  Fiber        *fiberarray;



  dx = 1.0; 
  dy = 1.0; 
  dz = 1.0;/*already dimensionless*/
  
  dx1 = 1.0/dx;
  dy1 = 1.0/dy;
  dz1 = 1.0/dz;

  // double  PI   = 3.14159265358979;
  double  c_x  = PI/(2.0*dx);
  double  c_y  = PI/(2.0*dy);
  double  c_z  = PI/(2.0*dz);
  double  c_64 = 1.0/(6.4e1);
  double   rx  = 0.0; 
  double   ry  = 0.0; 
  double   rz  = 0.0;

  Fluidnode *nodes;
  int        BI, BJ, BK, cube_size, cube_idx, node_idx;
  int        num_cubes_x, num_cubes_y, num_cubes_z;
  int        li, lj, lk;

  istart = 0;  
  jstart = 0; 
  kstart = 0; 
  s1=0; 
  s2=0; 
  s3=0;

  total_fibers_row  = gv->fiber_shape->sheets[0].num_rows;
  total_fibers_clmn = gv->fiber_shape->sheets[0].num_cols;
  fiberarray        = gv->fiber_shape->sheets[0].fibers;
  
  cube_size         = gv->cube_size;
  num_cubes_x       = gv->num_cubes_x;
  num_cubes_y       = gv->num_cubes_y;
  num_cubes_z       = gv->num_cubes_z;

  for (jf = 0; jf < total_fibers_row; jf++) { //along z -axis
    /*find the begining index of the computational 4x4 square for each fiber pt.
     * using transpose in order to be consistent to the rest of vectors */
    if(fiber2thread(jf, total_fibers_row, total_threads) ==tid){
    for (kf = 0; kf < total_fibers_clmn; kf++) {//along y-axis
      s1=0; 
      s2=0; 
      s3=0;
      
      istart = floor(fiberarray[jf].nodes[kf].x*dx1-2)+1;
      jstart = floor(fiberarray[jf].nodes[kf].y*dy1-2)+1;
      kstart = floor(fiberarray[jf].nodes[kf].z*dz1-2)+1;
     // printf("istart :%d, jstart :%d, kstart:%d\n", istart, jstart, kstart);
      for (i2 = 0; i2 <= 3; i2++){
      gi = istart+i2;
      rx = dx * gi - fiberarray[jf].nodes[kf].x;/*dimensionless*/
     for (j2 = 0; j2 <= 3; j2++){
        gj=jstart+j2;
        ry=dy * gj - fiberarray[jf].nodes[kf].y;
      for (k2 = 0; k2 <= 3; k2++){
          gk  = kstart + k2;
          rz = dz * gk -fiberarray[jf].nodes[kf].z;
          //printf("fiberarray[jf].nodes[kf].z:%.12f\n",fiberarray[jf].nodes[kf].z);
              if( gi < 0 || gi >= gv->fluid_grid->x_dim) {
          fprintf(stderr, "gi out of bound: %d, x bound=%d!\n", gi, gv->fluid_grid->x_dim); exit(1);
        }
              if( gj < 0 || gj >= gv->fluid_grid->y_dim) {
          fprintf(stderr, "gj out of bound: %d, y bound=%d!\n", gj, gv->fluid_grid->y_dim); exit(1);
        }
              if( gk < 0 || gk >= gv->fluid_grid->z_dim) {
          fprintf(stderr, "gk out of bound: %d, z bound=%d!\n", gk, gv->fluid_grid->z_dim); exit(1);
        }
        li= gi%cube_size; lj = gj%cube_size; lk = gk%cube_size;
        BI = gi/cube_size; BJ = gj/cube_size; BK = gk/cube_size;
        /*li= istart%cube_size; lj = jstart%cube_size; lk = kstart%cube_size;
        BI = istart/cube_size; BJ = jstart/cube_size; BK = kstart/cube_size;*/
        /*if(li <= cube_size-1)BI = gi/cube_size; 
         else{
          BI = BI + 1;
          li = 0;
        }
        if(lj <= cube_size-1)BJ = gj/cube_size;
         else{
          BJ = BJ +1;
          lj = 0;
        }
        if(lk <= cube_size-1)BK = gk/cube_size;
         else{
          BK = BK +1;
          lk =0;
         } */
      cube_idx = BI *num_cubes_y *num_cubes_z + BJ *num_cubes_z +BK;
      nodes    = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
     // node_idx = (gi%cube_size) *cube_size*cube_size +(gj%cube_size) *cube_size +gk%cube_size;
      node_idx = (li) *cube_size*cube_size +(lj) *cube_size +lk;
      /*Eqn 20 for all x, y and z direction...PN*/
      //printf("#####SVALUES### :%.12f, %.12f, %.12f \n", s1,s2,s3);
      //printf("For Cube <%d,%d,%d> and node <%d,%d,%d> \n", BI, BJ, BK, li, lj, lk);
      //printf("velocities at nodeidx and cube_idx %.12f,%.12f,%.12f \n",nodes[node_idx].vel_x,nodes[node_idx].vel_y,nodes[node_idx].vel_z);
      //printf("rx:%.12f, ry:%.12f,rz:%.12f\n", rx, ry, rz);
      
      s1 += nodes[node_idx].vel_x*(1.0+cos(c_x*rx))*(1.0+cos(c_y*ry))*(1.0+cos(c_z*rz))*c_64;/*factors in nondimensionalizing h^(-3) and dxdydz cancel each other*/
      /*notice the difference for spreading and for interpolation here, cf is used for spreading*/
      s2 += nodes[node_idx].vel_y*(1.0+cos(c_x*rx))*(1.0+cos(c_y*ry))*(1.0+cos(c_z*rz))*c_64;
      s3 += nodes[node_idx].vel_z*(1.0+cos(c_x*rx))*(1.0+cos(c_y*ry))*(1.0+cos(c_z*rz))*c_64;
      
     /* if(li==lookup_fluid_start_x%cube_size && lj== lookup_fluid_start_y%cube_size && lk ==lookup_fluid_start_z%cube_size 
      && BI ==lookup_fluid_start_x/cube_size && BJ == lookup_fluid_start_y/cube_size && BK ==lookup_fluid_start_z/cube_size ){
       // printf("rx:%.12f, ry:%.12f,rz:%.12f\n", rx, ry, rz);
        //printf("cx:%.12f, cy:%.12f,cz:%.12f\n", c_x, c_y, c_z);
      }*/
      /*Eqn 20 for all x, y and z direction ends...PN*/
 // printf("\n");
      }
     }
    } /*end of the 3 loops for i2 j2 k2*/
 /*  printf("Before moving \n");
   print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
     printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);
     printf("SVALUES :%.12f, %.12f, %.12f \n", s1,s2,s3);*/
     
      fiberarray[jf].nodes[kf].x += gv->dt * s1;
      fiberarray[jf].nodes[kf].y += gv->dt * s2;
      fiberarray[jf].nodes[kf].z += gv->dt * s3;
    
      /* printf("After moving \n");
     printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);
     printf("SVALUES :%.12f, %.12f, %.12f \n", s1,s2,s3);*/
//     }//if ends
    }/*end of all points for a fiber kf ends*/
   }// if fiber2thread ends
  } //jf ends i.e fiber along row
  #ifdef DEBUG_PRINT
  printf("**************  moveFibersheet   Exit**********\n");
  #endif //DEBUG_PRINT
}

void init_df1(GV gv) {
  Fluidnode     *nodes;
  int           ksi;
  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube

  cube_size = gv->cube_size;
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;
/*PTHREAD_Change*/
  for(BI =0; BI < num_cubes_x; ++BI)
   for(BJ =0; BJ < num_cubes_y; ++BJ)
    for(BK =0; BK < num_cubes_z; ++BK){
     cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
     nodes = gv->fluid_grid->sub_fluid_grid[cube_idx].nodes;
     starting_x = starting_y = starting_z = 0;
     stopping_x = stopping_y = stopping_z = cube_size-1;
     if(BI == 0) starting_x = 2;
      
     if(BI == num_cubes_x-1) stopping_x = cube_size-3;
      
     if(BJ == 0) starting_y = 2;
      
     if(BJ == num_cubes_y-1) stopping_y = cube_size-3;
     
     if(BK == 0) starting_z = 2;
     
     if(BK == num_cubes_z-1) stopping_z = cube_size-3;
      
     for(li=starting_x; li <= stopping_x; ++li)
      for(lj=starting_y; lj <= stopping_y; ++lj)
       for(lk=starting_z; lk <= stopping_z; ++lk){
        node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube.
        for (ksi=0; ksi<=18; ksi++) 
         nodes[node_idx].df1[ksi]=nodes[node_idx].dfeq[ksi];            
       }//local k
    }//for BK
   /*PTHREAD_Change*/
 //printf("****************init_df1 EXIT ***************\n");
}

void copy_inout_to_df2(LV lv){
  //printf("***********************Inside copy_inout_to_df2***********\n");
  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  
  Fluidgrid     *fluid_grid;
  Fluidnode     *nodes;
 
  int           ksi;
  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_y, starting_z, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube
  int           dim_z;

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;

  fluid_grid      = gv->fluid_grid;
  dim_z          = fluid_grid->z_dim;
     
  cube_size   = gv->cube_size; 
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

   //i = gv->ib;BJ=gv->jb-1; BJ<=gv->je+1;BK=gv->kb-1; BK<=gv->ke+1;
  li  = gv->ib;BI=0;
   for (BJ=0; BJ< num_cubes_y; BJ++)
    for (BK=0; BK< num_cubes_x; BK++){
      if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
      cube_idx = BI *num_cubes_y*num_cubes_z + BJ *num_cubes_z + BK;
      nodes    = fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_y = starting_z = 0;
      stopping_y = stopping_z = cube_size-1;
     if(BJ == 0) starting_y = gv->jb-1;
      
     if(BJ == num_cubes_y-1) stopping_y = cube_size-2;//gv->je+1
     
     if(BK == 0) starting_z = gv->kb-1;
      
     if(BK == num_cubes_z-1) stopping_z = cube_size-2;//ke +1
      
     for(lj=starting_y; lj <= stopping_y; ++lj)
      for(lk=starting_z; lk <= stopping_z; ++lk){
        node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. li =gv->ib
         for (ksi=0; ksi<=18; ksi++)
                nodes[node_idx].df2[ksi] = gv->fluid_grid->inlet->nodes[((BJ*cube_size+lj)*dim_z +BK*cube_size+lk)].df_inout[0][ksi];//is this correct?
       }//lk
     }//if cube2thread ends
   }

  
  //i=gv->ib-1;
  
    li = gv->ib-1; BI =0;
   for (BJ=0; BJ< num_cubes_y; BJ++)
    for (BK=0; BK< num_cubes_x; BK++){
      if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
      cube_idx = BI *num_cubes_y*num_cubes_z + BJ *num_cubes_z + BK;
      nodes    = fluid_grid->sub_fluid_grid[cube_idx].nodes;
      starting_y = starting_z = 0;
      stopping_y = stopping_z = cube_size-1;
     if(BJ == 0) starting_y = gv->jb-1;
      
     if(BJ == num_cubes_y-1) stopping_y = cube_size-2;//gv->je+1
     
     if(BK == 0) starting_z = gv->kb-1;
      
     if(BK == num_cubes_z-1) stopping_z = cube_size-2;//ke+1
      
    for(lj=starting_y; lj <= stopping_y; ++lj)
      for(lk=starting_z; lk <= stopping_z; ++lk){
        node_idx = li * cube_size * cube_size + lj * cube_size + lk; //local node index inside a cube. li =gv->ib-1
        for (ksi=0; ksi<=18; ksi++)       
        nodes[node_idx].df2[ksi] = gv->fluid_grid->inlet->nodes[((BJ*cube_size+lj)*dim_z +BK*cube_size+lk)].df_inout[0][ksi];
       }
     }//if cube2thread ends
  }
  
  
  /* outlet */
  //i=gv->ie;
  li = cube_size-3; BI = num_cubes_x-1;
   for (BJ=0; BJ< num_cubes_y; BJ++)
    for (BK=0; BK< num_cubes_x; BK++){
      if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
      cube_idx = BI *num_cubes_y*num_cubes_z + BJ *num_cubes_z + BK;
      nodes    = fluid_grid->sub_fluid_grid[cube_idx].nodes;
     starting_y = starting_z = 0;
     stopping_y = stopping_z = cube_size-1;
    if(BJ == 0) starting_y = gv->jb-1;
    
    if(BJ == num_cubes_y-1) stopping_y = cube_size-2;//gv->je+1
     
    if(BK == 0) starting_z = gv->kb-1;
     
    if(BK == num_cubes_z-1) stopping_z = cube_size-2;//ke +1
     
    for(lj=starting_y; lj <= stopping_y; ++lj)
      for(lk=starting_z; lk <= stopping_z; ++lk){
        node_idx = li * cube_size * cube_size + lj * cube_size + lk;//local node index inside a cube. li =gv->ie
        for (ksi=0; ksi<=18; ksi++)
        nodes[node_idx].df2[ksi] = gv->fluid_grid->outlet->nodes[((BJ*cube_size + lj)*dim_z +BK*cube_size+lk)].df_inout[1][ksi];
       }
     }//if cube2thread ends
    }   
  
  //i=gv->ie+1;
  li = cube_size-2; BI = num_cubes_x-1;
   for (BJ=0; BJ< num_cubes_y; BJ++)
    for (BK=0; BK< num_cubes_x; BK++){
      if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
      cube_idx = BI *num_cubes_y*num_cubes_z + BJ *num_cubes_z + BK;
      nodes    = fluid_grid->sub_fluid_grid[cube_idx].nodes;
     starting_y = starting_z = 0;
     stopping_y = stopping_z = cube_size-1;
    if(BJ == 0) starting_y = gv->jb-1;
    
    if(BJ == num_cubes_y-1) stopping_y = cube_size-2;//gv->je+1
     
    if(BK == 0) starting_z = gv->kb-1;
     
    if(BK == num_cubes_z-1) stopping_z = cube_size-2;//ke +1
    
    for(lj=starting_y; lj <= stopping_y; ++lj)
      for(lk=starting_z; lk <= stopping_z; ++lk){
        node_idx = li * cube_size * cube_size + lj * cube_size + lk;//local node index inside a cube. li =gv->ie+1
        for (ksi=0; ksi<=18; ksi++)
        nodes[node_idx].df2[ksi]   = gv->fluid_grid->outlet->nodes[((BJ*cube_size+lj)*dim_z +BK*cube_size+lk)].df_inout[1][ksi];
      }
     }// if cube2thread ends
    }   

   //printf("*********************** copy_inout_to_df2 Exit ***********\n");
   
}

void replace_old_DF(LV lv){
//printf("Inside replace_oldDF\n");
  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  
  Fluidgrid     *fluid_grid;
  Fluidnode     *nodes;
 
  int           ksi;
  int           BI, BJ, BK; //to identify the Sub grids
  int           cube_size, num_cubes_x, num_cubes_y, num_cubes_z, cube_idx;
  int           starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;//To identify buffer zone
  int           li, lj, lk, node_idx;//local access point inside cube

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;

  fluid_grid      = gv->fluid_grid;
     
  cube_size   = gv->cube_size; 
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;
  
 /* replacing the old d.f. values by the newly computed ones */
 for(BI =0; BI < num_cubes_x; ++BI)
   for(BJ =0; BJ < num_cubes_y; ++BJ)
    for(BK =0; BK < num_cubes_z; ++BK){
      if(cube2thread(BI, BJ, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) ==tid){
     cube_idx = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK;
     nodes = fluid_grid->sub_fluid_grid[cube_idx].nodes;
     starting_x = starting_y = starting_z = 0;
     stopping_x = stopping_y = stopping_z = cube_size-1;
     if(BI == 0) starting_x = 1;//ib-1
     
     if(BI == num_cubes_x-1) stopping_x = cube_size-2;//ie+1
     
     if(BJ == 0) starting_y = 1;//jb-1
      
     if(BJ == num_cubes_y-1) stopping_y = cube_size-2;//je+1
      
     if(BK == 0) starting_z = 1;//kb-1
      
     if(BK == num_cubes_z-1) stopping_z = cube_size-2;//ke+1
      
    for(li=starting_x; li <= stopping_x; ++li)
      for(lj=starting_y; lj <= stopping_y; ++lj)
       for(lk=starting_z; lk <= stopping_z; ++lk){
        node_idx = li * cube_size * cube_size + lj * cube_size + lk;
        for (ksi=0;ksi<=18;ksi++)
            nodes[node_idx].df1[ksi]=nodes[node_idx].df2[ksi];
        
      }//lk
     }
    }//BK

//printf("replace_oldDF EXIT\n");
}

void periodicBC(LV lv){
   /* use periodic boundary condition on y and z directions,as a result of this, we are actually
   *            doing simulation on a array of micro-channels */
  /*along y-direction,j=1 & n2-1 */
  //printf("***************Inside periodicBC********\n");
  int tid;
  GV gv = lv->gv;   
  tid   = lv->tid;
  
  int           ksi;
  Fluidgrid     *fluid_grid;
  Fluidnode     *nodes_first, *nodes_last;
  int           BI, BJ, BJ_first, BJ_last, BK, BK_first, BK_last;
  int           num_cubes_x, num_cubes_y, num_cubes_z, cube_size;
  int           li, lj, lk;
  int           cube_idx_first, cube_idx_last, node_idx_first, node_idx_last;
  int           starting_x, starting_y, starting_z, stopping_x, stopping_y, stopping_z;

  int P = gv->P;
  int Q = gv->Q;
  int R = gv->R;
  //pthread_mutex_t lock_Fluid = gv->lock_Fluid;
  
  fluid_grid   = gv->fluid_grid;
     
  cube_size   = gv->cube_size; 
  num_cubes_x = gv->num_cubes_x;
  num_cubes_y = gv->num_cubes_y;
  num_cubes_z = gv->num_cubes_z;

  BJ_first = 0; BJ_last = num_cubes_y-1;
 for(BI =0; BI < num_cubes_x; ++BI)
    for(BK =0; BK < num_cubes_z; ++BK){
      if(cube2thread(BI, BJ_first, BK, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R)==tid){
     cube_idx_first = BI * num_cubes_y * num_cubes_z + BJ_first * num_cubes_z + BK;
     nodes_first = fluid_grid->sub_fluid_grid[cube_idx_first].nodes;
    
    
     cube_idx_last  = BI * num_cubes_y * num_cubes_z + BJ_last * num_cubes_z + BK;     
     nodes_last  = fluid_grid->sub_fluid_grid[cube_idx_last].nodes;
   

     starting_x = starting_z = 0;
     stopping_x = stopping_z = cube_size-1;
     if(BI == 0) starting_x = 3;//ib+1
      
     if(BI == num_cubes_x-1) stopping_x = cube_size-4;//ie-1
      
     if(BK == 0) starting_z = gv->kb;//kb-1
     
     if(BK == num_cubes_z-1) stopping_z = cube_size-3;//ke
     
    for(li=starting_x; li <= stopping_x; ++li)
      for(lk=starting_z; lk <= stopping_z; ++lk){
     
        node_idx_first = li * cube_size * cube_size + (gv->jb-1) * cube_size + lk;
     
      
        node_idx_last = li * cube_size * cube_size + (cube_size-2) * cube_size + lk;//je+1
     
      for (ksi=0; ksi<=18; ksi++){

  nodes_first[node_idx_first].df2[ksi] = nodes_last[li*cube_size*cube_size + (cube_size-3)*cube_size + lk].df2[ksi];//gv->je*dim_z +k
  
 //pthread_mutex_lock(&(gv->lock_Fluid));
  nodes_last[node_idx_last].df2[ksi]   = nodes_first[li*cube_size*cube_size + (gv->jb)*cube_size + lk].df2[ksi];//gv->jb*dim_z +k
  //pthread_mutex_unlock(&(gv->lock_Fluid));
     }
    }
   }
  }
  //pthread_barrier_wait(&(gv->barr));
  //May need barrier

  /* along z-direction, z=1 & n3-1 */

        BK_first = 0; BK_last = num_cubes_z-1;
 for(BI =0; BI < num_cubes_x; ++BI)
    for(BJ =0; BJ < num_cubes_y; ++BJ){
      if(cube2thread(BI, BJ, BK_first, num_cubes_x, num_cubes_y, num_cubes_z, P, Q,  R) == tid){
     cube_idx_first = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK_first;
     cube_idx_last  = BI * num_cubes_y * num_cubes_z + BJ * num_cubes_z + BK_last;
     nodes_first = fluid_grid->sub_fluid_grid[cube_idx_first].nodes;
     nodes_last  = fluid_grid->sub_fluid_grid[cube_idx_last].nodes;
     starting_x = starting_y = 0;
     stopping_x = stopping_y = cube_size-1;
     if(BI == 0) starting_x = 3;//ib+1
      
     if(BI == num_cubes_x-1) stopping_x = cube_size-4;//ie-1
      
     if(BJ == 0) starting_y = gv->jb;//jb
      
     if(BJ == num_cubes_z-1) stopping_y = cube_size-3;//ke
      
      
    for(li=starting_x; li <= stopping_x; ++li)
      for(lj=starting_y; lj <= stopping_y; ++lj){
        node_idx_first = li * cube_size * cube_size + lj * cube_size + gv->kb-1;
        node_idx_last = li * cube_size * cube_size + lj * cube_size + cube_size-2;//ke+1
      for (ksi=0; ksi<=18; ksi++){
  nodes_first[node_idx_first].df2[ksi] = nodes_last[li*cube_size*cube_size + lj*cube_size + cube_size-3].df2[ksi];//j*dim_z +gv->ke
  //lock : So that BK_first thread may not write to node_last
  //pthread_mutex_lock(&(gv->lock_Fluid));
  nodes_last[node_idx_last].df2[ksi]   = nodes_first[li*cube_size*cube_size + lj*cube_size + gv->kb].df2[ksi];//j*dim_z +gv->kb
  //pthread_mutex_unlock(&(gv->lock_Fluid));
  //unlcok
      }
     }
    }
  }
  
  
 // printf("***************periodicBC Exit********\n");
}

//TODO: define compute_womega_(GV gv) here???
 
void* do_thread(void* v){
  LV lv= (LV) v;
  GV gv = lv->gv;
  int barrcheck;
  int tid   = lv->tid;
  char filename[80];
  int my_rank;
     
  #ifdef DEBUG_PRINT    
  printf("Inside do_thread :Started by Threadid: %d\n",lv->tid);
  #endif //DEBUG_PRINT
   while (gv->time <= gv->TIME_STOP) {
    #ifdef DEBUG_PRINT
     printf("\n\nStart time step %d ...\n", gv->time);
     #endif //DEBUG_PRINT
//     print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z);
//     print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);
     
     /* Pthread code starts...allocating as many threads as no of fibers*/
  /* Pthread code ends*/
     compute_bendingforce(lv);
     printf("Tid%d: After compute_bending_force\n", tid);
     /*printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);*/
     //print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
     //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

  //Error Chek#5
     pthread_barrier_wait(&(gv->barr));
     /*barrcheck = pthread_barrier_wait(&barr);
    if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD)
    {
        fprintf(stderr,"Could not wait on barrier\n");
        exit(1);
    }*/
#if 0
    if (tid == 0){
      my_rank = 1;
      sprintf(filename, "Fiber%d_bending_force%d.dat", my_rank, gv->time);
      save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
    }
    pthread_barrier_wait(&(gv->barr));
#endif

     compute_stretchingforce(lv);
     printf("Tid%d: After compute_streching_force.\n", tid);
     /*printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);
     print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);*/
     //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

     //Error Chek#6
     pthread_barrier_wait(&(gv->barr));
     /*barrcheck = pthread_barrier_wait(&barr);
    if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD)
    {
        fprintf(stderr,"Could not wait on barrier\n");
        exit(1);
    }*/
#if 0
    if (tid == 0){
      my_rank = 1;
      sprintf(filename, "Fiber%d_streching_force%d.dat", my_rank, gv->time);
      save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
    }
    pthread_barrier_wait(&(gv->barr));
#endif

     compute_elasticforce(lv);
     printf("Tid%d: B4 get_influence_domain\n", tid);
     /*printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);*/
     //print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z);
     //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

     //Error Chek#5
     pthread_barrier_wait(&(gv->barr));
    /* barrcheck = pthread_barrier_wait(&barr);
    if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD)
    {
        fprintf(stderr,"Could not wait on barrier\n");
        exit(1);
    }*/
#if 0
    if (tid == 0){
      my_rank = 1;
      sprintf(filename, "Fiber%d_elasitc_force%d.dat", my_rank, gv->time);
      save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
    }
    pthread_barrier_wait(&(gv->barr));
#endif

     get_influentialdomain_fluid_grid_and_SpreadForce(lv);
     #ifdef DEBUG_PRINT
     printf("After get_influence_domain\n");
     #endif //DEBUG_PRINT
     /*printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);*/
    // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z,gv->cube_size);
     //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

      //Error Chek#5
     pthread_barrier_wait(&(gv->barr));
    /* barrcheck = pthread_barrier_wait(&barr);
    if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD)
    {
        fprintf(stderr,"Could not wait on barrier\n");
        exit(1);
    }*/

#if 0
    if (tid == 0){
      my_rank = 0;
      sprintf(filename, "Fluid%d_get_SpreadForce_step%d.dat", my_rank, gv->time);
      save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);

      my_rank = 1;
      sprintf(filename, "Fiber%d_SpreadForce_step%d.dat", my_rank, gv->time);
      save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
    }
#endif

    compute_eqlbrmdistrfuncDF1(lv);
    #ifdef DEBUG_PRINT
    printf("After compute DF1 \n");
    #endif //DEBUG_PRINT
 /*    printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);*/
   //  print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z,gv->cube_size);
     //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

      //Error Chek#5
    pthread_barrier_wait(&(gv->barr));
    /* barrcheck = pthread_barrier_wait(&barr);
    if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD){
        fprintf(stderr,"Could not wait on barrier\n");
        exit(1);
    }*/

#if 0 //Verify results
    if (tid == 0){
      my_rank = 0;
      sprintf(filename, "Fluid%d_compute_DF1_step%d.dat", my_rank, gv->time);
      save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
    }
    pthread_barrier_wait(&(gv->barr));
#endif

    stream_distrfunc(lv);
    #ifdef DEBUG_PRINT
    printf("After streaming\n");
    #endif //DEBUG_PRINT

/*     printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);*/
    // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z,gv->cube_size);
     //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

      //Error Chek#5
    pthread_barrier_wait(&(gv->barr));
    /* barrcheck = pthread_barrier_wait(&barr);
    if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD){
        fprintf(stderr,"Could not wait on barrier\n");
        exit(1);
    }*/

#ifdef SAVE //Verify results
    if (tid == 0){
      my_rank = 0;
      sprintf(filename, "Fluid%d_streaming_step%d.dat", my_rank, gv->time);
      save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
    }
    pthread_barrier_wait(&(gv->barr));
#endif

    bounceback_rigidwalls(lv);
    #ifdef DEBUG_PRINT
    printf("After bounceback\n");
    #endif //DEBUG_PRINT
    /* printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);*/
     //print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z,gv->cube_size);
     //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

      //Error Chek#5
     pthread_barrier_wait(&(gv->barr));
    /* barrcheck = pthread_barrier_wait(&barr);
    if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD)
    {
        fprintf(stderr,"Could not wait on barrier\n");
        exit(1);
    }*/
     
    compute_rho_and_u(lv);
    #ifdef DEBUG_PRINT
    printf("After compute rho and u \n");
    #endif //DEBUG_PRINT
     /*printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);
     print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z,gv->cube_size);*/
     //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);
    //print_fluid_sub_grid(gv, 19,10,19, 21,12,12,gv->cube_size);

      //Error Chek#5
    pthread_barrier_wait(&(gv->barr));
    // barrcheck = pthread_barrier_wait(&barr);
    // if (barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD){
    //   fprintf(stderr,"Could not wait on barrier\n");
    //   exit(1);
    // }

// // VERIFY    
//      //compute_womega(gv) here
//      moveFiberSheet(lv);
//      #ifdef DEBUG_PRINT
//      printf("After moving fibersheet \n");
//      #endif //DEBUG_PRINT
//      /*printf("Printing for Corner Points(z,y) : 0,0 \n");
//      print_fiber_sub_grid(gv,0, 0, 0, 0);
//      printf("Printing for Corner Points(z,y) : 51,0 \n");
//      print_fiber_sub_grid(gv,0, 51, 0, 51);
//      printf("Printing for Corner Points(z,y) : 0,51 \n");
//      print_fiber_sub_grid(gv,51, 0, 51, 0);
//      printf("Printing for Corner Points (z,y): 51,51 \n");
//      print_fiber_sub_grid(gv,51, 51, 51, 51);
//      print_fluid_sub_grid(gv, 19,10,19, 21,12,12,gv->cube_size);
//      printf("LOOKUp INDICES START\n");*/
//      //print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z,gv->cube_size);
//      //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

//       //Error Chek#5
//      pthread_barrier_wait(&(gv->barr));
//     /* barrcheck = pthread_barrier_wait(&barr);
//     if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD)
//     {
//         fprintf(stderr,"Could not wait on barrier\n");
//         exit(1);
//     }*/
     
//      copy_inout_to_df2(lv);
//      #ifdef DEBUG_PRINT
//      printf("After  copy_inout_to_df2 \n");
//      #endif //DEBUG_PRINT
//      /*printf("Printing for Corner Points(z,y) : 0,0 \n");
//      print_fiber_sub_grid(gv,0, 0, 0, 0);
//      printf("Printing for Corner Points(z,y) : 51,0 \n");
//      print_fiber_sub_grid(gv,0, 51, 0, 51);
//      printf("Printing for Corner Points(z,y) : 0,51 \n");
//      print_fiber_sub_grid(gv,51, 0, 51, 0);
//      printf("Printing for Corner Points (z,y): 51,51 \n");
//      print_fiber_sub_grid(gv,51, 51, 51, 51);*/
//     // print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z,gv->cube_size);
//      //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

//       //Error Chek#5
//      pthread_barrier_wait(&(gv->barr));
//     /* barrcheck = pthread_barrier_wait(&barr);
//     if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD)
//     {
//         fprintf(stderr,"Could not wait on barrier\n");
//         exit(1);
//     }*/
    
//      replace_old_DF(lv);
//      #ifdef DEBUG_PRINT
//      printf("After  replace_old_DF \n");
//      #endif //DEBUG_PRINT
//      /*printf("Printing for Corner Points(z,y) : 0,0 \n");
//      print_fiber_sub_grid(gv,0, 0, 0, 0);
//      printf("Printing for Corner Points(z,y) : 51,0 \n");
//      print_fiber_sub_grid(gv,0, 51, 0, 51);
//      printf("Printing for Corner Points(z,y) : 0,51 \n");
//      print_fiber_sub_grid(gv,51, 0, 51, 0);
//      printf("Printing for Corner Points (z,y): 51,51 \n");
//      print_fiber_sub_grid(gv,51, 51, 51, 51);*/
//      //print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z,gv->cube_size);
//      //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

//       //Error Chek#5
//      pthread_barrier_wait(&(gv->barr));
//     /* barrcheck = pthread_barrier_wait(&barr);
//     if(barrcheck != 0 && barrcheck != PTHREAD_BARRIER_SERIAL_THREAD)
//     {
//         fprintf(stderr,"Could not wait on barrier\n");
//         exit(1);
//     }*/

//      // periodicBC(lv);
//      #ifdef DEBUG_PRINT
//      printf("After PeriodicBC\n");
//      #endif //DEBUG_PRINT
//      /*printf("Printing for Corner Points(z,y) : 0,0 \n");
//      print_fiber_sub_grid(gv,0, 0, 0, 0);
//      printf("Printing for Corner Points(z,y) : 51,0 \n");
//      print_fiber_sub_grid(gv,0, 51, 0, 51);
//      printf("Printing for Corner Points(z,y) : 0,51 \n");
//      print_fiber_sub_grid(gv,51, 0, 51, 0);
//      printf("Printing for Corner Points (z,y): 51,51 \n");
//      print_fiber_sub_grid(gv,51, 51, 51, 51);*/
//      //print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z, gv->cube_size);
//      //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);

// //VERIFY

    pthread_barrier_wait(&(gv->barr));

    #if 1
    if(tid == 0)
     printf("End of time step %d\n", gv->time);
    #endif //DEBUG_PRINT
     
     //if(tid==0)
    // pthread_mutex_lock(&(gv->lock_Fluid));
     if(tid==0)
     gv->time  += gv->dt;
     //if(gv->time == gv->TIME_STOP)
     // return NULL;
     pthread_barrier_wait(&(gv->barr));
     //pthread_mutex_lock(&(gv->lock_Fluid));
     
    //print_fluid_sub_grid(gv,lookup_fluid_start_x, lookup_fluid_start_y, lookup_fluid_start_z, lookup_fluid_end_x, lookup_fluid_end_y, lookup_fluid_end_z);
     //print_fiber_sub_grid(gv,lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z);
   }//while ends
   return NULL;
}

int fiber2thread(int fiber_row_i, int num_fibers, int num_threads){
  int tid;
  int Block_size;
  Block_size = num_fibers/num_threads;
  if( num_fibers%num_threads !=0  && Block_size * num_threads < num_fibers)
  {
    Block_size = Block_size +1;
  }
  tid        = fiber_row_i/Block_size;/*ith fiber row.*/
 return tid;
 /*return 0;for 1 thread*/
}

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

void needs_argument(int i, int argc, const char *flag) {
  if (i+1 >= argc) {
    fprintf(stderr, "error: Flag \"%s\" requires an argument\n", flag);
    abort();
  }
}

int main(int argc, char* argv[]) {
  LV lvs;
  GV gv;
  gv = (GV) calloc (1, sizeof(*gv));

  // Parse command line
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-steps")) {
      needs_argument(i, argc, "-steps");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-steps %d\" must be > 0\n", value);
        abort();
      }
      gv->TIME_STOP = value;
    }

    if (!strcmp(argv[i], "-fluid_grid_xyz")) {
      gv->fluid_grid = (Fluidgrid*) calloc(1, sizeof(Fluidgrid));
      // gv->fluid_grid->sub_fluid_grid = (Sub_Fluidgrid*) calloc (1, sizeof(Sub_Fluidgrid));

      needs_argument(i, argc, "-fluid_grid_x");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fluid_grid_x %d\" must be > 0\n", value);
        abort();
      }
      gv->fluid_grid->x_dim = value;

      needs_argument(i, argc, "-fluid_grid_y");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fluid_grid_y %d\" must be > 0\n", value);
        abort();
      }
      gv->fluid_grid->y_dim = value;

      needs_argument(i, argc, "-fluid_grid_z");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-fluid_grid_z %d\" must be > 0\n", value);
        abort();
      }
      gv->fluid_grid->z_dim = value;
    }

    if (!strcmp(argv[i], "-cube_size")) {
      needs_argument(i, argc, "-cube_size");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-cube_size %d\" must be > 0\n", value);
        abort();
      }
      gv->cube_size = value;
    }

    if (!strcmp(argv[i], "-thread_per_task_xyz")) {
      needs_argument(i, argc, "-thread_per_task_x");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-thread_per_task_x %d\" must be > 0\n", value);
        abort();
      }
      gv->P = value;

      needs_argument(i, argc, "-thread_per_task_y");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-thread_per_task_y %d\" must be > 0\n", value);
        abort();
      }
      gv->Q = value;

      needs_argument(i, argc, "-thread_per_task_z");
      value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-thread_per_task_z %d\" must be > 0\n", value);
        abort();
      }
      gv->R = value;
    }

    if (!strcmp(argv[i], "-num_fibersht")) {
      needs_argument(i, argc, "-num_fibersht");
      int value = atoi(argv[++i]);
      if (value <= 0) {
        fprintf(stderr, "error: Invalid flag \"-num_fibersht %d\" must be > 0\n", value);
        abort();
      }
      gv->fiber_shape = (Fibershape*) calloc (1, sizeof(Fibershape));
      gv->fiber_shape->num_sheets = value;
      gv->fiber_shape->sheets = (Fibersheet*) calloc (value, sizeof(Fibersheet));
    }

    if (!strcmp(argv[i], "-fibersht_width_height")) {
      for(int j = 0; j < gv->fiber_shape->num_sheets; j++){
        needs_argument(i, argc, "-fibersht_width");
        double value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_width %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].width = value;

        needs_argument(i, argc, "-fibersht_height");
        value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_height %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].height = value;
      }
    }

    if (!strcmp(argv[i], "-fibersht_row_clmn")) {
      for(int j = 0; j < gv->fiber_shape->num_sheets; j++){
        needs_argument(i, argc, "-fibersht_row");
        int value = atoi(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_row %d\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].num_rows = value;

        needs_argument(i, argc, "-fibersht_clmn");
        value = atoi(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_clmn %d\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].num_cols = value;
      }
    }

    if (!strcmp(argv[i], "-fibersht_xyz_0")) {
      for(int j = 0; j < gv->fiber_shape->num_sheets; j++){
        needs_argument(i, argc, "-fibersht_x0");
        double value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_x0 %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].x_orig = value;

        needs_argument(i, argc, "-fibersht_y0");
        value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_y0 %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].y_orig = value;

        needs_argument(i, argc, "-fibersht_z0");
        value = atof(argv[++i]);
        if (value <= 0) {
          fprintf(stderr, "error: Invalid flag \"-fibersht_z0 %f\" must be > 0\n", value);
          abort();
        }
        gv->fiber_shape->sheets[j].z_orig = value;
      }
    }
  }

  printf("***********IB Simulation using Pthreads cube v 4.1 starts************\n");
  int i;

  double fibersheet_w     = gv->fiber_shape->sheets[0].width;
  double fibersheet_h     = gv->fiber_shape->sheets[0].height;
  int total_fibers_row    = gv->fiber_shape->sheets[0].num_rows;  /* no of fibres along height */
  int total_fibers_clmn   = gv->fiber_shape->sheets[0].num_cols;  /* no of fibres along width  or column should be 512 for 1024:*/

  int fluidgrid_z         = gv->fluid_grid->x_dim; // along z
  int fluidgrid_y         = gv->fluid_grid->y_dim; // along y
  int fluidgrid_x         = gv->fluid_grid->z_dim; // along x
  double fibersheet_xo    = gv->fiber_shape->sheets[0].x_orig; //initial position of the fiber sheet chosen somewher in the middle of grid 
  double fibersheet_yo    = gv->fiber_shape->sheets[0].y_orig;  
  double fibersheet_zo    = gv->fiber_shape->sheets[0].z_orig;

  int k_cubedim           = gv->cube_size;
  //int total_threads       = atoi(argv[12]);

  int P                   = gv->P;
  int Q                   = gv->Q;
  int R                   = gv->R;

  int total_threads       = P*Q*R;
  printf("*****INPUT*********\n");
  printf("fibersheet_w:%f , fibersheet_h:%f, fibers_row:%d, fibers_clmn:%d \n ", fibersheet_w, fibersheet_w, total_fibers_row, total_fibers_row);
  printf("elem_z:%d , elem_y :%d, elem_x: %d\n ", fluidgrid_z, fluidgrid_y, fluidgrid_x);
  printf("fibersheet_xo:%f , fibersheet_yo :%f, fibersheet_zo: %f\n ", fibersheet_xo, fibersheet_yo, fibersheet_zo);
  printf("P :%d, Q:%d, R:%d \n TuningFactor:cubedim::%d \n",P,Q,R, k_cubedim);

  printf("***********total_threads = P*Q*R =%d  ************\n\n",total_threads);
  //Error Chek#2
  if (total_fibers_row%2 !=0 ||  total_threads%2 !=0){
    fprintf(stderr,"Check #Fiber_row ==:%d and total threads ==:%d :Should be multiple of 2", total_fibers_row, total_threads);
    exit(1);
  }

  if ((fluidgrid_z/k_cubedim) % R !=0){
    fprintf(stderr,"Check Fluid_elem_z== :%d and  R==:%d :Should be divisible \n", fluidgrid_z, R);
    exit(1);
  }

  if ((fluidgrid_y/k_cubedim) % Q !=0){
    fprintf(stderr,"Check Fluid_elem_y== :%d and  Q==:%d :Should be divisible \n", fluidgrid_y, Q);
    exit(1);
  }

  if ((fluidgrid_x/k_cubedim) % P !=0){
    fprintf(stderr,"Check Fluid_elem_x== :%d and  P==:%d :Should be divisible \n", fluidgrid_x, P);
    exit(1);
  }

   Fibershape* fiber_shape = gen_fiber_shape(gv, fibersheet_w, fibersheet_h, 
					     total_fibers_clmn, total_fibers_row,
					     fibersheet_xo, fibersheet_yo, fibersheet_zo);

  Fluidgrid*  fluid_grid = gen_fluid_grid(gv, fluidgrid_x, fluidgrid_y, fluidgrid_z, k_cubedim);
   
  init_gv(gv, fiber_shape, fluid_grid, k_cubedim);

  gv->P   = P;
  gv->Q   = Q;
  gv->R   = R;
  
  //allocating memory for mutex
   gv->lock_Fluid = (pthread_mutex_t*)malloc(sizeof(*gv->lock_Fluid)*total_threads);
  
  //Error Chek#3
   //initiallise mutex
  for(i=0; i <total_threads; i++){
    if (pthread_mutex_init(&gv->lock_Fluid[i], NULL)){
      fprintf(stderr,"Unable to initialize a mutex\n");
      exit(1);
    }
  }
  

  /*for (i = 0; i < MUTNUM; i++)
    pthread_mutex_init(&locks[i], NULL);*/

  //imitiallising barrier
  pthread_barrier_t barr;
  //Error Chek#4
  //pthread_barrier_init(&barr, NULL, total_threads+1);
  if (pthread_barrier_init(&barr, NULL, total_threads)){
    fprintf(stderr,"Could not create a barrier\n");
    exit(1);
  }
  gv->barr  = barr;

   //print_fluid_sub_grid(gv,0, 0, 0, 127, 127, 127,gv->cube_size);
  // print_fluid_cube(gv, 4,0,0, 5,5,5, gv->cube_size);
  gv->total_threads = total_threads;

   /*lookup_fluid_start_x = 30; lookup_fluid_start_y = gv->jb; lookup_fluid_start_z = 30;   
   lookup_fluid_end_x = 30; lookup_fluid_end_y = gv->je; lookup_fluid_end_z = 30;
   lookup_fiber_start_y = 0; lookup_fiber_start_z = 0; 
   lookup_fiber_end_y = 20;  lookup_fiber_end_z =20;*/


  printf("after generating the fiber shape:\n");
   //print_fiber_sub_grid(gv, 0, 0, 20, 20);   
#if 0
     printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);
#endif   
 
   //may need to move this methods to init
  init_eqlbrmdistrfuncDF0(gv); //changed as per serial code for calculation of eql distr func df0 which is stored as df1 in sequential
  init_df1(gv);   
  init_df_inout(gv);

#if 0 //PASS init
  char filename[80];
  int my_rank = 0;
  sprintf(filename, "Fluid%d_init.dat", my_rank);
  save_fluid_sub_grid(gv, 0, 0, 0, gv->fluid_grid->x_dim - 1, gv->fluid_grid->y_dim - 1, gv->fluid_grid->z_dim - 1, filename);
  printf("Pass save Fluid%d_init.dat\n", my_rank);

  my_rank = 1;
  sprintf(filename, "Fiber%d_init.dat", my_rank);
  save_fiber_sub_grid(gv, 0, 0, gv->fiber_shape->sheets[0].num_rows - 1, gv->fiber_shape->sheets[0].num_cols - 1, filename);
  printf("Pass save Fiber%d_init.dat\n", my_rank);
#endif   

   double startTime = get_cur_time();
   //do_thread(gv);
   pthread_t *threads;
   pthread_attr_t *attrs;
   void           *retval;
   threads  = (pthread_t*) malloc(sizeof(pthread_t)*total_threads);
   attrs = (pthread_attr_t*) malloc(sizeof(pthread_attr_t)*total_threads);
   lvs = (LV)malloc(sizeof(*lvs) * total_threads);
   
   gv->threads = threads; 
   for(i=0; i<total_threads; ++i){
    lvs[i].tid =i;
    lvs[i].gv = gv;
    if(pthread_create(threads +i, NULL, do_thread, lvs+i)){
     perror("Thread creation do_thread\n");
             printf("Error in Thread creation LBM_IB simulation********\n");
             exit(1);
    }
   }
       // pthread_barrier_wait(&(gv->barr));
  for(i = 0; i < total_threads; i++) {
    pthread_join(threads[i], &retval);
    #ifdef DEBUG_PRINT
    printf("Thread %d is finished\n", i);
    #endif //DEBUG_PRINT
  } 
  printf("***********IB Simulation using Pthreads cube v 4.1 ends************\n\n\n");
#if 0  
     printf("Printing for Corner Points(z,y) : 0,0 \n");
     print_fiber_sub_grid(gv,0, 0, 0, 0);
     printf("Printing for Corner Points(z,y) : 51,0 \n");
     print_fiber_sub_grid(gv,0, 51, 0, 51);
     printf("Printing for Corner Points(z,y) : 0,51 \n");
     print_fiber_sub_grid(gv,51, 0, 51, 0);
     printf("Printing for Corner Points (z,y): 51,51 \n");
     print_fiber_sub_grid(gv,51, 51, 51, 51);
#endif
  
   double endTime = get_cur_time();
   free(gv->fluid_grid->sub_fluid_grid);
   free(gv->fluid_grid);
   free(gv);
   pthread_mutex_destroy(&gv->lock_Fluid[total_threads]);
   
   printf("TOTAL TIME TAKEN IN Seconds:%f\n", endTime -startTime);
   return 0; 
}