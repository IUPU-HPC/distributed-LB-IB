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

#ifndef DO_THREAD_H
#define DO_THREAD_H

#include <cstdio>
#include <cstdlib>
#include <pthread.h>
#include <cmath>
#include <mpi.h>

/* DF0 === g0== DFEQ = g0  computed every time step, at t==0 DFEQ
*  DF1 == newely computed Distribution in a fluid node
*  DF2 == Destination buffer to store DF1 in neighbour fluid nodes
*  DF0 does not exist
*/

/* MPI
* In move Fiber sheet function, we can use the same influential domain for later versions
* Fiber reads velocity, rho from the fluid points it spreads the force to.
* Partial s1s(vel x, y and z) can be sent from fluid machine to fiber machines and computed in fiber machine
*/

/*
* ASSUMPTIONS: Cube Size >=4
* There are no irregular cubes
* The boundary cubes will be always be first and last cubes for eg BI = BJ = BK = 0
* The number of threads n = P*Q*R where n = 2^i P,Q,R = 2^j fluid_elem_x= fluid_elem_y, fluid_elem_z = 2^k :FOR SIMPLE Distribution
* Using Block Distribution on Threads for both 1D (FIber to thread:fiber2thread ) and 3D (Fluid_Sub_cube to thread :cube2thread)
*/


#define DEBUG_PRINT
#define lookup_fluid_start_x 20
#define lookup_fluid_start_y 11
#define lookup_fluid_start_z 20
#define lookup_fluid_end_x 20
#define lookup_fluid_end_y 11
#define lookup_fluid_end_z 20
//int lookup_fiber_start_y, lookup_fiber_start_z, lookup_fiber_end_y, lookup_fiber_end_z ;//for Debug purpose
//pthread_mutex_t mutex_EF_forFluid;
//pthread_mutex_t mutex_EF_forVel;

// related to direction
#define X_transfer_1D 1
#define Y_transfer_1D 5
#define Z_transfer_1D 3

// just random value
#define X_transfer_2D 1
#define Y_transfer_2D 2
#define Z_transfer_2D 3

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

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
  Fibernode* nodes; // 8 bytes, pointing to an array of fiber nodes
  int num_nodes; // 4 bytes, how many grid points on a fiber
} Fiber;

/* one fiber sheet */
typedef struct fiber_sheet_t {
  Fiber* fibers; // i.e., an array of fibers
  double width, height; // TODO: NOT NEEDED. floating point values, width is the dimension in z, height is the dimension in y.
  long num_cols, num_rows; /* number of fibres along width or column should be 512 for 1024:*/
  /* bottom left corner: located at <min_y, min_z> */
  double x_orig; // starting point for fibersheet x0 = 20
  double y_orig; // starting point for fibersheet y0 = 21.5
  double z_orig; // starting point for fibersheer z0 = 11.5
} Fibersheet; // initial configuration of the sheet is given by users

/* a set of fiber sheets compose a general IB structure*/
typedef struct fiber_shape_t {
  Fibersheet* sheets;
  int num_sheets; //right now, only one rectangle sheet. Will be more general later.
} Fibershape; // IBShape


/* Defining a fluid grid */

/* A single fluid node */
typedef struct fluid_node_t {
  double vel_x, vel_y, vel_z; // u,v and w given initally as 1 and later computed using rho from eqn 12
  double rho, pressure; // add pressure, a scalar given initially as 1.0 and for next time step computed from df2[] using eqn 11
  double df1[19]; // old value, corresponds to 19 different values of distribution function along 19 different discrete velocities ξj given on FIg1 computed from  df0[] using eq 10
  double df2[19]; // new value, after streaming
  double dfeq[19]; /*???*/            // added since getting errroneous value for timestep>1  That Is, df0 //TODO: no needed
  double df_inout[2][19]; /*???*/    // df_inout[0] for inlet and df_inout[1] for outlet    //TODO: Too Expensive!!
  double elastic_force_x, elastic_force_y, elastic_force_z;  // bending + stretching force, which is the little f in equation 19.
} Fluidnode;

/* A 2-D plane (or surface) of fluid nodes */
typedef struct fluid_surface_t {
  Fluidnode* nodes;  // Pointing to a 2-D array of fluid nodes; nodes or 2DarrayOfNodes???
                     // We assume 1st surface is inlet, last one is outlet. Surfaces between is the actual fluid grid.
} Fluidsurface;

/*PTHREAD_Change*/
/* A sub fluid grid adding up together to form the entire fluid grid
The bigger fluid grid composed of smaller cubes of sub fluid grid */
typedef struct sub_fluid_grid_t {
  long sub_coordx_start, sub_coordx_end;  // To identify the cube along X :accessed by FluidGrid->SubFluidGrid->Surfaces->sub_x_dim
  long sub_coordy_start, sub_coordy_end;  // To identify the cube along Y direction (or #columns) accessed by FluidGrid->SubFluidGrid->Surfaces->sub_y_dim
  long sub_coordz_start, sub_coordz_end;  // To identify the cube along Z direction (or #rows)    accessed by FluidGrid->SubFluidGrid->Surfaces->sub_z_dim
  Fluidnode* nodes; // element
  bool isboundary;  // To check if this subgrid is part of inlet, outlet or boundary to handle LBM boundary cases
   //=k input from the user, shluld be dim_x,y , z :Some cube might not be k*K*K
  long grid_dim_x;
  long grid_dim_y;
  long grid_dim_z;
  //int subgridID;//To identify individual subgrid identified by identified by globalI/cube_size, globalJ/cube_size, globalK/cube_size
} Sub_Fluidgrid;
/*PTHREAD_Change*/

/* A fluid grid is simply a stack of fluid surfaces, where surface itself is a 2D matrix */

typedef struct fluid_grid_t {
  long x_dim; //Surface dimension along X direction (How many surfaces a long the X direction)
  long y_dim; //Surface dimension along Y direction (or #columns)
  long z_dim; //Surface dimension along Z direction (or #rows)
  long num_cubes_x, num_cubes_y, num_cubes_z;
  Fluidsurface* inlet;  //with constant values*/
  Fluidsurface* outlet; //with constant values*/
  /*Added following for Pthread version*//*PTHREAD_Change*/
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
  long dt, time, time1, timesteps, TIME_WR, TIME_WR1, N_WR;
  long ib, ie, jb, je, kb, ke; // For Fluid Grid's actual computation part

/*PTHREAD_Change*/
  //thread
  int tx, ty, tz; //number of threads = tx*ty*tz within each fluid task
  int threads_per_task;
  int cube_size; // Tuning factor k having dimension of the subgrid

  pthread_t *threads;
  //pthread_mutex_t lock_Fluid;
  /*Lock for every thread for optimisation*/
  // pthread_mutex_t *lock_Fluid;
  pthread_barrier_t barr; // to put a bariier after every routine and also in some parts of Fiber's force compuatation
/*PTHREAD_Change*/

/*MPI Changes*/
  //task
  int taskid; //my_rank is the rank assigned by COMM_WORLS to diff machines..to be stored in GV.
  int rank[3], size[3];
  int total_tasks;// Total number of tasks
  int num_fluid_task_x, num_fluid_task_y, num_fluid_task_z; //number of fluid task along x, y, z
  int num_fluid_tasks;
  int rankCoord[3]; //x, y, z
  MPI_Comm cartcomm;

  // Fiber <--> Fluid influential domain
  char** ifd_bufpool;          //fiber task pack message in bufferpool, and then send message to fluid tasks 
  long* ifd_bufpool_msg_size;  //each message size
  long* ifd_last_pos;          //Track last position of each fluid process's buffer
  char* ifd_recv_buf;
  int ifd_recv_count;          //Track number of char in received message
  pthread_mutex_t* lock_ifd_fluid_task_msg;
  int* influenced_macs;        //Fluid taskid of influential domain
  int num_influenced_macs;
  int ifd_max_bufsize;

  //streaming
  char* stream_msg[19];
  pthread_mutex_t lock_stream_msg[19];
  int stream_last_pos[19];
  char* stream_recv_buf;
  long stream_recv_max_bufsize;
  // int streamdir[19];
  int streamDest[19]; //used in MPI_SendRecv when streaming
  int streamSrc[19];
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

/* Debug help function */
void print_fiber_sub_grid(GV gv, int start_y, int start_z,
 int end_y,   int end_z);
void print_fluid_cube(GV gv, int BI_start, int BJ_start, int BK_start,
  int BI_end, int BJ_end, int BK_end,
  int cube_size);

/* Debug help function */
/* Requires following 2 more Debug Help Functions
* Cube Distribution Function: Map a cube to Threadid 0 to N-1
* Distribtion is Arbitary Surface or  Block Cyclic
* Distribuion is Block based where every thread is working on a block
* # of threads = P*Q*R--
*/

/* Mapping function*/
/*cube2thread not used in MPI version*/
int cube2thread(int BI, int BJ, int BK, int num_cubes_x, int num_cubes_y, int num_cubes_z, int P, int Q, int R);
/*For MPI*/
int cube2task(long BI, long BJ, long BK, GV gv);
int cube2thread_and_task(int BI, int BJ, int BK, GV gv, int *dst_task);
int global2task(long x, long y, long z, GV gv);

/*Fiber Distribution Function: Map a fiber row to each thread
  fiber_row -ith row*/
int fiber2thread(int fiber_row, int num_fibers, int num_threads);

/* Initiallize the global GV *//*PTHREAD_Change*/
void init_gv_constant(GV gv);
void init_gv(GV gv); //TODO: pass all necessary information from users to initialize GV.

/* Create a fiber shape (user defined arbitray shape)*/
void gen_fiber_sheet(Fibersheet* sheet);

/* Create a 3D fluid grid */ /*PTHREAD_Change*/
void gen_fluid_grid(Fluidgrid *fluid_grid, int cube_size, int taskid, GV gv);


/* Methods for IB computations */
/*Local to fiber machine*/
void compute_bendingforce(LV lv);//eqn 18

/*Local to fiber machine*/
void compute_stretchingforce(LV lv);//eqn 16

/*Local to fiber machine*/
void compute_elasticforce(LV lv);//F for eqn 19

/*One way Message Passing involvng both fluid and fiber machines.
Spreading force from fiber to influenced fluid nodes*/
void fiber_SpreadForce(LV); // eqn 19
void fluid_get_SpreadForce(LV); // eqn 19


/* Methods for LBM computations */
/*No Communication
Finction is used to compute dfeq but the purpose is to initiallize DF1, because in 1st iteration DF1 is empty*/
double computeEquilibrium(int iPop, double rho,
                          double ux, double uy, double uz, double uSqr);
void init_eqlbrmdistrfuncDF0(GV);//eqn13 for df0 or dfeq, df0 does not exist

/*No Communication*/
void init_df_inout(GV);           //used by inlet, outlet, before the simuation starts.

void init_df1(GV);

/*Local to fluid machine*/
void compute_eqlbrmdistrfuncDF1(LV); //eqn10 for df1

/*One way Message Passing from fluid to fluid machine.
*/
void stream_distrfunc(GV, LV);           //eqn10 for df2

/*Local to fluid machine*/
void bounceback_rigidwalls(LV);      //after streaming, special handling of boundary fluid nodes

/*Local to fluid machine*/
void compute_rho_and_u(LV);          //eqn 11 and eqn 12 for rho and velocity

/*Two Way message Passing involving fiber to fluid
1st way fiber to fluid for influential domain
2nd way influenced fluid machine sends sum of velocity (s1, s2, s3)of fluid nodes back to fiber machine*/
void fluid_SpreadVelocity(LV lv);            //eqn 22 to move fiber sheet, also  where eqn 20 is implemented.
void fiber_get_SpreadVelocity(LV lv);

/*Local to fluid machine*/
void copy_inout_to_df2(LV);          //computed in every time step

/*Local to fluid machine*/
void periodicBC(LV);

/*Local to fluid machine*/
void replace_old_DF(LV);// it copies DF2 to DF1

/*TODO: move everything above to a separate header file */

void* do_thread(void* v); //IB-LBM simulation moved to do_thread called via ptheread_create
int check_1d(int x, int y, int z, int* dir);
int check_2d(int x, int y, int z, int* dir);

static inline double get_cur_time();

#endif //DO_THREAD_H