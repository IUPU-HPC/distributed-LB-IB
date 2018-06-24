/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>
#include <stdbool.h>
#include <mpi.h>

/*DF0 === g0== DFEQ = g0  computed every time step, at t==0 DFEQ
  DF1 == newely computed Distribution in a fluid node
  DF2 == Destination buffer to store DF1 in neighbour fluid nodes
  DF0 does not exist*/

/*MPI

To compile the code for OpenMPI run the following before make on BigredII
1)module swap PrgEnv-cray PrgEnv-gnu
2)module load openmpi/ccm/gnu/1.7.2
make

In move Fiber sheet function we can use the same influential domain for later versions
Fiber reads velocity, rho  from the fluid points it spreads the force to.
Partial s1s(vel x, y and z) can be sent from fluid machine to fiber machines and computed in fiber machine


*/
/*
ASSUMPTIONS: Cube Size >=4
There are no irregular cubes
The boundary cubes will be always be first and last cubes for eg BI=BJ =BK =0

The number of threads n = P*Q*R where n = 2^i P,Q,R = 2^j fluid_elem_x= fluid_elem_y, fluid_elem_z = 2^k :FOR SIMPLE Distribution

Using Block Distribution on Threads for both 1D (FIber to thread:fiber2thread ) and 3D (Fluid_Sub_cube to thread :cube2thread)

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


//------------------------------------------------
/* Defining immersed flexible structures consisting of fibers */
//------------------------------------------------
/* A node on a fiber */
typedef struct fiber_node_t {
  double  x, y, z;
  double  bend_force_x, bend_force_y, bend_force_z;
  double  stretch_force_x, stretch_force_y, stretch_force_z;
  double  elastic_force_x, elastic_force_y, elastic_force_z;// Bendingforce + Stretchingforce
} Fibernode;

/* a single fiber consisting of a number of fiber nodes */
typedef struct fiber_t {
  Fibernode* nodes;/* 8 bytes*/       /*pointing to an array of fiber nodes*/
  int        num_nodes; /* 4 bytes*/  /*how many grid points on a fiber*/
} Fiber;

/* one fiber sheet */
typedef struct fiber_sheet_t {
  Fiber*     fibers;           // i.e., an array of fibers
  double     width,    height; // Todo: NOT NEEDED.floating point values, width is the dimension in z, height is the dimension in y.
  int        num_cols, num_rows;
  /* bottom left corner: located at <min_y, min_z> */
  double     x_orig;           // starting point for fibersheet x0 = 20
  double     y_orig;           // starting point for fibersheet y0 = 21.5
  double     z_orig;           // starting point for fibersheer z0 = 11.5
 } Fibersheet; /* initial configuration of the sheet is given by users!!!*/

/* a set of fiber sheets compose a general IB structure*/
  typedef struct fiber_shape_t {
   Fibersheet *sheets;
   int        num_sheets; //right now, only one rectangle sheet. Will be more general later.
} Fibershape;// -->IBShape


//------------------------------------------------
/* Defining a fluid grid */
//------------------------------------------------
/* A single fluid node */
typedef struct fluid_node_t {
  double vel_x, vel_y, vel_z;       // u,v and w given initally as 1 and later computed using rho from eqn 12
  double rho, pressure;             /*add pressure*/  // a scalar given initially as 1.0 and for next time step computed from df2[] using eqn 11
  double df1[19];                   // old value, corresponds to 19 different values of distribution function along 19 different discrete velocities ξj given on FIg1 computed from  df0[] using eq 10
  double df2[19];                   // new value, after streaming
  double dfeq[19]; /*???*/          // added since getting errroneous value for timestep>1  That Is, df0
  double df_inout[2][19]; /*???*/   // df_inout[0] for inlet and df_inout[1] for outlet    //Todo: Too Expensive!!
  double elastic_force_x, elastic_force_y, elastic_force_z;  // bending +stretching force, which is the little f in equation 19.
} Fluidnode;

/* A 2-D plane (or surface) of fluid nodes */
typedef struct fluid_surface_t {
  Fluidnode* nodes;     //Pointing to a 2-D array of fluid nodes; nodes or 2DarrayOfNodes???
                        //We assume 1st surface is inlet, last one is outlet. Surfaces between is the actual fluid grid.
} Fluidsurface;

/*PTHREAD_Change*/
/* A sub fluid grid adding up together to form the entire fluid grid
The bigger fluid grid composed of smaller cubes of sub fluid grid */
typedef struct sub_fluid_grid_t {
  int            sub_coordx_start, sub_coordx_end;  //To identify the cube across X :accessed by FluidGrid->SubFluidGrid->Surfaces->sub_x_dim
  int            sub_coordy_start, sub_coordy_end;  //To identify the cube along Y direction (or #columns) accessed by FluidGrid->SubFluidGrid->Surfaces->sub_y_dim
  int            sub_coordz_start, sub_coordz_end;  //To identify the cube along Z direction (or #rows)    accessed by FluidGrid->SubFluidGrid->Surfaces->sub_z_dim
  Fluidnode*     nodes;//element
  bool           isboundary; //To check if this subgrid is part of inlet, outlet or boundary to handle LBM boundary cases
   //=k input from the user, shluld be dim_x, y , z :Some cube might not be k*K*K
  int            grid_dim_x;
  int            grid_dim_y;
  int            grid_dim_z;
//  int            subgridID;//To identify individual subgrid identified by identified by globalI/cube_size, globalJ/cube_size, globalK/cube_size
} Sub_Fluidgrid;
/*PTHREAD_Change*/

/* A fluid grid is simply a stack of fluid surfaces, where surface itself is a 2D matrix */

typedef struct fluid_grid_t {
  int            x_dim;    //Surface dimension along X direction (How many surfaces a long the X direction)
  int            y_dim;    //Surface dimension along Y direction (or #columns)
  int            z_dim;    //Surface dimension along Z direction (or #rows)
  Fluidsurface*  inlet;    //with constant values*/
  Fluidsurface*  outlet;   //with constant values*/
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
  double     tau, nu_l, u_l, rho_l, L_l, M_l, g_l, Ks_l, Kb_l, Kst_l;
  double     Re, cs_l, Ma, Kn, Kshat, Ksthat, Kbhat, Mhat, Fr;
  int        dt, time, time1, TIME_STOP, TIME_WR, TIME_WR1, N_WR;
  double     c[19][3];                   //stores 19 different velocity directions for ksi
  int        ib, ie, je, jb, ke, kb;     //For Fluid Grid's actual computation part

  /*PTHREAD_Change*/
  pthread_t *threads;
  //pthread_mutex_t lock_Fluid;
  /*Lock for every thread for optimisation*/
  // pthread_mutex_t *lock_Fluid;

  int total_threads;
  int cube_size; //Tuning factor k having dimension of the subgrid
  int num_cubes_x, num_cubes_y, num_cubes_z;

  int P,Q,R; //number of threads = P*Q*R :: Used for thread distribtuion

  int Px,Py,Pz; //machine dimension along x, y, z

  pthread_barrier_t barr;// to put a bariier after every routine and also in some parts of Fiber's force compuatation
  /*PTHREAD_Change*/
  /*MPI Changes*/
  int my_rank;//my_rank is the rank assigned by COMM_WORLS to diff machines..to be stored in GV.]
  int my_rank_x, my_rank_y, my_rank_z;
  int num_macs;// Total number of processes or machines
  int dest_mac[19];

  //ifd
  char** ifd_bufpool;          //fiber process send bufferpool
  int* ifd_bufpool_msg_size;   //each buffer size
  int* ifd_last_pos;              //Track last position of each fluid process's buffer
  char* ifd_recv_buf;
  int ifd_recv_count;             //Track number of char in received message
  pthread_mutex_t* lock_ifd_fluid_mac_msg;
  int* influenced_macs;         //machines rank of influential domain
  int num_influenced_macs;
  int ifd_max_bufsize;

  //streaming
  char** stream_msg;
  pthread_mutex_t* lock_stream_msg;
  int* stream_last_pos;
  char* stream_recv_buf;
  int stream_recv_max_bufsize;
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
double get_cur_time();
/* Requires following 2 more Debug Help Functions

Cube Distribution Function: Map a cube to Threadid 0 to N-1
Distribtion is Arbitary Surface or  Block Cyclic
Distribuion is Block based where every thread is working on a block
 # of threads = P*Q*R--

 */
/*cube2thread not used in MPI version*/
int cube2thread(int BI, int BJ, int BK, int num_cubes_x, int num_cubes_y, int num_cubes_z, int P, int Q, int R);
/*For MPI*/
int cube2thread_and_machine(int BI, int BJ, int BK, GV gv, int *mac_rank);

/*Fiber Distribution Function: Map a fiber row to each thread
  fiber_row -ith row*/
int fiber2thread(int fiber_row, int num_fibers, int num_threads);

/* Initiallize the global GV *//*PTHREAD_Change*/
void  init_gv(GV, Fibershape*, Fluidgrid*, int cube_size); //TODO: pass all necessary information from users to initialize GV.


/* Create a fiber shape (user defined arbitray shape)*/
Fibershape*  gen_fiber_shape(double width, double height, int num_cols, int num_rowsi, double x0, double y0, double z0, GV gv);

/* Create a 3D fluid grid */ /*PTHREAD_Change*/
Fluidgrid*  gen_fluid_grid(int elem_x, int elem_y, int elem_z, int sub_fluid_grid_dim_k, GV gv);


/* Methods for IB computations */
/*Local to fiber machine*/
void compute_bendingforce(LV);//eqn 18

/*Local to fiber machine*/
void compute_stretchingforce(LV);//eqn 16

/*Local to fiber machine*/
void compute_elasticforce(LV);//F for eqn 19

/*One way Message Passing involvng both fluid and fiber machines.
Spreading force from fiber to influenced fluid nodes*/
void fiber_SpreadForce(LV); // eqn 19
void fluid_get_SpreadForce(LV); // eqn 19


/* Methods for LBM computations */
/*No Communication
Finction is used to compute dfeq but the purpose is to initiallize DF1, because in 1st iteration DF1 is empty*/
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

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
