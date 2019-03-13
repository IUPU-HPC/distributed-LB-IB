#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #define MAX_ERROR 1.1e-3
#define MAX_ERROR 1e-6

/* A node on a fiber */
typedef struct fiber_node_t {
  int i, j;
  double x, y, z;
  double bend_force_x, bend_force_y, bend_force_z;
  double stretch_force_x, stretch_force_y, stretch_force_z;
  double elastic_force_x, elastic_force_y, elastic_force_z; // Bendingforce + Stretchingforce
} Fibernode;

typedef struct fiber_node_t_2 {
  int i, j;
  double x, y, z;
  double bend_force_x, bend_force_y, bend_force_z;
  double stretch_force_x, stretch_force_y, stretch_force_z;
  double elastic_force_x, elastic_force_y, elastic_force_z; // Bendingforce + Stretchingforce
  int ifd_msg_pos[4][4][4];
} Fibernode2;

int main(int argc, char *argv[]){

  int time, TIME_STOP = 200;
  int dump = 10;
  int rank1 = 1; //cube_thread code output
  int rank2 = 1; //mpi code output

  char filename[80];
  char* buffer0 = NULL;
  Fibernode  *buffer1;
  Fibernode2 *buffer2;
  int start_y, start_z, end_y, end_z;
  int num_row, num_col;

  num_row = 52;
  num_col = 52;

  start_y = start_z = 0;
  end_y = num_row - 1;
  end_z = num_col - 1;

  int total_points = num_row * num_col;

#if OG
  buffer0 = (char*) malloc((sizeof(int)*2 + sizeof(double)*3)*total_points);
#endif

  buffer1 = (Fibernode *) malloc(sizeof(Fibernode)*total_points);
  buffer2 = (Fibernode2 *) malloc(sizeof(Fibernode2)*total_points);

  int i, j;
  FILE* iFile = NULL;
  Fibernode *node1;
  Fibernode2 *node2;

  for(time=1; time<=TIME_STOP; time++){

  	if(time % dump == 0){

  	printf("Timestep = %d\n", time);

#ifdef OG
    sprintf(filename, "/home/qoofyk/ScalableIB/ParallelCode/fpL_kbp005_re150_dump_step%d.bin", time);
    iFile = fopen(filename, "rb");
    if (iFile == NULL){
      perror("OG: Failed_open: ");
      return 1;
    }
    fread(buffer0, (sizeof(int)*2 + sizeof(double)*3)*total_points, 1, iFile);
    fclose(iFile);
#endif

  sprintf(filename, "../cube_pthread/Fiber%d_dump_step%d.bin", rank1, time);
	// sprintf(filename, "/pylon5/ac561jp/qoofyk/sharedmem-LBIB/4829945/Fiber%d_dump_step%d.bin", rank1, time);
  // sprintf(filename, "/pylon5/ac561jp/qoofyk/distributed-LB-IB/4830042/Fiber%d_dump_step%d.bin", rank1, time);
	iFile = fopen(filename, "rb");
	if (iFile == NULL){
	  perror("SM: Failed_open: ");
      return 1;
	}
	fread(buffer1, sizeof(Fibernode)*total_points, 1, iFile);
	fclose(iFile);

  sprintf(filename, "../mpi/Fiber%d_dump_step%d.bin", rank1, time);
	// sprintf(filename, "/pylon5/ac561jp/qoofyk/distributed-LB-IB/4840369/Fiber%d_dump_step%d.bin", rank2, time);
	iFile = fopen(filename, "rb");
	if (iFile == NULL){
	  perror("DS: Failed_open: ");
      return 1;
	}
	fread(buffer2, sizeof(Fibernode2)*total_points, 1, iFile);
	fclose(iFile);


	double relative_forward_error_x, relative_forward_error_y, relative_forward_error_z;
  int count = 0;
	for (i = start_y; i <= end_y; ++i) {
    for (j = start_z; j <= end_z; ++j) {

      node1 = buffer1 + i * num_col + j;
      node2 = buffer2 + i * num_col + j;

#if 0
#ifdef OG
    // single core
    printf("SC (%d,%d):(%.24f,%.24f,%.24f)\n",
      ((int*)(buffer0 + count))[0], ((int*)(buffer0 + count))[1],
      ((double *)(buffer0 + count + sizeof(int)*2))[0],
      ((double *)(buffer0 + count + sizeof(int)*2))[1],
      ((double *)(buffer0 + count + sizeof(int)*2))[2]);
    count += sizeof(int)*2 + sizeof(double)*3;
#endif
	  printf("SM (%d,%d):(%.24f,%.24f,%.24f)\n",
	  	node1->i, node1->j,
	  	node1->x, node1->y, node1->z);

	  printf("MPI(%d,%d):(%.24f,%.24f,%.24f)\n",
	  	node2->i, node2->j,
	  	node2->x, node2->y, node2->z);
#endif

	  //perform forward error analysis
#if 1
	  relative_forward_error_x = abs(node2->x - node1->x) / node1->x;
	  relative_forward_error_y = abs(node2->y - node1->y) / node1->y;
	  relative_forward_error_z = abs(node2->z - node1->z) / node1->z;

	  if (relative_forward_error_x > MAX_ERROR || relative_forward_error_y > MAX_ERROR
	  	|| relative_forward_error_z > MAX_ERROR){
	    printf("step=%d, FWE(%d,%d): (%.8f, %.8f, %.8f)\n",
	  	time, i, j,
	  	relative_forward_error_x, relative_forward_error_y, relative_forward_error_z);
	  }

#endif
      }
    }
	}
  }

  return 0;
}
