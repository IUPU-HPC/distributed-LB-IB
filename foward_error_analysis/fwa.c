#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ERROR 1e-4

/* A node on a fiber */
typedef struct fiber_node_t {
  int i, j;
  double x, y, z;
  double bend_force_x, bend_force_y, bend_force_z;
  double stretch_force_x, stretch_force_y, stretch_force_z;
  double elastic_force_x, elastic_force_y, elastic_force_z; // Bendingforce + Stretchingforce
} Fibernode;

int main(int argc, char *argv[]){

  int time, TIME_STOP=100;
  int dump = 10;

  char filename[80];
  Fibernode  *buffer1, *buffer2;
  int start_y, start_z, end_y, end_z;
  int num_row, num_col;

  num_row = 52; 
  num_col = 52;

  start_y = start_z = 0;
  end_y = num_row - 1;
  end_z = num_col - 1;

  int total_points = num_row * num_col;

  buffer1 = (Fibernode *) malloc(sizeof(Fibernode)*total_points);
  buffer2 = (Fibernode *) malloc(sizeof(Fibernode)*total_points);


  int my_rank, i, j;
  FILE* iFile = NULL;
  Fibernode *node1, *node2;

  for(time=1; time<=TIME_STOP; time++){

  	if(time % dump == 0){

  	printf("Timestep = %d\n", time);

  	my_rank = 1;
	sprintf(filename, "../cube_pthread/Fiber%d_dump_step%d.bin", my_rank, time);
	iFile = fopen(filename, "rb");
	if (iFile == NULL){
	  perror("Failed_open: ");
      return 1;
	}
	fread(buffer1, sizeof(Fibernode)*total_points, 1, iFile);
	fclose(iFile);

	my_rank = 1;
	sprintf(filename, "../mpi/Fiber%d_dump_step%d.bin", my_rank, time);
	iFile = fopen(filename, "rb");
	if (iFile == NULL){
	  perror("Failed_open: ");
      return 1;
	}
	fread(buffer2, sizeof(Fibernode)*total_points, 1, iFile);
	fclose(iFile);

	double relative_forward_error_x, relative_forward_error_y, relative_forward_error_z;
	for (i = start_y; i <= end_y; ++i) {
    for (j = start_z; j <= end_z; ++j) {

      node1 = buffer1 + i * num_col + j;
      node2 = buffer1 + i * num_col + j;

#if 0
	  printf("SM(%d,%d):(%f,%f,%f)\n", 
	  	node1->i, node1->j, 
	  	node1->x, node1->y, node1->z);

	  printf("MPI(%d,%d):(%f,%f,%f)\n", 
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
	    printf("FWA (%d,%d): (%f, %f, %f)\n", 
	  	i, j,
	  	relative_forward_error_x, relative_forward_error_y, relative_forward_error_z);
	  }
	  
#endif	  	
      }
    }
	}
  }

  return 0;  
}