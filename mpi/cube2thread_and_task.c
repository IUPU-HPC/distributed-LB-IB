/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"
// #include <math.h>

/*
BI, BJ, BK: Thread Block index in whole fluid domain
*/
int cube2thread_and_task(int BI, int BJ, int BK, GV gv, int *dst_task){
	int tid, i, j, k;
	int num_cubes_x, num_cubes_y, num_cubes_z, tx, ty, tz;
	int cubes_per_task_x, cubes_per_task_y, cubes_per_task_z;
	int num_fluid_task_x, num_fluid_task_y, num_fluid_task_z;

	num_cubes_x = gv->fluid_grid->num_cubes_x;
	num_cubes_y = gv->fluid_grid->num_cubes_y;
	num_cubes_z = gv->fluid_grid->num_cubes_z;
	tx          = gv->tx; //number of fluid thread per fluid task along x
	ty          = gv->ty;
	tz          = gv->tz;
	num_fluid_task_x = gv->num_fluid_task_x;
	num_fluid_task_y = gv->num_fluid_task_y;
	num_fluid_task_z = gv->num_fluid_task_z;

	cubes_per_task_x = num_cubes_x / num_fluid_task_x; // along x: how many cubes in each task
	cubes_per_task_y = num_cubes_y / num_fluid_task_y; // along y: how many cubes in each task
	cubes_per_task_z = num_cubes_z / num_fluid_task_z; // along y: how many cubes in each task

	int Block_size_x  =  cubes_per_task_x / tx;
 	int Block_size_y  =  cubes_per_task_y / ty;
 	int Block_size_z  =  cubes_per_task_z / tz;

	// cyclic fluid thread division
	// i = (BI % (cubes_per_task_x)) % tx;
	// j = (BJ % (cubes_per_task_y)) % ty;
	// k = (BK % (cubes_per_task_z)) % tz;

	// consecutive fluid thread division
	// i = (BI) >> ((int)(log2(cubes_per_task_x/tx)));
	// j = (BJ) >> ((int)(log2(cubes_per_task_y/ty)));
	// k = (BK) >> ((int)(log2(cubes_per_task_z/tz)));

	i = (BI % cubes_per_task_x) / Block_size_x;
	j = (BJ % cubes_per_task_y) / Block_size_y;
	k = (BK % cubes_per_task_z) / Block_size_z;

	// return thread id
	tid  = i * ty * tz + j * tz + k;

	// task id -- 3D distribution mapping
	// rank_z = BK / cubes_per_task_z;
	// rank_y = BJ / cubes_per_task_y;
	// rank_x = BI / cubes_per_task_x;
	*dst_task = (BI / cubes_per_task_x) * num_fluid_task_y * num_fluid_task_z + (BJ / cubes_per_task_y) * num_fluid_task_z + BK / cubes_per_task_z;

	return tid;
}

int cube2task(int BI, int BJ, int BK, GV gv){
	int num_cubes_x, num_cubes_y, num_cubes_z;
	int cubes_per_task_x, cubes_per_task_y, cubes_per_task_z;
	int num_fluid_task_x, num_fluid_task_y, num_fluid_task_z;

	num_cubes_x = gv->fluid_grid->num_cubes_x;
	num_cubes_y = gv->fluid_grid->num_cubes_y;
	num_cubes_z = gv->fluid_grid->num_cubes_z;

	num_fluid_task_x = gv->num_fluid_task_x;
	num_fluid_task_y = gv->num_fluid_task_y;
	num_fluid_task_z = gv->num_fluid_task_z;

	// printf("num_cubes_x %ld, num_cubes_y %ld, num_cubes_z %ld\n", num_cubes_x, num_cubes_y, num_cubes_z);
	// printf("num_fluid_task_x %d, num_fluid_task_y %d, num_fluid_task_z %d\n", num_fluid_task_x, num_fluid_task_y, num_fluid_task_z);
	// fflush(stdout);

	cubes_per_task_x = num_cubes_x / num_fluid_task_x; // along x: how many cubes in each task
	cubes_per_task_y = num_cubes_y / num_fluid_task_y; // along y: how many cubes in each task
	cubes_per_task_z = num_cubes_z / num_fluid_task_z; // along y: how many cubes in each task

	// printf("cubes_per_task_x %d, cubes_per_task_y %d, cubes_per_task_z %d\n", cubes_per_task_x, cubes_per_task_y, cubes_per_task_z);
	// fflush(stdout);

	// task -- 3D distribution mapping
	// rank_z = BK / cubes_per_task_z;
	// rank_y = BJ / cubes_per_task_y;
	// rank_x = BI / cubes_per_task_x;
	int dst_task = (BI / cubes_per_task_x) * num_fluid_task_y * num_fluid_task_z + (BJ / cubes_per_task_y) * num_fluid_task_z + BK / cubes_per_task_z;

	return dst_task;
}

int global2task(int X, int Y, int Z, GV gv){
	int num_cubes_x, num_cubes_y, num_cubes_z;
	int cubes_per_task_x, cubes_per_task_y, cubes_per_task_z;
	int  num_fluid_task_x, num_fluid_task_y, num_fluid_task_z;
	int BI, BJ, BK;
	int cube_size = gv->cube_size;

	BI = X / cube_size;
    BJ = Y / cube_size;
    BK = Z / cube_size;

	num_cubes_x = gv->fluid_grid->num_cubes_x;
	num_cubes_y = gv->fluid_grid->num_cubes_y;
	num_cubes_z = gv->fluid_grid->num_cubes_z;

	num_fluid_task_x = gv->num_fluid_task_x;
	num_fluid_task_y = gv->num_fluid_task_y;
	num_fluid_task_z = gv->num_fluid_task_z;

	// printf("num_cubes_x %ld, num_cubes_y %ld, num_cubes_z %ld\n", num_cubes_x, num_cubes_y, num_cubes_z);
	// printf("num_fluid_task_x %d, num_fluid_task_y %d, num_fluid_task_z %d\n", num_fluid_task_x, num_fluid_task_y, num_fluid_task_z);
	// fflush(stdout);

	cubes_per_task_x = num_cubes_x / num_fluid_task_x; // along x: how many cubes in each task
	cubes_per_task_y = num_cubes_y / num_fluid_task_y; // along y: how many cubes in each task
	cubes_per_task_z = num_cubes_z / num_fluid_task_z; // along y: how many cubes in each task

	// printf("cubes_per_task_x %d, cubes_per_task_y %d, cubes_per_task_z %d\n", cubes_per_task_x, cubes_per_task_y, cubes_per_task_z);
	// fflush(stdout);

	// task -- 3D distribution mapping
	// rank_z = BK / cubes_per_task_z;
	// rank_y = BJ / cubes_per_task_y;
	// rank_x = BI / cubes_per_task_x;
	int dst_task = (BI / cubes_per_task_x) * num_fluid_task_y * num_fluid_task_z + (BJ / cubes_per_task_y) * num_fluid_task_z + BK / cubes_per_task_z;

	return dst_task;
}