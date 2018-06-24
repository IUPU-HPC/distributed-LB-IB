/*  -- Scalable LB-IB --
    Indiana University Purdue University Indianapolis, USA

    @file

    @date

    @author Yuankun Fu
*/

#include "do_thread.h"

/*
BI, BJ, BK: Thread Block index in whole fluid domain
*/
int cube2thread_and_machine(int BI, int BJ, int BK, GV gv, int *dest_mac){
	int tid, i, j, k;
	int num_cubes_x, num_cubes_y, num_cubes_z, P, Q, R;
	int num_cubes_per_proc_x, num_cubes_per_proc_y, num_cubes_per_proc_z;
	int num_macs, Px, Py, Pz;
	int segment_size_x, segment_size_y, segment_size_z;
	int temp_dest_mac;

	num_cubes_x = gv->num_cubes_x;
	num_cubes_y = gv->num_cubes_y;
	num_cubes_z = gv->num_cubes_z;
	P           = gv->P;						//along x, num of fluid thread per fluid process
	Q           = gv->Q;
	R           = gv->R;
	num_macs = gv->num_macs;
	Px           = gv->Px;
	Py           = gv->Py;
	Pz           = gv->Pz;

	// num_thread_blocks_x  =  P*Px;
	// num_thread_blocks_y  =  Q*Py;
	// num_thread_blocks_z  =  R*Pz;

	num_cubes_per_proc_x = num_cubes_x / Px; // along x: how many cubes in each process
	num_cubes_per_proc_y = num_cubes_y / Py; // along y: how many cubes in each process
	num_cubes_per_proc_z = num_cubes_z / Pz; // along y: how many cubes in each process

	// use % to cyclic to the first thread block
	i = (BI % (num_cubes_per_proc_x)) / P;
	j = (BJ % (num_cubes_per_proc_y)) / Q;
	k = (BK % (num_cubes_per_proc_z)) / R;

	// return thread id
	tid  = i*Q*R + j*R + k;


	// machine or process rank -- 3D distribution mapping
	// rank_z = BK / num_cubes_per_proc_z;
	// rank_y = BJ / num_cubes_per_proc_y;
	// rank_x = BI / num_cubes_per_proc_x;
	*dest_mac = (BI / num_cubes_per_proc_x) * Py * Pz + (BJ / num_cubes_per_proc_y) * Pz + BK / num_cubes_per_proc_z;


	return tid;
}
