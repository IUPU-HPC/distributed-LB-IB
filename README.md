# The Parallel LBM-IB Library for Manycore Clusters

## cube_sequential (completed)
This is the sequential LBM-IB code with cube algorithm design.


## cube_pthread (completed)
This is the shared-memory pthread LBM-IB code with cube algorithm design.
This code can run on 1 machine with many cores.
Each core run as 1 pthread.


## mpi_pthread
This is the distributed memory MPI/pthread version of LBM-IB algorithm (ongoing)

To compile the code for OpenMPI run the following before make on BigredII
1)module swap PrgEnv-cray PrgEnv-gnu
2)module load openmpi/ccm/gnu/1.7.2
make