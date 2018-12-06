####################################################
# common commands for all experiments

date

#prepare output result directory
PBS_RESULTDIR=${SCRATCH_DIR}/results
# OUTPUT_DIR=${SCRATCH_DIR}/distributed-lb-ib/${directory}
mkdir -pv ${PBS_RESULTDIR}
# lfs setstripe --stripe-size 1m --stripe-count ${tune_stripe_count} ${PBS_RESULTDIR}

#generate hostfile
HOST_DIR=$SCRATCH_DIR/hosts
mkdir -pv $HOST_DIR
rm -f $HOST_DIR/hostfile*
#all tasks run the following command
srun -o $HOST_DIR/hostfile-dup hostname
cat $HOST_DIR/hostfile-dup | sort | uniq | sed "s/$/:${nproc_per_mac}/" >$HOST_DIR/hostfile-all
# # cd ${PBS_RESULTDIR}

#SET TOTAL MPI PROC
export SLURM_NTASKS=$total_proc

#MPI_Init_thread(Multiple level)
export MV2_ENABLE_AFFINITY=0
# export MV2 SHOW CPU BINDING=1

#Turn on Debug Info on Bridges
# export PGI_ACC_NOTIFY=3

# print SLURM Environment variables
echo "SLURM_NNODES=$SLURM_NNODES" #Number of nodes allocated to job
# echo "SLURM_SUBMIT_HOST=$SLURM_SUBMIT_HOST" #Names of nodes allocated to job

# find number of threads for OpenMP
# find number of MPI tasks per node
echo "SLURM_TASKS_PER_NODE=$SLURM_TASKS_PER_NODE"
# TPN=`echo $SLURM_TASKS_PER_NODE | cut -d '(' -f 1`
# echo "TPN=$TPN"
# # find number of CPU cores per node
echo "SLURM_JOB_CPUS_PER_NODE=$SLURM_JOB_CPUS_PER_NODE"
PPN=`echo $SLURM_JOB_CPUS_PER_NODE | cut -d '(' -f 1`
echo "PPN=$PPN"
# THREADS=$(( PPN / TPN ))
# export OMP_NUM_THREADS=$THREADS

# echo "TPN=$TPN, PPN=$PPN, THREADS=$THREADS, lp=$lp, SLURM_NTASKS=$SLURM_NTASKS"

#
export OMP_NUM_THREADS=$(( PPN / nproc_per_mac ))
# export OMP_NUM_THREADS=2
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS, SLURM_NTASKS=$SLURM_NTASKS"

LAUNCHER="mpirun_rsh"
my_run="$LAUNCHER -export -hostfile $HOST_DIR/hostfile-all -np $total_proc $BIN"
echo "MV2_CPU_BINDING_POLICY=$MV2_CPU_BINDING_POLICY"
env|grep '^OMP'
env|grep '^MV2'
# LAUNCHER="mpirun"
# my_run_exp1="$LAUNCHER -np $total_proc $BIN"

###################################################################################

# mpirun -np $total_proc ${BIN} -steps $timestep \
# -fluid_grid_xyz ${fluid_grid_x} ${fluid_grid_y} ${fluid_grid_z} \
# -cube_size $cube_size \
# -fluid_task_xyz ${fluid_task_x} ${fluid_task_y} ${fluid_task_z} \
# -thread_per_task_xyz ${thread_per_task_x} ${thread_per_task_y} ${thread_per_task_z} \
# -num_fibersht ${num_fibersht} \
# -fibersht_width_height ${fiber_width} ${fiber_height} \
# -fibersht_row_clmn ${fiber_row} ${fiber_clmn} \
# -fibersht_xyz_0 ${fiber_x0} ${fiber_y0} ${fiber_z0}

$my_run -steps $timestep \
-fluid_grid_xyz ${fluid_grid_x} ${fluid_grid_y} ${fluid_grid_z} \
-cube_size $cube_size \
-fluid_task_xyz ${fluid_task_x} ${fluid_task_y} ${fluid_task_z} \
-thread_per_task_xyz ${thread_per_task_x} ${thread_per_task_y} ${thread_per_task_z} \
-num_fibersht ${num_fibersht} \
-fibersht_width_height ${fiber_width} ${fiber_height} \
-fibersht_row_clmn ${fiber_row} ${fiber_clmn} \
-fibersht_xyz_0 ${fiber_x0} ${fiber_y0} ${fiber_z0}