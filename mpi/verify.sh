timestep=100
dump=10

# #fluid
# fluid_grid_z=32
# fluid_grid_y=32
# fluid_grid_x=32
# fluid_task_z=1
# fluid_task_y=1
# fluid_task_x=1
# cube_size=4
# thread_per_task_x=2
# thread_per_task_y=2
# thread_per_task_z=2

# #fiber
# num_fibersht=1
# fiber_width=2
# fiber_height=2
# fiber_row=20
# fiber_clmn=20
# fiber_x0=3
# fiber_y0=3
# fiber_z0=3

# #fluid
# fluid_grid_z=64
# fluid_grid_y=64
# fluid_grid_x=64
# fluid_task_y=4
# fluid_task_z=4
# fluid_task_x=4
# cube_size=4
# thread_per_task_x=2
# thread_per_task_y=2
# thread_per_task_z=2

# #fiber
# num_fibersht=1
# fiber_width=20
# fiber_height=20
# fiber_row=52
# fiber_clmn=52
# fiber_x0=20
# fiber_y0=21.5
# fiber_z0=11.5

#fluid
fluid_grid_z=128
fluid_grid_y=128
fluid_grid_x=128
fluid_task_y=1
fluid_task_z=1
fluid_task_x=1
cube_size=4
thread_per_task_x=2
thread_per_task_y=2
thread_per_task_z=2

#fiber
num_fibersht=1
fiber_width=20
fiber_height=20
fiber_row=52
fiber_clmn=52
fiber_x0=20
fiber_y0=21.5
fiber_z0=11.5

num_proc=$(( fluid_task_x * fluid_task_y * fluid_task_z + num_fibersht))

export OMP_NUM_THREADS=28

mpirun -np $num_proc ./distributed-lb-ib -steps $timestep -dump $dump \
-fluid_grid_xyz ${fluid_grid_x} ${fluid_grid_y} ${fluid_grid_z} \
-cube_size $cube_size \
-fluid_task_xyz ${fluid_task_x} ${fluid_task_y} ${fluid_task_z} \
-thread_per_task_xyz ${thread_per_task_x} ${thread_per_task_y} ${thread_per_task_z} \
-num_fibersht ${num_fibersht} \
-fibersht_width_height ${fiber_width} ${fiber_height} \
-fibersht_row_clmn ${fiber_row} ${fiber_clmn} \
-fibersht_xyz_0 ${fiber_x0} ${fiber_y0} ${fiber_z0}
