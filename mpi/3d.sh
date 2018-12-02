timestep=1

#fluid
fluid_grid_z=192
fluid_grid_y=192
fluid_grid_x=192
fluid_task_z=3
fluid_task_y=3
fluid_task_x=3
cube_size=2
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

# #fluid
# fluid_grid_z=24
# fluid_grid_y=24
# fluid_grid_x=24
# fluid_task_z=3
# fluid_task_y=3
# fluid_task_x=3
# cube_size=2
# thread_per_task_x=2
# thread_per_task_y=2
# thread_per_task_z=2

# #fiber
# num_fibersht=1
# fiber_width=5
# fiber_height=5
# fiber_row=6
# fiber_clmn=6
# fiber_x0=5
# fiber_y0=5
# fiber_z0=2.5

num_proc=$(( fluid_task_x * fluid_task_y * fluid_task_z + num_fibersht))

mpirun -np $num_proc ./distributed-lb-ib -steps $timestep \
-fluid_grid_xyz ${fluid_grid_x} ${fluid_grid_y} ${fluid_grid_z} \
-cube_size $cube_size \
-fluid_task_xyz ${fluid_task_x} ${fluid_task_y} ${fluid_task_z} \
-thread_per_task_xyz ${thread_per_task_x} ${thread_per_task_y} ${thread_per_task_z} \
-num_fibersht ${num_fibersht} \
-fibersht_width_height ${fiber_width} ${fiber_height} \
-fibersht_row_clmn ${fiber_row} ${fiber_clmn} \
-fibersht_xyz_0 ${fiber_x0} ${fiber_y0} ${fiber_z0}