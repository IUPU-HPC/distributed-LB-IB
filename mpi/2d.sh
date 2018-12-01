timestep=1

#fluid
fluid_grid_z=16
fluid_grid_y=64
fluid_grid_x=64
fluid_task_z=1
fluid_task_y=4
fluid_task_x=4
cube_size=2
thread_per_task_x=2
thread_per_task_y=2
thread_per_task_z=2

#fiber
num_fibersht=1
fiber_width=5
fiber_height=5
fiber_row=13
fiber_clmn=13
fiber_x0=10
fiber_y0=10
fiber_z0=5

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