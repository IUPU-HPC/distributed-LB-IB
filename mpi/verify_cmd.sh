timestep=6

#fluid
fluid_grid_z=128
fluid_grid_y=128
fluid_grid_x=128
fluid_task_x=2
fluid_task_y=2
fluid_task_z=2
cube_size=2
thread_per_task_x=2
thread_per_task_y=2
thread_per_task_z=2

#fiber
num_fibersht=1

num_proc=$(( fluid_task_x * fluid_task_y * fluid_task_z + num_fibersht))


mpirun -np $num_proc ./distributed-lb-ib -steps $timestep -fluid_grid_xyz 128 128 128 -cube_size $cube_size -fluid_task_xyz 2 2 2 -thread_per_task_xyz ${thread_per_task_x} ${thread_per_task_y} ${thread_per_task_z} -num_fibersht 1 -num-fibersht_width_height 20 20 -fibersht_row_clmn 52 52 -fibersht_xyz_0 20 21.5 11.5