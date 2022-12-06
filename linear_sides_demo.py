from psim import pre_builts
from psim.matplotlib_plotter import MatplotLibPlotter

mp = MatplotLibPlotter()

length = 1000

cell_x_dim = 10
ydim = 100
num_cells = 50
thigh = 330
tlow = 270
tavg = (thigh + tlow) / 2
spec = 1
step_interval = 8
start_time = 0
duration = 0

# 0 -> steady state 1 -> periodic 2 -> transient
sim_types = [0,1,2] 
names = ['ss', 'per', 'trans']
for sim_type, name in zip(sim_types,names):  
    step_interval = 0 if sim_type == 0 else 4
    if (name == 'trans'):
        start_time = 0.1 
        duration = 0.15
    b = pre_builts.simple_linear_sides(length//cell_x_dim, thigh, tlow, tavg, tavg, 
                                       cell_x_dim, ydim, spec, 
                                       sim_type=sim_type, step_interval=step_interval,
                                       start_time=start_time, duration=duration)
    b.setMeasurements(1000)
    b.setSimTime(0.5)
    b.setNumPhonons(10000000)
    
    filename = f'demo_json/linear_sides_demo_{name}.json'
    b.export(filename)     

names = ['ss', 'per', 'trans']
for name in names:
    plot_file = f'linear_sides_demo_{name}.json'
    json_filepath = f'demo_json/{plot_file}'      
    results_filepath = f'demo_results/{name}_{plot_file[:-4]}txt'
    
    if (name != 'ss'):
        results_filepath = f'demo_results/per_{plot_file[:-4]}txt'
        mp.animate(json_filepath, results_filepath)
    else:
        mp.plot(json_filepath, results_filepath, square_axes=False)
    
    
    
    
    