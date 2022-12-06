from psim import pre_builts


# Simple linear prebuilt construction 
num_cells = 20
t_high = 310
t_low = 290
t_init = 300
# Set t_eq to 0 to run a non deviational simulation 
# Initial t_eq should be in interval [t_low, t_high], otherwise it may take
# a very long time for the simulation to converge
t_eq = 300 
cell_x_len = 50 
cell_y_len = 200 

# Default is a steady-state simulation with fully specular surfaces
b = pre_builts.simple_linear(num_cells, t_high, t_low, t_init, t_eq, cell_x_len, cell_y_len)
b.setSimTime(10) # nanoseconds
b.setMeasurements(1000) # number of measurements over simulation time (50 measurements/ns here)
b.setNumPhonons(5000000) # more phonons -> longer to converge but less noise up to a point

filename = 'linear_demo.json'
json_filepath = f'demo_json/{filename}'
results_filepath = f'demo_results/ss_{filename[:-4]}txt'
b.export(json_filepath)


from psim.matplotlib_plotter import MatplotLibPlotter
mp = MatplotLibPlotter()


# To visualize the system, the C++ code should catch invalid geometries
mp.visualize([json_filepath], square_axes=False)

# Once the JSON file has been processed we can plot the results
# ss_{filename}.txt -> simulation results
mp.plot(filename, results_filepath, square_axes=False) # Not so useful for this geometry

# Plotly linear plotter expects the midpoint of each cell as an array
axis1 = [i for i in range(cell_x_len//2, cell_x_len*num_cells,cell_x_len)] 
axis2 = [i for i in range(0,cell_x_len*num_cells,cell_x_len)]
labels = ["good", "bad"]

from psim.plotly_plotter import PlotlyPlotter
pp = PlotlyPlotter()

pp.plot_linear([axis1, axis2], [results_filepath, results_filepath], main_title='main_title',
                temp_title='temp title here', x_flux_title='x_flux_title_here', y_flux_title='y_flux_title_here',
                legend_labels=labels, normalize_axis=True, outfile='linear_demo')