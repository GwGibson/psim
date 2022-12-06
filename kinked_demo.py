from psim import pre_builts
from psim.matplotlib_plotter import MatplotLibPlotter

mp = MatplotLibPlotter()

# Kinked wires
t_high = 310
t_low = 290
t_init = 300
t_eq = 300
num_layers = 21 # wire is comprised of 21 mini wires (more detail in visualization)
l = 120 # The physical dimensions of the kinked wire are based off l (length)
spec = 1 # perfectly specular

#angles = [i for i in range(0,76,5)] # Can iteratively produce configs
angles = [35] # for demo purposes
for angle in angles:

    b = pre_builts.layered_kinked(
        angle, l, spec, t_high, t_low, t_init, t_eq, num_layers) 
    b.setMeasurements(5000)
    b.setSimTime(10)
    b.setNumPhonons(20000000)
    
    filename = f'demo_json/kinked_demo_{l}_{angle}_spec.json'
    b.export(filename)
    

plot_file = f'kinked_demo_{l}_35_spec.json'
json_filepath = f'demo_json/{plot_file}'
results_filepath = f'demo_results/ss_{plot_file[:-4]}txt'
mp.plot(plot_file, results_filepath, square_axes=False)

# Adjusting the flux scaling
mp.plot(plot_file, results_filepath, square_axes=False,
        fxmax=2.5e9, fxmin=0, fymax=1.45e9, fymin=-1.45e9)
