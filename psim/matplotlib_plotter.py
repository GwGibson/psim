# -*- coding: utf-8 -*-

# This needs tons of work
from .plotting_tools import deserialize
from .plotting_tools import parse_ss_data
from .plotting_tools import parse_periodic_data
from .plotting_tools import UpdatablePatchCollection
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, FancyArrow
from matplotlib import animation

from typing import List
import numpy as np
import matplotlib.pyplot as plt
import copy
import sys
        
class MatplotLibPlotter:
    
    def mp_flux(self, json_file, ss_file, outfile):
        min_max = lambda x1, x2, x3: (min(x1,x2,x3), max(x1,x2,x3))
        midpoints = []
        
        md = deserialize(json_file)
        x_max = y_max = sys.float_info.min
        x_min = y_min = sys.float_info.max
        s_id = 0
        for sensor, triangles in md.cell_dict.items():
            next_id = sensor.id
            if (s_id != next_id):
                x_avg = (x_max + x_min)/2
                y_avg = (y_max + y_min)/2
                midpoints.append((x_avg, y_avg))
                s_id = next_id
                x_max = y_max = sys.float_info.min
                x_min = y_min = sys.float_info.max
            
            for t in triangles:
                (xl, xh) = min_max(t.p1.x, t.p2.x, t.p3.x)
                (yl, yh) = min_max(t.p1.y, t.p2.y, t.p3.y)
                
                if (xh > x_max):
                    x_max = xh
                if (yh > y_max):
                    y_max = yh  
                if (xl < x_min):
                    x_min = xl
                if (yl < y_min):
                    y_min = yl  
        
        ss_data = parse_ss_data(ss_file)    
        [temps, xfs, yfs] = ss_data.getPrimaryData()
        with open(outfile, 'w') as f:            
            for (mp, xf, yf) in zip(midpoints, xfs, yfs):
                f.write(f'{mp[0]} {mp[1]} {xf} {yf}\n')
           
    def plot(self, json_file: str, ss_file: str = None, square_axes=True,
             outfile: str = None, flux_scale = 1, **kwargs):
        
        options = {
            'fxmax': None,
            'fxmin': None,
            'fymax': None,
            'fymin': None
        }
        options.update(kwargs)   
        
        md = deserialize(json_file)

        t1 = md.cell_dict[next(iter(md.cell_dict))][0]
        xmax = t1.p1.x
        xmin = t1.p1.x
        ymax = t1.p1.y
        ymin = t1.p1.y

        title = json_file
        patches = []
        tcolors = []
        xcolors = []
        ycolors = []

        if (ss_file != None):  # Use simulation data
            ss_data = parse_ss_data(ss_file)
            title = ss_data.title
            [temps, xfs, yfs] = ss_data.getPrimaryData()
            for temp, x_flux, y_flux, triangles in zip(temps, xfs, yfs, md.cell_dict.values()):
                for t in triangles:
                    tcolors.append(temp)
                    xcolors.append(x_flux/flux_scale)
                    ycolors.append(y_flux/flux_scale)

                    p1 = (t.p1.x, t.p1.y)
                    p2 = (t.p2.x, t.p2.y)
                    p3 = (t.p3.x, t.p3.y)
                    patches.append(Polygon([p1, p2, p3]))

                    xvals = (t.p1.x, t.p2.x, t.p3.x)
                    yvals = (t.p1.y, t.p2.y, t.p3.y)

                    xmax = max(max(xvals), xmax)
                    xmin = min(min(xvals), xmin)
                    ymax = max(max(yvals), ymax)
                    ymin = min(min(yvals), ymin)

        else:  # Depict initial geometry prior to simulation
            for sensor, triangles in md.cell_dict.items():
                for t in triangles:
                    tcolors.append(sensor.t_init)
                    p1 = (t.p1.x, t.p1.y)
                    p2 = (t.p2.x, t.p2.y)
                    p3 = (t.p3.x, t.p3.y)
                    patches.append(Polygon([p1, p2, p3]))

                    xvals = (t.p1.x, t.p2.x, t.p3.x)
                    yvals = (t.p1.y, t.p2.y, t.p3.y)

                    xmax = max(max(xvals), xmax)
                    xmin = min(min(xvals), xmin)
                    ymax = max(max(yvals), ymax)
                    ymin = min(min(yvals), ymin)

            xcolors = [0] * len(tcolors)
            ycolors = [0] * len(tcolors)

        # Common
        xavg = (xmax + xmin) / 2 * 0.1
        yavg = (ymax + ymin) / 2 * 0.1
        xmax += xavg
        xmin -= xavg
        ymax += yavg
        ymin -= yavg
        
        if (square_axes):
            if (xmax > ymax):
                ymax = xmax
            else:
                xmax = ymax
        

        
        
        fig = plt.figure(figsize=(22, 18), dpi=200, constrained_layout=True)
        plt.rc('font', size=16)
        

        tax = fig.add_subplot(2, 1, 1)
        xax = fig.add_subplot(2, 2, 3)
        yax = fig.add_subplot(2, 2, 4)

        fxmin = min(xcolors) 
        fxmax = max(xcolors)
        
        fymin = min(ycolors)
        fymax = max(ycolors)
        
        # fxmin = min(-fxmax, fxmin)
        
        if (flux_scale != 1):
            fxmax = 1
            fxmin = -1
            fymax = 1
            fymin = -1
            
            
        if options['fxmax'] != None:
            fxmax = options['fxmax']
        if options['fxmin'] != None:
            fxmin = options['fxmin']
        if options['fymax'] != None:
            fymax = options['fymax']
        if options['fymin'] != None:
            fymin = options['fymin']
        #print(f'{fmin} {fmax}')
        xp = PatchCollection(patches, cmap='viridis', alpha=.9, antialiased=True,
                             snap=True)
        xp.set_array(np.array(xcolors))
        xp.set_clim([fxmin, fxmax])
        xax.add_collection(xp)

        yp = PatchCollection(patches, cmap='viridis', alpha=.9, antialiased=True,
                             snap=True)
        yp.set_array(np.array(ycolors))
        yp.set_clim([fymin, fymax])
        yax.add_collection(yp)

        for s in md.surfaces:
            dx = s.p2.x - s.p1.x
            dy = s.p2.y - s.p1.y
            patches.append(FancyArrow(s.p1.x, s.p1.y, dx, dy, head_width=0,
                                      width=np.sqrt(dx*dx + dy*dy)*0.10))
            tcolors.append(s.temp)

        tmax = int(max(tcolors))
        tmin = int(min(tcolors))

        tp = PatchCollection(patches, cmap='plasma', alpha=.9, antialiased=True,
                             snap=True)
        tp.set_array(np.array(tcolors))
        tp.set_clim([tmin, tmax])
        tax.add_collection(tp)

        xlabel = 'Length (nm)'
        ylabel = 'Width (nm)'

        fig.colorbar(tp, ax=tax, ticks=[int(x)
                                        for x in np.linspace(tmin, tmax, 3)])
        fig.colorbar(xp, ax=xax, ticks=[int(x)
                                        for x in np.linspace(fxmin, fxmax, 2)])
        fig.colorbar(yp, ax=yax, ticks=[int(x)
                                        for x in np.linspace(fymin, fymax, 2)])
        tax.axes.set_xlim([xmin, xmax])
        tax.axes.set_ylim([ymin, ymax])
        tax.set_title('Temperature (K)')
        tax.set_ylabel(ylabel)
        xax.axes.set_xlim([xmin, xmax])
        xax.axes.set_ylim([ymin, ymax])
        yax.axes.set_xlim([xmin, xmax])
        yax.axes.set_ylim([ymin, ymax])
        if (flux_scale == 1):
            xax.set_title('X-Flux (W/m^2)')
            yax.set_title('Y-Flux (W/m^2)')
        else:
            xax.set_title('X-Flux (Scaled)')
            yax.set_title('Y-Flux (Scaled)')
        xax.set_xlabel(xlabel, fontsize=16)
        xax.set_ylabel(ylabel, fontsize=16)
        yax.set_xlabel(xlabel, fontsize=16)
        fig.suptitle(title, fontweight="bold")
        # fig.suptitle('Full Simulation', fontweight="bold")
        fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0.2,
                                        wspace=0.2)

        if (outfile):
            plt.savefig("{}.png".format(outfile), dpi=200)
        else:
            plt.show()
        plt.close()


    def animate(self, json_file: str, per_file: str, **kwargs):

        
        options = {
            'fxmax': None,
            'fxmin': None,
            'fymax': None,
            'fymin': None
        }

        options.update(kwargs)  

        def scale(mult, temps, xfs, yfs):
            ts = []
            xs = []
            ys = []
            for m, t, x, y in zip(mult, temps, xfs, yfs):
                t = [t]*m
                x = [x]*m
                y = [y]*m
                ts.extend(t)
                xs.extend(x)
                ys.extend(y)
            return (ts, xs, ys)

        md = deserialize(json_file)

        t1 = md.cell_dict[next(iter(md.cell_dict))][0]
        xmax = t1.p1.x
        xmin = t1.p1.x
        ymax = t1.p1.y
        ymin = t1.p1.y

        title = json_file
        patches = []

        per_data = parse_periodic_data(per_file)
        title = per_data.title

        for triangles in md.cell_dict.values():
            for t in triangles:

                p1 = (t.p1.x, t.p1.y)
                p2 = (t.p2.x, t.p2.y)
                p3 = (t.p3.x, t.p3.y)
                patches.append(Polygon([p1, p2, p3]))

                xvals = (t.p1.x, t.p2.x, t.p3.x)
                yvals = (t.p1.y, t.p2.y, t.p3.y)

                xmax = max(max(xvals), xmax)
                xmin = min(min(xvals), xmin)
                ymax = max(max(yvals), ymax)
                ymin = min(min(yvals), ymin)

        xavg = (xmax + xmin) / 2 * 0.1
        yavg = (ymax + ymin) / 2 * 0.1
        xmax += xavg
        xmin -= xavg
        ymax += yavg
        ymin -= yavg
        fig = plt.figure(figsize=(22, 18), dpi=200, constrained_layout=True)
        tax = fig.add_subplot(2, 1, 1)
        xax = fig.add_subplot(2, 2, 3)
        yax = fig.add_subplot(2, 2, 4)

        [temps, xfs, yfs] = per_data.getPrimaryData()
        steps = per_data.measurement_steps
        step_time = md.sim_time / len(steps) / (steps[1] - steps[0]) 

        mult = [len(triangles) for triangles in md.cell_dict.values()]
        [ts, xs, ys] = scale(mult, temps[0], xfs[0], yfs[0])
        fmin = min(min(xfs+yfs, key=lambda x: min(x)))
        fmax = max(max(xfs+yfs, key=lambda x: max(x)))
        fmin = min(-fmax, fmin)

        fxmax = fmax
        fxmin = fmin
        fymax = fmax
        fymin = fmin
               

        if options['fxmax'] != None:
            fxmax = options['fxmax']
        if options['fxmin'] != None:
            fxmin = options['fxmin']
        if options['fymax'] != None:
            fymax = options['fymax']
        if options['fymin'] != None:
            fymin = options['fymin']
            
        xp = PatchCollection(patches, cmap='viridis', alpha=.9, antialiased=True,
                             snap=True)
        xp.set_array(np.array(xs))
        xp.set_clim([fxmin, fxmax])
        xax.add_collection(xp)

        yp = PatchCollection(patches, cmap='viridis', alpha=.9, antialiased=True,
                             snap=True)
        yp.set_array(np.array(ys))
        yp.set_clim([fymin, fymax])
        yax.add_collection(yp)
        
        s_temps = [s.temp for s in md.surfaces]
                
        tmax = int(max(max(temps+[s_temps], key=lambda x: max(x))))
        tmin = int(min(min(temps+[s_temps], key=lambda x: min(x))))
        

        tp = UpdatablePatchCollection(patches, cmap='plasma', alpha=.9, 
                                      antialiased=True, snap=True)
        tp.set_array(np.array(ts))
        tp.set_clim([tmin, tmax])
        tax.add_collection(tp)

        xlabel = 'Length (nm)'
        ylabel = 'Width (nm)'

        fig.colorbar(tp, ax=tax, ticks=[int(x)
                                        for x in np.linspace(tmin, tmax, 2)])
        fig.colorbar(xp, ax=xax, ticks=[int(x)
                                        for x in np.linspace(fxmin, fxmax, 2)])
        fig.colorbar(yp, ax=yax, ticks=[int(x)
                                        for x in np.linspace(fymin, fymax, 2)])
        tax.axes.set_xlim([xmin, xmax])
        tax.axes.set_ylim([ymin, ymax])
        tax.set_title('Temperature (K)')
        tax.set_ylabel(ylabel)
        xax.axes.set_xlim([xmin, xmax])
        xax.axes.set_ylim([ymin, ymax])
        yax.axes.set_xlim([xmin, xmax])
        yax.axes.set_ylim([ymin, ymax])
        xax.set_title('X-Flux (W/m^2)')
        xax.set_xlabel(xlabel)
        xax.set_ylabel(ylabel)
        yax.set_title('Y-Flux (W/m^2)')
        yax.set_xlabel(xlabel)
        fig.suptitle(title, fontweight="bold")
        fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0.2,
                                        wspace=0.2)

        def animate(frame):
            print("frame: {}".format(frame))
            [ts, xs, ys] = scale(mult, temps[frame], xfs[frame], yfs[frame])
            tp.patches = copy.deepcopy(patches)
            # Display emitting surfaces if the are currently emitting phonons
            step = steps[frame]
            cur_time = step * step_time
            for s in md.surfaces:
                if (cur_time < s.start_time or cur_time > s.start_time + s.duration):
                    pass
                else:
                     dx = s.p2.x - s.p1.x
                     dy = s.p2.y - s.p1.y
                     tp.patches.append(FancyArrow(s.p1.x, s.p1.y, dx, dy, head_width=0, 
                                          width=np.sqrt(dx*dx + dy*dy)*0.01))
                     ts.append(s.temp)
            tp.set_array(np.array(ts))
            xp.set_array(np.array(xs))
            yp.set_array(np.array(ys))
            fig.suptitle(title + '\nStep: {}   Time: {} ns'.format(steps[frame],
                         round(cur_time, 3)), fontweight="bold")

        ani = animation.FuncAnimation(
            fig, animate, frames=range(0, len(steps)), blit=False)

        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(
            artist='Graham Gibson', title='Phonons'), bitrate=1800)
        ani.save('{}.mp4'.format(json_file[:-5]), writer=writer)
        

    def visualize(self, geo_files: List[str], square_axes=True, outfile: str=None, 
                  ss_files: List[str]=None, scale_flux=1):
        
        # Find x & y axes limits (max value for constant values)
        def limits(geo_files, ss_files):
            patches = []
            ts = []
            xs = []
            ys = []
            titles = []
            for i, file in enumerate(geo_files):
                ss_data = parse_ss_data(ss_files[i]) if ss_files else None
                titles.append(file)
                file_patches = []
                temps = []
                md = deserialize(file)
                for j, (sensor, triangles) in enumerate(md.cell_dict.items()):
                    for t in triangles:
                        p1 = (t.p1.x, t.p1.y)
                        p2 = (t.p2.x, t.p2.y)
                        p3 = (t.p3.x, t.p3.y)
                        xvals = (t.p1.x, t.p2.x, t.p3.x)
                        yvals = (t.p1.y, t.p2.y, t.p3.y)
                        file_patches.append(Polygon([p1, p2, p3]))
                        if ss_data:
                            temps.append(ss_data.temps[j]) 
                        else:
                            temps.append(sensor.t_init)
                        xs.append(max(xvals))
                        xs.append(min(xvals))
                        ys.append(max(yvals))
                        ys.append(min(yvals))    
                for s in md.surfaces:
                    dx = s.p2.x - s.p1.x
                    dy = s.p2.y - s.p1.y
                    file_patches.append(FancyArrow(s.p1.x, s.p1.y, dx, dy, head_width=0,
                                                  width=np.sqrt(dx*dx + dy*dy)*0.04))
                    temps.append(s.temp)      
                ts.append(temps)
                # print(f'{file} {max(temps)} {min(temps)}')
                patches.append(file_patches)
            return (patches, ts, max(xs), min(xs), max(ys), min(ys), titles)     
        
        def draw_fig(patches, temps, title): 
            tax = fig.add_subplot(1, 1, 1)
            tp = PatchCollection(patches, cmap='plasma', alpha=.9, antialiased=True,
                             snap=True)
            tp.set_array(np.array(temps))
            # tmax = max(temps)
            # tmin = -tmax
            tp.set_clim([tmin, tmax])
            tax.add_collection(tp)
            tp.set_array(np.array(temps))
            tp.set_clim([tmin, tmax])
            tax.add_collection(tp)
            fig.colorbar(tp, ax=tax, ticks=[int(x)
                                            for x in np.linspace(tmin, tmax, 3)])
           
            tax.axes.set_xlim([xmin, xmax])
            tax.axes.set_ylim([ymin, ymax])
            tax.set_title('Temperature (K)')
            tax.set_xlabel(xlabel)
            tax.set_ylabel(ylabel)
            tax.set_aspect('equal')
            fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0.2,
                                            wspace=0.2)
            fig.suptitle(title, fontweight="bold")            
        
      
        def animate(frame):
            fig.clear()
            draw_fig(patches[frame], temps[frame], titles[frame])
        

        (patches, temps, xmax, xmin, ymax, ymin, titles) = limits(geo_files, ss_files)

        tmin = min([min(t) for t in temps])
        tmax = max([max(t) for t in temps])
        # tmin = -tmax

        # Common
        xavg = (xmax + xmin) / 2 * 0.1
        yavg = (ymax + ymin) / 2 * 0.1
        xmax += xavg
        xmin -= xavg
        ymax += yavg
        ymin -= yavg
        
        if (square_axes):
            if (xmax > ymax):
                ymax = xmax
            else:
                xmax = ymax
        
        fig = plt.figure(dpi=1200)
        
        xlabel = 'Length (nm)'
        ylabel = 'Width (nm)'
        frames = [i for i in range(0, len(geo_files)) for _ in range(4)]
        

        
        if len(geo_files) > 1 and outfile:
            ani = animation.FuncAnimation(
                fig, animate, frames=frames, interval=1000, blit=False)
    
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=10, metadata=dict(
                artist='Graham Gibson', title='Phonons'), bitrate=1800)
            ani.save(f'{outfile}.mp4', writer=writer, dpi=600)
        else:
            draw_fig(patches[0], temps[0], titles[0])
            plt.show()

