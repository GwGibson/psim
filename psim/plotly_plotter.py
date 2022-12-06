# -*- coding: utf-8 -*-

# This needs tons of work
from .plotting_tools import parse_ss_data
from .plotting_tools import parse_periodic_data
from .plotting_tools import parse_avg_flux
from .plotting_tools import parse_avg_kink_flux
from plotly.subplots import make_subplots
from plotly.validators.scatter.marker import SymbolValidator
from random import randint

from typing import List
import numpy as np
import plotly.graph_objects as go
import os


class PlotlyPlotter:

    def plot_linear(self, x_axes: List[List[float]], ss_files: List[str], **kwargs):
          
        options = {
            'legend_labels': None,
            'main_title': None,
            'temp_title': None,
            'x_flux_title': None,
            'y_flux_title': None,
			'x_label': None,
            'tlabel': "Temperature (K)",
            'fxlabel': 'X Flux (W/m\u00B2)',
            'fylabel': 'Y Flux (W/m\u00B2)',
            'normalize_axis': False,
            'outfile': None,
            'legend': None
        }
       
        options.update(kwargs)          
    
        """ Plots steady state as a linear plot. The x_axis must be entered
        manually. The points in this list should correspond to the center location
        of each cell."""
		
        # Main figure
        fig = make_subplots(rows=5, cols=2,
                            specs=[[{"rowspan": 3, "colspan": 2}, None],
                                   [None, None],
                                   [None, None],
                                   [{"rowspan": 2, "colspan": 1}, {
                                    "rowspan": 2, "colspan": 1}],
                                   [None, None]],
                            shared_yaxes=False,
                            shared_xaxes=True,
                            vertical_spacing=0.2,
                            subplot_titles=(options['temp_title'], options['x_flux_title'],
                                            options['y_flux_title']))
        t_max = None
        t_min = None
        x_max = x_axes[0][0]
        x_min = x_axes[0][0]
        f_max = 0
        f_min = 0
        x_label = options['x_label']
		
        lines = ('solid', 'dot', 'dash', 'longdash') * 6
        lines = ('solid', 'dot') * 6
        colors = ('rgba(37, 0, 251, 1)', 'rgba(251, 125, 0, 1)', 'rgba(16, 251, 0, 1)', 'rgba(251, 0, 0, 1)') 
        fcs = ('rgba(37, 0, 251, 0.2)', 'rgba(251, 125, 0, 0.2)', 'rgba(16, 251, 0, 0.2)', 'rgba(251, 0, 0, 0.2)') 
        colors = ('rgba(37, 0, 251, 1)', 'rgba(251, 125, 0, 1)') * 6
        fcs = ('rgba(37, 0, 251, 0.2)', 'rgba(251, 125, 0, 0.2)') * 6
        
        for i, (x_axis, file) in enumerate(zip(x_axes, ss_files)):
            data = parse_ss_data(file)
            if len(x_axis) != len(data.temps):
                raise ValueError('The x-axis length ({}) does not match ss data '
                                 '({}) from file: {}.'
                                 .format(len(x_axis), len(data.temps), file))

            if (x_axis[-1] > x_max):
                x_max = x_axis[-1]
            if (x_axis[0] < x_min):
                x_min = x_axis[0]
                
            # Get cell offsets (temp measurement at center of cell)
            start_range = x_min - (x_axis[1] - x_axis[0]) / 2.
            end_range = x_max + (x_axis[-1] - x_axis[-2]) / 2.
            
            if (options['normalize_axis']):             
                x_axis = np.array([i/end_range for i in x_axis])
                start_range = 0
                end_range = 1

            [temps, temps_std] = data.getTempData()
            [xf, xf_std] = data.getXFluxData()
            [yf, yf_std] = data.getYFluxData()

            x_rev = x_axis[::-1]
            [y_upper, y_lower] = self.__get_fill_lines(temps, temps_std)

            color = 'rgba(' + str(randint(0, 255)) + ',' + str(randint(0, 255)) + \
                ',' + str(randint(0, 255)) + ','

            if (t_max == None and t_min == None):
                t_max = max(temps) + 5
                t_min = min(temps) - 5

            if (max(temps) + 5 > t_max):
                t_max = max(temps) + 5
            if (min(temps) - 5 < t_min):
                t_min = min(temps) - 5

            ft_max = max(max(xf), max(yf))
            ft_min = min(min(xf), min(yf))
            if ft_max > f_max:
                f_max = ft_max
            if ft_min < f_min:
                f_min = ft_min
            max_dev = max(max(xf_std), max(yf_std))
            f_max += max_dev
            f_min -= max_dev
            f_min = min(-f_max, f_min)
            f_max = max(-f_min, f_max)
            

            
            name = options['legend_labels'][i] if options['legend_labels'] else os.path.basename(file)
            lc = color + str(1) + ')'
            lc = colors[i]

            ''' Start Plots '''

            # Temperature Plot
            # This is the std filled lines
            fig.add_trace(go.Scatter(
                x=x_axis.tolist()+x_rev.tolist(),
                y=y_upper+y_lower,
                fill='toself',
                fillcolor=fcs[i],
                # fillcolor=color + str(0.2) + ')',
                line_color='rgba(255,255,255,0)',
                name=name,
                legendgroup=f'group{i}',
                showlegend=False,
            ), row=1, col=1)

            # Main Line
            fig.add_trace(go.Scatter(
                x=x_axis, y=temps,
                showlegend=True,
                legendgroup=f'group{i}',
                line_color=lc,
                line=dict(dash=lines[i]),
                name=name
            ), row=1, col=1)

            # # x-flux
            [y_upper, y_lower] = self.__get_fill_lines(xf, xf_std)
            # This is the std filled lines
            fig.add_trace(go.Scatter(
                x=x_axis.tolist()+x_rev.tolist(),
                y=y_upper+y_lower,
                fill='toself',
                fillcolor=fcs[i],
                # fillcolor=color + str(0.2) + ')',
                line_color='rgba(255,255,255,0)',
                name=name,
                legendgroup=f'group{i}',
                showlegend=False,
            ), row=4, col=1)

            # Main Line
            fig.add_trace(go.Scatter(
                x=x_axis, y=xf,
                line_color=lc,
                name=name,
                showlegend=False,
                legendgroup=f'group{i}',
                line=dict(dash=lines[i]),
            ), row=4, col=1)

            fig.update_traces(mode='lines')

            fig.update_yaxes(dict(
                showexponent='last',
                exponentformat='e'
            ), row=4, col=1)


            fig.update_xaxes(title = x_label,
                range=[start_range, end_range], row=4, col=1)
            fig.update_yaxes(title = options['fxlabel'],
                range=[f_min, f_max], row=4, col=1)

            # y-flux
            [y_upper, y_lower] = self.__get_fill_lines(yf, yf_std)
            # This is the std filled lines
            fig.add_trace(go.Scatter(
                x=x_axis.tolist()+x_rev.tolist(),
                y=y_upper+y_lower,
                fill='toself',
                fillcolor=fcs[i],
                # fillcolor=color + str(0.2) + ')',
                line_color='rgba(255,255,255,0)',
                name=name,
                legendgroup=f'group{i}',
                showlegend=False,
            ), row=4, col=2)

            # Main Line
            fig.add_trace(go.Scatter(
                x=x_axis, y=yf,
                line_color=lc,
                name=name,
                showlegend=False,
                legendgroup=f'group{i}',
                line=dict(dash=lines[i]),
            ), row=4, col=2)

            fig.update_traces(mode='lines')

            fig.update_yaxes(dict(
                showexponent='last',
                exponentformat='e'
            ), row=4, col=2)

            fig.update_xaxes(title = x_label,
                range=[start_range, end_range], row=4, col=2)
            fig.update_yaxes(title = options['fylabel'],
                             range=[f_min, f_max], row=4, col=2)

            fig.update_layout(yaxis_range=[t_min, t_max],
                              xaxis_range=[start_range, end_range],
                              xaxis_title = x_label,
                              yaxis_title = options['tlabel'],
                              title_text = options['main_title'],
                              template="simple_white",
                              font=dict(
                                  size=28,
                              )
                        )
            
            fig.update_traces(mode='lines')

        # x = np.array([i for i in range(0,5001,500)])
        # y = np.array([500, 465, 422, 400, 375, 350, 325, 305, 282, 272, 250])
        # new_x = np.linspace(min(x), max(x), 5000)
        # coefs = np.polyfit(x,y,3)
        # new_line = np.polyval(coefs, new_x)
        # print(new_line)
        
        # l = 5001
        # x = np.array([i for i in range(0, l, 10)])
        # y = 310 - 20/l * x  
        # y1 = 302 - 4/l * x
        
        
        
        # fig.add_trace(go.Scatter(
        #     x=x, y=y,
        #     line_color="black",
        #     showlegend=True,
        #     name="Linear",
        # ), row=1, col=1)
        
        # fig.add_trace(go.Scatter(
        #     x=x, y=y1,
        #     line_color="black",
        #     showlegend=True,
        #     name="Linear",
        # ), row=1, col=1)

        # fig.add_trace(go.Scatter(
        #     x=x, y=[149 / 5000e-9 * 4]*len(y),
        #     line_color="black",
        #     showlegend=True,
        #     name="Bulk Flux (X Flux Plot)",
        # ), row=4, col=1)

        # fig.add_trace(go.Scatter(
        #     x=new_x, y = new_line,
        #     line_color="black",
        #     showlegend=True,
        #     name="Heat Diffusion Equation",
        # ), row=1, col=1)     

        # fig.add_trace(go.Scatter(
        #     x=x_axis, y = [5557585859]*len(x_axis),
        #     line_color="black",
        #     showlegend=True,
        #     name="Flux @ 375 K",
        # ), row=4, col=1) 

        # fig.add_trace(go.Scatter(
        #     x=x, y = [17.075]*(len(x_axis)+2),
        #     line_color="black",
        #     showlegend=True,
        #     name="Stefan Boltzmann Law (17.075 K)",
        # ), row=1, col=1) 
        
        # fig.add_trace(go.Scatter(
        #     x=x, y = [36.029]*(len(x_axis)+2),
        #     line_color="black",
        #     showlegend=False,
        #     name="Stefan Boltzmann Law (36.028 K)",
        # ), row=1, col=1) 
        
        # fig.update_traces(mode='lines')
  
        fig.update_layout(
            yaxis = dict(
                showexponent = 'all',
                exponentformat = 'SI'
            )
        )
        if (options['legend']):
            fig.update_layout(legend=options['legend'])  
        fig.update_xaxes(nticks=5)    
        for trace in fig['data']: 
            if(trace['name'] == 'a'): trace['showlegend'] = False
            
        if options['outfile']:
            with open(options['outfile'] + '.html', 'w') as f:
                f.write(fig.to_html(include_plotlyjs='cdn'))
        else:
            fig.show()


    def plot_kink_flux(self, x_axis, ss_files: List[List[str]], outfile,
                       **kwargs):
                
        options = {
            'paths': None,
            'title': None,
            'xlabel': "Kink Angle (Degrees)",
            'ylabel': 'Flux (W/m\u00B2)',
            'legend_labels': None,
            'legend': None
            }
       
        options.update(kwargs)
        symbols = SymbolValidator().values
        legend = True if options['legend_labels'] else False
        
        fig = go.Figure()
        fig.update_layout(title=options['title'], showlegend=legend,
                  yaxis_zeroline=False, xaxis_zeroline=False,
                  template="simple_white")
        fig.update_xaxes(title=options['xlabel'])
        fig.update_yaxes(title=options['ylabel'])
        fig.update_layout(showlegend=True, yaxis={'tickformat': '.1e'})
        
        flux = []
        flux_std = []
        
        for i, filenames in enumerate(ss_files):
            path = options['paths'][i] if options['paths'] else None
            tx = []
            txs = []

            # Should be filename in filenames
            for j, file in enumerate(filenames):
                if (path):
                    file = path + file
                data = parse_avg_kink_flux(file, x_axis[j])
                tx.append(data.x_flux)
                txs.append(data.x_flux_std)

            flux.append(tx)
            flux_std.append(txs)  
            
        
        for i, (x, x_std) in enumerate(zip(flux, flux_std)):
            for j in range(len(x)):
                print(f"{x_axis[j]} {x[j]} {x_std[j]}")
            
            
            label = options['legend_labels'][i] if legend else None
            fig.add_trace(go.Scatter(x=x_axis, y=x, mode='markers', 
                                     name=label, marker = dict(size = 10, symbol = symbols[(i+2)*3]),
                                     customdata = x_std,
                                     hovertemplate='%{y:.1e} \u00B1 %{customdata:.1e} W/m\u00B2',
                                     error_y = dict(
                type='data', 
                array=x_std,
                visible=True,
            ) ))                                     
            fig.update_xaxes(tickvals=x_axis)
                
        if (options['legend']):
            fig.update_layout(legend=options['legend'])    
        fig.show()  
        if outfile:
            with open(outfile + '.html', 'w') as f:
                f.write(fig.to_html(include_plotlyjs='cdn'))

    def plot_avg_flux(self, x_axis: tuple,
                      ss_files: List[List[str]], plot_x: bool,
                      plot_y: bool, outfile: str = None, **kwargs):
        options = {
            'paths': None,
            'title': None,
            'xlabel': None,
            'ylabel': 'Flux (W/m\u00B2)',
            'legend_labels': None,
            'xlens': None,
            'ylens': None,
            'y_format': '.2e',
            'tx_diff': 1,
            'ty_diff': 1,
            'legend': None}

        options.update(kwargs)
        SF = 1e-9
        symbols = SymbolValidator().values
        legend = True if options['legend_labels'] else False
        
        fig = go.Figure()
        fig.update_layout(title=options['title'], showlegend=legend,
                  yaxis_zeroline=True, xaxis_zeroline=True,
                  template="simple_white")
        fig.update_xaxes(title=options['xlabel'])
        fig.update_yaxes(title=options['ylabel'])
        
        fig.update_layout(showlegend=True, yaxis={'tickformat': options['y_format']})
        x_flux = []
        x_flux_std = []
        y_flux = []
        y_flux_std = []

        for i, filenames in enumerate(ss_files):
            path = options['paths'][i] if options['paths'] else None
            xls = options['xlens'][i] if options['xlens'] else None
            yls = options['ylens'][i] if options['ylens'] else None
            tx = []
            txs = []
            ty = []
            tys = []
            # Should be filename in filenames
            for j, file in enumerate(filenames):
                if (path):
                    file = path + file
                data = parse_avg_flux(file)
                tx.append(data.x_flux)
                txs.append(data.x_flux_std)
                if (xls):
                    length = xls[j] * SF
                    tx[-1] *= length / options['tx_diff']
                    txs[-1] *= length / options['tx_diff']
                ty.append(data.y_flux)
                tys.append(data.y_flux_std)
                if (yls):
                    length = yls[j] * SF
                    ty[-1] *= length / options['ty_diff']
                    tys[-1] *= length / options['ty_diff']
            x_flux.append(tx)
            x_flux_std.append(txs)
            y_flux.append(ty)
            y_flux_std.append(tys)
            
        legend_pos = 0
        if (plot_x):
            for i, (x, x_std) in enumerate(zip(x_flux, x_flux_std)):
                label = options['legend_labels'][legend_pos] if legend else None
                fig.add_trace(go.Scatter(x=x_axis, y=x, mode='markers',
                                         name=label, marker = dict(size = 10, symbol = symbols[(i+2)*3]),
                                         error_y = dict(type='data', 
                                                        array=x_std,
                                                        visible=False,
                                                        )                                            
                                         ))
                legend_pos += 1
                fig.update_xaxes(tickvals=x_axis)
                
        if (plot_y):
            for i, (y, y_std) in enumerate(zip(y_flux, y_flux_std)):
                label = options['legend_labels'][legend_pos] if legend else None
                fig.add_trace(go.Scatter(x=x_axis, y=y, mode='markers', 
                                         name=label,
                                         error_y = dict(type='data', 
                                                        array=y_std,
                                                        visible=False)                                            
                                         )) 
                legend_pos += 1
                
                
        fig.add_trace(go.Scatter(x=x_axis, y=[60.31918149003619]*9, #149.42 #102.035
                            mode='lines',
                            line_color='black',
                            name='&#954;<sub>bulk@300</sub>'))
        
        fig.add_trace(go.Scatter(x=x_axis, y=[43.32872330356556]*9, #149.42 #102.035
                    mode='lines',
                    line_color='gray',
                    name='&#954;<sub>bulk@400</sub>'))
        
        fig.update_layout(
            xaxis = dict(
                tickmode = 'array',
                tickvals = [0.1, 0.5, 1, 2, 3, 4, 5, 6],
                )
            )
        
        fig.update_xaxes(range=[0, x_axis[-1]+x_axis[0]/2], tickfont = dict(size = 20))
        fig.update_xaxes(title_font=dict(size=30))
        fig.update_yaxes(range=[0, 100], title_font=dict(size=30))
        fig.update_layout(legend = dict(font = dict(size = 30, color = "black")))
        
        if (options['legend']):
            fig.update_layout(legend=options['legend'])
        if outfile:
            self.__fig_to_file(outfile, fig)
        else:
            fig.show()

    def plot_u_kink(self, filename:str, x_axis: List[float], y_axis: List[float],
                      t_min: float, t_max: float):
        

        def arrange_data(ss_values: List[float]) -> List[List[float]]:
            arranged_values = []   
            for i in range(2):
                row_temps = [None]*(len(x_axis)-1)
                for i in range(3):
                    row_temps[i] = ss_values[i]
                    row_temps[-1-i] = ss_values[-1-i]
                del ss_values[-3:]
                del ss_values[:3]
                arranged_values.append(row_temps) 
            
            for i in range(2):
                row_temps = [None]*(len(x_axis)-1)
                for i in range(2):
                    row_temps[i+1] = ss_values[i]
                    row_temps[-2-i] = ss_values[-1-i]                   
                del ss_values[-2:]
                del ss_values[:2]
                arranged_values.append(row_temps)

            for i in range(2):
                row_temps = [None]*(len(x_axis)-1)
                for i in range(2):
                    row_temps[i+1] = ss_values[i]
                    row_temps[-2-i] = ss_values[-1-i] 
                    row_temps[i+3] = ss_values[i+4]
                del ss_values[-2:]
                del ss_values[:2]
                arranged_values.append(row_temps)                

            return arranged_values                    

        data = parse_ss_data(filename)
        f_min, f_max = self.__get_flux_min_max(data.x_flux, data.y_flux)
                                               
        
        arranged_data = (arrange_data(data.temps), arrange_data(data.temps_std), 
                         arrange_data(data.x_flux), arrange_data(data.x_flux_std),
                         arrange_data(data.y_flux), arrange_data(data.y_flux_std))
        
        
        fig=self.__common_plot(arranged_data, data.title, x_axis, y_axis, t_min, t_max,
                           f_min, f_max) 
        
        self.__fig_to_file(filename, fig)  
   
    def plot_wafer(self, filename, axis, t_min, t_max):
       
        def arrange_data(ss_values: List[float]) -> List[List[float]]:
            arranged_values = []
            row_length = len(axis)-1
            for i in range(row_length):
                row_temps = []
                for i in range(row_length):
                    row_temps.append(ss_values[i])
                del ss_values[:row_length]
                arranged_values.append(row_temps) 
            
            
            return arranged_values                    

        data = parse_ss_data(filename)
        f_min, f_max = self.__get_flux_min_max(data.x_flux, data.y_flux)
        f_min, f_max =(-9e10, 9e10)
        
        arranged_data = (arrange_data(data.temps), arrange_data(data.temps_std), 
                         arrange_data(data.x_flux), arrange_data(data.x_flux_std),
                         arrange_data(data.y_flux), arrange_data(data.y_flux_std))
        
        fig = self.__common_plot(arranged_data, data.title, axis, axis, t_min, t_max,
                           f_min, f_max)
        
        fig.update_layout(
                  font=dict(
                      size=24,
                  )
            )
 

        
        fig.update_layout(
            xaxis = dict(
                tickmode = 'linear',
                dtick = 400e-9
            ),
            yaxis = dict(
                tickmode = 'linear',
                dtick = 400e-9
            )
        )
        
        fig.update_xaxes(dtick=400e-9, row=3, col=1)
        fig.update_xaxes(dtick=400e-9, row=3, col=2)
        fig.update_yaxes(dtick=400e-9, row=3, col=1)
        
        fig.update_annotations(font_size=28)
        
        fig.update_layout(
            font = dict(size=26),
            title={
        'text': "Germanium Deviational Simulation",
        'y':0.99,
        'x':0.45,
        'xanchor': 'center',
        'yanchor': 'top'})

        self.__fig_to_file(filename, fig)

    def animate_linear(self, filename, axis, t_min, t_max, sim_time=0):
        data = parse_periodic_data(filename)
        f_min, f_max = self.__get_flux_min_max(data.x_flux, data.y_flux, True)
        
        n = len(axis)-1
        
        
        
        print(data.temps)
        # self.__animate(filename, axis, axis, t_min, t_max, f_min, f_max, sim_time,
        #        data.title, data.measurement_steps, temp, x_flux, y_flux)
        
        
    def animate_wafer(self, filename, axis, t_min, t_max, sim_time=0):
        data = parse_periodic_data(filename)
        f_min, f_max = self.__get_flux_min_max(data.x_flux, data.y_flux, True)
        
        n = len(axis)-1
        temp = [[lst[i:i+n] for i in range(0, n*n, n)] for lst in data.temps]
        x_flux = [[lst[i:i+n] for i in range(0, n*n, n)] for lst in data.x_flux]
        y_flux = [[lst[i:i+n] for i in range(0, n*n, n)] for lst in data.y_flux]
        
        self.__animate(filename, axis, axis, t_min, t_max, f_min, f_max, sim_time,
                       data.title, data.measurement_steps, temp, x_flux, y_flux)
    
    # TODO: Finish
    def animate_u_kink(self, filename, x_axis, y_axis, t_min, t_max):
        data = parse_periodic_data(filename)
        f_min, f_max = self.__get_flux_min_max(data.x_flux, data.y_flux, True)


    def __animate(self, filename, x_axis, y_axis, t_min, t_max, f_min, f_max, sim_time, 
                  title, measurement_steps, temp, x_flux, y_flux):
        t_fig = self.__animate_data(temp, title+' - [Temperature]', x_axis, y_axis, t_min, 
                                  t_max, measurement_steps, None, sim_time)
        xf_fig = self.__animate_data(x_flux, title+' - [X-Flux]', x_axis, y_axis, f_min, 
                                  f_max, measurement_steps, 'Viridis', sim_time)
        yf_fig = self.__animate_data(y_flux, title+' - [Y-Flux]', x_axis, y_axis, f_min, 
                                  f_max, measurement_steps, 'Viridis', sim_time)        
        
        outfile = os.path.basename(filename)[:-4]
        with open('{}_temp.html'.format(outfile), 'w') as f:
            f.write(t_fig.to_html(include_plotlyjs='cdn'))
        with open('{}_xflux.html'.format(outfile), 'w') as f:
            f.write(xf_fig.to_html(include_plotlyjs='cdn'))
        with open('{}_yflux.html'.format(outfile), 'w') as f:
            f.write(yf_fig.to_html(include_plotlyjs='cdn'))        

    def __animate_data(self, data, title, x_axis, y_axis, d_min, d_max,
                       measurement_steps, scheme, sim_time=0):
        step_time = sim_time / len(measurement_steps) / measurement_steps[1]
        def get_label(frame):
            if sim_time == 0:
                return f'Step: {measurement_steps[frame]}'
            else:
                return f'{round(measurement_steps[frame]*step_time,3)} ns'        
     
        # Create figure
        fig = go.Figure()

        # Add traces, one for each slider step
        for step_data in data:
            fig.add_trace(
                go.Heatmap(
                    visible=False,
                    z=step_data,
                    x=x_axis,
                    y=y_axis,
                    name='',
                    zmin = d_min,
                    zmax = d_max,
                    colorscale=scheme,
                    zsmooth = 'best',
                    hoverongaps = False
                )
            )
    
        # Make 0th trace visible
        fig.data[0].visible = True
        
        # Create and add slider
        steps = []
        for i in range(len(fig.data)):
            step = dict(
                method="update",
                label = get_label(i),
                args=[{"visible": [False] * len(fig.data)}],  # layout attribute
            )
            step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
            steps.append(step)

        sliders_dict = {
            'active': 0,
            'yanchor': 'top',
            'xanchor': 'left',
            'currentvalue': {
                'font': {'size': 20},
                'visible': True,
                'xanchor': 'right'
            },
            'transition': {'duration': 300, 'easing': 'cubic-in-out'},
            'pad': {'b': 10, 't': 50},
            'len': 1,
            'x': 0,
            'y': 0,
            'steps': steps
        }
        
        fig.update_layout(
            sliders=[sliders_dict],
            title=title
        )
    
        return fig 

    def __common_plot(self, arranged_data, title, x_axis, y_axis, t_min, t_max,
                      f_min, f_max):
        temps, std_temps, x_flux, x_std, y_flux, y_std = arranged_data    
        
        fig = make_subplots(rows=3, cols=2, 
            specs=[[{"rowspan": 2, "colspan": 2}, None],
                  [None, None],
                  [{}, {}]],
            shared_yaxes=True,
            shared_xaxes=True,
            vertical_spacing=0.05,
            print_grid=False,
            subplot_titles=("Temperature","X Flux", "Y Flux"),         
            )    
        
        t_hovertemplate = '%{z:.2f} \u00B1 %{customdata:.2f} K'
        f_hovertemplate = '%{z:.1e} \u00B1 %{customdata:.1e} W/m\u00B2'
        
        fig.add_trace(go.Heatmap(
                  z=x_flux,
                  x=x_axis,
                  y=y_axis,
                  zmin = f_min,
                  zmax = f_max,
                  customdata = x_std,
                  name='',
                  colorscale='Viridis',
                  zsmooth = 'best',
                  hoverongaps = False,
                  hovertemplate=f_hovertemplate,
              ), row=3, col=1)
        
        fig.update_traces(showscale=False)
        
        fig.add_trace(go.Heatmap(
                           z=temps,
                           x=x_axis,
                           y=y_axis,
                           zmin = t_min,
                           zmax = t_max,
                           hoverongaps = False,
                           name='',
                           customdata = std_temps,
                           hovertemplate=t_hovertemplate,
                           type = 'heatmap',
                           zsmooth='best',
                           colorbar=dict(
                              len=.65, x= 1.020, y=.67,
                              tickmode="array",
                              tickvals=[t_min, t_max],
                              ticktext=["{} K".format(int(t_min)), "{} K".format(int(t_max))],
                              ticks="outside"       
                           )          
          ), row=1, col=1)
        
         
        fig.add_trace(go.Heatmap(
                  z=y_flux,
                  x=x_axis,
                  y=y_axis,
                  zmin = f_min,
                  zmax = f_max,
                  customdata = y_std,
                  name='',
                  colorscale='Viridis',
                  zsmooth = 'best',
                  hoverongaps = False,
                  hovertemplate=f_hovertemplate,
                  colorbar=dict(
                  len=.3, x=1.020, y=.150,
                  tickmode="array",
                  tickvals=[f_min, f_max],
                  ticktext=["{:0.1e}W/m\u00B2".format(f_min), "{:0.1e}W/m\u00B2".format(f_max)],
                  ticks="outside"       
                  )         
              ), row=3, col=2)
        
        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=False)
        fig.update_layout(title = title)
        
        return fig
      

    def __get_flux_min_max(self, x_flux, y_flux, animate=False):
        if animate:
           f_max = max(np.max(x_flux), np.max(y_flux))
           f_min = min(np.min(x_flux), np.min(y_flux))
        else:
            f_max = max(max(x_flux), max(y_flux))
            f_min = min(min(x_flux), min(y_flux))
       
        return(min(-f_max, f_min), max(-f_min, f_max))
    
    def __fig_to_file(self, filename, fig):
        fig.show()
        outfile = os.path.basename(filename)[:-4]
        with open('{}.html'.format(outfile), 'w') as f:
            f.write(fig.to_html(include_plotlyjs='cdn'))
            
    def __get_fill_lines(self, main_data: List[float], std_data: List[float]):
        y_upper = [sum(i) for i in zip(main_data, std_data)]
        y_lower = [sum(i) for i in zip(main_data, [-i for i in std_data])]
        return (y_upper, y_lower[::-1])