# -*- coding: utf-8 -*-

from matplotlib.collections import PatchCollection
from .builder_tools import Model
from dataclasses import dataclass

import numpy as np
import statistics
import json


@dataclass
class AvgFluxData:
    title: str
    x_flux: float
    x_flux_std: float
    y_flux: float
    y_flux_std: float


class ParsedData:
    def __init__(self, title: str):
        self.title = title
        self.temps = []
        self.temps_std = []
        self.x_flux = []
        self.x_flux_std = []
        self.y_flux = []
        self.y_flux_std = []
        self.measurement_steps = []  # periodic results only

    def addTemp(self, temp: float, temp_std: float):
        self.temps.append(float(temp))
        self.temps_std.append(float(temp_std))

    def addXFlux(self, x_flux: float, x_flux_std: float):
        self.x_flux.append(float(x_flux))
        self.x_flux_std.append(float(x_flux_std))

    def addYFlux(self, y_flux: float, y_flux_std: float):
        self.y_flux.append(float(y_flux))
        self.y_flux_std.append(float(y_flux_std))

    def getTempData(self):
        return (self.temps, self.temps_std)

    def getXFluxData(self):
        return (self.x_flux, self.x_flux_std)

    def getYFluxData(self):
        return (self.y_flux, self.y_flux_std)

    def getPrimaryData(self):
        return(self.temps, self.x_flux, self.y_flux)


class ModelData:
    def __init__(self, model: Model):
        self.sim_time = model.settings.sim_time
        self.surfaces = model.emit_surfaces
        self.cell_dict = dict(zip(model.sensors,
                                  [[] for i in range(len(model.sensors))]))
        for cell in model.cells:
            self.cell_dict[cell.sensorID].append(cell.triangle)
            
class UpdatablePatchCollection(PatchCollection):
    def __init__(self, patches, *args, **kwargs):
        self.patches = patches
        PatchCollection.__init__(self, patches, *args, **kwargs)

    def get_paths(self):
        self.set_paths(self.patches)
        return self._paths     
    
def deserialize(json_file: str) -> ModelData:
    with open(json_file, 'r') as f:
        model = Model.from_json(json.load(f))
        return ModelData(model)

def parse_ss_data(ss_file: str) -> ParsedData:
    with open(ss_file, 'r') as f:
        ss_data = f.read().splitlines()
        data = ParsedData(ss_data[0])
        ss_data = ss_data[1:]

        for line in ss_data:
            line = line.split()
            data.addTemp(line[0], line[1])
            data.addXFlux(line[2], line[3])
            data.addYFlux(line[4], line[5])

        return data

def parse_avg_flux(ss_file: str) -> AvgFluxData:
    data = parse_ss_data(ss_file)
    x_flux = statistics.mean(data.x_flux)
    x_flux_std = statistics.stdev(data.x_flux)
    y_flux = statistics.mean(data.y_flux)
    y_flux_std = statistics.stdev(data.y_flux)
    return AvgFluxData(data.title, x_flux, x_flux_std, y_flux, y_flux_std)

# Hard coded for specific kink geometries -> should fix this
def parse_avg_kink_flux(ss_file: str, angle: float) -> AvgFluxData:
    
    def remove_transitions(vals: list):
        if (angle > 30):
            del vals[3:5]
            del vals[6:8]
            del vals[9:11]
        elif (angle != 0):
            del vals[3]
            del vals[6]
            del vals[9]
        return vals
    
    def scale_vals(x_vals: list, y_vals: list) -> list:     
        ret_vals = []
        for i, (fx, fy) in enumerate(zip(x_vals, y_vals)):
            r_angle = np.radians(angle)
            if (i < 3 or i > 8):
                ret_vals.append(fx)
            else:
                ret_vals.append(fx*np.cos(r_angle) + np.abs(fy*np.sin(r_angle)))
        return ret_vals


    data = parse_ss_data(ss_file)
    
    t_flux = scale_vals(remove_transitions(data.x_flux), remove_transitions(data.y_flux))
    t_std = scale_vals(remove_transitions(data.x_flux_std), remove_transitions(data.y_flux_std))
    flux = statistics.mean(t_flux)
    flux_std = statistics.mean(t_std)
    
    return AvgFluxData(data.title, flux, flux_std, 0, 0)

def parse_periodic_data(per_file: str) -> ParsedData:
    with open(per_file, 'r') as f:
        per_data = f.read().splitlines()
        data = ParsedData(per_data[0])
        per_data = per_data[1:]

        num_cells = int(per_data[1])
        for i in range(0, len(per_data), num_cells+2):
            step_data = per_data[i:i+num_cells+2]
            data.measurement_steps.append(int(step_data[0]))
            step_data = step_data[2:]
            data.temps.append([float(l[0])
                              for l in [line.split() for line in step_data]])
            data.temps_std = None
            data.x_flux.append([float(l[1])
                               for l in [line.split() for line in step_data]])
            data.x_flux_std = None
            data.y_flux.append([float(l[2])
                               for l in [line.split() for line in step_data]])
            data.y_flux_std = None

    return data    
