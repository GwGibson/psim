# -*- coding: utf-8 -*-

from typing import List
from enum import Enum
import json
import math
import warnings



class DispersionData:
    """Holds the 3 coefficients needed for quadratic fits of the
    dispersion data for a given material and the maximum frequency from
    all the branches."""
    def __init__(self, la_data: tuple,  max_freq_la: float,
                 ta_data: tuple, max_freq_ta: float):
        """         
        :param la_data: Coefficients for the LA branch dispersion curve 
        starting with the K^2 coefficient.
        ::param max_freq_la: The maximum frequency of the LA branch
        :param ta_data: Coefficients for the TA branch dispersion curve 
        starting with the K^2 coefficient.
        ::param max_freq_ta: The maximum frequency of the TA branch
        """  
        self.la_data = la_data
        self.max_freq_la = max_freq_la
        self.ta_data = ta_data
        self.max_freq_ta = max_freq_ta

    @classmethod
    def from_json(cls, data):
        return cls(**data)


class RelaxationData:
    def __init__(self, b_l: float, b_tn: float, b_tu: float, b_i: float, w: float):
        self.b_l = b_l
        self.b_tn = b_tn
        self.b_tu = b_tu
        self.b_i = b_i
        self.w = w

    @classmethod
    def from_json(cls, data):
        return cls(**data)


class Material:
    def __init__(self, name: str, d_data: DispersionData,
                 r_data: RelaxationData):
        self.name = name
        self.d_data = d_data
        self.r_data = r_data

    def __eq__(self, other):
        return True if self.name == other.name else False

    def __hash__(self):
        return hash(self.name)

    @classmethod
    def from_json(cls, data):
        name = data['name']
        d_data = DispersionData.from_json(data['d_data'])
        r_data = RelaxationData.from_json(data['r_data'])
        return cls(name, d_data, r_data)


class Sensor:
    def __init__(self, id: int, material: str, t_init: float):
        self.id = id
        self.material = material
        self.t_init = t_init

    def __eq__(self, other):
        if self.id == other or self.id == other.id:
            return True
        return False

    def __hash__(self):
        return hash(self.id)

    @classmethod
    def from_json(cls, data):
        return cls(**data)

class Point:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return True if self.x == other.x and self.y == other.y else False

    def __hash__(self):
        return hash((self.x, self.y))

    @classmethod
    def from_json(cls, data):
        return cls(**data)

class Line:
    def __init__(self, p1: Point, p2: Point):
        self.p1 = p1
        self.p2 = p2
        
    def midpoint(self):
        return Point((self.p1.x + self.p2.x)/2, (self.p1.y + self.p2.y)/2)

class Triangle:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        
    def midpoint(self):
        l1 = Line(self.p1, self.p2)
        l1_mid = l1.midpoint()
        
        x = self.p3.x + 2/3*(l1_mid.x - self.p3.x)
        y = self.p3.y + 2/3*(l1_mid.y - self.p3.y)
        
        return Point(x,y)
    

    def __eq__(self, other):
        if (self.p1 == other.p1 and self.p2 == other.p2
                and self.p3 == other.p3):
            return True
        else:
            return False

    def __hash__(self):
        return hash((self.p1, self.p2, self.p3))

    @classmethod
    def from_json(cls, data):
        return Triangle(Point.from_json(data['p1']),
                        Point.from_json(data['p2']),
                        Point.from_json(data['p3']))


class Cell:
    def __init__(self, triangle: Triangle, sensorID: int, specularity: float):
        self.triangle = triangle
        self.sensorID = sensorID
        self.specularity = specularity

    def __eq__(self, other):
        return True if self.triangle == other.triangle else False

    def __hash__(self):
        return hash(self.triangle)

    @classmethod
    def from_json(cls, data):
        triangle = Triangle.from_json(data['triangle'])
        sensorID = data['sensorID']
        specularity = data['specularity']
        return cls(triangle, sensorID, specularity)


class EmitSurface:
    def __init__(self, p1: Point, p2: Point, temp: float, 
                 duration: float, start_time: float):
        self.p1 = p1
        self.p2 = p2
        self.temp = temp
        self.duration = duration
        self.start_time = start_time
        self.length = (p2.y - p1.y)**2 + (p1.x - p2.x)**2

    def __hash__(self):
        return hash((self.p1, self.p2))

    @classmethod
    def from_json(cls, data):
        p1 = Point.from_json(data['p1'])
        p2 = Point.from_json(data['p2'])
        return cls(p1, p2, data['temp'], data['duration'], data['start_time'])


class SimulationSettings:
    def __init__(self, num_measurements: int, sim_time: float, num_phonons: int, 
                 t_eq: float, sim_type: int, step_interval: int, phasor_sim: bool):    
        self.num_measurements = num_measurements
        self.sim_time = sim_time
        self.num_phonons = num_phonons
        self.t_eq = t_eq
        self.sim_type = sim_type
        self.step_interval = step_interval
        self.phasor_sim = phasor_sim

    @classmethod
    def from_json(cls, data):
        return cls(**data)


class Model:
    def __init__(self, settings: SimulationSettings, materials: List[Material],
                 sensors: List[Sensor], cells: list, emit_surfaces: List[EmitSurface]):
        self.settings = settings
        self.materials = materials
        self.sensors = sensors
        self.cells = cells
        self.emit_surfaces = emit_surfaces

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, indent=4)

    @classmethod
    def from_json(cls, data):
        settings = SimulationSettings.from_json(data['settings'])
        materials = list(map(Material.from_json, data['materials']))
        sensors = list(map(Sensor.from_json, data['sensors']))
        cells = list(map(Cell.from_json, data['cells']))
        surfaces = list(map(EmitSurface.from_json, data['emit_surfaces']))
        return cls(settings, materials, sensors, cells, surfaces)


class ModelBuilder:
    """Takes model components and serializes them to JSON format.
    Sensor ID's start at 0 and automatically increment by 1 each time
    a new sensor is added. """

    def __init__(self):
        self.num_measurements = 0
        self.sim_time = 0
        self.num_phonons = 0
        self.t_eq = 0
        self.sim_type = 0
        self.step_interval = 0
        self.phasor_sim = False
        self.materials = []
        self.sensors = []
        self.cells = []
        self.surfaces = []
        
        self.__sensor_id = 0
              
    def setSimType(self, sim_type: int, step_interval: int=0):
        """Set the simulation type. Default is SteadyState simulation.
        Periodic -> use argument 1 or "Periodic"
        Transient -> use argument 2 or "Transient"
        Steady State -> default
        
        :param sim_type: Type of simulation
        :param step_interval: For periodic and transient simulations only. 
         The distance between steps for which measurements are recorded.
        """        
        self.step_interval = step_interval
        if (sim_type == 1 or sim_type == "Periodic"):
            self.sim_type = 1
        elif (sim_type == 2 or sim_type == "Transient"):
            self.sim_type = 2
            
        if (self.sim_type == 0 and step_interval > 0):
            raise ValueError('Steady state simulation should have a '
                             'step interval of 0.')
        if ( (self.sim_type == 1 or self.sim_type == 2) and step_interval == 0):
            raise ValueError('Transient and periodic simulations cannot have '
                             'a step interval of 0.')            

    def setPhasorSim(self):
        """Turn emitting surfaces into phasors -> ejected phonons all have
        the same velocity and direction vector perfectly oriented to the 
        normal vector of the surface. """ 
        
        if (self.sim_type == 0):
            warnings.warn("STEADY_STATE simulation with phasor engaged!")
        self.phasor_sim = True
           
    def setMeasurements(self, num_measurements: int):
        """Set the number of measurements to take throughout the simulation.

        Smaller number of measurements will increase
        performace but decrease accuracy. 1000 measurements seems adequate
        for simulations in the 10 - 50ns length.
        """
        if (num_measurements <= 10):
            raise ValueError('The number of measurements must be a '
                             'positive integer greater that 10.')
        self.num_measurements = num_measurements

    def setSimTime(self, sim_time: float):
        """Set the total simulation time.

        The units are in nanoseconds (1e-9). A sim time of 1 = 1ns.
        """
        if (sim_time < 0.):
            raise ValueError('The simulation time must be positive.')
        self.sim_time = sim_time

    def setNumPhonons(self, num_phonons: int):
        """Set the number of phonons to simulate.

        A higher number of phonons will decrease the effective energy per 
        phonon and give results with less variance but increase runtime. 
        This is similar to the weighting factor used in the M&M and 
        Lacroix papers.
        """
        if (num_phonons < 0):
            raise ValueError('The number of phonons must be positive.')
        self.num_phonons = num_phonons

    def setTeq(self, t_eq: float):
        """Set the equilibrium or linearization temperature of the system.

        This should only be used when the temperature difference of the system
        is small (<50) and the temperatures are high (>200). Setting this to
        0 will run a simulation where every phonon is simulated and not just
        the 'useful' ones.
        """
        if (t_eq < 0.):
            raise ValueError('The equilibrium temperature cannot be below 0.')
        self.t_eq = t_eq
       
    def addRawMaterial(self, name: str, d_data: DispersionData,
                    r_data: RelaxationData):
        """Add a material to the model.

        Dispersion and relaxation data must be manually specified at this time.
        Materials must have different names.
        """
        self.materials.append(Material(name, d_data, r_data))
      
    def addMaterial(self, material: Material):
        'Add a pre build material to the model.'
        self.materials.append(material)       
        

    def addSensor(self, material_name: str, t_init: float) -> int:
        """Add a sensor to the model.

        A sensor if linked to a material via the material name and sensor id's 
        must be unique.
        """
        if (self.__sensor_id < 0):
            raise ValueError('Please use an unsigned int for the sensor ID.')
        if (t_init < 0):
            raise ValueError('The initial temperature cannot be below 0.')
        self.sensors.append(Sensor(self.__sensor_id, material_name, t_init))
        self.__sensor_id += 1
        return self.__sensor_id - 1

    def addTriangleCell(self, p1: tuple, p2: tuple, p3: tuple, sensorID: int,
                        specularity: float):
        """Add a trianglular cell to the model.

        There is only minor error checking at this time. The user should 
        verify the geometry of the model by plotting it via the JSON output. 
        The C++ code will catch invalid geometries. Units are nm.
        """
        if (specularity > 1. or specularity < 0.):
            raise ValueError('Specularity must be in the range [0,1].')
        if (p1 == p2 or p2 == p3 or p3 == p1):
            raise ValueError('Cell points must all be different.')
        if (sensorID > self.__sensor_id):
            raise ValueError(f'Invalid sensorID: {sensorID}.')
            
        self.cells.append(Cell(Triangle(Point(p1[0], p1[1]), Point(p2[0], p2[1]),
                                        Point(p3[0], p3[1])), sensorID, specularity))

    def addRectangularCell(self, p1: tuple, p2: tuple, sensorID: int,
                           specularity: float):
        """Add a rectangular cell to the system.

        The lower left and top right points are represented by p1 and p2.
        Units are nm
        """
               
        self.addTriangleCell(p1, (p1[0], p2[1]), (p2[0], p1[1]), sensorID,
                             specularity)
        self.addTriangleCell(p2, (p2[0], p1[1]), (p1[0], p2[1]), sensorID,
                             specularity)
        
    def addPolygonCell(self, p1: tuple, p2: tuple, p3: tuple, p4: tuple,
                             sensorID: int, specularity: float):
        """Adds a 4 sided polygon.

        Two triangles are created. one from p1, p2, p3 
        & another from p2, p3, p4.
        """
                
        self.addTriangleCell(p1, p2, p3, sensorID, specularity)
        self.addTriangleCell(p2, p3, p4, sensorID, specularity)

    def addQuadrilateralCell(self, p1: tuple, dx: float, dy: float,
                             angle: float, sensorID: int, specularity: float):
        """Add a quadilateral cell to the system - 4 sided polygon.

        Can cause loss of precision
        p1: anchor point
        dx: width
        dy: height
        angle: rotation in degrees counterclockwise about p1
        """
        
        p2 = (p1[0] + dx*math.cos(math.radians(angle)),
              p1[1] + dx*math.sin(math.radians(angle)))
        p3 = (p1[0] + dy*math.cos(math.radians(angle+90)),
              p1[1] + dy*math.sin(math.radians(angle+90)))
        p4 = (p2[0] + p3[0] - p1[0], p2[1] + p3[1] - p1[1])
                
        self.addTriangleCell(p1, p2, p3, sensorID, specularity)
        self.addTriangleCell(p4, p2, p3, sensorID, specularity)

    def addSurface(self, p1: tuple, p2: tuple, temp: float, 
                   duration: float=0., start_time: float=0.):
        """Add an emitting surface to the model.

        Boundary and transition surfaces are automatically created 
        where applicable.
        """
        if (p1 == p2):
            raise ValueError('A surface cannot be formed from identical points.')
        if (temp < 0.):
            raise ValueError('Surface temperature cannot be below 0.')
        if (start_time < 0. or duration < 0.):
            raise ValueError('Start time and/or duration cannot be below 0.')      
        self.surfaces.append(EmitSurface(Point(p1[0], p1[1]),
                                         Point(p2[0], p2[1]), temp, 
                                         duration, start_time))

    def __unique(self, lst: list):
        seen = set()
        return not any(i in seen or seen.add(i) for i in lst)

    def __verifyMaterials(self):
        if not self.__unique(self.materials):
            raise ValueError("Duplicate material name in materials list.")

    def __verifySensors(self):
        if not self.__unique(self.sensors):
            raise ValueError("Duplicate sensorID in sensors list.")

        for sensor in self.sensors:
            s_mat_valid = False
            for material in self.materials:
                if (sensor.material == material.name):
                    s_mat_valid = True
                    break
            if (not s_mat_valid):
                raise ValueError('Material {} from sensor {} does not exist'
                                 ' in the materials list.'
                                 .format(sensor.material, sensor.ID))
            s_id_valid = False
            for cell in self.cells:
                if sensor.id == cell.sensorID:
                    s_id_valid = True
                    break
            if (not s_id_valid):
                raise ValueError('No cells are linked to Sensor {}'
                                 .format(sensor.id))

    def __verifyCells(self):
        if not self.__unique(self.cells):
            raise ValueError("Duplicate cell in cells list.")

        for cell in self.cells:
            valid = False
            for sensor in self.sensors:
                if cell.sensorID == sensor.id:
                    valid = True
                    break
            if (not valid):
                raise ValueError('A cell has a sensorID of {} but no '
                                 'sensor with this ID exists.'
                                 .format(cell.sensorID))

    def __verifySurfaces(self):
        if not self.__unique(self.surfaces):
            raise ValueError("Duplicate surface in surfaces list.")
        for surface in self.surfaces:
            if surface.duration + surface.start_time > self.sim_time:
                raise ValueError("Emit surface duration cannot exceed sim time")    
            if surface.duration == 0.:
                surface.duration = self.sim_time
        # Sort surfaces here to increase performance in c++ code
        self.surfaces.sort(reverse=True, key=lambda x: x.length)
        
    def export(self, filepath: str):
        """Serialize the model to JSON format."""
        if (self.num_measurements == 0):
            raise ValueError('Number of measurements not set.')
        if (self.sim_time == 0.):
            raise ValueError('Simulation time not set.')
        if (self.num_phonons == 0):
            raise ValueError('Number of phonons is not set.')

        self.__verifyMaterials()
        self.__verifySensors()
        self.__verifyCells()
        self.__verifySurfaces()

        ss = SimulationSettings(self.num_measurements, 
                                self.sim_time, self.num_phonons, self.t_eq,
                                self.sim_type, self.step_interval, self.phasor_sim)
        model = Model(ss, self.materials, self.sensors, self.cells,
                      self.surfaces)
        
        with open(filepath, 'w') as f:
            f.write(model.toJSON())

