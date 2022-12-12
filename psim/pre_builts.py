# -*- coding: utf-8 -*-
import numpy as np
from .builder_tools import ModelBuilder
from .builder_tools import DispersionData
from .builder_tools import RelaxationData
from .builder_tools import Material


# Data from Jean2014 AIP Appendix
# Silicon data - B_i = 1.2e-45
sd_data = DispersionData((-2.22e-7, 9.26e3, 0.), 7.63916048e13,
                        (-2.28e-7, 5.24e3, 0.), 3.0100793072e13)
sr_data = RelaxationData(1.3e-24, 9.0e-13, 1.9e-18, 0, 2.42e13)


# Holland
sd_data_h = DispersionData((-2.01e-7, 9.01e3, 0.), 7.63916048e13,
                        (-2.26e-7, 5.23e3, 0.), 3.0100793072e13)
sr_data_h = RelaxationData(2.0e-24, 9.3e-13, 5.5e-18, 0., 2.417e13)

# Wong 2011
sd_data_w = DispersionData((-2.234e-7, 9.28e3, 0.), 7.63916048e13,
                        (-2.278e-7, 5.24e3, 0.), 3.0100793072e13)
sr_data_w = RelaxationData(2.0e-24, 9.3e-13, 1.7e-18, 0., 2.417e13)

# Germanium data - B_i = 24.0e-45
gd_data = DispersionData((-1.50e-7, 5.63e3, 0.), 4.45236386e13,
                        (-1.13e-7, 2.60e3, 0.), 1.4937724175e13)
gr_data = RelaxationData(2.3e-24, 30.0e-13, 1.5e-18, 0., 1.23e13)

silicon = Material("Silicon", sd_data, sr_data)
germanium = Material("Germanium", gd_data, gr_data) # not used


def simple_linear(num_cells: int, t_high: float, t_low: float, t_init: float, t_eq: float,
                  x_base: float, y_base: float, spec: int=1, sim_type: int=0,
                  step_interval: int=0) -> ModelBuilder:
    b = ModelBuilder() 
    # Specify general model settings
    b.setSimType(sim_type, step_interval)
    b.t_eq = t_eq
    # The model will be comprised of a single material
    b.addMaterial(silicon)
    
    avg_temp = (t_high + t_low) / 2
    # Create the cell component
    for i in range(num_cells):
        # Sensor must be added first
        s_id = b.addSensor(silicon.name, avg_temp)
        lower_left_point = (i*x_base, 0.)
        top_right_point = ((i+1)*x_base, y_base)
        # The cell must be attached to a sensor
        b.addRectangularCell(lower_left_point, top_right_point, s_id, spec)
    # Add left side surface
    b.addSurface((0., 0.), (0., y_base), t_high)
    # Add right side surface
    b.addSurface((num_cells*x_base, 0.), (num_cells*x_base, y_base), t_low)

    return b

def simple_linear_sides(num_cells: int, t_high: float, t_low: float, t_init: float, t_eq: float,
                  x_base: float, y: float, spec: int=1, **kwargs) -> ModelBuilder:
    
    options = {
        'sim_type': 0,
        'step_interval': 0,
        'start_time': 0,
        'duration': 0
    }
    options.update(kwargs)  
    
    start_time = options['start_time']
    duration = options['duration']
            
    b = ModelBuilder() 
    # Specify general model settings
    b.setSimType(options['sim_type'], options['step_interval'])
    b.t_eq = t_eq
    # The model will be comprised of a single material
    b.addMaterial(silicon)
    
    avg_temp = (t_high + t_low) / 2

    system_length = num_cells * x_base    
    surface_length = system_length * 0.1
    num_surfaces = int(surface_length/x_base)
    
    num_cells_y = int(y/x_base)
    
    # Create the cell component
    for i in range(num_cells):
        for j in range(num_cells_y):
            # Sensor must be added first
            s_id = b.addSensor(silicon.name, avg_temp)
            lower_left_point = (i*x_base, j*x_base)
            top_right_point = ((i+1)*x_base, (j+1)*x_base)
            # The cell must be attached to a sensor
            b.addRectangularCell(lower_left_point, top_right_point, s_id, spec)
        
    top_left_x = system_length * 0.3   
    top_right_x = system_length * 0.6
    for i in range(num_surfaces):
        
        if (duration != 0):
            # Add top left side surface
            b.addSurface((top_left_x, y), (top_left_x+x_base, y), t_high, 
                         duration, start_time+duration)
            # Add bot left side surface
            b.addSurface((top_left_x, 0.), (top_left_x+x_base, 0), t_high, 
                         duration, start_time)
        else:
            # Add top left side surface
            b.addSurface((top_left_x, y), (top_left_x+x_base, y), t_high)
            # Add bot left side surface
            b.addSurface((top_left_x, 0.), (top_left_x+x_base, 0), t_high)
        
        top_left_x += x_base
                
        if (options['duration'] != 0):
            # Add top right side surface
            b.addSurface((top_right_x, y), (top_right_x+x_base, y), t_low,
                         duration, start_time)
            # Add bot right side surface
            b.addSurface((top_right_x, 0.), (top_right_x+x_base, 0.), t_low,
                         duration, start_time+duration)
        else:
            # Add top right side surface
            b.addSurface((top_right_x, y), (top_right_x+x_base, y), t_low)
            # Add bot right side surface
            b.addSurface((top_right_x, 0.), (top_right_x+x_base, 0.), t_low)
    
        top_right_x += x_base
        
    for i in range(num_cells_y):          
        # Add left side surface
        b.addSurface((0., i*x_base), (0., (i+1)*x_base), avg_temp)
        # Add right side surface
        b.addSurface((num_cells*x_base, i*x_base), (num_cells*x_base, (i+1)*x_base), avg_temp)
    return b


def u_kink(l: int, t_high: float, t_low: float, t_eq: float, 
           spec: float) -> ModelBuilder:
    t_init = t_eq    
    bx = l/3
    by = bx*2
    b = ModelBuilder()
    b.addMaterial(silicon)
    b.setTeq(t_eq)
    
    def add_cells(num_cells: int, x: float, y: float) -> tuple:
        for i in range(num_cells):
            s_id = b.addSensor(silicon.name, t_init)
            b.addRectangularCell((x,y), (x+bx, y+bx), s_id, spec)
            x += bx
    
    # Horizontal segment
    for i in range(2):
       add_cells(3, 0, i*bx)
 
    # Vertical ascending segment
    for i in range(4):
        add_cells(2, bx, by+(i*bx))
        
    # Top horizontal segment    
    for i in range(2):
        add_cells(2, l, 2*l-((i+1)*bx))
        
    # Vertical descending segment
    for i in range(4):
        add_cells(2, l+by, 2*l-((i+1)*bx))
    
    # Horizontal segment
    for i in range(2):
        add_cells(3, l+by, by-((i+1)*bx))
    
     
    # Add left side surface
    b.addSurface((0., 0.), (0., bx), t_high)
    b.addSurface((0., bx), (0., by), t_high)
    # Add right side surface
    b.addSurface((2*l+by, 0), (2*l+by, bx), t_low)
    b.addSurface((2*l+by, bx), (2*l+by, by), t_low)
    return b    


def wafer2D(cell_dim: float, num_cells: int, t_init: float, spec: float,
            sim_time: float, num_measurements: int, num_phonons: int, t_eq: float,
            interval: int=None, **kwargs) -> ModelBuilder:
    
    options = {
        'left': None,
        'top': None,
        'right': None,
        'bot': None,
        }

    options.update(kwargs)
    
    def get_points(row: int, col: int):
        p1_x = cell_dim * row
        p1_y = cell_dim * col
        p2_x = p1_x + cell_dim
        p2_y = p1_y + cell_dim  
        return ( (p1_x, p1_y), (p2_x, p2_y) )
    
    def add_surface(surface_details: tuple, c1: int, cx: bool):
        if (surface_details):
            temp, (c2, c3), duration, start_time = surface_details
            for i in range(c2, c3, cell_dim):
                if (cx):
                    p1 = (c1, i)
                    p2 = (c1, i + cell_dim)
                else:
                    p1 = (i, c1)
                    p2 = (i+cell_dim, c1)
                if duration:
                    b.addSurface(p1, p2, temp, duration, start_time)
                    b.setSimType(2, interval)
                else:
                    b.addSurface(p1, p2, temp)

    b = ModelBuilder()
    b.addMaterial(silicon)
    b.setSimTime(sim_time)
    b.setMeasurements(num_measurements)
    b.setNumPhonons(num_phonons)
    b.setTeq(t_eq)
    if (interval):
        b.setSimType(1, interval)
        
    # Add sensors and geometric cells
    for col in range(num_cells):
        for row in range(num_cells):
            s_id = b.addSensor(silicon.name, t_init)
            p1, p2 = get_points(row, col)           
            b.addRectangularCell(p1, p2, s_id, spec)
    
    # Add emitting surfaces
    add_surface(options['left'], 0, True)
    add_surface(options['top'], num_cells*cell_dim, False)
    add_surface(options['right'], num_cells*cell_dim, True)
    add_surface(options['bot'], 0, False)
    
    return b


# TODO: clean this up
def layered_kinked(angle: float, l: int, spec: float, t_high: float,
                 t_low: float, t_init: float, t_eq: float, num_rows: int,
                 interval: int=0) -> ModelBuilder:
    
    if (angle < 0 or angle > 90):
        raise ValueError('Angle must be in the range [0,90].')
    
    # trying to get more control over the cell sizes
    # doesn't really work -> should start this from scratch if time permits
    straight_factor = 1 # reduces
    kink_factor = 1
    knee_factor = 1
    top_factor = 1
    
    num_segs = num_rows
    r_angle = np.radians(angle)
    sin = np.sin(r_angle)
    cos = np.cos(r_angle)
    tan = np.tan(r_angle)
    tan2 = np.tan(r_angle/2)
    
    
    bx = l/num_segs
    by = 2 * l / 3
    dy = by/num_rows
    
    kl = 2
    
    top_tri_length = by * tan 
    kink_length = kl*l - by/2*tan - by/2*tan2
    if (top_tri_length/2 - by/2*tan2 > 2*l):
        raise ValueError('Kink angle is not possible under these conditions.')

    b = ModelBuilder()
    material_name = silicon.name
    b.setTeq(t_eq)
    b.full_simulation = False if t_eq > 0 else True
    b.addMaterial(silicon)
    if (interval):
        b.setSimType(1, interval)
    
    x_dist = by * np.tan(r_angle/2)
    hor_seg = num_segs*bx - (x_dist - by/2*np.tan(r_angle/2))
    for row in range(num_rows):
        temp = t_init
        x = 0
        y = dy * row
        x_dist = (by-y) * tan2
        x_dist_top = (by-y-dy) * tan2

        
        # Build first straight segment
        for _ in range(num_segs):
            if (_%straight_factor == 0):
                s_id = b.addSensor(material_name, temp)
            b.addRectangularCell((x, y), (x+hor_seg/num_segs, y+dy), s_id, spec)
            x += hor_seg/num_segs       
        
          # Add first joining segment
        if (angle > 0 and angle < 360):
            xtf = x_dist_top+cos*x_dist_top
            xbf = x_dist+cos*x_dist
            ytf = sin*x_dist_top
            ybf = sin*x_dist
            
            t_dist = (x_dist - x_dist_top)
            s_id = b.addSensor(material_name, temp)

            # Joining Segment
            if (row == num_rows-1):
                b.addTriangleCell((x,y), (x,y+dy), (x+x_dist,y), s_id, spec)  
                x_bot = x+xbf
                y_bot = y+ybf
                #s_id = b.addSensor(material_name, temp)
                b.addTriangleCell((x_bot,y_bot), (x,y+dy), (x+x_dist,y), s_id, spec)
            else:
                segs = (num_rows - row - 1)
                dx = (x_dist_top)/segs
                for s in range(segs):
                    b.addRectangularCell((x+(s*dx),y), (x+(s+1)*dx,y+dy), s_id, spec)
                    if (s%knee_factor == 0):
                        s_id = b.addSensor(material_name, temp)
                         
                b.addTriangleCell((x+x_dist_top,y), (x+x_dist, y), (x+x_dist_top,y+dy), s_id, spec)
                # b.addPolygonCell((x,y), (x,y+dy), (x+x_dist, y), (x+x_dist_top,y+dy), s_id, spec)
                # s_id = b.addSensor(material_name, temp)
                nx = x+x_dist+t_dist*cos
                ny = y+t_dist*sin
                b.addTriangleCell((nx,ny),(x+x_dist, y),  (x+x_dist_top,y+dy), s_id, spec)
                
                for s in range(segs):
                    if (s%knee_factor == 0):
                        s_id = b.addSensor(material_name, temp)
                    b.addQuadrilateralCell((nx+(s*dx*cos),ny+(s*dx*sin)), dx, dy, angle, s_id, spec)
                
                # b.addPolygonCell((x+xbf,y+ybf), (x+xtf,y+dy+ytf), (x+x_dist, y), (x+x_dist_top,y+dy), s_id, spec)
            x = x+xbf
            y = y+ybf
        
        # Kink
        inc_dist = kink_length/num_segs
        for _ in range(num_segs):
            if (_%kink_factor == 0):
                s_id = b.addSensor(material_name, temp)
            b.addQuadrilateralCell(
                (x, y), inc_dist, dy, angle, s_id, spec)
            x += cos*inc_dist
            y += sin*inc_dist                
        
          # Build top joining segment
        if (angle > 0 and angle < 360):
            
            # s_id = b.addSensor(material_name, temp-row)
            if (row == 0):
               s_id = b.addSensor(material_name, temp-row)
               b.addTriangleCell((x,y), (x,y+dy/cos), (x-dy*sin,y+dy*cos), s_id, spec)    
               s_id = b.addSensor(material_name, temp)
               b.addTriangleCell((x,y), (x,y+dy/cos), (x+dy*sin,y+dy*cos), s_id, spec) 
               x += sin*dy
               y += cos*dy
            else:
                xf = dy*row
                yf = sin*xf*tan
                dist = xf*tan/row
                for r in range(row):
                    if (r%top_factor==0):
                        s_id = b.addSensor(material_name, temp) 
                    b.addQuadrilateralCell((x+(r*dist*cos),y+(r*dist*sin)), dist, dy, angle, s_id, spec)
                s_id = b.addSensor(material_name, temp)
                nx = x+xf*sin
                ny = y+yf
                b.addTriangleCell((nx-sin*dy,ny+cos*dy), (nx,ny), (nx, ny+dy/cos), s_id, spec)
                s_id = b.addSensor(material_name, temp)
                b.addTriangleCell((nx+sin*dy,ny+cos*dy), (nx,ny), (nx, ny+dy/cos), s_id, spec)
                for r in range(row):
                    if (r%top_factor==0):
                        s_id = b.addSensor(material_name, temp) 
                    b.addQuadrilateralCell((nx+(r*dist*cos),ny-(r*dist*sin)), dist, dy, 360-angle, s_id, spec)
                # b.addPolygonCell((x,y), (x-sin*dy,y+cos*dy), 
                #                  (x+xf, y+yf), 
                #                  (x+xf, y+yf+dy/cos), s_id, spec)
                # s_id = b.addSensor(material_name, temp-row-2)
                # b.addPolygonCell((x+2*xf,y), (x+2*xf+sin*dy,y+cos*dy), 
                #                   (x+xf, y+yf), 
                #                   (x+xf, y+yf+dy/cos), s_id, spec)
                x = x+2*xf*sin+sin*dy
                y = y+cos*dy
        else:
            x += sin*dy
            y += cos*dy 
        # Build second kink segment
        for _ in range(num_segs):
            if (_%kink_factor == 0):
                s_id = b.addSensor(material_name, temp)
            b.addQuadrilateralCell(
                (x, y), inc_dist, -dy, 360-angle, s_id, spec)
            x += cos*inc_dist
            y -= sin*inc_dist  
            
       # Add second joining segment
        if (angle > 0 and angle < 360):            
              
            if (row == num_rows-1):
                s_id = b.addSensor(material_name, temp)
                b.addTriangleCell((x,y), (x-sin*dy,y-cos*dy), (x-x_dist,y-dy), s_id, spec)  
                x_bot = x+xbf
                y_bot = y-ybf
                #s_id = b.addSensor(material_name, temp)
                b.addTriangleCell((x,y), (x,y-dy), (x-x_dist,y-dy), s_id, spec)
            else:
                segs = (num_rows - row - 1)
                dx = (x_dist_top)/segs
                nx = x-sin*dy
                ny = y-cos*dy
                for s in range(segs):
                    if (s%knee_factor == 0):
                        s_id = b.addSensor(material_name, temp)
                    b.addQuadrilateralCell((nx + (s*dx*cos), ny - (s*dx*sin)), dx, dy, 360-angle, s_id, spec)
                s_id = b.addSensor(material_name, temp)
                nx = x+cos*x_dist_top
                ny = y-ytf
                b.addTriangleCell((nx-sin*dy, ny-cos*dy), (nx, ny), (x+xtf-x_dist,ny-dy), s_id, spec)
                
                # b.addPolygonCell((x,y), (x-sin*dy,y-cos*dy), 
                #                  (x+cos*x_dist_top, y-ytf), (x+xtf-x_dist,y-dy-ytf), s_id, spec)

                # s_id = b.addSensor(material_name, temp)
                b.addTriangleCell((x+cos*x_dist_top, ny-dy), (nx, ny), (x+xtf-x_dist,ny-dy), s_id, spec)
                for s in range(segs):
                    if (s%knee_factor == 0):
                        s_id = b.addSensor(material_name, temp)
                    b.addRectangularCell((nx+(s*dx),ny-dy), (nx+(s+1)*dx,ny), s_id, spec)
                    
                
                # b.addPolygonCell((x+xtf,y-ytf), (x+xtf,y-dy-ytf), 
                #                  (x+cos*x_dist_top, y-ytf), (x+xtf-x_dist,y-dy-ytf), s_id, spec)
            x = x+xbf
            y = y+ybf
        # Build second straight segment
        x = l+by/2*np.tan(r_angle/2)+l*kl*cos*2
        y = dy * row
        for _ in range(num_segs):
            if (_%straight_factor==0):
                s_id = b.addSensor(material_name, temp)
            b.addRectangularCell((x, y), (x+hor_seg/num_segs, y+dy), s_id, spec)
            x += hor_seg/num_segs
            
        b.addSurface((0., y), (0., y+dy), t_high)
        b.addSurface((x, y), (x, y+dy), t_low)
        
        # if row%2==0:    
        #     b.addSurface((0., y), (0., y+dy), t_high)
        #     b.addSurface((x, y), (x, y+dy), t_high)
        # else:
        #     b.addSurface((0., y), (0., y+dy), t_low)
        #     b.addSurface((x, y), (x, y+dy), t_low)
    return b

def waferInternalHeat() -> ModelBuilder:
    b = ModelBuilder();
    
    
    pass
    
