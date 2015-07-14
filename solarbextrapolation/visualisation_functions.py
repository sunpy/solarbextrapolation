# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 18:45:05 2015

@author: alex_
"""

# Temp

from classes import *
import numpy as np
import os
import sunpy.map as mp
from astropy import units as u
import math

# Visulisation
import yt
import mayavi.sources
from mayavi import mlab
from mayavi_seed_streamlines import SeedStreamline
from mayavi.tools.sources import vector_field
from tvtk.tools import visual


def visualise(aMap3D, **kwargs):
    """
    Basic function for visualising a vector field from an extrapolator.
    """
    # Optional parameters
    boo_debug = kwargs.get('debug', False)
    np_seeds = kwargs.get('seeds', None)
    boundary = kwargs.get('boundary', None)
    qu_origin = kwargs.get('origin', u.Quantity([ 0.0, 0.0, 0.0] * u.meter))
    print "qu_origin: " + str(qu_origin)
    
    mayavi_unit_length = kwargs.get('unit_length', 1.0 * u.Mm)
    units_boundary = kwargs.get('boundary_units', mayavi_unit_length)
    units_volume = kwargs.get('volume_units', mayavi_unit_length)
    show_boundary_axes = kwargs.get('show_boundary_axes', True)
    show_volume_axes = kwargs.get('show_volume_axes', True)
    
    # Setup the arc to length equivilence
    obs_distance = aMap3D.dsun - aMap3D.rsun_meters
    radian_length = [ (u.radian, u.meter, lambda x: obs_distance * x, lambda x: x / obs_distance) ]
    
    # Slice (scale) the fields to make the vectors usable in mayavi.
    int_slice_scale = 1
    npm_3d_sliced   = aMap3D.data[::int_slice_scale,::int_slice_scale,::int_slice_scale,:]
    
    # Plot the main vector field (volume).
    fig = mlab.figure()
    
    # Make 3D coords for ever point in the 3D grid.
    print "xobsrange: " + str(aMap3D.xobsrange.to(u.meter, equivalencies=radian_length))
    x_range = aMap3D.xobsrange.to(u.meter, equivalencies=radian_length) - qu_origin[0]
    print "x_range: " + str(x_range)
    y_range = aMap3D.yobsrange.to(u.meter, equivalencies=radian_length) - qu_origin[1]
    z_range = aMap3D.zrange.to(u.meter, equivalencies=radian_length) - qu_origin[2]
    x_range_scaled = (x_range/mayavi_unit_length).decompose().value
    y_range_scaled = (y_range/mayavi_unit_length).decompose().value
    z_range_scaled = (z_range/mayavi_unit_length).decompose().value
    if boo_debug: print '\n\n' + str(aMap3D.xrange) + '\n' + str(aMap3D.xrange)
    if boo_debug: print str(x_range_scaled) + '\n' + str(x_range_scaled) + '\n\n'
    X, Y, Z = np.mgrid[x_range_scaled[0]:x_range_scaled[1]:npm_3d_sliced.shape[0]*1j,
                       y_range_scaled[0]:y_range_scaled[1]:npm_3d_sliced.shape[1]*1j,
                       z_range_scaled[0]:z_range_scaled[1]:npm_3d_sliced.shape[2]*1j]
    vec_field = vector_field(X, Y, Z, npm_3d_sliced[:,:,:,0], npm_3d_sliced[:,:,:,1], npm_3d_sliced[:,:,:,2],
                             name='Magnetic Vector Field', figure=fig)
    vec_field_mag = mlab.pipeline.extract_vector_norm(vec_field, name="Magnetic Field Magnitude")
    
    
    # Place a small outline around the data cube
    mlab.outline()
    
    if show_volume_axes:
        # Label axes
        axes = mlab.axes()
        x_range_scaled = (x_range/units_volume).decompose()
        y_range_scaled = (y_range/units_volume).decompose()
        z_range_scaled = (z_range/units_volume).decompose()
        axes.axes.ranges = np.array([ x_range_scaled[0],  x_range_scaled[1],  y_range_scaled[0],  y_range_scaled[1],  z_range_scaled[0],  z_range_scaled[1]])
        str_units = ' (' + str(units_volume) + ')'
        axes.axes.x_label = 'Solar X ' + str_units
        axes.axes.y_label = 'Solar Y' + str_units
        axes.axes.z_label = 'Z' + str_units
    
    # Plot streamlines
    if np_seeds is None:
        tup_seed_grid = (5, 5) # default (x, y) divisions
        x_inc = (x_range_scaled[1] - x_range_scaled[0]) / (tup_seed_grid[0] + 1.0)
        y_inc = (y_range_scaled[1] - y_range_scaled[0]) / (tup_seed_grid[1] + 1.0)
        
        # Avoid placing seeds in ghost cells
        z_pos = 1 * x_inc

        # List for storing the seeds as [x,y,z] in MayaVi coordinates.        
        lis_seeds = []
        for i in range(0, tup_seed_grid[0]):
            x_pos = x_range_scaled[0] + (i + 1) * x_inc
            for j in (0,1,2,3,4):#range(0, tup_seed_grid[1]):
                y_pos = y_range_scaled[0] + (j + 1) * y_inc     
                lis_seed = [ x_pos.value, y_pos.value, z_pos.value ]
                lis_seeds.append(lis_seed)
        np_seeds = np.array(lis_seeds)
        """
        tupSeedGrid = (5, 5) # (x, y) divisions, if a file isn't declared.

        tupCentrePixel = [ int(aMap3D.data.shape[0] / 2.0), int(aMap3D.data.shape[0] / 2.0)]
        while math.isnan(aMap3D.data[tupCentrePixel[0],tupCentrePixel[1],intSeedZ, 0]):
            intSeedZ += 1
            if boo_debug: print 'intSeedZ: ' + str(intSeedZ)
        floZ = intSeedZ
        
        intIncX = ((aMap3D.xrange[1] * 1.0 - aMap3D.xrange[0]) / (tupSeedGrid[1] + 1.0))
        intIncY = ((aMap3D.zrange[1] * 1.0 - aMap3D.zrange[0]) / (tupSeedGrid[0] + 1.0))
        lisSeeds = []
        for i in range(1, tupSeedGrid[0] + 1):
            for j in range(1, tupSeedGrid[1] + 1):
                lisSeed = [ i * intIncX + aMap3D.xrange[0], j * intIncY + aMap3D.yrange[0], floZ ]
                lisSeed = [ lisSeed[0].value, lisSeed[1].value, lisSeed[2] ]
                lisSeeds.append(lisSeed)
                if boo_debug: print lisSeed
        """
        np_seeds = np.array(lis_seeds)

    # Plot the seed points
    points = mlab.points3d(np_seeds[:,0], np_seeds[:,1], np_seeds[:,2])
    # Make the points smaller
    points.glyph.glyph.scale_factor = 10.0 #mayavi_scale
    # Make the points blue
    points.actor.property.color = (0.2,0,1)
    # Create the custom streamline object
    slines = SeedStreamline(seed_points=np_seeds)
    # Add the streamline object to the plot and make it use the magentic field data,
    # by adding it as a child of the field we created earlier.
    # We add it to the magnitude field (which is in itself a child of bfield)
    # so that it picks up the scalar values and colours the lines.
    vec_field_mag.add_child(slines)
    
    # Make the lines colour coordinated
    slines.module_manager.scalar_lut_manager.lut_mode = 'winter'    
    slines.stream_tracer.integration_direction = 'both'
    
    
    # Add the boundary data 2D map
    if boundary:
        x_range = boundary.xrange.to(u.meter, equivalencies=radian_length) - qu_origin[0]
        if boo_debug: print '\nboundary: x_range: ' + str(x_range)
        y_range = boundary.yrange.to(u.meter, equivalencies=radian_length) - qu_origin[1]
        x_range_scaled = (x_range/mayavi_unit_length).decompose().value
        if boo_debug: print '\nboundary: x_range_scaled: ' + str(x_range_scaled)
        y_range_scaled = (y_range/mayavi_unit_length).decompose().value
        #x_range_scaled = boundary.xrange.to(x_range/mayavi_unit_length, equivalencies=radian_length).decompose().value
        #y_range_scaled = boundary.yrange.to(y_range/mayavi_unit_length, equivalencies=radian_length).decompose().value
        print '\n\n' + str(aMap2D.xrange) + '\n' + str(mayavi_unit_length) + '\n\n'
        # Create explicit points in 3D space
        X, Y = np.mgrid[x_range_scaled[0]:x_range_scaled[1]:aMap2D.data.shape[0]*1j,
                        y_range_scaled[0]:y_range_scaled[1]:aMap2D.data.shape[1]*1j]
        # Plot and add to the current figure
        img_boundary = mlab.pipeline.array2d_source(X, Y, aMap2D.data, figure=fig)
        img_boundary = mlab.pipeline.image_actor(img_boundary, figure = fig)
        
        # Color the image according to the data ############ Issue HERE
        #mayavi_ct = aMap2D.cmap(range(255))
        #img_boundary.module_manager.scalar_lut_manager.lut.table = mayavi_ct
        
        # Legend details
        img_boundary.module_manager.scalar_lut_manager.show_legend = True #module_manager2.scalar_lut_manager.show_legend = True
        img_boundary.module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.1,  0.1 ])

        # Place a small outline around the data cube
        mlab.outline()
        if show_boundary_axes:
            axes = mlab.axes()
            x_range = boundary.xrange.to(units_boundary.unit, equivalencies=radian_length)
            if boo_debug: print '\nboundary: x_range: ' + str(x_range)
            y_range = boundary.yrange.to(units_boundary.unit, equivalencies=radian_length)
            x_range_scaled = (x_range/units_boundary).decompose().value
            if boo_debug: print '\nboundary: x_range_scaled: ' + str(x_range_scaled)
            y_range_scaled = (y_range/units_boundary).decompose().value
            axes.axes.ranges = np.array([ x_range_scaled[0],  x_range_scaled[1],  y_range_scaled[0],  y_range_scaled[1],  0,  0])
            str_units = ' (' + str(units_boundary) + ')'
            axes.axes.x_label = 'Solar X (' + unit_label(units_boundary) + ')'
            axes.axes.y_label = 'Solar Y (' + unit_label(units_boundary) + ')'
    # And open mayavi
    mlab.show()

def unit_label(quantity):
    if quantity.value == 1.0:
        return str(quantity.unit)
    return str(quantity)

def angle_to_length(map, arc):
    """
    Approximate a length on the surface from the arc length.
    Uses the small angle approximation.
    """    
    length = ((map.dsun - map.rsun_meters) * arc.to(u.radian))
    return length.to(u.m, equivalencies=u.dimensionless_angles())


if __name__ == '__main__':
    # 2D boundary map
    aMap2D = mp.Map('C:/fits/2011-02-14__20-35-25__01_hmi.fits').submap([-10,10] * u.arcsec, [-10,10] * u.arcsec)
    aMap2D_cropped = aMap2D.submap([-8,8] * u.arcsec, [-8,8] * u.arcsec)

    aMap2D = mp.Map('F://Dropbox//Study//2014-2015//QM Project-Dissertation//Python//data_fits//2011-02-14__20-35-25__01_hmi.fits')
    aMap2D_cropped = aMap2D.submap([50,300] * u.arcsec, [-350,-100] * u.arcsec)
    aMap2D = aMap2D.submap([0,350] * u.arcsec, [-400,-50] * u.arcsec)
    
    # Using the dimensions of the boundary data to get the volume size
    dim = [ aMap2D_cropped.data.shape[0], aMap2D_cropped.data.shape[1] ]
    
    # Observational ranges are from the boundary data
    x_obs_range = aMap2D_cropped.xrange
    y_obs_range = aMap2D_cropped.yrange
    
    # Setup the arc to length equivilence
    obs_distance = aMap2D_cropped.dsun - aMap2D_cropped.rsun_meters
    radian_length = [ (u.radian, u.meter, lambda x: obs_distance * x, lambda x: x / obs_distance) ]
    
    # MayaVi scale and origin details
    scale_length = 1.0 * u.Mm
    origin_x = aMap2D_cropped.xrange[0].to(u.meter, equivalencies=radian_length)
    origin_y = aMap2D_cropped.yrange[0].to(u.meter, equivalencies=radian_length)
    origin = u.Quantity([ origin_x, origin_y, 0.0 * u.meter ])
    #print "origin: " + str(origin)

    
    # Getting the ranges. x/y from the boundary, z is defined manually.    
    
    #x_range = scaled_length(x_obs_range, scale_length, aMap2D_cropped) * scale_length.unit
    #y_range = scaled_length(y_obs_range, scale_length, aMap2D_cropped) * scale_length.unit
    #z_range = [0.0, 10.0] * u.Mm
    #x_range=[0.0, 20.0]*u.Mm
    #y_range=[0.0, 20.0]*u.Mm
    #z_range=[0.0, 200.0]*u.Mm

    x_range = aMap2D_cropped.xrange.to(u.meter, equivalencies=radian_length)
    #print 'x_range: ' + str(x_range)
    x_range = u.Quantity([ x_range[0] - origin[0], x_range[1] - origin[0] ])
    #print 'x_range: ' + str(x_range)
    y_range = aMap2D_cropped.yrange.to(u.meter, equivalencies=radian_length)
    #print 'y_range: ' + str(y_range)
    y_range = u.Quantity([ y_range[0] - origin[1], y_range[1] - origin[1] ])
    #print 'y_range: ' + str(y_range)
    z_range = x_range

    # Vector field as a 3D map
    #a4DArray = np.random.rand(dim[0],dim[1],dim[0],3)
    a4DArray = np.load('F://Dropbox//Study/2014-2015//QM Project-Dissertation//Python//data_extrapolations//2011-02-14__20-35-25__02_Bxyz.npy')
    aMetaDict = { 'file': 'test SunPy Map object'}
    aMap3D = Map3D(a4DArray, aMetaDict, xrange=x_range, yrange=y_range, zrange=z_range, xobsrange=x_obs_range, yobsrange=y_obs_range)
    
    
    #visualise(aMap3D, boundary=aMap2D, scale=1.0*u.Mm, boundary_units=1.0*u.arcsec, show_volume_axes=False)
    seeds = np.array([[4,4,2], [-4,4,2], [4,-4,2], [-4,-4,2], [2,2,2], [-2,2,2], [2,-2,2], [-2,-2,2]])
    seeds = None
    visualise(aMap3D, boundary=aMap2D, seeds=seeds, scale=1.0*u.Mm, boundary_units=1.0*u.arcsec, show_boundary_axes=False, show_volume_axes=True)#, origin=origin)
    """
    # Vector field as a 3D map
    a4DArray = np.random.rand(dim[0],dim[1],dim[0],3)
    '''
    for x in range(0,10):
        for y in range(0,10):
            a4DArray[x,y,0,0] = float('nan')
    '''
    
    
    

    
    x_range=[0.0, 20.0]*u.Mm
    y_range=[0.0, 20.0]*u.Mm
    z_range=[0.0, 20.0]*u.Mm
    x_obs_range = [-15,-15] * u.arcsec
    y_obs_range = [-15,-15] * u.arcsec

    aMap3D = Map3D(a4DArray, aMetaDict, xrange=x_range, yrange=y_range, zrange=z_range)
    #aMap3D = Map3D(a4DArray, aMetaDict)
    # Note: consider scipy.interpolate.griddata for 3D subsampling
    
    # Boundary data as a sunpy map
    #aMap2D = mp.Map('C://git//solarextrapolation//solarextrapolation//data//example_data_(100x100)__01_hmi.fits')
    #aMap2D = mp.Map('F://Dropbox//Study//2014-2015//QM Project-Dissertation//Python//data_fits//2011-02-14__20-35-25__01_hmi.fits').submap([-22,22], [-22,22])
    aMap2D = mp.Map('C:/fits/2011-02-14__20-35-25__01_hmi.fits').submap([-22,22] * u.arcsec, [-22,22] * u.arcsec)
    
    # Plotthe lot
    visualise(aMap3D, boundary=aMap2D, seeds=np.array([[5,5,3]]), scale=1000*u.kilometer, boundaryscale=1*u.arcsec)
    #visualise(aMap3D, boundary=aMap2D)
    """