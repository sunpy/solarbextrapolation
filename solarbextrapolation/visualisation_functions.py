# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 18:45:05 2015

@author: alex_
"""

# Temp
import numpy as np
from classes import *
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
    np_seeds = kwargs.get('seeds', None)
    boundary = kwargs.get('boundary', None)
    
    mayavi_unit_length = kwargs.get('unit_length', 1.0 * u.Mm)
    units_boundary = kwargs.get('boundary_units', mayavi_unit_length)
    units_volume = kwargs.get('volume_units', mayavi_unit_length)
    show_boundary_axes = kwargs.get('show_boundary_axes', True)
    show_volume_axes = kwargs.get('show_volume_axes', True)
    
    # Setup the equivilence
    obs_distance = aMap3D.dsun - aMap3D.rsun_meters
    radian_length = [ (u.radian, u.meter, lambda x: obs_distance * x, lambda x: x / obs_distance) ]
    
    # Slice (scale) the fields make the vectors usable in mayavi.
    int_slice_scale = 1
    npm_3d_sliced   = aMap3D.data[::int_slice_scale,::int_slice_scale,::int_slice_scale,:]
    
    # Plot the main vector field (volume).
    fig = mlab.figure()
    
    # Make 3D coords for ever point in the 3D grid.
    x_range = aMap3D.xobsrange.to(u.meter, equivalencies=radian_length)
    y_range = aMap3D.yobsrange.to(u.meter, equivalencies=radian_length)
    z_range = aMap3D.zrange.to(u.meter, equivalencies=radian_length)
    x_range_scaled = (x_range/mayavi_unit_length).decompose().value
    y_range_scaled = (y_range/mayavi_unit_length).decompose().value
    z_range_scaled = (z_range/mayavi_unit_length).decompose().value
    print '\n\n' + str(aMap3D.xrange) + '\n' + str(aMap3D.xrange)
    print str(x_range_scaled) + '\n' + str(x_range_scaled) + '\n\n'
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
        tupSeedGrid = (4, 4) # (x, y) divisions, if a file isn't declared.
        # Avoid placing seeds in ghost cells
        intSeedZ = 3
        tupCentrePixel = [ int(aMap3D.data.shape[0] / 2.0), int(aMap3D.data.shape[0] / 2.0)]
        while math.isnan(aMap3D.data[tupCentrePixel[0],tupCentrePixel[1],intSeedZ, 0]):
            intSeedZ += 1
            print 'intSeedZ: ' + str(intSeedZ)
        floZ = intSeedZ
        
        intIncX = ((aMap3D.xrange[1] * 1.0 - aMap3D.xrange[0]) / (tupSeedGrid[1] + 1.0))
        intIncY = ((aMap3D.zrange[1] * 1.0 - aMap3D.zrange[0]) / (tupSeedGrid[0] + 1.0))
        lisSeeds = []
        for i in range(1, tupSeedGrid[0] + 1):
            for j in range(1, tupSeedGrid[1] + 1):
                lisSeed = [ i * intIncX + aMap3D.xrange[0], j * intIncY + aMap3D.yrange[0], floZ ]
                lisSeeds.append(lisSeed)
                print lisSeed
        np_seeds = np.array(lisSeeds)

    # Plot the seed points
    points = mlab.points3d(np_seeds[:,0], np_seeds[:,1], np_seeds[:,2])
    # Make the points smaller
    points.glyph.glyph.scale_factor = 1.0 #mayavi_scale
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
        x_range = boundary.xrange.to(u.meter, equivalencies=radian_length)
        y_range = boundary.yrange.to(u.meter, equivalencies=radian_length)
        x_range_scaled = (x_range/mayavi_unit_length).decompose().value
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
            x_range_scaled = (x_range/units_boundary).decompose().value
            y_range_scaled = (y_range/units_boundary).decompose().value
            axes.axes.ranges = np.array([ x_range_scaled[0],  x_range_scaled[1],  y_range_scaled[0],  y_range_scaled[1],  0,  0])
            str_units = ' (' + str(units_boundary) + ')'
            axes.axes.x_label = 'Solar X ' + str_units
            axes.axes.y_label = 'Solar Y' + str_units
    # And open mayavi
    mlab.show()

def scaled_length(length, scale_length, map):
    if length.unit.is_equivalent(u.radian):
        #print '\nscaled_length(' + str(length) + ', ' + str(scale_length) + ', map)'
        distance = (map.dsun - map.rsun_meters)
        #print 'distance: ' + str(distance)
        length = (distance * length.to(u.radian)).to(u.m, equivalencies=u.dimensionless_angles())
        #print 'length: ' + str(length)
        #print 'output: ' + str((length / scale_length).decompose().value)
    return (length / scale_length).decompose().value


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
    
    # Using the dimensions of the boundary data to get the volume size
    dim = [ aMap2D_cropped.data.shape[0], aMap2D_cropped.data.shape[1] ]
    
    # Observational ranges are from the boundary data
    x_obs_range = aMap2D_cropped.xrange
    y_obs_range = aMap2D_cropped.yrange
    
    # Getting the ranges. x/y from the boundary, z is defined manually.    
    scale_length = 1.0 * u.Mm
    x_range = scaled_length(x_obs_range, scale_length, aMap2D_cropped) * scale_length.unit
    y_range = scaled_length(y_obs_range, scale_length, aMap2D_cropped) * scale_length.unit
    z_range = [0.0, 10.0] * u.Mm

    # Vector field as a 3D map
    a4DArray = np.random.rand(dim[0],dim[1],dim[0],3)
    aMetaDict = { 'file': 'test SunPy Map object'}
    aMap3D = Map3D(a4DArray, aMetaDict, xrange=x_range, yrange=y_range, zrange=z_range, xobsrange=x_obs_range, yobsrange=y_obs_range)
    
    visualise(aMap3D, boundary=aMap2D, seeds=np.array([[5,5,3]]), scale=1.0*u.Mm, boundary_units=1*u.arcsec, show_volume_axes=False)
    
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