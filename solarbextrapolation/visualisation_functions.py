# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 18:45:05 2015

@author: alex_
"""

# Temp
import numpy as np
from classes import *

# Visulisation
import yt
import mayavi.sources
from mayavi import mlab
from mayavi_seed_streamlines import SeedStreamline
from tvtk.tools import visual


if __name__ == '__main__':
    # Vector field as a 3D map
    a4DArray = np.random.rand(5,5,5)
    aMetaDict = { 'file': 'test SunPy Map object'}
    aMap3D = Map3D(a4DArray, aMetaDict)
    # Note: consider scipy.interpolate.griddata for 3D subsampling
    
    # Boundary data as a sunpy map
    aMap2D = sunpy.map.Map('C://git//solarextrapolation//solarextrapolation//data//example_data_(100x100)__01_hmi.fits')
    
    
    
    
    # Plot the main vector field (volume).
    fig = mlab.figure()
    
    # Scale from astronomical units to something more usable with mayavi
    mayavi_scale = 10**(-6)

    # Slice (scale) the fields make the vectors usable in mayavi.
    int_slice_scale = 20
    npm_3d_sliced   = aMap3D.data[::int_scale,::int_scale,::int_scale,:]
    # Make 3D coords for ever point in the 3D grid.
    X, Y, Z = np.mgrid[aMap3D.xrange[0]*mayavi_scale:aMap3D.xrange[1]*mayavi_scale:npm_3d.shape[0]*1j,
                       aMap3D.yrange[0]*mayavi_scale:aMap3D.yrange[1]*mayavi_scale:npm_3d.shape[1]*1j,
                       aMap3D.zrange[0]*mayavi_scale:aMap3D.zrange[1]*mayavi_scale:npm_3d.shape[2]*1j]
    vec_field = vector_field(X, Y, Z, npm_3d_sliced[:,:,:,0], npm_3d_sliced[:,:,:,1], npm_3d_sliced[:,:,:,2],
                             name='Magnetic Vector Field', figure=fig)  
    vec_field_mag = mlab.pipeline.extract_vector_norm(vec_field, name="Magnetic Field Magnitude")
    
    
    
    
    # For the seeds, either load in the csv file OR use a simple grid.
    npSeeds = np.array([])
    tupSeedGrid = (10, 10) # (x, y) divisions, if a file isn't declared.
    intSeedZ = 2000000 * mayavi_scale
    str_seeds_path = 'file'
    if os.path.isfile(str_seeds_path):
        npSeeds = np.genfromtxt(str_seeds_path, dtype=int, delimiter=',')
    else:
        # Consider seed points.
        intIncX = ((tup_boundaries[0][1] * 1.0 / tup_boundaries[0][0] ) / (tupSeedGrid[0] + 1.0))
        intIncX = ((tup_boundaries[0][1] * 1.0 - tup_boundaries[0][0] ) / (tupSeedGrid[0] + 1.0)) * mayavi_scale
        print 'intIncX: ' + str(intIncX)
        intIncY = ((tup_boundaries[1][1] * 1.0 / tup_boundaries[1][0] ) / (tupSeedGrid[1] + 1.0) )
        intIncY = ((tup_boundaries[1][1] * 1.0 - tup_boundaries[1][0] ) / (tupSeedGrid[1] + 1.0) ) * mayavi_scale
        lisSeeds = []    
        for i in range(1, tupSeedGrid[0] + 1):
            for j in range(1, tupSeedGrid[1] + 1):
                lisSeed = [ i * intIncX + tup_boundaries[0][0] * mayavi_scale, j * intIncY + tup_boundaries[1][0] * mayavi_scale, intSeedZ ]
                lisSeeds.append(lisSeed)
                print lisSeed
        npSeeds = np.array(lisSeeds)
        #npSeeds = np.array([[-30.0,0.0, 2.0]])#[[-81.0,-131.0, 2.0],[-81.0,-113.0, 2.0],[-81.0,-95.0, 2.0],[-81.0,-77.0, 2.0],[-81.0,-59.0, 2.0],[-81.0,-41.0, 2.0],[-81.0,-23.0, 2.0],[-81.0,-5.0, 2.0],[-81.0,13.0, 2.0],[-81.0,31.0, 2.0],[-81.0,49.0, 2.0],[-81.0,67.0, 2.0],[-81.0,85.0, 2.0],[-81.0,103.0, 2.0],[-81.0,121.0, 2.0],[-63.0,-131.0, 2.0],[-63.0,-113.0, 2.0],[-63.0,-95.0, 2.0],[-63.0,-77.0, 2.0],[-63.0,-59.0, 2.0],[-63.0,-41.0, 2.0],[-63.0,-23.0, 2.0],[-63.0,-5.0, 2.0],[-63.0,13.0, 2.0],[-63.0,31.0, 2.0],[-63.0,49.0, 2.0],[-63.0,67.0, 2.0],[-63.0,85.0, 2.0],[-63.0,103.0, 2.0],[-63.0,121.0, 2.0],[-45.0,-131.0, 2.0],[-45.0,-113.0, 2.0],[-45.0,-95.0, 2.0],[-45.0,-77.0, 2.0],[-45.0,-59.0, 2.0],[-45.0,-41.0, 2.0],[-45.0,-23.0, 2.0],[-45.0,-5.0, 2.0],[-45.0,13.0, 2.0],[-45.0,31.0, 2.0],[-45.0,49.0, 2.0],[-45.0,67.0, 2.0],[-45.0,85.0, 2.0],[-45.0,103.0, 2.0],[-45.0,121.0, 2.0],[-27.0,-131.0, 2.0],[-27.0,-113.0, 2.0],[-27.0,-95.0, 2.0],[-27.0,-77.0, 2.0],[-27.0,-59.0, 2.0],[-27.0,-41.0, 2.0],[-27.0,-23.0, 2.0],[-27.0,-5.0, 2.0],[-27.0,13.0, 2.0],[-27.0,31.0, 2.0],[-27.0,49.0, 2.0],[-27.0,67.0, 2.0],[-27.0,85.0, 2.0],[-27.0,103.0, 2.0],[-27.0,121.0, 2.0],[-9.0,-131.0, 2.0],[-9.0,-113.0, 2.0],[-9.0,-95.0, 2.0],[-9.0,-77.0, 2.0],[-9.0,-59.0, 2.0],[-9.0,-41.0, 2.0],[-9.0,-23.0, 2.0],[-9.0,-5.0, 2.0],[-9.0,13.0, 2.0],[-9.0,31.0, 2.0],[-9.0,49.0, 2.0],[-9.0,67.0, 2.0],[-9.0,85.0, 2.0],[-9.0,103.0, 2.0],[-9.0,121.0, 2.0],[9.0,-131.0, 2.0],[9.0,-113.0, 2.0],[9.0,-95.0, 2.0],[9.0,-77.0, 2.0],[9.0,-59.0, 2.0],[9.0,-41.0, 2.0],[9.0,-23.0, 2.0],[9.0,-5.0, 2.0],[9.0,13.0, 2.0],[9.0,31.0, 2.0],[9.0,49.0, 2.0],[9.0,67.0, 2.0],[9.0,85.0, 2.0],[9.0,103.0, 2.0],[9.0,121.0, 2.0],[27.0,-131.0, 2.0],[27.0,-113.0, 2.0],[27.0,-95.0, 2.0],[27.0,-77.0, 2.0],[27.0,-59.0, 2.0],[27.0,-41.0, 2.0],[27.0,-23.0, 2.0],[27.0,-5.0, 2.0],[27.0,13.0, 2.0],[27.0,31.0, 2.0],[27.0,49.0, 2.0],[27.0,67.0, 2.0],[27.0,85.0, 2.0],[27.0,103.0, 2.0],[27.0,121.0, 2.0],[45.0,-131.0, 2.0],[45.0,-113.0, 2.0],[45.0,-95.0, 2.0],[45.0,-77.0, 2.0],[45.0,-59.0, 2.0],[45.0,-41.0, 2.0],[45.0,-23.0, 2.0],[45.0,-5.0, 2.0],[45.0,13.0, 2.0],[45.0,31.0, 2.0],[45.0,49.0, 2.0],[45.0,67.0, 2.0],[45.0,85.0, 2.0],[45.0,103.0, 2.0],[45.0,121.0, 2.0],[63.0,-131.0, 2.0],[63.0,-113.0, 2.0],[63.0,-95.0, 2.0],[63.0,-77.0, 2.0],[63.0,-59.0, 2.0],[63.0,-41.0, 2.0],[63.0,-23.0, 2.0],[63.0,-5.0, 2.0],[63.0,13.0, 2.0],[63.0,31.0, 2.0],[63.0,49.0, 2.0],[63.0,67.0, 2.0],[63.0,85.0, 2.0],[63.0,103.0, 2.0],[63.0,121.0, 2.0]])
        
    # Plot the seed points
    points = mlab.points3d(npSeeds[:,0], npSeeds[:,1], npSeeds[:,2])
    # Make the points smaller
    points.glyph.glyph.scale_factor = 2000000*mayavi_scale
    # Make the points blue
    points.actor.property.color = (0.2,0,1)
    # Create the custom streamline object
    slines = SeedStreamline(seed_points=npSeeds)
    # Add the streamline object to the plot and make it use the magentic field data,
    # by adding it as a child of the field we created earlier.
    # We add it to the magnitude field (which is in itself a child of bfield)
    # so that it picks up the scalar values and colours the lines.
    bmag.add_child(slines)
    
    # Make the lines a cool colour
    slines.module_manager.scalar_lut_manager.lut_mode = 'winter'    
    
    
    # Plot the boundary data
    '''
    # Make 3D coords for ever point in the 3D grid.
    X, Y = np.mgrid[aMap2D.xrange[0]*mayavi_scale:aMap2D.xrange[1]*mayavi_scale:npm_3d.shape[0]*1j,
                    aMap2D.yrange[0]*mayavi_scale:aMap2D.yrange[1]*mayavi_scale:npm_3d.shape[1]*1j]
    img_source_hmi_xyz = mlab.pipeline.array2d_source(X, Y, map_hmi_data_cropped_surface.data, figure=fig_hmi_xyz)
    img_hmi_xyz = mlab.pipeline.image_actor(img_source_hmi_xyz, figure = fig_hmi_xyz)
    img_hmi_xyz.module_manager.scalar_lut_manager.show_legend = True #module_manager2.scalar_lut_manager.show_legend = True
    img_hmi_xyz.module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.1,  0.1 ])
    '''