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
    
    # Slice (scale) the fields make the vectors usable in mayavi.
    int_slice_scale = 1
    npm_3d_sliced   = aMap3D.data[::int_slice_scale,::int_slice_scale,::int_slice_scale,:]
    
    # Plot the main vector field (volume).
    fig = mlab.figure()

    # Scale from astronomical units to something more usable with mayavi
    mayavi_scale = 1#10**(-6)

    # Make 3D coords for ever point in the 3D grid.
    X, Y, Z = np.mgrid[aMap3D.xrange[0]*mayavi_scale:aMap3D.xrange[1]*mayavi_scale:npm_3d_sliced.shape[0]*1j,
                       aMap3D.yrange[0]*mayavi_scale:aMap3D.yrange[1]*mayavi_scale:npm_3d_sliced.shape[1]*1j,
                       aMap3D.zrange[0]*mayavi_scale:aMap3D.zrange[1]*mayavi_scale:npm_3d_sliced.shape[2]*1j]
    vec_field = vector_field(X, Y, Z, npm_3d_sliced[:,:,:,0], npm_3d_sliced[:,:,:,1], npm_3d_sliced[:,:,:,2],
                             name='Magnetic Vector Field', figure=fig)
    vec_field_mag = mlab.pipeline.extract_vector_norm(vec_field, name="Magnetic Field Magnitude")
    
    # Place a small outline around the data cube
    mlab.outline()
    axes = mlab.axes()
    
    # Plot streamlines
    if np_seeds is None:
        tupSeedGrid = (4, 4) # (x, y) divisions, if a file isn't declared.
        # Avoid placing seeds in ghost cells
        intSeedZ = 20# * mayavi_scale
        tupCentrePixel = [ int(aMap3D.data.shape[0] / 2.0), int(aMap3D.data.shape[0] / 2.0)]
        while math.isnan(aMap3D.data[tupCentrePixel[0],tupCentrePixel[1],intSeedZ, 0]):
            intSeedZ += 1
            print 'intSeedZ: ' + str(intSeedZ)
        floZ = intSeedZ# * mayavi_scale
        
        intIncX = ((aMap3D.xrange[1] * 1.0 - aMap3D.xrange[0]) / (tupSeedGrid[1] + 1.0)) * mayavi_scale
        intIncY = ((aMap3D.zrange[1] * 1.0 - aMap3D.zrange[0]) / (tupSeedGrid[0] + 1.0)) * mayavi_scale
        lisSeeds = []
        for i in range(1, tupSeedGrid[0] + 1):
            for j in range(1, tupSeedGrid[1] + 1):
                lisSeed = [ i * intIncX + aMap3D.xrange[0] * mayavi_scale, j * intIncY + aMap3D.yrange[0] * mayavi_scale, floZ ]
                lisSeeds.append(lisSeed)
                print lisSeed
        np_seeds = np.array(lisSeeds)

    # Plot the seed points
    points = mlab.points3d(np_seeds[:,0], np_seeds[:,1], np_seeds[:,2])
    # Make the points smaller
    points.glyph.glyph.scale_factor = mayavi_scale
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


    if boundary:
        # CReate explicit points in 3D space
        X, Y = np.mgrid[aMap2D.xrange[0]*mayavi_scale:aMap2D.xrange[1]*mayavi_scale:aMap2D.data.shape[0]*1j,
                        aMap2D.yrange[0]*mayavi_scale:aMap2D.yrange[1]*mayavi_scale:aMap2D.data.shape[1]*1j]
        # Plot and add to the current figure
        img_boundary = mlab.pipeline.array2d_source(X, Y, aMap2D.data, figure=fig)
        img_boundary = mlab.pipeline.image_actor(img_boundary, figure = fig)
        
        # Color the image according to the data ############ Issue HERE
        #mayavi_ct = aMap2D.cmap(range(255))
        #img_boundary.module_manager.scalar_lut_manager.lut.table = mayavi_ct
        
        # Legend details
        img_boundary.module_manager.scalar_lut_manager.show_legend = True #module_manager2.scalar_lut_manager.show_legend = True
        img_boundary.module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.1,  0.1 ])
    
    
    # And open mayavi
    mlab.show()

if __name__ == '__main__':
    # Vector field as a 3D map
    a4DArray = np.random.rand(10,10,10,3)
    for x in range(0,10):
        for y in range(0,10):
            a4DArray[x,y,0,0] = float('nan')
    aMetaDict = { 'file': 'test SunPy Map object'}
    aMap3D = Map3D(a4DArray, aMetaDict)
    # Note: consider scipy.interpolate.griddata for 3D subsampling
    
    # Boundary data as a sunpy map
    #aMap2D = mp.Map('C://git//solarextrapolation//solarextrapolation//data//example_data_(100x100)__01_hmi.fits')
    #aMap2D = mp.Map('F://Dropbox//Study//2014-2015//QM Project-Dissertation//Python//data_fits//2011-02-14__20-35-25__01_hmi.fits').submap([-22,22], [-22,22])
    aMap2D = mp.Map('C:/fits/2011-02-14__20-35-25__01_hmi.fits').submap([-22,22] * u.arcsec, [-22,22] * u.arcsec)
    
    # Plotthe lot
    #visualise(aMap3D, boundary=aMap2D, seeds=np.array([[0,0,1]]))
    visualise(aMap3D, boundary=aMap2D)