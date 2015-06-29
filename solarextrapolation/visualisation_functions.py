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
    map.extent
    
    
    
    
    # Plot the main vector field (volume).
    fig = mlab.figure()
    
    # Scale from astronomical units to something more usable with mayavi
    mayavi_scale = 10**(-6)

    # Slice (scale) the fields make the vectors usable in mayavi.
    int_slice_scale = 20
    npm_3d_cropped     = npm_B_q[::int_scale,::int_scale,::int_scale,:]
    # Make 3D coords for ever point in the 3D grid.
    X, Y, Z = np.mgrid[tup_boundaries[0][0]*mayavi_scale:tup_boundaries[0][1]*mayavi_scale:npm_B_TD_cropped.shape[0]*1j,
                       tup_boundaries[1][0]*mayavi_scale:tup_boundaries[1][1]*mayavi_scale:npm_B_TD_cropped.shape[1]*1j,
                       tup_boundaries[2][0]*mayavi_scale:tup_boundaries[2][1]*mayavi_scale:npm_B_TD_cropped.shape[2]*1j]
    
    
    
    
    # Plot the boundary data

    # Make 3D coords for ever point in the 3D grid.
    X, Y, Z = np.mgrid[tup_boundaries[0][0]*mayavi_scale:tup_boundaries[0][1]*mayavi_scale:npm_B_TD_cropped.shape[0]*1j,
                       tup_boundaries[1][0]*mayavi_scale:tup_boundaries[1][1]*mayavi_scale:npm_B_TD_cropped.shape[1]*1j,
                       tup_boundaries[2][0]*mayavi_scale:tup_boundaries[2][1]*mayavi_scale:npm_B_TD_cropped.shape[2]*1j]