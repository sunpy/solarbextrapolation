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
    a4DArray = np.random.rand(5,5,5)
    aMetaDict = { 'file': 'test SunPy Map object'}
    aMap3D = Map3D(a4DArray, aMetaDict)
    # Note: consider scipy.interpolate.griddata for 3D subsampling
    
    fig = mlab.figure()
    