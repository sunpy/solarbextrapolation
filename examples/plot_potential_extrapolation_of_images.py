# -*- coding: utf-8 -*-
"""
=====================================
Potential Extrapolation from an Image
=====================================

This example is to demonstrate using the potential extrapolator on an image.
It was built for a bit of fun.
"""

from __future__ import print_function

# General imports
from astropy import units as u
from scipy import misc
from mayavi import mlab
import numpy as np

# Module imports
from solarbextrapolation.potential_field_extrapolator import PotentialExtrapolator
from solarbextrapolation.example_data_generator import generate_example_data, dummyDataToMap
from solarbextrapolation.visualisation_functions import visualise

# The input parameters:
arr_grid_shape = [ 50, 50, 50 ]         # [ y-size, x-size ]
xrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
yrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
zrange = u.Quantity([ 0,     20.0 ] * u.arcsec)

# Manual Pole Details
arrA0 = [ u.Quantity([ 25, 25 ] * u.percent), 10.0 * u.percent,  0.2 * u.T ]
arrA1 = [ u.Quantity([ 75, 75 ] * u.percent), 10.0 * u.percent, -0.2 * u.T ]

# Generate the data and make into a map
arr_data = generate_example_data(arr_grid_shape[0:2], xrange, yrange, arrA0, arrA1)
print('\n' + str(type(arr_data.dtype)))
print(str(arr_data.shape))
print(str(arr_data.dtype) + '\n')
arr_image = (misc.imread('sunpy_powered_50x50.png')[...,:3] - 127.5) / 1270.0
arr_data = np.zeros(arr_image.shape[:2])
for i in range(0,arr_data.shape[0]):     # Row/Y
    for j in range(0,arr_data.shape[1]): # Column/X
        arr_data[i,j] = ((arr_image[i,j,0] + arr_image[i,j,1] + arr_image[i,j,2]) / 3.0)
print('\n' + str(type(arr_data.dtype)))
print(str(arr_data.shape))
print(str(arr_data.dtype) + '\n')

map_boundary = dummyDataToMap(arr_data, xrange, yrange)

# Use potential extrapolator to generate field
aPotExt = PotentialExtrapolator(map_boundary, zshape=arr_grid_shape[2], zrange=zrange)
aMap3D = aPotExt.extrapolate(enable_numba=True)
print('extrapolator_duration:' + str(aMap3D.meta['extrapolator_duration']))


# Visualise
visualise(aMap3D,
          boundary=map_boundary,
          volume_units=[1.0*u.arcsec, 1.0*u.arcsec, 1.0*u.Mm],
          show_boundary_axes=False,
          boundary_units=[1.0*u.arcsec, 1.0*u.arcsec],
          show_volume_axes=True,
          debug=False)
mlab.show()
