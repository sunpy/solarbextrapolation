# -*- coding: utf-8 -*-
"""
==========================
Example Data Extrapolation
==========================

Example of creating a basic example dataset, extrapolating using the potential
extrapolator and visualising.
"""

# General imports
from astropy import units as u
from mayavi import mlab

# Module imports
from solarbextrapolation.classes import *
from solarbextrapolation.potential_field_extrapolator import *
from solarbextrapolation.utilities import *
from solarbextrapolation.example_data_generator import *
from solarbextrapolation.visualisation_functions import *

# The input parameters:
arr_grid_shape = [ 20, 22, 20 ]         # [ y-size, x-size ]
xrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
yrange = u.Quantity([ -11.0, 11.0 ] * u.arcsec)
zrange = u.Quantity([ 0,     20.0 ] * u.arcsec)

# Manual Pole Details
arrA0 = [ u.Quantity([ 25, 25 ] * u.percent), 10.0 * u.percent,  0.2 * u.T ]
arrA1 = [ u.Quantity([ 75, 75 ] * u.percent), 10.0 * u.percent, -0.2 * u.T ]

# Generate the data and make into a map
arr_data = generate_example_data(arr_grid_shape[0:2], xrange, yrange, arrA0, arrA1)
map_boundary = dummyDataToMap(arr_data, xrange, yrange)

# Use potential extrapolator to generate field
aPotExt = PotentialExtrapolator(map_boundary, zshape=arr_grid_shape[2], zrange=zrange)
aMap3D  = aPotExt.extrapolate(enable_numba=True)
#print '\nextrapolation duration: ' + str(np.round(aMap3D.meta['extrapolator_duration'],3)) + ' s\n'

# Visualise
fig = visualise(aMap3D,
                boundary=map_boundary,
                volume_units=[1.0*u.arcsec, 1.0*u.arcsec, 1.0*u.Mm],
                show_boundary_axes=False,
                boundary_units=[1.0*u.arcsec, 1.0*u.arcsec],
                show_volume_axes=True,
                debug=False)

#mlab.show()
