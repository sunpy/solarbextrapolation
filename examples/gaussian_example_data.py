# -*- coding: utf-8 -*-
"""
=========================================
Generating Example Gaussian Boundary Data
=========================================

In this example you will be generating some example data and extrapolate this
using the basic potential extrapolator.
"""

##############################################################################
# You can start by importing the necessary module components.

# Module imports
from solarbextrapolation.example_data_generator import generate_example_data, dummyDataToMap

##############################################################################
# You also need the ability to convert astropyunits.
import astropy.units as u

##############################################################################
# You need to define the parameters of the eare, includsing the x and y ranges
# as astropy quantities with angular or distance units and the grid shape.

# Input parameters:
arr_grid_shape = [ 20, 22 ]         # [ y-size, x-size ]
qua_xrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
qua_yrange = u.Quantity([ -11.0, 11.0 ] * u.arcsec)

##############################################################################
# The generated data will consist of a 2D space with 2 Gaussian spots, one
# positive and one negative, on a background of 0.0.
# solarbextrapolation.example_data_generator provides many ways to achieve this,
# including letting it randomly generate the position, magnitude and size of
# each spot/pole.

# To randomly generate 2 poles simply don't add any pole parameters:
arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange)
# Note: each time you run this pole positions/magnitudes will change.

##############################################################################
# We can now convert this into a a sunpy map object:
aMap = dummyDataToMap(arr_Data, qua_xrange, qua_yrange)

##############################################################################
# We can see this map using peek:
aMap.peek()

##############################################################################
# To manually position poles, simply build lists of parameters for each pole.
# It's often easiest to use percentage units for location/size, wheer we compare
# to the maps region.
# arrA0 = [ Position, size, Max Magnitude ]
arrA0 = [ u.Quantity([ 25, 25 ] * u.percent), 10.0 * u.percent,  0.2 * u.T ]
arrA1 = [ u.Quantity([ 75, 75 ] * u.percent), 10.0 * u.percent, -0.2 * u.T ]

# To generate and view:
arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange, arrA0, arrA1)
aMap = dummyDataToMap(arr_Data, qua_xrange, qua_yrange)
aMap.peek()

##############################################################################
# But absolute positioning using the map range units is also possible
arrA2 = [ u.Quantity([ -6,  6 ] * u.arcsec), 2 * u.arcsec, -0.2 * u.T ]
arrA3 = [ u.Quantity([  6, -7 ] * u.arcsec), 2 * u.arcsec,  0.2 * u.T ]

# To generate and view:
arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange, arrA2, arrA3)
aMap = dummyDataToMap(arr_Data, qua_xrange, qua_yrange)
aMap.peek()

##############################################################################
# You can add as many poles as you want:
arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange, arrA0, arrA1, arrA2, arrA3)
aMap = dummyDataToMap(arr_Data, qua_xrange, qua_yrange)
aMap.peek()

##############################################################################
# And being a map you can use all the normal SunPy functions, such as saving
# the map using aMap.save(filepath).
