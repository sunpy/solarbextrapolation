# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 17:53:45 2015

@author: alex_

Example of creating a basic example dataset, extrapolating using the potential
extrapolator and visualising.
"""

from astropy import units as u

# Module Imports
from classes import *
from potential_field_extrapolator import *
from utilities import *
from example_data_generator import *
from visualisation_functions import *

if __name__ == '__main__':
    # Generate an example map
    # The input parameters:
    arr_grid_shape = [ 20, 22, 20 ]         # [ y-size, x-size ]
    xrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
    yrange = u.Quantity([ -11.0, 11.0 ] * u.arcsec)
    zrange = u.Quantity([ 0,     20.0 ] * u.arcsec)
    
    # Manual Pole Details
    arrA0 = [ u.Quantity([ 25, 25 ] * u.percent), 10.0 * u.percent,  0.2 * u.T ]
    arrA1 = [ u.Quantity([ 75, 75 ] * u.percent), 10.0 * u.percent, -0.2 * u.T ]

    # Generate the data and save to a map
    arr_data = generate_example_data(arr_grid_shape[0:2], xrange, yrange, arrA0, arrA1)
    map_boundary = dummyDataToMap(arr_data, xrange, yrange)
    
    # Do the extrapolation
    aPotExt = PotentialExtrapolator(map_boundary, zshape=arr_grid_shape[2], zrange=zrange)
    aMap3D = aPotExt.extrapolate()
    
    # Visualise this
    visualise(aMap3D,
              boundary=map_boundary,
              volume_units=[1.0*u.arcsec, 1.0*u.arcsec, 1.0*u.Mm],
              show_boundary_axes=False,
              boundary_units=[1.0*u.arcsec, 1.0*u.arcsec],
              show_volume_axes=True,
              debug=False)