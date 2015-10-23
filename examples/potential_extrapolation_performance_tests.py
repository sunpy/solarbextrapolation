# -*- coding: utf-8 -*-
"""
=====================================
Testing Perfomance of an Extrapolator
=====================================

This tests the speed of it all.
"""

# General imports
from astropy import units as u
from astropy.table import Table
import numpy as np

# Module imports
from solarbextrapolation.extrapolators import PotentialExtrapolator
from solarbextrapolation.example_data_generator import generate_example_data, dummyDataToMap

# The input parameters:
lis_grid_shapes = [ [ 100, 100, 100 ] ]#, [ 20, 20, 20 ]]#, [ 30, 30, 30 ]]#, [ 100, 100, 100 ]]#[ 10, 10, 10 ],[ 50, 50, 50 ], [ 100, 100, 100 ], [ 200, 200, 200 ] ]
xrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
yrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
zrange = u.Quantity([ 0,     20.0 ] * u.arcsec)

# Manual Pole Details
arrA0 = [ u.Quantity([ 25, 25 ] * u.percent), 10.0 * u.percent,  0.2 * u.T ]
arrA1 = [ u.Quantity([ 75, 75 ] * u.percent), 10.0 * u.percent, -0.2 * u.T ]
arrA2 = [ u.Quantity([ 25, 75 ] * u.percent), 10.0 * u.percent,  0.1 * u.T ]
arrA3 = [ u.Quantity([ 75, 25 ] * u.percent), 10.0 * u.percent, -0.1 * u.T ]

# Generate the datasets and maps
lis_maps = []
lis_extrapolators = []

# A table for storing the data
t = Table(names=('grid size', 'time (min)', 'time (ave)', 'time (std)'), meta={'name': 'times tables'}, dtype=('S24', 'f8', 'f8', 'f8'))
t['time (min)'].unit = u.s
t['time (ave)'].unit = u.s
t['time (std)'].unit = u.s

# A list of all the necessary parameters to perform an extrapolation
lis_datasets = []
for shape in lis_grid_shapes:
    lis_datasets.append([ str(shape), shape[2], zrange,
                          dummyDataToMap(generate_example_data(shape[0:2], xrange, yrange, arrA0, arrA1, arrA2, arrA3), xrange, yrange) ])
int_trials = 1 # The times to repeat each extrapolation.

# Iterate through the extrapolations
for extrapolation in lis_datasets:
    # Setup the extrapolator and table
    aPotExt = PotentialExtrapolator(extrapolation[3], zshape=extrapolation[1], zrange=extrapolation[2])

    # List to store the trial
    lis_times = []

    # Run the extrapolation without numba for each dataset (map and ranges).
    for i in range(0, int_trials):
        aMap3D = aPotExt.extrapolate(enable_numba=False)
        lis_times.append(aMap3D.meta['extrapolator_duration'])
    t.add_row([extrapolation[0], np.round(np.min(lis_times), 2), np.round(np.average(lis_times), 2), np.round(np.std(lis_times), 2)])

    # List to store the trial
    lis_times = []

    # Run the extrapolation with numba for each dataset (map and ranges).
    for i in range(0, int_trials):
        aMap3D = aPotExt.extrapolate(enable_numba=True)
        lis_times.append(aMap3D.meta['extrapolator_duration'])
    t.add_row(['(numba)'+extrapolation[0], np.round(np.min(lis_times), 2), np.round(np.average(lis_times), 2), np.round(np.std(lis_times), 2)])

# Show the data table
print t
print '\n\n'


