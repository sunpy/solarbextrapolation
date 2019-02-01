# -*- coding: utf-8 -*-
"""
=====================================
Extrapolator Perfomance Testing
=====================================

In this example you will be running the potential field extrapolator both with
numba enabled and disabled over a number of datasets and tabulating the results
into an astropy table.

Note: if you don't have conda numba installed the code should work but the
results should not show any speed difference.
"""

##############################################################################
# You can start by importing the necessary module components.

# Module imports
from solarbextrapolation.extrapolators import PotentialExtrapolator
from solarbextrapolation.example_data_generator import generate_example_data, dummyDataToMap

##############################################################################
# You also need the ability to convert astropyunits, use numpy arrays and
# astropy tables.

# General imports
from astropy import units as u
from astropy.table import Table
import numpy as np

##############################################################################
# You are going to create a series of volume grids with a given shape and then
# attribute arbitrary units to it's axes.
lis_grid_shapes = [ [ 50, 50, 50 ] ]
xrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
yrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
zrange = u.Quantity([ 0,     20.0 ] * u.arcsec)

##############################################################################
# Note that you could easily choose any grid dimensions:
# e.g. [ [ 100, 100, 200 ] ]
# or add more then one grid shape in the list:
# e.g. [ [ 10, 10, 10 ],[ 50, 50, 50 ], [ 100, 100, 100 ] ]
# to make the test more grid-size agnostic, but this will notably increase
# runtime.

##############################################################################
# You will create an example dataset using Gaussian spots, as show in the
# Generating Example Data example. We use the parameters:

# Manual Pole Details
arrA0 = [ u.Quantity([ 25, 25 ] * u.percent), 10.0 * u.percent,  0.2 * u.T ]
arrA1 = [ u.Quantity([ 75, 75 ] * u.percent), 10.0 * u.percent, -0.2 * u.T ]
arrA2 = [ u.Quantity([ 25, 75 ] * u.percent), 10.0 * u.percent,  0.1 * u.T ]
arrA3 = [ u.Quantity([ 75, 25 ] * u.percent), 10.0 * u.percent, -0.1 * u.T ]

# Generate the datasets and maps
# lis_maps = []
# lis_extrapolators = []

##############################################################################
# You will create an astropy table to store the runtimes of the extrapolations.

# A table for storing the data
t = Table(names=('grid size', 'time (min)', 'time (ave)', 'time (std)'), meta={'name': 'times tables'}, dtype=('S24', 'f8', 'f8', 'f8'))
t['time (min)'].unit = u.s
t['time (ave)'].unit = u.s
t['time (std)'].unit = u.s

##############################################################################
# You will store all the datasets to test with in a list.
# In this case the datasets will simply be the various generated example
# boundary data maps for the list of grid sizes, which is simply one example.
lis_datasets = []
for shape in lis_grid_shapes:
    lis_datasets.append([ str(shape), shape[2], zrange,
                          dummyDataToMap(generate_example_data(shape[0:2], xrange, yrange, arrA0, arrA1, arrA2, arrA3), xrange, yrange) ])

##############################################################################
# You may wish to run each test more than once, so you can use a parameter to
# autimate this.
int_trials = 1 # The times to repeat each extrapolation.

##############################################################################
# You iterate through the extrapolations on each dataset, adding teh runtime to
# the table.
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

##############################################################################
# You can now see the results in the table.
print t
