# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 15:38:51 2015

Function for creating dummy boundary map datas for use with extrapolator
routines.

@author: alex_
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import random
import sunpy.map as mp
import re
from astropy import units as u
from datetime import datetime

# Function to generate the grid with Gaussian points.
# Arguments are:
#    - in_arr_area: 2 tuple for the x and y dimensions.
#    - *argv: manual parameters for all the spots. Optional: defaults to 2 random spots.
def generate_example_data(shape, xrange, yrange, *argv):
    """
    A function to generate a 2D numpy.array of example data for testing
    extrapolation code.
    The result is a mid-value region with a number of gausian spots with
    positive/negative values.
    The gausians can be specifially defined, or randomly generated.

    Parameters
    ----------

    shape : list
        A list of the axis grid sizes, (nx,ny).

    xrange : astropy.units.Quantity
        The xrange for the returned dataset.

    yrange : astropy.units.Quantity
        The yrange for the returned dataset.

    *argv : int or list, optional
        Either given the integer number of the number of poles to randomly
        generate, which defaults to 2.
        Otherwise, the user can put in lists of parameters that define a pole.
        Each list contains:
        position : astropy.units.Quantity
            both x and y coordinates as physical or percentage units
        sigma : astropy.units.Quantity
            spot size as physical or percentage units
        max : astropy.units.Quantity
            the maximum spot intensity


    """
    # If the list is empty then create random data.
    arr_args = []
    if not argv: # If no particle parameters or numbers were given.
        arr_args = [2] # [ random.randrange(1, 6) ]
    else:
        arr_args = list(argv)
    arr_poles = []

    # If we are only given the number, then generate randomly.
    if isinstance( arr_args[0], ( int, long ) ):
        for pole in range(0, arr_args[0]):
            # random parameters in percentage
            sigma = random.uniform(2, 15) * u.percent
            x_pos = random.uniform(2.0 * sigma.value, 100.0 - 2.0 * sigma.value)
            y_pos = random.uniform(2.0 * sigma.value, 100.0 - 2.0 * sigma.value)
            An_max = random.uniform(0.1, 0.2) * ((float(pole % 2) * 2.0) - 1) * u.T # Alternate pos/neg

            arrPole = [ u.Quantity([x_pos, y_pos] * u.percent), sigma, An_max ]
            arr_poles.append(arrPole)

    else:
        # We are given the hard-coded parameters, so use them.
        arr_poles = arr_args

    # Build the empty data array
    arr_data = np.zeros((shape[1], shape[0]))

    # Grid pixel shape
    qua_pixel = u.Quantity([ ( xrange[1] - xrange[0] ) / shape[0], ( yrange[1] - yrange[0] ) / shape[1] ])

    # Convert percentage positions/sigmas to physical units (units from ranges)
    for pole in range(0, len(arr_poles)):
        if arr_poles[pole][0].unit is u.percent:
            position = u.Quantity([ (arr_poles[pole][0][0].value / 100.0) * (xrange[1] - xrange[0]) + xrange[0],
                                    (arr_poles[pole][0][1].value / 100.0) * (yrange[1] - yrange[0]) + yrange[0] ])
            arr_poles[pole] = [ position, arr_poles[pole][1], arr_poles[pole][2] ]
        if arr_poles[pole][1].unit is u.percent:
            sigma = (arr_poles[pole][1].value / 100.0) * (xrange[1] - xrange[0])
            arr_poles[pole] = [ arr_poles[pole][0], sigma, arr_poles[pole][2] ]


    # Iterate through the 2D array/matrix.
    for i in range(0,shape[0]):     # Row/Y
        for j in range(0,shape[1]): # Column/X
            # The current position
            floXPrime = i * qua_pixel[0]
            floYPrime = j * qua_pixel[1]

            # A variable to store the sum of the magnetic fields for this point.
            flo_value = 0.0

            # Add all the contributions.
            for tupPole in arr_poles:
                # A0 (positive) and A1 (negative) parameters
                An_max   = tupPole[2].value
                An_x     = tupPole[0][0]
                An_y     = tupPole[0][1]
                An_Dx    = floXPrime - An_x + xrange[0]
                An_Dy    = floYPrime - An_y + yrange[0]
                An_DxSqu = np.power(An_Dx.value, 2.0)
                An_DySqu = np.power(An_Dy.value, 2.0)
                An_Sigma = tupPole[1].value

                # So this contibution is calculated and added.
                flo_An_cont = An_max * math.exp( - ( (An_DxSqu + An_DySqu) / (2 * np.power(An_Sigma, 2.0)) ))
                flo_value += flo_An_cont

            # Now add this to the data array.
            arr_data[j][i] = flo_value

    # Now return the 2D numpy array.
    return arr_data

# A function that creates a dummy header and saves the input as a fits file.
def dummyDataToMap(data, xrange, yrange, **kwargs):
    """
    Basic function for taking generated data and returning a valid sunpy.map.
    """
    # The kwargs
    dic_user_def_meta = kwargs.get('meta', {})

    # Create a header dictionary.
    dicHeader = {
                  't_obs':    datetime.now().isoformat(),
                  'bunit':    'Tesla', #'Gauss',
                  'bitpix':   64, #re.search('\\d+', 'float64')[0],#64, # Automatic
                  'naxis':    2,  # Automatic
                  'naxis1':   data.shape[1],  # Automatic
                  'naxis2':   data.shape[0],  # Automatic
                  'cdelt1':   (xrange[1].value - xrange[0].value) / data.shape[1],  # 0.504295,
                  'cdelt2':   (yrange[1].value - yrange[0].value) / data.shape[0],
                  'cunit1':   str(xrange.unit), #'arcsec',
                  'cunit2':   str(yrange.unit), #'arcsec',
                  'crpix1':   data.shape[1] / 2.0 + 0.5, # central x-pixel.
                  'crpix2':   data.shape[0] / 2.0 + 0.5, # cnetral y-pixel.
                  'rsun_ref': 696000000,
                  'dsun_ref': 149597870691,
                  'datamax':  data.max(),
                  'datamin':  data.min(),
                  'datavals': data.shape[0] * data.shape[1],
                  'CRVAL1':   (xrange[0].value + xrange[1].value)/2.0, #0.000000,
                  'CRVAL2':   (yrange[0].value + yrange[1].value)/2.0
    }

    # Add the user defined meta entries
    for key, value in dic_user_def_meta.iteritems():
        dicHeader[key] = value
        #print str(key) + ': ' + str(value)

    # Create and return a sunpy map from the data
    return mp.Map((data, dicHeader))

if __name__ == '__main__':
    # Generate an example map
    # The input parameters:
    arr_grid_shape = [ 20, 22 ]         # [ y-size, x-size ]
    qua_xrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
    qua_yrange = u.Quantity([ -11.0, 11.0 ] * u.arcsec)

    # Manual Pole Details
    #arrA0 = [ u.Quantity([ 1.0, 1.0 ] * u.arcsec), 2.0 * u.arcsec, 0.2 * u.T ]
    arrA0 = [ u.Quantity([ 25, 25 ] * u.percent), 10.0 * u.percent,  0.2 * u.T ]
    arrA1 = [ u.Quantity([ 75, 75 ] * u.percent), 10.0 * u.percent, -0.2 * u.T ]

    # Generate the data and save to a map
    arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange, arrA0, arrA1)#, arrA0, arrA1)
    #arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange)#, arrA0, arrA1)
    aMap = dummyDataToMap(arr_Data, qua_xrange, qua_yrange)
    aMap.save('C://fits//temp6.fits')
