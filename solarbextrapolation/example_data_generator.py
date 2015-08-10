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

# Function to generate the grid with Gaussian points.
# Arguments are:
#    - in_arr_area: 2 tuple for the x and y dimensions. 
#    - *argv: manual parameters for all the spots. Optional: defaults to 2 random spots.
def generate_example_data(shape, ranges, debug, *argv):
    """
    A function to generate a 2D numpy.array of example data for testing
    extrapolation code.
    The result is a mid-value region with a number of gausian spots with
    positive/negative values.
    The gausians can be specifially defined, or randomly generated.    
    """
    # If the list is empty then create random data.
    arrArgs = []
    if not argv: # If no particle parameters or numbers were given.
        arrArgs = [2] # [ random.randrange(1, 6) ]
    else:
        arrArgs = list(argv)
    arrPoles = []
    
    #If we are only given the number, then generate randomly.
    if isinstance( arrArgs[0], ( int, long ) ):
        for pole in range(0, arrArgs[0]):
            arrPole = [ 0.0, 0.0, 0.0, 0.0 ]
            arrPole[2] = random.uniform(int(shape[0]*0.05), int(shape[0]*0.1))                # Gaussian Diameter
            arrPole[0] = random.uniform(int(arrPole[2] + 1 + arrPole[2]), shape[0] - int(arrPole[2] + 1 + arrPole[2]))    # i/y-pos
            arrPole[1] = random.uniform(int(arrPole[2] + 1 + arrPole[2]), shape[0] - int(arrPole[2] + 1 + arrPole[2]))   # j/x-pos
            floPolarity = ((float(pole % 2) * 2.0) - 1)      # Alternate pos/neg
            arrPole[3] = floPolarity * random.uniform(1000, 2000)  # Max value.
            arrPoles.append(arrPole)
            
            # Print out the parameters for debugging.
            if debug:
                print str(pole) + ': radius: ' + str(arrPole[2]) + ', pos: (' + str(arrPole[0]) + ', ' + str(arrPole[1]) + '), max: ' + str(arrPole[3])
    else:
        # We are given the hard-coded parameters, so use them.
        arrPoles = arrArgs

    # Print out the parameters for the points for debugging.
    for pole in range(0, len(arrArgs)):
        if debug:
            arrPole = arrArgs[pole]
            print str(pole) + ': radius: ' + str(arrPole[2]) + ', pos: (' + str(arrPole[0]) + ', ' + str(arrPole[1]) + '), max: ' + str(arrPole[3])

    # Build the empty data array
    arrData = np.zeros((shape[0], shape[1]))
    
    # Iterate through the 2D array/matrix.
    for i in range(0,shape[0]):     # Row/Y
        for j in range(0,shape[1]): # Column/X
            # The current position           
            floXPrime = i
            floYPrime = j
            
            # A variable to store the sum of the magnetic fields for this point.
            floValue = 0.0
            
            # Add all the contributions.
            for tupPole in arrPoles:
                # A0 (positive) and A1 (negative) parameters
                An = tupPole[3]
                An_x = tupPole[0]
                An_y = tupPole[1]
                An_Dx = floXPrime - An_x
                An_Dy = floYPrime - An_y
                An_DxSqu = math.pow(An_Dx, 2)
                An_DySqu = math.pow(An_Dy, 2)
                floSigma_n = tupPole[2]
                
                # So this contibution is calculated and added.
                floAnCont = An * math.exp( - ( (An_DxSqu + An_DySqu) / (2 * math.pow(floSigma_n, 2)) ))
                floValue += floAnCont
            
            # Now add this to the data array.
            arrData[i][j] = floValue
            
    # Now return the 2D numpy array.
    return arrData
    
# A function that creates a dummy header and saves the input as a fits file.
def dummyDataToMap(data, ranges):
    """
    Basic function for taking generated data and retuning a valid map.
    """
    # Create a header dictionary.
    dicHeader = {
                  'bunit':  'Gauss',
                  'bitpix':   64, #re.search('\\d+', 'float64')[0],#64, # Automatic        
                  'naxis':    2,  # Automatic
                  'naxis1':   8,  # Automatic
                  'naxis2':   inArrData.shape[0],  # Automatic
                  'cdelt1':   (ranges[1].value - ranges[0].value) / data.shape[0],  # 0.504295,
                  'cdelt2':   (ranges[3].value - ranges[2].value) / data.shape[1],
                  'cunit1':   str(ranges.unit), #'arcsec',
                  'cunit2':   str(ranges.unit), #'arcsec',
                  'crpix1':   data.shape[1] / 2.0 + 0.5, # central x-pixel.
                  'crpix2':   data.shape[0] / 2.0 + 0.5, # cnetral y-pixel.
                  'rsun_ref': 696000000,
                  'dsun_ref': 149597870691,
                  'datamax':  inArrData.max(),
                  'datamin':  inArrData.min(),
                  'CRVAL1':   (ranges[0].value + ranges[1].value)/2.0, #0.000000,
                  'CRVAL2':   (ranges[2].value + ranges[3].value)/2.0
    }
    
    
    # Create and return a sunpy map from the data
    return mp.Map((data, dicHeader))
    
if __name__ == '__main__':
    # Generate an example map
    
    # The input parameters:
    booRandom = True
    arr_grid_shape = [ 20, 22 ]         # [ y-size, x-size ]
    qua_grid_ranges = u.Quantity([ -10.0, 10.0, -11.0, 11.0 ] * u.arcsec)
    
    # Manual Pole Details
    arrA0 = [ 10, 18, 8.0, 0 ] # [A1_x, A1_y, sigma_0 (diameter), A1_stength]
    arrA1 = [ arr_grid_shape[0] - 10, arr_grid_shape[1] - 10, 4.0, 0 ] # [A1_x, A1_y, sigma_0 (diameter), A1_stength]
    arrA2 = [ 15, 10, 2.0, 2.0 ] # [A1_x, A1_y, sigma_0 (diameter), A1_stength]
    arrA3 = [ 5, 5, 2.0, -2.0 ] # [A1_x, A1_y, sigma_0 (diameter), A1_stength]

    # Generate the data and save to a map
    inArrData = generate_example_data(arr_grid_shape, qua_grid_ranges, False)#, arrA0, arrA1)
    aMap = dummyDataToMap(inArrData, qua_grid_ranges)
    aMap.save('C://fits//temp5.fits')
    
    