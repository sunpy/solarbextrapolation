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


# Function to generate the grid with Gaussian points.
# Arguments are:
#    - in_arr_area: 2 tuple for the x and y dimensions. 
#    - *argv: manual parameters for all the spots. Optional: defaults to 2 random spots.
def generate_example_data(in_arr_area, debug, *argv):
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
            arrPole[2] = random.uniform(int(in_arr_area[0]*0.05), int(in_arr_area[0]*0.1))                # Gaussian Diameter
            arrPole[0] = random.uniform(int(arrPole[2] + 1 + arrPole[2]), in_arr_area[0] - int(arrPole[2] + 1 + arrPole[2]))    # i/y-pos
            arrPole[1] = random.uniform(int(arrPole[2] + 1 + arrPole[2]), in_arr_area[0] - int(arrPole[2] + 1 + arrPole[2]))   # j/x-pos
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
    arrData = np.zeros((in_arr_area[0], in_arr_area[1]))
    
    # Iterate through the 2D array/matrix.
    for i in range(0,in_arr_area[0]):     # Row/Y
        for j in range(0,in_arr_area[1]): # Column/X
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
def dummyDataToMap(data):
    """
    Basic function for taking generated data and retuning a valid map.
    """
    # Create a header dictionary.
    dicHeader = {
                  'bunit':  'Gauss',
                  'bitpix':   64, # Automatic        
                  'naxis':    2,  # Automatic
                  'naxis1':   8,  # Automatic
                  'naxis2':   inArrData.shape[0],  # Automatic
                  'cdelt1':   0.504295,
                  'cdelt2':   0.504295,
                  'cunit1':   'arcsec',
                  'cunit2':   'arcsec',
                  'crpix1':   data.shape[1] / 2.0 + 0.5, # central x-pixel.
                  'crpix2':   data.shape[0] / 2.0 + 0.5, # cnetral y-pixel.
                  'rsun_ref': 696000000,
                  'dsun_ref': 149597870691,
                  'datamax':  inArrData.max(),
                  'datamin':  inArrData.min(),
                  'CRVAL1':   0.000000,
                  'CRVAL2':   0.000000
    }
    
    
    # Create and return a sunpy map from the data
    return mp.Map((data, dicHeader))
    
if __name__ == '__main__':
    # Generate an example map
    
    # The input parameters:
    booRandom = True
    arrArea = [ 40, 50 ]         # [ y-size, x-size ]
    arrA0 = [ 10, 18, 8.0, 0 ] # [A1_x, A1_y, sigma_0 (diameter), A1_stength]
    arrA1 = [ arrArea[0] - 15, arrArea[1] - 15, 4.0, 0 ] # [A1_x, A1_y, sigma_0 (diameter), A1_stength]
    arrA2 = [ 15, 10, 2.0, 2.0 ] # [A1_x, A1_y, sigma_0 (diameter), A1_stength]
    arrA3 = [ 5, 5, 2.0, -2.0 ] # [A1_x, A1_y, sigma_0 (diameter), A1_stength]
    
    # Generate the data and save to a map
    inArrData = generate_example_data(arrArea, False)#, arrA0, arrA1)
    aMap = dummyDataToMap(inArrData)
    aMap.save('C://fits//temp2.fits')
    
    