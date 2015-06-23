# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 19:49:49 2014

@author: Alex
"""

import numpy as np
import time


# A function that takes the 3D magnetic scalar potential field and returns the 3D vector field.
def determineVec_FORTRAN(phi, D = 1, debug = False):
    # Time this function.
    start = time.time()
    
    # Create an empty 3D matrix from the output.
    # ATM, for simplicity, I make the same size as the potential field, though the outer 2 layers are all 0.0.
    tupVolShape = phi.shape
    npmVecSpace = np.zeros((tupVolShape[0], tupVolShape[1], tupVolShape[2], 3)) # in Order XYZC (C = component directions)
    
    # For each cell we use data from 2 in each direction, this means we need to reduce the volume by 2 in eaach direction.
    for k in range(2, tupVolShape[2]-2):          # Z - Only done first so I can say when an XY slice has been rendered.
        for j in range(2, tupVolShape[1]-2):      # Y
            for i in range(2, tupVolShape[0]-2):  # X               
                # 
                npmVecSpace[i,j,k,0]=-(phi[i-2,j,k]-8.0*phi[i-1,j,k]+8.0*phi[i+1,j,k]-phi[i+2,j,k])/(12.0*D)
                npmVecSpace[i,j,k,1]=-(phi[i,j-2,k]-8.0*phi[i,j-1,k]+8.0*phi[i,j+1,k]-phi[i,j+2,k])/(12.0*D)
                npmVecSpace[i,j,k,2]=-(phi[i,j,k-2]-8.0*phi[i,j,k-1]+8.0*phi[i,j,k+1]-phi[i,j,k+2])/(12.0*D)
        if debug and k%int(tupVolShape[2]*0.1) is 0:
            print '(Bx,By,Bz) calculated for layer ' + str(k) + '.'
            
    # Print the time ellapsed if debug mode.
    if debug:
        finish = time.time()
        print 'determineVec_FORTRAN time ' + str(tupVolShape) + ': ' + str(time.strftime("%H:%M:%S", time.gmtime(finish - start)))
    return npmVecSpace
    
############
#
#  Examples
#
#######
if __name__ == '__main__':
    a3DArray = np.random.rand(5,5,5)
    # Make this array structured:
    delta = 0.5
    for k in range(0,a3DArray.shape[0]):
        for j in range(0,a3DArray.shape[1]):
            for i in range(0,a3DArray.shape[2]):
                # Increase in all directions
                a3DArray[i,j,k] = i * delta + j * delta + k * delta
                # Increase in 2D
                #a3DArray[i,j,k] = i * delta + j * delta
                # Increase in 1D
                #a3DArray[i,j,k] = i * delta
    
    # Uncomment below to make it a random number volume again
    #a3DArray = np.random.rand(5,5,5)

    # Print the scalar volume array
    print '\n\n'
    print a3DArray
    
    # Print the centre cell vector from my function
    print '\n\n'
    a4DArray = determineVec_FORTRAN(a3DArray)    
    print a4DArray[2,2,2]
    
    # Print the centre cell vector from the numpy gradient function
    b4DArray = np.gradient(a3DArray)
    print '[ ' + str(b4DArray[0][2,2,2]) + ', ' + str(b4DArray[1][2,2,2]) + ', ' + str(b4DArray[2][2,2,2]) + ']'
    
