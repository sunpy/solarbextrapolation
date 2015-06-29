# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 13:01:54 2015

@author: Alex
"""
# Import the SolarBExtrapolation API
from classes import *
import numpy as np

__all__ = ['PotentialExtrapolator']

class PotentialExtrapolator(Extrapolators):
    """
    | This is a greens function for extrapolating the potential (scalar) field
    | above a given magnetogram.
    | Exquations are from the following book:
    |     Title:      Physics of the Solar Corona
    |     Author:     T. J. M. Boyd and J. J. Sanderson
    |     Publisher:  Springer Books and Praxis Publishing
    |     ISBN:       978-3-540-30766-2
    | See chapter 5 on potential fields.
    """
    def __init__(self, map_magnetogram, **kwargs):
        super(PotentialExtrapolator, self).__init__(map_magnetogram, **kwargs)
        self.meta['extrapolator_routine'] = 'Potential Field Extrapolator'

    def _extrapolation(self):
        """
        Override the primary execution method from the extrapolation class.
        The process is to extrapolate the potential (scalar) field (phi) and
        then use numerical differentiation (gradient) to find the vector field
        (Bxyz).
        """
        phi = self._extrapolate_phi()
        Bxyz = self._determine_vec(phi, D = 1, debug = False)

        return Map3D(Bxyz, self.meta)

    # Greens function.
    def _Gn_5_2_26(self, inR, inRPrime):
        """
        Continious Greens Function
        Treats the magnetic field silimarly to the electrostatic potential and
        uses the assumption that all magnetic potential is from the Sun below.
        At any one point in space (:math:`\mathbf{r}`), we have a contribution 
        from the point on the boundary map (:math:`\mathbf{r}^\prime`) given by:
        .. math::
        G_n(\mathbf{r}, \mathbf{r}^\prime) = \frac{1}{2\pi|\mathbf{r} - \mathbf{r}^\prime|}
        Which can be used to find the total field at a point as a result of
        integrating over the whole boundary:
        
        """
        floModDr = np.linalg.norm(inR - inRPrime)
        floOut = 1.0 / (2.0 * np.pi * floModDr)
        return floOut
    
    # Greens function.
    def _Gn_5_2_29(self, i, j, k, i_prime, j_prime, d, d_com):    
        d_i = i - i_prime
        d_j = j - j_prime
        d_k = k - d_com # d / np.sqrt(2.0 * np.pi)
        floModDr = np.sqrt(d_i * d_i + d_j * d_j + d_k * d_k)
        
        #print 'floModDr: ' + str(floModDr)
        floOut = 1.0 / (2.0 * np.pi * floModDr)
        #print 'floOut: ' + str(floOut)
        return floOut
    
    
    # A function to extrapolate the magnetic field above the given boundary.
    def _phi_extrapolation_python(self, boundary, d):    
        # Volume size.
        M = boundary.shape[0]
        N = boundary.shape[1]
        Z = self.z
        
        # Derived parameters
        d = 1.0
        d_squ = np.power(d,2)
        d_com = d / np.sqrt(2.0 * np.pi)
        
        # Create the empty numpy volume array.
        D = np.empty((M, N, Z), dtype=np.float)
        
        # Itherate though the 3D space.
        for i in range(0, M):
            for j in range(0, N):
                for k in range(0, Z):
                    # Variable holding running total for the contributions to point.
                    point_phi_sum = 0.0
                    # Iterate through the boundary data.
                    for i_prime in range(0,M):
                        for j_prime in range(0,N):                        
                            # Find the components for this contribution product.
                            B_n = boundary[i_prime, j_prime]
                            
                            G_n = self._Gn_5_2_29(i, j, k, i_prime, j_prime, d, d_com)
    
                            # Add the contribution.
                            point_phi_sum += B_n * G_n * d_squ
                    # Now add this to the 3D grid.
                    D[i, j, k] = point_phi_sum
        return D
    
    def _extrapolate_phi(self, debug=False):
        """
        A function to extrapolate the magnetic field above the given boundary.
        Assumes the input B-field boundary data is near normal (the image must
        be near the centre of the HMI data).
        P183 (5.2.28)
        """
        if debug:
            print "extrapolatePhi({},{},{})".format(self.map_boundary_data.data.shape, inZ, debug)
        
        # Parameters
        arr_boundary = self.map_boundary_data.data
        d = 1.0        
    
        # Now return the volume.
        return self._phi_extrapolation_python(arr_boundary, d) #np.empty((1, 1, 1), dtype=np.float)


    def _determine_vec(self, phi, D = 1, debug = False):        
        """
        Create an empty 3D matrix from the output.
        ATM, for simplicity, I make the same size as the potential field, though the outer 2 layers are all 0.0.
        """
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

        return npmVecSpace
        
if __name__ == '__main__':
    aMap2D = sunpy.map.Map('C://git/solarextrapolation/solarextrapolation/data/example_data_(10x10)__01_hmi.fits')
    aPotExt = PotentialExtrapolator(aMap2D, filepath='C://git/solarextrapolation/solarextrapolation/3Dmap.m3d')
    aMap3D = aPotExt.extrapolate()