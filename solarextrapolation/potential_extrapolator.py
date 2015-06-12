# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 13:01:54 2015

@author: Alex
"""
# Import the SolarBExtrapolation API
from classes import *

class PotentialExtrapolator(Extrapolators):
    def __init__(self, map_magnetogram, **kwargs):
        super(PotentialExtrapolator, self).__init__(map_magnetogram, **kwargs)
        self.meta['extrapolator_routine'] = 'Potential Field Extrapolator'

    def _extrapolation(self):
        # Adding in custom parameters to the meta
        phi = self._extrapolate_phi()
        Bxyz = self._determine_vec(phi, D = 1, debug = False)

        return Map3D(( Bxyz, self.meta ))

    # Greens function.
    def _Gn_5_2_26(self, inR, inRPrime):    
        floModDr = np.linalg.norm(inR - inRPrime)
        floOut = 1.0 / ( 2.0 * math.pi * floModDr)
        return floOut
    
    # Greens function.
    def _Gn_5_2_29(self, i, j, k, i_prime, j_prime, d, d_com):    
        d_i = i - i_prime
        d_j = j - j_prime
        d_k = k - d_com # d / np.sqrt(2.0 * math.pi)
        floModDr = np.sqrt( d_i * d_i + d_j * d_j + d_k * d_k)
        
        #print 'floModDr: ' + str(floModDr)
        floOut = 1.0 / ( 2.0 * math.pi * floModDr)
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
        d_squ = math.pow(d,2)
        d_com = d / math.sqrt(2.0 * math.pi)
        
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
    
    # A function to extrapolate the magnetic field above the given boundary.
    # Assumes the input B-field boundary data is near normal (the image must be
    # near the centre of the HMI data).
    # P183 (5.2.28)
    def _extrapolate_phi(self, debug=False):
        if debug:
            print 'extrapolatePhi(' + str(self.map_boundary_data.data.shape) + ', ' + str(inZ) + ', ' + str(debug) + ')'
        
        # Parameters
        arr_boundary = self.map_boundary_data.data
        d = 1.0        
    
        # refuse to do an extrapolation with any dimension > 300.
        if arr_boundary.shape[0] > 300 or arr_boundary.shape[0] > 300 or arr_boundary.shape[0] > 300:
            print 'extrapolatePhi(' + str(arr_boundary.shape) + ', ' + str(self.z) + ', ' + str(debug) + ')'
            print 'Dimensions to large.'
            return np.array([])
        
        # Now return the volume.
        return self._phi_extrapolation_python(arr_boundary, d) #np.empty((1, 1, 1), dtype=np.float)


    def _determine_vec(self, phi, D = 1, debug = False):
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
            print '_determine_vec time ' + str(tupVolShape) + ': ' + str(time.strftime("%H:%M:%S", time.gmtime(finish - start)))
        return npmVecSpace
        
if __name__ == '__main__':
    aMap2D = sunpy.map.Map('C:/Users/Alex/solarextrapolation/solarextrapolation/data/example_data_(10x10)__01_hmi.fits')
    aPotExt = PotentialExtrapolator(aMap2D, filepath='C://Users/Alex/solarextrapolation/solarextrapolation/3Dmap.m3d')
    aMap3D = aPotExt.extrapolate()