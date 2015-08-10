# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 13:01:54 2015

@author: Alex
"""

import numpy as np
import sunpy.map as mp

# Module Imports
from classes import *
from utilities import *
from example_data_generator import *
from visualisation_functions import *


__all__ = ['PotentialExtrapolator']

class PotentialExtrapolator(Extrapolators):
    """
    This is a greens function for extrapolating the potential (scalar) field
    above a given magnetogram.
    Equations are from the following book:
    
    |     Title:      Physics of the Solar Corona
    |     Author:     T. J. M. Boyd and J. J. Sanderson
    |     Publisher:  Springer Books and Praxis Publishing
    |     ISBN:       978-3-540-30766-2

    See chapter 5 on potential fields.
    Which references to the paper Takashi Sakurai 1982:
        http://adsabs.harvard.edu/full/1982SoPh...76..301S
    
    """
    def __init__(self, map_magnetogram, **kwargs):
        super(PotentialExtrapolator, self).__init__(map_magnetogram, **kwargs)
        self.meta['extrapolator_routine'] = 'Potential Field Extrapolator'
        
        # Convert the map to SI units. (Add to extrapolator class API???)
        self.map_boundary_data = si_this_map(self.map_boundary_data)
        
        # More specific parameters (Add to extrapolator class API???)
        self.Dx = (self.xrange[1] - self.xrange[0]) / self.shape[0]
        self.Dy = (self.yrange[1] - self.yrange[0]) / self.shape[1]
        self.Dz = (self.zrange[1] - self.zrange[0]) / self.shape[2]
        

    def _extrapolation(self):
        """
        Override the primary execution method from the extrapolation class.
        The process is to extrapolate the potential (scalar) field (phi) and
        then use numerical differentiation (gradient) to find the vector field
        (Bxyz).
        """
        phi = self._extrapolate_phi()
        Bxyz = self._determine_vec(phi, D = 1, debug = False)

        return Map3D(Bxyz, self.meta, xrange=self.xrange, yrange=self.yrange, zrange=self.zrange)

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
    
    def _Gn_5_2_29(self, x, y, z, xP, yP, DxDy_val, z_submerge):#i, j, k, i_prime, j_prime, d, d_com):  
        """
        Disctete Greens Function
        Extends _Gn_5_2_26 by taking the starting position of each magnetic
        monopole as 1/root(2 pi) z grid cells below the surface. (as described
        in Sakurai 1982)
        """
        d_i = x - xP
        d_j = y - yP
        d_k = z - z_submerge
        floModDr = np.sqrt(d_i * d_i + d_j * d_j + d_k * d_k)
        
        floOut = 1.0 / (2.0 * np.pi * floModDr)
        return floOut
    
    
    def _phi_extrapolation_python(self, boundary, d):
        """
        Function to extrapolate the magnetic field above the given boundary data.
        """
        
        # Derived parameters
        Dx_val = self.Dx.value
        Dy_val = self.Dy.value
        Dz_val = self.Dz.value
        DxDy_val = Dx_val * Dy_val
        
        # From Sakurai 1982 P306, we submerge the monopole
        z_submerge = Dz_val / np.sqrt(2.0 * np.pi)
        
        # Create the empty numpy volume array.
        D = np.empty((self.shape[0], self.shape[1], self.shape[2]), dtype=np.float)
        
        # Iterate though the 3D space.
        for i in range(0, self.shape[0]):
            for j in range(0, self.shape[1]):
                for k in range(0, self.shape[2]):
                    # Position of point in 3D space
                    x = i * Dx_val
                    y = j * Dy_val
                    z = k * Dz_val
                    
                    # Variable holding running total for the contributions to point.
                    point_phi_sum = 0.0
                    
                    # Iterate through the boundary data.
                    for i_prime in range(0, self.shape[0]):
                        for j_prime in range(0, self.shape[1]):
                            # Position of contributing point on 2D boundary
                            xP = i_prime * Dx_val
                            yP = j_prime * Dy_val
                            
                            # Find the components for this contribution product
                            B_n = self.map_boundary_data.data[i_prime, j_prime]
                            G_n = self._Gn_5_2_29(x, y, z, xP, yP, DxDy_val, z_submerge)
    
                            # Add the contributions
                            point_phi_sum += B_n * G_n * DxDy_val
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
    #aMap2D = sunpy.map.Map('C://git/solarextrapolation/solarextrapolation/data/example_data_(10x10)__01_hmi.fits')
    str_folder = 'C://fits//'
    str_dataset = 'temp5'
    #str_dataset = '2011-02-14__20-35-25__01_hmi'
    str_extrapolation = str_dataset + '_3Dmap.m3d'
    str_boundary = str_dataset + '.fits'
    
    
    aMap2D = mp.Map(str_folder + str_boundary)
    
    if not os.path.isfile(str_folder+str_saved):
        aPotExt = PotentialExtrapolator(aMap2D, filepath=str_folder+str_saved, zshape=50, zrange=u.Quantity([0, 15] * u.Mm))
        aMap3D = aPotExt.extrapolate()
    aMap3D = Map3D.load(str_folder + str_extrapolation)
    #print '\n\n'
    #print aMap3D.xrange
    #print aMap3D.yrange
    #print aMap3D.zrange
    
    
    # Visualise this
    visualise(aMap3D, boundary=aMap2D, scale=1.0*u.Mm, boundary_unit=1.0*u.arcsec, show_boundary_axes=False, show_volume_axes=True, debug=False)
    
    