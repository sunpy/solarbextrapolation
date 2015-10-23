# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 13:01:54 2015

@author: Alex
"""

import numpy as np
import sunpy.map as mp
import astropy.units as u

# Module Imports
#from classes import *
#from solarbextrapolation.utilities import *
from solarbextrapolation.extrapolators import Extrapolators
from solarbextrapolation.utilities import si_this_map
from solarbextrapolation.map3dclasses import Map3D
#from solarbextrapolation.visualisation_functions import visualise


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

    def _extrapolation(self, enable_numba=True, **kwargs):
        """
        Override the primary execution method from the extrapolation class.
        The process is to extrapolate the potential (scalar) field (phi) and
        then use numerical differentiation (gradient) to find the vector field
        (Bxyz).
        """

        if enable_numba:
            # Test that numba and the numba'ed extrpolator can be imported
            try:
                import numba
                from potential_field_extrapolator_numba import phi_extrapolation_numba
            except ImportError:
                enable_numba = False

        phi = self._extrapolate_phi(enable_numba, **kwargs)

        if enable_numba:
            from numba.decorators import autojit
            determine_vec = autojit(self._determine_vec)
        else:
            determine_vec = self._determine_vec

        npmVecSpace = np.zeros(list(phi.shape)+[3]) # in Order XYZC (C = component directions)
        Bxyz = determine_vec(phi, 1, npmVecSpace)

        return Map3D(Bxyz, self.meta, xrange=self.xrange, yrange=self.yrange, zrange=self.zrange)

    def _extrapolate_phi(self, enable_numba, debug=False, **kwargs):
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

        print(enable_numba)
        if enable_numba:
            from .potential_field_extrapolator_numba import phi_extrapolation_numba as phi_extrapolation
        else:
            from .potential_field_extrapolator_python import phi_extrapolation_python as phi_extrapolation

        return phi_extrapolation(arr_boundary, self.shape, self.Dx.value, self.Dy.value, self.Dz.value)


    # Make this a static method so it is more efficient to numba
    @staticmethod
    def _determine_vec(phi, D, npmVecSpace):
        """
        Create an empty 3D matrix from the output.
        ATM, for simplicity, I make the same size as the potential field, though the outer 2 layers are all 0.0.
        """
        tupVolShape = npmVecSpace.shape

        # For each cell we use data from 2 in each direction, this means we need to reduce the volume by 2 in eaach direction.
        for k in range(2, tupVolShape[2]-2):          # Z - Only done first so I can say when an XY slice has been rendered.
            for j in range(2, tupVolShape[1]-2):      # Y
                for i in range(2, tupVolShape[0]-2):  # X
                    npmVecSpace[i,j,k,0]=-(phi[i-2,j,k]-8.0*phi[i-1,j,k]+8.0*phi[i+1,j,k]-phi[i+2,j,k])/(12.0*D)
                    npmVecSpace[i,j,k,1]=-(phi[i,j-2,k]-8.0*phi[i,j-1,k]+8.0*phi[i,j+1,k]-phi[i,j+2,k])/(12.0*D)
                    npmVecSpace[i,j,k,2]=-(phi[i,j,k-2]-8.0*phi[i,j,k-1]+8.0*phi[i,j,k+1]-phi[i,j,k+2])/(12.0*D)

        return npmVecSpace

