# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 12:59:44 2015

@author: alex_
"""

import numpy as np


def Gn_5_2_29(x, y, z, xP, yP, DxDy_val, z_submerge):
    """
    Discrete Greens Function
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


def phi_extrapolation_python(boundary, shape, Dx, Dy, Dz):
    """
    Function to extrapolate the scalar magnetic field above the given boundary
    data.
    This implementation runs in python and so is very slow for larger datasets.

    Parameters
    ----------
    boundary : array-like
        Magnetogram boundary data
    shape : array-like
        Dimensions of of the extrapolated volume, (nx,ny,nz)
    Dx : `float`
        Spacing in x-direction, in units of the boundary map
    Dy : `float`
        Spacing in y-direction, in units of the boundary map
    Dz : `float`
        Spacing in z-direction, in chosen units
    """

    # Derived parameters
    DxDy = Dx * Dy

    # From Sakurai 1982 P306, we submerge the monopole
    z_submerge = Dz / np.sqrt(2.0 * np.pi)

    # Create the empty numpy volume array.
    D = np.empty((shape[1], shape[0], shape[2]), dtype=np.float)

    i_prime, j_prime = np.indices((shape[1], shape[0]))
    xP = i_prime * Dx
    yP = j_prime * Dy

    # Iterate though the 3D space.
    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            for k in range(0, shape[2]):
                # Position of point in 3D space
                x = i * Dx
                y = j * Dy
                z = k * Dz

                # Variable holding running total for the contributions to point.
                point_phi_sum = 0.0

                G_n = Gn_5_2_29(x, y, z, xP, yP, DxDy, z_submerge)

                # Now add this to the 3D grid.
                D[j, i, k] = np.sum(boundary * G_n * DxDy)
    return D
