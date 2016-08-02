# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 12:59:44 2015

@author: alex_
"""

import numpy as np

# Imports for numba JIT compilation
from numba import jit


def phi_extrapolation_numba(boundary, shape, Dx, Dy, Dz):
    """
    Function to extrapolate the scalar magnetic field above the given boundary
    data.
    This implementation runs in python and so is very slow for larger datasets.
    """

    # Create the empty numpy volume array.
    D = np.empty((shape[1], shape[0], shape[2]), dtype=np.float)

    D = outer_loop(D, Dx, Dy, Dz, boundary)

    return D


@jit(nopython=True)
def outer_loop(D, Dx, Dy, Dz, boundary):
    shape = D.shape
    # From Sakurai 1982 P306, we submerge the monopole
    z_submerge = Dz / np.sqrt(2.0 * np.pi)
    # Iterate though the 3D space.
    for i in range(0, shape[1]):
        for j in range(0, shape[0]):
            for k in range(0, shape[2]):
                # Position of point in 3D space
                x = i * Dx
                y = j * Dy
                z = k * Dz

                # Now add this to the 3D grid.
                D[j, i, k] = inner_loop(shape, Dx, Dy, x, y, z, boundary, z_submerge)
    return D


@jit(nopython=True)
def inner_loop(shape, Dx, Dy, x, y, z, boundary, z_submerge):
    DxDy = Dx * Dy
    # Variable holding running total for the contributions to point.
    point_phi_sum = 0.0
    # Iterate through the boundary data.
    for i_prime in range(0, shape[1]):
        for j_prime in range(0, shape[0]):
            # Position of contributing point on 2D boundary
            xP = i_prime * Dx
            yP = j_prime * Dy

            # Find the components for this contribution product
            B_n = boundary[j_prime, i_prime]
            G_n = Gn_5_2_29(x, y, z, xP, yP, DxDy, z_submerge)

            # Add the contributions
            point_phi_sum += B_n * G_n * DxDy

    return point_phi_sum


@jit(nopython=True)
def Gn_5_2_29(x, y, z, xP, yP, DxDy_val, z_submerge):
    """
    Discrete Greens Function
    Extends _Gn_5_2_26 by taking the starting position of each magnetic
    monopole as 1/root(2 pi) z grid cells below the surface. (as described
    in Sakurai 1982)
    This implementation runs using Anaconda numba JIT compilation to speed up
    the process.
    """
    d_i = x - xP
    d_j = y - yP
    d_k = z - z_submerge
    floModDr = np.sqrt(d_i * d_i + d_j * d_j + d_k * d_k)

    floOut = 1.0 / (2.0 * np.pi * floModDr)
    return floOut
