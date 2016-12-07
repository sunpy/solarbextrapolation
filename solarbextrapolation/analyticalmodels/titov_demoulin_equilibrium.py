# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 00:45:30 2015

This code is intended to implement an analytical solution for a flux loop.
The Titov Demoulin Equilibrium is from the paper:
    Basic topology of twisted magnetic configurations in solar flares
        V.S. Titov 1 and P. Démoulin
            1999

The magnetic field over the region is created as the sum of 3 components:
  B_I:     Field over circular area a due to current I
  B_q:     Field created by the active region spots +q and -q that lie on I_0
  B_theta: Field created by flow I_0
With:
  I_0: the line of sub-surface current going allong the x-axis.
With parameters:
  d: depth of I_0 below the photosphere.
  L: distances of q_plus and q_minus from the centre of I_0.
  R: radii of the large curcle of current I about I_0.

Useful numpy functions:
np.add(x,y), np.add(x,y): element-wise addition/subtraction.

@author: alex_
"""

import numpy as np
import math
import scipy
from scipy import interpolate
from astropy import units as u


# Imports for numba JIT compilation
from numba import double
from numba.decorators import jit, autojit

# Module Imports
###from solarbextrapolation import *
from classes import *
from utilities import *
from example_data_generator import *
from visualisation_functions import *

# Universal values.
from scipy.constants import  mu_0, pi, au
# mu_0 = 1.256637061 * 10^-6 # T m A^-1


# Loop Configuration Parameters
qua_TD_L = 50.0*10**6 * u.m  # m  # Distance of +q/-q from centre of volume.
qua_TD_d = 50.0*10**6 * u.m  # m  # Depth of I_0 below photosphere.
qua_TD_R = 85.0*10**6 * u.m # m  # Radius of circle of current I (makes area a)
qua_TD_q = 100.0 * 10**12 * u.m * u.m * u.T  # T m^2  # ABS Charge of +q and -q.
qua_TD_a = 31.0*10**6 * u.m  # m      # Radius of uniform current I.

qua_TD_I_0 = -7.0*10**12 * u.A # A # - 7.0 TA #
flo_TD_li = 1.0 / 2.0 # for uniform distribution of current over toroidal flux tube.

# Convert all these into SI units
flo_TD_L   = qua_TD_L.to(u.m).value
flo_TD_d   = qua_TD_d.to(u.m).value
flo_TD_R   = qua_TD_R.to(u.m).value
flo_TD_q   = qua_TD_q.value # This doesn't convert into SI units
flo_TD_a   = qua_TD_a.to(u.m).value
flo_TD_I_0 = qua_TD_I_0.to(u.A).value
flo_TD_I   = 8.0 * pi * flo_TD_q * flo_TD_L * flo_TD_R * (flo_TD_R**2.0 + flo_TD_L**2)**(-3.0/2.0) / ( mu_0 * (np.log(8.0 * flo_TD_R / flo_TD_a) - (3.0/2.0) + (flo_TD_li / 2.0)) )


"""
TD_L = 50.0*10**6  # m  # Distance of +q/-q from centre of volume.
TD_d = 50.0*10**6  # m  # Depth of I_0 below photosphere.
TD_R = 85.0*10**6 # m  # Radius of circle of current I (makes area a)
TD_q = 100.0*10**12   # T m^2  # ABS Charge of +q and -q.
TD_a = 31.0*10**6  # m      # Radius of uniform current I.

TD_I_0 = -7.0*10**12 # A # - 7.0 TA #
TD_li = 1.0/2.0 # for uniform distribution of current over toroidal flux tube.
TD_I = 8.0 * pi * TD_q * TD_L * TD_R * (TD_R**2.0 + TD_L**2)**(-3.0/2.0) / ( mu_0 * (np.log(8.0 * TD_R / TD_a) - (3.0/2.0) + (TD_li / 2.0)) ) # 11 # 11000 GA # Equilibrium


# Unit strings, for filenames
TD_L_unit, TD_d_unit, TD_R_unit, TD_a_unit = 'm', 'm', 'm', 'm'
TD_q_unit = 'Tm2'
TD_I_0_unit = 'A'
"""


###############################################################################
#                                                                             #
#                              B_theta components                             #
#                                                                             #
###############################################################################


# The Heaviside chi(X) function
# Returns:
#   chi(X) = 1  if  X > 0
#   chi(X) = 1 otherwise
#   else  returns 1 and prints an error message in the terminal.
def chi(X, debug = 0):
    if X > 0.0:
        return 1.0;
    else:
        return 0.0;

# Safe Heaviside chi(X) function
# Returns:
#   chi_safe(X, val) = val  if  X > 0
#   chi_safe(X, val) = 0 otherwise
# This is designed to correct for errors where you have X <= 0 and val = NaN,
# which causes:
#   chi(X) * val = NaN
# Rather then the prefered:
#   chi(X) * val = 0.0
def chi_safe(X, val, debug = 0):
    out = 0.0
    if X > 0.0:
        out = val

    # Debugging
    if debug > 0:
        print('chi_safe(' + str(X) + '): ' + str(out))
    # Output
    return out

# r_perpendicular
def r_perpendicular_scalar(y, z, d, debug = 0):
    # Output
    out = ( y**2.0 + (z + d)**2.0 )**(0.5)

    # Debugging
    if debug > 0:
        print('r_perpendicular_scalar: ' + str(out))
    # Output
    return out

#
def rho(x, y, z, R, d, debug = 0):
    # Output
    out = ( x**2.0 + (r_perpendicular_scalar(y, z, d, debug - 2) - R)**2.0 )**(0.5)

    # Debugging
    if debug > 0:
        print('rho: ' + str(out))
    # Output
    return out

# theta_hat function
def theta_hat(y, z, d, debug = 0):
    # r_perp
    r_perp = r_perpendicular_scalar(y, z, d, debug - 2)

    # Components
    theta_hat_x = 0.0
    theta_hat_y = - ( ( z + d ) / r_perp )
    theta_hat_z = y / r_perp

    # Output
    out = np.array([ theta_hat_x, theta_hat_y, theta_hat_z ])

    # Debugging
    if debug > 0:
        print('theta_hat: ' + str(out))
    # Output
    return out


# Note, independent of x???
def B_theta(x, y, z, a, d, R, I, I_0, debug = 0):
    # Find rho
    rho_val = rho(x, y, z, R, d, debug - 2)

    # The parts.
    part_1  = ( mu_0 * I_0 ) / (2.0 * pi )
    part_2a = R**(-2.0)
    part_2b = (2.0 ) / ( a**2.0 )
    part_2c = (I**2.0) / (I_0**2.0)
    part_2d = 1.0 - ((rho_val**2.0)/(a**2.0))
    part_2  = (part_2a + chi_safe(a - rho_val, part_2b * part_2c * part_2d))**(0.5)
    part_3  = (y**2.0 + (z + d)**2.0)**(-0.5)
    part_4  = R**(-1.0)

    # Now put it together.
    scalar = part_1 * (part_2 + part_3 - part_4)
    vector = theta_hat(y, z, d)

    # Output
    out = np.array([scalar * vector[0], scalar * vector[1], scalar * vector[2]])

    # Debugging
    if debug > 0:
        print('B_theta: ' + str(out))
    if debug > 1:
        print ('  B_theta  part_1: ' + str(part_1))
        print ('  B_theta  part_2: ' + str(part_2))
        print ('  B_theta    part_2a: ' + str(part_2a))
        print ('  B_theta    part_2b: ' + str(part_2b))
        print ('  B_theta    part_2c: ' + str(part_2c))
        print ('  B_theta    part_2d: ' + str(part_2d))
        print ('  B_theta  part_3: ' + str(part_3))
        print ('  B_theta  part_4: ' + str(part_4))
        print ('  B_theta  scalar: ' + str(scalar))
        print ('  B_theta  vector: ' + str(vector) + '\n')
    # Output
    return out


###############################################################################
#                                                                             #
#                                B_q components                               #
#                                                                             #
###############################################################################

# sign is either +1.0 or -1.0
def r_plusminus(x, y, z, L, d, sign, debug = 0):
    # Components
    r_x = x + - sign * L
    r_y = y * 1.0
    r_z = z + d

    # Output B_q vector.
    out = np.array([ r_x, r_y, r_z ])

    # Debugging
    if debug > 0:
        print('r_plusminus(..., ' + str(sign) + '): ' + str(out))
    # Output
    return out

# Returns the B-q vector (numpy array) at given x, y, z.
def B_q(x, y, z, L, d, q, debug = 0):
    # Getting the r+- vectors
    r_plus = r_plusminus(x, y, z, L, d, 1.0, debug - 2)
    r_minus = r_plusminus(x, y, z, L, d, -1.0, debug - 2)

    # Get the modulus of these
    mod_r_plus = (r_plus[0]**2.0 + r_plus[1]**2.0 + r_plus[2]**2.0)**(0.5)
    mod_r_minus = (r_minus[0]**2.0 + r_minus[1]**2.0 + r_minus[2]**2.0)**(0.5)

    # Get the two fractions form (20)
    frac_r_plus = r_plus / (mod_r_plus**3.0)
    frac_r_minus = r_minus / (mod_r_minus**3.0)

    # Output B_q vector.
    out = q * np.subtract(frac_r_plus, frac_r_minus)

    # Debugging
    if debug > 0:
        print('B_q: ' + str(out))
    if debug > 1:
        print('  B_q  r_plus: ' + str(r_plus))
        print('  B_q  r_minus: ' + str(r_minus))
        print('  B_q  mod_r_plus: ' + str(mod_r_plus))
        print('  B_q  mod_r_minus: ' + str(mod_r_minus))
        print('  B_q  frac_r_plus: ' + str(frac_r_plus))
        print('  B_q  frac_r_minus: ' + str(frac_r_minus) + '\n')
    # Output
    return out




###############################################################################
#                                                                             #
#                                B_I components                               #
#                                                                             #
###############################################################################

# Note, we use the scipy complete elliptic integrals of the 1st and 2nd kind
#   scipy.special.ellipk(flo_m) # first kind
#   scipy.special.ellipe(flo_m) # second kind


# r_perpendicular as a vector
def r_perpendicular_vector(y, z, d, debug = 0):
    # Output
    out = np.array([ 0.0, (z + d)*1.0, y*1.0])

    # Debugging
    if debug > 0:
        print('r_perpendicular_vector: ' + str(out))
    # Output
    return out


def k_func(x, y, z, d, R, debug = 0):
    # Parameters
    r_perp = r_perpendicular_scalar(y, z, d, debug - 2)

    # fraction inside root
    frac = (r_perp * R)/((r_perp + R)**2.0 + x**2.0)

    # Output
    out = 2 * (frac)**(0.5)

    # Debugging
    if debug > 0:
        print('k_a_func: ' + str(out))
    if debug > 1:
        print('  k_func  r_perp: ' + str(r_perp))
        print('  k_func  frac: ' + str(frac))
    # Output
    return out


def k_a_func(y, z, d, R, a, debug = 0):
    # Parameters
    r_perp = r_perpendicular_scalar(y, z, d, debug - 2)

    # fraction inside root
    frac = (r_perp * R)/(4.0 * r_perp * R + a**2.0)

    # Output
    out = 2 * (frac)**(0.5)

    # Debugging
    if debug > 0:
        print('k_a_func: ' + str(out))
    if debug > 1:
        print('  k_a_func  r_perp: ' + str(r_perp))
        print('  k_a_func  frac: ' + str(frac))
    # Output
    return out


def A_of_k(k, debug = 0):
    # Output
    out = (k**(-1.0))*((2.0 - k**2.0) * scipy.special.ellipk(k))

    # Debugging
    if debug > 0:
        print('A_of_k(' + str(k) + '): ' + str(out))
    # Output
    return out

def A_prime_of_k(k, debug = 0):
    # Building the parts
    numerator_1 = (2.0 - k**2.0) * scipy.special.ellipe(k)
    numerator_2 = 2.0 * (1.0 - k**2.0) * scipy.special.ellipk(k)
    denominator = k**2.0 * (1.0 - k**2.0)

    # Output
    out = (numerator_1 - numerator_2) / denominator

    # Debugging
    if debug > 0:
        print ('A_prime_of_k: ' + str(out))
    if debug > 1:
        print('  A_prime_of_k  numerator_1: ' + str(numerator_1))
        print('  A_prime_of_k  numerator_2: ' + str(numerator_2))
        print('  A_prime_of_k  denominator: ' + str(denominator))
    # Output
    return out

def A_tilde_I_in(k, k_a, r_perp, R, I, debug = 0):
    # Building the parts
    part_1 = (mu_0 * I)/(2.0 * pi)
    part_2 = ((R)/(r_perp))**(0.5)
    part_3 = A_of_k(k_a, debug - 2) + A_prime_of_k(k_a, debug - 2) * (k - k_a)

    # Output
    out = part_1 * part_2 * part_3

    # Debugging
    if debug > 0:
        print('A_tilde_I_in: ' + str(out))
    if debug > 1:
        print('  A_tilde_I_in  part_1: ' + str(part_1))
        print('  A_tilde_I_in  part_2: ' + str(part_2))
        print('  A_tilde_I_in  part_3: ' + str(part_3))
    # Output
    return out

def A_I_ex(k, r_perp, R, I, debug = 0):
    # Building the parts
    part_1 = (mu_0 * I)/(2.0 * pi)
    part_2 = ((R)/(r_perp))**(0.5)
    part_3 = A_of_k(k)

    # Output
    out = part_1 * part_2 * part_3

    # Debugging
    if debug > 0:
        print('A_I_ex: ' + str(out))
    if debug > 1:
        print('  A_I_ex  part_1: ' + str(part_1))
        print('  A_I_ex  part_2: ' + str(part_2))
        print('  A_I_ex  part_3: ' + str(part_3))
    # Output
    return out

def A_I(x, y, z, R, a, d, I, debug = 0):
    # Values
    rho_val = rho(x, y, z, R, d, debug - 2)
    k = k_func(x, y, z, d, R, debug - 2)
    k_a = k_a_func(y, z, d, R, a, debug - 2)
    r_perp = r_perpendicular_scalar(y, z, d, debug - 2)

    # Parts
    part_1 = chi_safe(a - rho_val, A_tilde_I_in(k, k_a, r_perp, R, I))
    part_2 = chi_safe(rho_val - a, A_I_ex(k, r_perp, R, I))
    out = part_1 + part_2

    # Debugging
    if debug > 0:
        print('A_I: ' + str(out))
    if debug > 1:
        print('  A_I  part_1: ' + str(part_1))
        print('  A_I  part_2: ' + str(part_2))
    # Output
    return out

# function to return an interpolator object for A_I_from_r_perp.
def interpolate_A_I_from_r_perp(R, a, d, I, r_perp_max, resolution = 100000, debug = 0):
    # parameters to pass in:
    dr_perp = 1.0 * r_perp_max / resolution

    # 1D array of vectors for A_I(r_perp)
    npm_A_I = np.zeros((resolution, 2))

    # If we lock x = 0, z = -d then we know r_perp = y
    x = 0.0
    z = - d
    dy = dr_perp
    for i in range(1, resolution): # Can't start at 0.
        y = i * dy
        npm_A_I[i][0] = y
        npm_A_I[i][1] = A_I(x, y, z, R, a, d, I)
    # Make the first row very close to 0.
    npm_A_I[0][0] = 0.001
    npm_A_I[i][1] = A_I(x, npm_A_I[0][0], z, R, a, d, I)

    # Make/return the interpolation object.
    interpolator = scipy.interpolate.interp1d(npm_A_I[:,0], npm_A_I[:,1], kind='linear', fill_value=0.0, bounds_error=False)
    return interpolator


# dr_perp should be notably smaller then the grid size in the original 3D space.
def dA_I_dr_perp(r_perp, dr_perp, R, a, d, I, interpolator, debug = 0):
    # Get my 2 values of A_I
    A_Ia = interpolator(r_perp - dr_perp)
    A_Ib = interpolator(r_perp + dr_perp)

    # Numerical differentiation
    out = (A_Ib - A_Ia) / (2.0 * dr_perp)

    # Debugging
    if debug > 0:
        print('dA_I_dr_perp: ' + str(out))
    if debug > 1:
        print('  dA_I_dr_perp  A_Ia: ' + str(A_Ia))
        print('  dA_I_dr_perp  A_Ib: ' + str(A_Ib))
    return out


# The numerical derivative of A_I
def dA_I_dx(x, y, z, R, a, d, I, Dx, debug = 0):
    A_I_x_minus_1 = A_I(x - Dx, y, z, R, a, d, I, debug - 2)
    A_I_x_plus_1  = A_I(x + Dx, y, z, R, a, d, I, debug - 2)

    # Out
    out = (A_I_x_plus_1 - A_I_x_minus_1)/(2.0*Dx)

    # Debugging
    if debug > 0:
        print('dA_I_dx: ' + str(out))
    if debug > 1:
        print('  dA_I_dx  A_I_x_minus_1: ' + str(A_I_x_minus_1))
        print('  dA_I_dx  A_I_x_plus_1: ' + str(A_I_x_plus_1))
    # Output
    return out


# The resulting B_I function.
def B_I(x, y, z, R, a, d, I, Dx, A_I_r_perp_interpolator, debug = 0):
    # Values
    A_I_val = A_I(x, y, z, R, a, d, I, debug - 2)
    dA_I_dx_val = dA_I_dx(x, y, z, R, a, d, I, Dx, debug - 2)
    r_perp = r_perpendicular_scalar(y, z, d, debug - 2)
    r_perp_vec = r_perpendicular_vector(y, z, d, debug - 2)

    # To get dA_I_dr_perp we use inperpolation to get A_I(r_perp).
    dA_I_dr_perp_val = dA_I_dr_perp(r_perp, Dx * 0.2, R, a, d, I, A_I_r_perp_interpolator, debug - 2)

    # Parts
    part_1 = - dA_I_dx_val * ( r_perp_vec / r_perp )
    #print '1: ' + str(part_1)
    part_2 = np.array([dA_I_dr_perp_val + ( A_I_val / r_perp ), 0, 0])
    #print '2: ' + str(part_2) + '\n'

    # Output
    out = np.add(part_1,part_2)

    # Debugging
    if debug > 0:
        print('B_I: ' + str(out))
    if debug > 1:
        print('  B_I  part_1: ' + str(part_1))
        print('  B_I  part_2: ' + str(part_2) + '\n')
    return out


if __name__ == '__main__':
    # User-specified parameters
    tup_shape = ( 20, 20, 20 )
    x_range = ( -80.0, 80 ) * u.Mm
    y_range = ( -80.0, 80 ) * u.Mm
    z_range =  ( 0.0, 120 ) * u.Mm

    # Derived parameters (make SI where applicable)
    x_0 = x_range[0].to(u.m).value
    Dx = (( x_range[1] - x_range[0] ) / ( tup_shape[0] * 1.0 )).to(u.m).value
    x_size = Dx * tup_shape[0]
    y_0 = y_range[0].to(u.m).value
    Dy = (( y_range[1] - y_range[0] ) / ( tup_shape[1] * 1.0 )).to(u.m).value
    y_size = Dy * tup_shape[1]
    z_0 = z_range[0].to(u.m).value
    Dz = (( z_range[1] - z_range[0] ) / ( tup_shape[2] * 1.0 )).to(u.m).value
    z_size = Dy * tup_shape[2]

    # For B_I field only, to save re-creating this interpolator for every cell.
    A_I_r_perp_interpolator = interpolate_A_I_from_r_perp(flo_TD_R, flo_TD_a, flo_TD_d, flo_TD_I, (x_size**2 + y_size**2 + z_size**2)**(0.5)*1.2, 10000)

    field = np.zeros( ( tup_shape[0], tup_shape[1], tup_shape[2], 3 ) )
    for i in range(0, tup_shape[0]):
        for j in range(0, tup_shape[1]):
            for k in range(0, tup_shape[2]):
                # Position of this point in space
                x_pos = x_0 + ( i + 0.5 ) * Dx
                y_pos = y_0 + ( j + 0.5 ) * Dy
                z_pos = z_0 + ( k + 0.5 ) * Dz

                #field[i,j,k] = B_theta(x_pos, y_pos, z_pos, flo_TD_a, flo_TD_d, flo_TD_R, flo_TD_I, flo_TD_I_0)
                #field[i,j,k] = B_q(x_pos, y_pos, z_pos, flo_TD_L, flo_TD_d, flo_TD_q)
                #field[i,j,k] = B_I(x_pos, y_pos, z_pos, flo_TD_R, flo_TD_a, flo_TD_d, flo_TD_I, Dx, A_I_r_perp_interpolator)
                field[i,j,k] = B_theta(x_pos, y_pos, z_pos, flo_TD_a, flo_TD_d, flo_TD_R, flo_TD_I, flo_TD_I_0) + B_q(x_pos, y_pos, z_pos, flo_TD_L, flo_TD_d, flo_TD_q) + B_I(x_pos, y_pos, z_pos, flo_TD_R, flo_TD_a, flo_TD_d, flo_TD_I, Dx, A_I_r_perp_interpolator)




    map_field = Map3D( field, {}, xrange=x_range, yrange=y_range, zrange=z_range )
    np_boundary_data = field[:,:,0,2].T
    dummyDataToMap(np_boundary_data, x_range, y_range)

    #dic_boundary_data = { 'datavals': np_boundary_data.data.shape[0]**2, 'dsun_obs': 147065396219.34,  }
    visualise(map_field, scale=1.0*u.Mm, show_volume_axes=True, debug=True)
