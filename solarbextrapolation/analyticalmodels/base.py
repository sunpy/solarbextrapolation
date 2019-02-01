# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 19:30:22 2015

@author: alex_
"""

# General Imports
import matplotlib as mpl
mpl.use('TkAgg') # Force mpl backend not to use qt. Else we have a conflict.
import numpy as np
#import pickle
import time
from datetime import datetime
#from collections import namedtuple
import warnings
import inspect
#from sunpy.sun._constants import physical_constants as con

# SunPy imports
import sunpy.map
from sunpy.sun import constants, sun
from sunpy.time import parse_time, is_time
from astropy.table import Table
import astropy.units as u
from mayavi import mlab

# Internal imports
#from solarbextrapolation.utilities import si_this_map
from solarbextrapolation.map3dclasses import Map3D

class AnalyticalModel(object):
    """
    Common class for the development of anylitical models of magnetic fields.
    Use the models to evaluate the accuracy of an extrapolation routine with
    the figures of merit.
    """
    def __init__(self, **kwargs):
        # Default grid shape and physical ranges for the volume the model covers.
        self.shape  = kwargs.get('shape', u.Quantity([5, 5, 5] * u.pixel)) # (x,y,z)
        self.xrange = kwargs.get('xrange', u.Quantity([-10, 10] * u.Mm))
        self.yrange = kwargs.get('yrange', u.Quantity([-10, 10] * u.Mm))
        self.yrange = kwargs.get('zrange', u.Quantity([0, 20] * u.Mm))

        # Metadata
        self.meta = {'ZNAXIS': 3, 'ZNAXIS1': self.shape[0].value, 'ZNAxXIS2': self.shape[0].value, 'ZNAXIS3': self.shape[0].value}
        self.meta['analytical_model_notes'] = kwargs.get('notes', '')
        self.meta['BUNIT'] = kwargs.get('bunit', u.T)
        # CRVALn, CDELTn and NAXIS (alreadu in meta) used for storing range in 2D fits files.
        self.filepath = kwargs.get('filepath', None)
        self.routine = kwargs.get('analytical_model_routine', type(self))

        # Default 3D magnetic field
        #X,Y,Z = np.zeros(self.shape.value), np.zeros(self.shape.value), np.zeros(self.shape.value)
        npField = np.zeros([3]+self.shape.value)
        self.field = Map3D(npField, self.meta)

        # Default magnetic field on boundary
        magnetogram = np.zeros(self.shape[0:2].value)
        magnetogram_header  = {'ZNAXIS': 2, 'ZNAXIS1': self.shape[0].value, 'ZNAXIS2': self.shape[1].value}
        self.magnetogram = sunpy.map.Map((magnetogram, magnetogram_header))

    def _generate_field(self, **kwargs):
        """
        The method for running a model to generate the field.
        This is the primary method to be edited in subclasses for specific
        model implementations.
        """
        # Model code goes here.
        arr_4d = np.zeros([self.map_boundary_data.data.shape[0], self.map_boundary_data.data.shape[1], 1, 3])

        # Turn the 4D array into a Map3D object.
        map_output = Map3D( arr_4d, self.meta, xrange=self.xrange, yrange=self.yrange, zrange=self.zrange, xobsrange=self.xrange, yobsrange=self.yrange )

        return map_output

    def generate(self, **kwargs):
        """
        Method to be called to calculate the vector field and return as a Map3D object.
        Times and saves the extrapolation where applicable.
        """
        # Record the time and duration of the extrapolation.
        dt_start = datetime.now()
        tim_start = time.time()
        arr_output = self._generate_field(**kwargs)
        tim_duration = time.time() - tim_start

        # Add the duration and time to the meta/header data.
        arr_output.meta['extrapolator_start_time']    = dt_start.isoformat()
        arr_output.meta['extrapolator_duration']      = tim_duration
        arr_output.meta['extrapolator_duration_unit'] = u.s

        # Save the Map3D if a filepath has been set. (to avoid loosing work)
        if self.filepath:
            arr_output.save(self.filepath)

        # Add the output map to the object and return.
        self.map = arr_output
        return arr_output

    def to_los_magnetogram(self, **kwargs):
        """
        Calculate the LoS vector field as a SunPy map and return.

        Generally this will require that you have run generate(self, ``**kwargs``)
        first, so in the base class this is checked, but it is not always the
        case as some models may allow this to be determined without calculating
        the full field.

        .. I'm not sure if this is a good default.
        """
        return self.magnetogram

    def to_vec_magnetogram(self, **kwargs):
        """
        Calculate the vector field as a SunPy map and return.

        Generally this will require that you have run ``generate(self, **kwargs)``
        first, so in the base class this is checked, but it is not always the
        case as some models may allow this to be determined without calculating
        the full field. ######### I'm not sure if this is a good default.
        """
        return self.magnetogram

if __name__ == '__main__':
    # User-specified parameters
    tup_shape = ( 20, 20, 20 )
    x_range   = ( -80.0,  80 ) * u.Mm
    y_range   = ( -80.0,  80 ) * u.Mm
    z_range   = (   0.0, 120 ) * u.Mm

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




    # Define the extrapolator as a child of the Extrapolators class
    class AnaOnes(AnalyticalModel):
        def __init__(self, **kwargs):
            super(AnaOnes, self).__init__(**kwargs)

        def _generate_field(self, **kwargs):
            # Adding in custom parameters to the metadata
            self.meta['analytical_model_routine'] = 'Ones Model'

            # Generate a trivial field and return (X,Y,Z,Vec)
            arr_4d = np.ones(self.shape.value.tolist() + [3])
            return Map3D( arr_4d, self.meta )


    # Setup an anylitical model
    xrange = u.Quantity([  50,  300] * u.arcsec)
    yrange = u.Quantity([-350, -100] * u.arcsec)
    zrange = u.Quantity([   0,  250] * u.arcsec)

    aAnaMod = AnaOnes()
    aMap3D = aAnaMod.generate()


    # Visualise the 3D vector field
    from solarbextrapolation.visualisation_functions import visualise
    """
    fig = visualise(aMap3D,
                    show_boundary_axes=False,
                    boundary_units=[1.0*u.arcsec, 1.0*u.arcsec],
                    show_volume_axes=True,
                    debug=False)
    """
    fig = visualise(aMap3D,
                    show_boundary_axes=False,
                    show_volume_axes=False,
                    debug=False)
    mlab.show()


    """
    # For B_I field only, to save re-creating this interpolator for every cell.
    A_I_r_perp_interpolator = interpolate_A_I_from_r_perp(flo_TD_R, flo_TD_a, flo_TD_d, flo_TD_I, (x_size**2 + y_size**2 + z_size**2)**(0.5)*1.2, 1000`0)

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
    """
