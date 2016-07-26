# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 19:17:31 2015

@author: alex_
"""

# General Imports
import matplotlib as mpl
mpl.use('TkAgg') # Force mpl backend not to use qt. Else we have a conflict.
import numpy as np
import pickle
import time
from datetime import datetime
from collections import namedtuple
import warnings
import inspect
#from sunpy.sun._constants import physical_constants as con

# SunPy imports
import sunpy.map
from sunpy.sun import constants, sun
from sunpy.time import parse_time, is_time
from astropy.table import Table
import astropy.units as u

# Internal imports
#from solarbextrapolation.utilities import *
from solarbextrapolation.map3dclasses import Map3D

class Extrapolators(object):
    """
    Common class for all 3D vector field extrapolation routines.
    Each routine, created by building a subclass, will have wildly varying
    capabilities and input arguments so this have been left intentionally
    minimal.
    The primary method to override is extrapolation(), the primary method to
    call is extrapolate() which will both call extrapolation() and save the
    result if a filepath argument is given.

    Parameters
    ----------

    map_magnetogram : `sunpy.map.GenericMap`
        The sunpy map containing the boundary magnetogram data.

    filepath : `string`
        The optional filepath for automatic saving of extrapolation results.

    notes : `string`
        The optional notes regarding thius run of the extrapolation routine.

    extrapolator_routine : `string`
        The name for the extrapolation routine.

    zshape : `int`
        The vertical grid size.

    xrange : `astropy.unit.Quantity`, optional
        The x edge to edge coordinates. If defined will manually scale the
        boundary data.

    yrange : `astropy.units.quantity.Quantity`, optional
        The y edge to edge coordinates. If defined will manually scale the
        boundary data.

    zrange : `astropy.unit.Quantity`
        The vertical edge to edge coordinates for the vertical range.

    notes : `string`
        User specified notes that will be added to the metadata.
    """

    def __init__(self, map_magnetogram, **kwargs):
        """
        Construct an extrapolator using the given 2D map.
        """
        self.map_boundary_data = map_magnetogram
        self.meta = { 'boundary_1_meta': self.map_boundary_data.meta }
        self.meta['extrapolator_notes'] = kwargs.get('notes', '')

        # Normalise the units to SI May possible be added here

        # Crop the boundary data if required.
        self.xrange = kwargs.get('xrange', self.map_boundary_data.xrange)
        self.yrange = kwargs.get('yrange', self.map_boundary_data.yrange)
        self.map_boundary_data = self.map_boundary_data.submap(self.xrange, self.yrange)
        self.xobsrange = self.map_boundary_data.xrange
        self.yobsrange = self.map_boundary_data.yrange

        #print '\n\nHelp for u:'
        #print 'help(u): ' + str(help(u))
        #print '\n\n'
        self.zrange = kwargs.get('zrange', u.Quantity([0.0, 1.0] * u.Mm) )
        self.shape = np.asarray([self.map_boundary_data.data.shape[0],
                      self.map_boundary_data.data.shape[1],
                      long(kwargs.get('zshape', 1L))])
        self.filepath = kwargs.get('filepath', None)
        self.routine = kwargs.get('extrapolator_routine', type(self))


    def _angle_to_length(self, arc, **kwargs):
        """
        Approximate a surface length from the observed arc length.
        Uses the small angle approximation.
        """
        r = self.map_boundary_data.dsun - self.map_boundary_data.rsun_meters
        length = (r * arc.to(u.radian))
        return length.to(u.m, equivalencies=u.dimensionless_angles())

    def _to_SI(self, **kwargs):
        """

        """
        # Scale the x/y ranges
        # Setup the equivilence
        obs_distance = self.map_boundary_data.dsun - self.map_boundary_data.rsun_meters
        radian_length = [ (u.radian, u.meter, lambda x: obs_distance * x, lambda x: x / obs_distance) ]

        # Extract the maps x/yrange values and convert to length units
        #x_range = self.map_boundary_data.xrange
        #x_range = ( decompose_ang_len(x_range[0], equivalencies=radian_length),
        #            decompose_ang_len(x_range[1], equivalencies=radian_length) )
        #x_range =
        #y_range = self.map_boundary_data.yrange
        """
        x_range = self.map_boundary_data.xrange.to(u.meter, equivalencies=radian_length)
        y_range = self.map_boundary_data.yrange.to(u.meter, equivalencies=radian_length)
        # Normalise to start at 0.0
        x_range = [self.map_boundary_data.xrange[0] - self.map_boundary_data.xrange[0],
                   self.map_boundary_data.xrange[1] - self.map_boundary_data.xrange[0]]
        y_range = [self.map_boundary_data.yrange[0] - self.map_boundary_data.yrange[0],
                   self.map_boundary_data.yrange[1] - self.map_boundary_data.yrange[0]]
        """
        # Scale the magnetic field units
        ori_bunit = u.Unit(self.map_boundary_data.meta.get('bunit', 'Tesla'))
        scale_factor = ori_bunit.to(u.T)
        self.map_boundary_data = self.map_boundary_data * scale_factor
        self.map_boundary_data.meta['bunit'] = 'Tesla'
        self.meta['boundary_1_meta']['bunit'] = 'Tesla'

    def _extrapolation(self, **kwargs):
        """
        The method for running an extrapolation routine.
        This is the primary method to be edited in subclasses for specific
        extrapolation routine implementations.
        """
        # Add some type checking, we want a map object, check for .unit attribute.
        # Extrapolation code goes here.
        arr_4d = np.zeros([self.map_boundary_data.data.shape[0], self.map_boundary_data.data.shape[1], 1, 3])

        # Calculate the ranges in each dimension in length units (meters)
        x_range = self._angle_to_length(self.xrange)
        y_range = self._angle_to_length(self.yrange)
        z_range = self.zrange

        # Turn the 4D array into a Map3D object.
        map_output = Map3D( arr_4d, self.meta, xrange=x_range, yrange=y_range, zrange=z_range, xobsrange=self.xobsrange, yobsrange=self.yobsrange )

        return map_output

    def extrapolate(self, **kwargs):
        """
        Method to be called to run the extrapolation.
        Times and saves the extrapolation where applicable.
        """
        # Record the time and duration of the extrapolation.
        dt_start = datetime.now()
        tim_start = time.time()
        arr_output = self._extrapolation(**kwargs)
        tim_duration = time.time() - tim_start

        # Add the duration and time to the meta/header data.
        arr_output.meta['extrapolator_start_time'] = dt_start.isoformat()
        arr_output.meta['extrapolator_duration'] = tim_duration
        arr_output.meta['extrapolator_duration_unit'] = u.s

        # Save the Map3D if a filepath has been set. (to avoid loosing work)
        if self.filepath:
            arr_output.save(self.filepath)
        return arr_output
