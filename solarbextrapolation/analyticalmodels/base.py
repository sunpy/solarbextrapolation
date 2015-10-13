# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 19:30:22 2015

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
        self.shape = kwargs.get('shape', u.Quantity([5, 5, 5] * u.pixel)) # (x,y,z)
        self.xrange = kwargs.get('xrange', u.Quantity([-10, 10] * u.Mm))
        self.yrange = kwargs.get('yrange', u.Quantity([-10, 10] * u.Mm))
        self.yrange = kwargs.get('zrange', u.Quantity([0, 20] * u.Mm))

        # Metadata 
        self.meta = {'ZNAXIS': 3, 'ZNAXIS1': self.shape[0].value, 'ZNAXIS2': self.shape[0].value, 'ZNAXIS3': self.shape[0].value}
        self.meta['analytical_model_notes'] = kwargs.get('notes', '')
        self.meta['BUNIT'] = kwargs.get('bunit', u.T)
        # CRVALn, CDELTn and NAXIS (alreadu in meta) used for storing range in 2D fits files.
        
        # Default 3D magnetic field
        X, Y, Z = np.zeros(self.shape.value), np.zeros(self.shape.value), np.zeros(self.shape.value)
        self.field = Map3D(X, Y, Z, self.meta)
        
        # Default magnetic field on boundary
        magnetogram = np.zeros(self.shape[0:2].value)
        magnetogram_header  = {'ZNAXIS': 2, 'ZNAXIS1': self.shape[0].value, 'ZNAXIS2': self.shape[1].value}
        self.magnetogram = sunpy.map.Map((magnetogram, magnetogram_header))

    def generate(self, **kwargs):
        """
        Calculate the vector field and return as a Map3D object.
        """
        return self.map

    def to_los_magnetogram(self, **kwargs):
        """
        Calculate the LoS vector field as a SunPy map and return.
        
        Generally this will require that you have run generate(self, **kwargs)
        first, so in the base class this is checked, but it is not always the
        case as some models may allow this to be determined without calculating
        the full field. ######### I'm not sure if this is a good default.
        """
        return self.magnetogram

    def to_vec_magnetogram(self, **kwargs):
        """
        Calculate the vector field as a SunPy map and return.
        
        Generally this will require that you have run generate(self, **kwargs)
        first, so in the base class this is checked, but it is not always the
        case as some models may allow this to be determined without calculating
        the full field. ######### I'm not sure if this is a good default.
        """
        return self.magnetogram