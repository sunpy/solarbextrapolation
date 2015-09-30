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
from utilities import *

class AnalyticalModel(object):
    """
    Common class for the development of anylitical models of magnetic fields.
    Use the models to evaluate the accuracy of an extrapolation routine with
    the figures of merit.
    """
    def __init__(self, **kwargs):
        dim = 16
        # 
        header = {'ZNAXIS': 3, 'ZNAXIS1': dim, 'ZNAXIS2': dim, 'ZNAXIS3': dim}
        X, Y, Z = np.zeros([dim, dim, dim]), np.zeros([dim, dim, dim]), np.zeros([dim, dim, dim])
        self.field = Map3D(X, Y, Z, header)
        
        magnetogram = np.zeros([dim, dim])
        magnetogram_header  = {'ZNAXIS': 2, 'ZNAXIS1': dim, 'ZNAXIS2': dim}
        self.magnetogram = sunpy.map.Map((magnetogram, magnetogram_header))

    def generate(self, **kwargs):
        """
        Calculate the vector field and return as a Map3D object.
        """
        return self.map

    def to_magnetogram(self, **kwargs):
        """
        Calculate the LoS vector field as a SunPy map and return.
        """
        return self.magnetogram