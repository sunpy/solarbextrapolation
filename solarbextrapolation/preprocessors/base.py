# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 19:18:58 2015

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

class Preprocessors(object):
    """
    A common class for all 2D pre-processing routines, tools used to pre-process
    the 2D sunpy map data for use in extrapolations.
    Usage can include basic filters for noise/contrast or algorythms to
    compensate for extrapolator assumptions, such as the force-free assumption
    that is assumed in many extrapolations, but isn't true in the photosphere
    where magnetogram observations are generally taken.
    
    Parameters
    ----------
    
    map_data : `sunpy.map.GenericMap`
        The sunpy map containing the data to be processed.
    filepath : `string`
        The optional filepath for automatic saving of preprocessed results.
    notes : `string`
        User specified notes that will be added to the metadata.
    """
    def __init__(self, map_data, **kwargs):
        """
        Method for creating a preprocessor object, using a sunpy map.
        """
        # Add some type checking, we want a map object, check for .unit attribute.
        self.map_input = map_data
        self.routine = kwargs.get('preprocessor_routine', type(self))
        self.meta = self.map_input.meta
        self.meta['preprocessor_notes'] = kwargs.get('notes', '')
        self.meta['preprocessor_routine'] = self.routine
        self.filepath = kwargs.get('filepath', None)

    def _preprocessor(self, **kwargs):
        """
        Method running the and returning a sunpy map.
        For tracability this should add entries into the metadata that
        include any parameters used for the given run.
        """
        map_output = sunpy.map.Map(self.map_input.data, self.meta)
        return map_output

    def preprocess(self, **kwargs):
        
        """
        Method to be called to run the preprocessor.
        Times the process and saves output where applicable.
        """
        dt_start = datetime.now()
        tim_start = time.time()
        map_output = self._preprocessor()
        tim_duration = time.time() - tim_start

        map_output.meta['preprocessor_start_time'] = dt_start.isoformat()
        map_output.meta['preprocessor_duration'] = tim_duration

        if self.filepath:
            map_output.save(self.filepath)
        return map_output