# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pytest
import tempfile
import pkg_resources

import numpy as np
import os
from astropy import units as u

import sunpy.map

#from ..classes import *
from solarbextrapolation.map3dclasses import Map3D#, Map3DCube
from solarbextrapolation.preprocessors import Preprocessors
from solarbextrapolation.extrapolators import Extrapolators, PotentialExtrapolator


# Tests for the map3D class

def test_create_scalar_Map3d():
     aNumpyArray = np.zeros((2,2,2))
     aMetaDict = { 'file': 'test Map3D object'}
     aMap3D = Map3D(aNumpyArray, aMetaDict)
     assert aMap3D.is_scalar

@pytest.fixture
def test_create_vector_Map3d():
     aNumpyArray = np.zeros((2,2,2,2))
     aMetaDict = { 'file': 'test Map3D object'}
     aMap3D = Map3D(aNumpyArray, aMetaDict)
     assert not is_scalar
     return aMap3D

@pytest.fixture
def text_save_Map3d(test_create_Map3d):
    afilename = tempfile.NamedTemporaryFile(suffix='np').name
    test_create_vector_Map3d.save(afilename)
    assert os.path.isfile(afilename)
    return afilename

def text_load_Map3d(text_save_Map3d):
    aMap3D = Map3D.load(text_save_Map3d)
    # Compare the returned data array
    assert (aMap3D.data == np.zeros((2,2,2,2))).all()


# Tests for the preprocessor class

@pytest.fixture
def test_create_preprocessor():
    aNumpyArray = np.zeros((2,2))
    aMetaDict = { 'file': 'test SunPy Map object'}
    aMap2D = sunpy.map.Map(aNumpyArray, aMetaDict)
    aPreprocessor = Preprocessors(aMap2D)
    return aPreprocessor

def test_preprocessor_preprocess_method(test_create_preprocessor):
    test_create_preprocessor.preprocess()

@pytest.fixture
def test_create_extrapolator():
    aNumpyArray = np.zeros((2,2))
    aMetaDict = { 'file': 'test SunPy Map object'}
    aMap2D = sunpy.map.Map(aNumpyArray, aMetaDict)
    aExtrapolator = Extrapolators(aMap2D)
    return aExtrapolator

def test_extrapolator_extrapolate_method(test_create_extrapolator):
    test_create_extrapolator.extrapolate()

def test_create_and_run_preprocessor_subclass():
    # Define the preprocessor as a child of the Preprocessors class
    class PreZeros(Preprocessors):
        def __init__(self, map_magnetogram):
            super(PreZeros, self).__init__(map_magnetogram)

        def _preprocessor(self):
            # Adding in custom parameters to the meta
            self.meta['preprocessor_routine'] = 'Zeros Preprocessor'

            # Creating the trivial zeros map of teh shape of the input map
            map_output = sunpy.map.Map((np.zeros(self.map_input.data.shape),
                                        self.meta))

            # Outputting the map.
            return map_output

    # Instansiate the new child class
    aNumpyArray = np.zeros((2,2))
    aMetaDict = { 'file': 'test Map object'}
    aMap2D = sunpy.map.Map(aNumpyArray, aMetaDict)
    aPrePro = PreZeros(aMap2D.submap([0,10]*u.arcsec, [0,10]*u.arcsec))
    aPreProData = aPrePro.preprocess()

	
# Tests for the extrapolator class

def test_create_and_run_extrapolator_subclass():
    # Define the extrapolator as a child of the Extrapolators class
    class ExtZeros(Extrapolators):
        def __init__(self, map_magnetogram, **kwargs):
            super(ExtZeros, self).__init__(map_magnetogram, **kwargs)
            # This is a test fixing hack: TODO:  FIXME
            self.z = 1

        def _extrapolation(self):
            # Adding in custom parameters to the meta
            self.meta['extrapolator_routine'] = 'Zeros Extrapolator'

            arr_4d = np.zeros([self.map_boundary_data.data.shape[0], self.map_boundary_data.data.shape[0], self.z, 3])
            return Map3D( arr_4d, self.meta )

    # Instansiate the new child class
    afilename = tempfile.NamedTemporaryFile(suffix='np').name
    aNumpyArray = np.zeros((2,2))
    aMetaDict = { 'file': 'test SunPy Map object'}
    aMap2D = sunpy.map.Map(aNumpyArray, aMetaDict)
    aExt = ExtZeros(aMap2D, filepath=afilename)
    aMap3D = aExt.extrapolate()
    assert os.path.isfile(afilename)

# Tests for the extrapolator.PotentialExtrapolator class

def test_potential_extrapolator_subclass():
    # Parameters for the extrapolator
    xrange = u.Quantity([50,    300] * u.arcsec)
    yrange = u.Quantity([-350, -100] * u.arcsec)
    zrange = u.Quantity([0,     250] * u.arcsec)
    shape = u.Quantity([5, 5, 5] * u.pixel)

    # Load HMI map (from fits file) then submap and resample.
    hmi_filename = pkg_resources.resource_filename('solarbextrapolation',
                                                   'data/sdo-hmi_2011-02-14_20-34-12.fits')
    map_boundary = sunpy.map.Map(hmi_filename)
    map_boundary = map_boundary.submap(xrange, yrange).resample(shape[0:2], method='linear')
    
    # Extrapolate using python native code
    aPotExt = PotentialExtrapolator(map_boundary, zshape=shape[2].value, zrange=zrange)
    aMap3D = aPotExt.extrapolate(enable_numba=False)
    
    # Extrapolate using numba
    aMap3D = aPotExt.extrapolate()


