# -*- coding: utf-8 -*-
"""
Created on Fri Jun 05 12:49:03 2015

@author: Alex
"""

import sunpy.map
import numpy as np
#import astropy.units as u

class Preprocessors(object):
    'Common class for all 2D pre-processing routines.'
    def __init__(self, map_magnetogram):
        # Add some type checking, we want a map object, check for .unit attribute.
        self.data_map = map_magnetogram

    def preprocess(self):
        return self.data_map


class Extrapolators(object):
    'Common class for all 3D vector field extrapolation routines.'
    
    def __init__(self, map_magnetogram, z=10):
        self.data_boundary_map = map_magnetogram
        self.z = z
    
    def extrapolate(self):
        # Add some type checking, we want a map object, check for .unit attribute.        
        # Extrapolation code goes here.
        arr_4d = np.zeros([0, 0, 0, 0])
        return Map3D(arr_4d, self.data_boundary_map.meta)

        
class AnalyticalModel(object):
    def __init__(self):
        dim = 16
        header = {'ZNAXIS': 3, 'ZNAXIS1': dim, 'ZNAXIS2': dim, 'ZNAXIS3': dim}
        X, Y, Z = np.zeros([dim, dim, dim]), np.zeros([dim, dim, dim]), np.zeros([dim, dim, dim])
        self.field = Map3D(X, Y, Z, header)
        magnetogram = np.zeros([dim, dim])
        magnetogram_header  = {'ZNAXIS': 2, 'ZNAXIS1': dim, 'ZNAXIS2': dim}
        self.magnetogram = sunpy.map.Map((magnetogram, magnetogram_header))
        
        
    def generate(self):
        # Extrapolate the vector field and return.
        return self.map
        
    def to_magnetogram(self):
        # Extrapolate the vector field and return.
        return self.magnetogram


class Map3D(object):
    def __init__(self, arr_data, dic_meta, vectors = False):
        self.meta = dic_meta
        self.data = arr_data
        self.is_vector = vectors
        self.is_scalar = not(vectors) # Redundant
        
    def is_vector(self):
        return self.is_vector

    def is_scalar(self):
        return self.is_scalar

    def meta(self):
        return self.meta
    
    def data(self):
        return self.data



from sunpy.util import expand_list


class Map3DCube:
    def __init__(self, *args, **kwargs):

        # Hack to get around Python 2.x not backporting PEP 3102.
        #sortby = kwargs.pop('sortby', 'date')
        #derotate = kwargs.pop('derotate', False)
        
        self.maps = expand_list(args)
        
        for m in self.maps:
            if not isinstance(m, Map3D):
                raise ValueError(
                           'CompositeMap expects pre-constructed map objects.')
    
    def __getitem__(self, key):
        """Overriding indexing operation.  If the key results in a single map,
        then a map object is returned.  This allows functions like enumerate to
        work.  Otherwise, a mapcube is returned."""

        if isinstance(self.maps[key],Map3D):
            return self.maps[key]
        else:
            return Map3DCube(self.maps[key])
            
    def __len__(self):
        """Return the number of maps in a mapcube."""
        return len(self.maps)
        
        
        
        
        
############
#
#  Examples
#
#######

class PreZeros(Preprocessors):

    def __init__(self, map_magnetogram):
        super(PreZeros, self).__init__(map_magnetogram)

    def preprocess(self):
        #return self.data_map        
        return sunpy.map.Map((np.zeros(self.data_map.data.shape),self.data_map.meta))

aMap = sunpy.map.Map('C://Users//Alex//Dropbox//Study//2014-2015//SoCiS//Coding//Python//Random//data//hmi.m_720s.2015.05.01_00-12-00_TAI.magnetogram.fits')
aPrePro = PreZeros(aMap.submap([0,10], [0,10]))
aPreProData = aPrePro.preprocess()


class ExtZeros(Extrapolators):
    def __init__(self, map_magnetogram, z=10):
        super(ExtZeros, self).__init__(map_magnetogram, z)

    def extrapolate(self):
        arr_4d = np.zeros([self.data_boundary_map.data.shape[0], self.data_boundary_map.data.shape[0], self.z, 3])
        return Map3D(arr_4d, self.data_boundary_map.meta, True)



aExt = ExtZeros(aPreProData, 15)
a3DMap = aExt.extrapolate()



#class Greens_Potential(Extrapolators):
#    super(Extrapolators, self).__init__(map_magnetogram)
