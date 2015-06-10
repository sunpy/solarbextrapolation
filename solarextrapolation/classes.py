# -*- coding: utf-8 -*-
"""
Created on Fri Jun 05 12:49:03 2015

@author: Alex
"""

import sunpy.map
import numpy as np
import pickle
import time
#import astropy.units as u

class Preprocessors(object):
    """
    A common class for all 2D pre-processing routines, tools used to pre-process
    the 2D sunpy map data for use in extrapolations.
    Usage can include basic filters for noise/contrast or algorythms to
    compensate for extrapolator assumptions, such as the force-free assumption
    that is assumed in many extrapolations, but isn't true in the photosphere
    where magnetiogram observations are generally taken.
    """
    def __init__(self, map_magnetogram):
        """Method for creating a preprocessor object, using a sunpy map.
        """
        # Add some type checking, we want a map object, check for .unit attribute.
        self.data_map = map_magnetogram

    def preprocess(self):
        """Method running the preprocessor and returning a sunpy map.
        """
        return self.data_map


class Extrapolators(object):
    """
    Common class for all 3D vector field extrapolation routines.
    Each routine, created by building a subclass, will have wildly variging
    capabilities and input arguments so this have been left intentionally
    minimal.
    The primary method to override is extrapolation(), the primary method to
    call is extrapolate() which will both call extrapolation() and save the
    result if a filepath argument is given.
    """

    def __init__(self, map_magnetogram, **kwargs):
        """Construct an extrapolator using the given 2D map.
        """
        self.data_boundary_map = map_magnetogram
        self.meta = { 'boundary_1_meta': map_magnetogram }
        self.z = 10
        if kwargs.get('z'):
            self.filepath = kwargs['z']
        self.filepath = ''
        if kwargs.get('filepath'):
            self.filepath = kwargs['filepath']

    def extrapolation(self):
        """
        The routine for running an extrapolation routine.
        This is the primary method to be edited in subclasses for specific
        extrapolation routine implementations.
        """
        # Add some type checking, we want a map object, check for .unit attribute.
        # Extrapolation code goes here.
        arr_4d = np.zeros([0, 0, 0, 0])
        self.meta['extrapolation_routine'] = 'default extrapolator class'
        return Map3D( arr_4d, self.meta, True )

    def extrapolate(self):
        """
        Method to be called to run the extrapolation.
        Times and saves the extrapolation where applicable.
        """
        tim_start = time.time()
        arr_output = self.extrapolation()
        tim_duration = time.time() - tim_start

        arr_output.meta['extrapolation_start_time'] = tim_start
        arr_output.meta['extrapolation_duration'] = tim_duration

        if self.filepath:
            arr_output.save(self.filepath)
        return arr_output



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
    """
    A basic data structure for holding a 3D numpy array of floats or 3-float
    vectors and metadata.
    The structure can be saved/loaded (using pickle ATM).
    """
    def __init__(self, arr_data = np.zeros([0,0,0]), dic_meta = {}, vectors = False):
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

# #### I/O routines #### #
    def save(self, filepath, filetype='auto', **kwargs):
        """Saves the Map3D object to a file.

        Currently uses Python pickle.
        https://docs.python.org/2/library/pickle.html
        In the future support will be added for saving to other formats.

        Parameters
        ----------
        filepath : string
            Location to save file to.

        filetype : string
            'auto' or any supported file extension
        """
        #io.write_file(filepath, self.data, self.meta, filetype=filetype,
        #              **kwargs)
        pickle.dump( self, open( filepath, "wb" ), **kwargs)

    def _read_file(self, fname, **kwargs):
        loaded = pickle.load( open( fname, "rb" ) )
        self.meta = loaded.meta
        self.data = loaded.data
        self.is_vector = loaded.is_vector
        self.is_scalar = loaded.is_scalar







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

class Map3DComparer(object):
    def __init__(self, map3D, *args, **kwargs):
        self.maps = expand_list(args)

        for m in self.maps:
            if not isinstance(m, Map3D):
                raise ValueError(
                         'CompositeMap expects pre-constructed map3D objects.')

        self.normalize_all = True
        if kwargs.get('normalize'):
            self.normalize = kwargs['normalize']
        self.normalize_to = 0
        if kwargs.get('normalize_to'):
            self.normalize = kwargs['normalize_to']

    def L_infin_norm(map_field, benchmark):
        # Placeholder for the maximum value.
        output = - 10.0**15

        # Iterate through the volume
        ni, nj, nk, D = map_field.shape
        for i in range(0, ni):
            for j in range(0, nj):
                for k in range(0, nk):
                    # Get the sum of the components
                    component_sum = 0.0
                    for component in map_field[i][j][k]:
                        component_sum += component

                    # If this is bigger then the current max value.
                    if output < component_sum:
                        output = component_sum



        # Output
        return output

    def compare_all(self):
        num_tests = 1
        num_maps = len(self.maps)
        arr_data = np.zeros([num_tests, num_maps])

        for map3D in self.maps:
            arr_data = map3D.data
            result = self.L_infin_norm(arr_data)



############
#
#  Examples
#
#######
if __name__ == '__main__':
    class PreZeros(Preprocessors):
        def __init__(self, map_magnetogram):
            super(PreZeros, self).__init__(map_magnetogram)

        def preprocess(self):
            #return self.data_map
            return sunpy.map.Map((np.zeros(self.data_map.data.shape),self.data_map.meta))

    aMap2D = sunpy.map.Map('C://Users//Alex//Dropbox//Study//2014-2015//SoCiS//Coding//Python//Random//data//hmi.m_720s.2015.05.01_00-12-00_TAI.magnetogram.fits')
    aPrePro = PreZeros(aMap2D.submap([0,10], [0,10]))
    aPreProData = aPrePro.preprocess()
    # aPreProData = aMap2D.submap([0,10], [0,10])


    class ExtZeros(Extrapolators):
        def __init__(self, map_magnetogram, **kwargs):
            super(ExtZeros, self).__init__(map_magnetogram, **kwargs)

        def extrapolation(self):
            arr_4d = np.zeros([self.data_boundary_map.data.shape[0], self.data_boundary_map.data.shape[0], self.z, 3])
            return Map3D(arr_4d, self.meta, True)

    aMap2D = sunpy.map.Map('C://Users//Alex//Dropbox//Study//2014-2015//SoCiS//Coding//Python//Random//data//hmi.m_720s.2015.05.01_00-12-00_TAI.magnetogram.fits')
    aPreProData = aMap2D.submap([0,10], [0,10])
    aExt = ExtZeros(aPreProData, filepath='C://Users/Alex/solarextrapolation/solarextrapolation/3Dmap.m3d')
    aMap3D = aExt.extrapolate()

    aExt = ExtZeros(aPreProData)




    #class Greens_Potential(Extrapolators):
    #    super(Extrapolators, self).__init__(map_magnetogram)
