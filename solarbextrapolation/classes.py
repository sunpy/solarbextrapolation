# -*- coding: utf-8 -*-
"""
Created on Fri Jun 05 12:49:03 2015

@author: Alex

This module was developed with funding provided by ESA Summer of Code in Space
(2015).
"""

import sunpy.map
import numpy as np
import pickle
import time
from datetime import datetime

from astropy.table import Table
#import astropy.units as u

__all__ = ['Preprocessors', 'Extrapolators', 'Map3D', 'Map3DCube']

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

    def _preprocessor(self):
        """
        Method running the and returning a sunpy map.
        For tracability this should add entries into the metadata that
        include any parameters used for the given run.
        """
        map_output = sunpy.map.Map(self.map_input.data, self.meta)
        return map_output

    def preprocess(self):
        
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
    z : `int`
        The vertical grid size.
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
        self.z = kwargs.get('z', 10)
        self.filepath = kwargs.get('filepath', None)
        self.routine = kwargs.get('extrapolator_routine', type(self))
        self.meta = { 'original_map_header': self.map_boundary_data.meta, 'extrapolator_routine': self.routine }

    def _extrapolation(self):
        """
        The routine for running an extrapolation routine.
        This is the primary method to be edited in subclasses for specific
        extrapolation routine implementations.
        """
        # Add some type checking, we want a map object, check for .unit attribute.
        # Extrapolation code goes here.
        arr_4d = np.zeros([0, 0, 0, 0])
        map_output = Map3D( arr_4d, self.meta )
        
        return map_output

    def extrapolate(self):
        """
        Method to be called to run the extrapolation.
        Times and saves the extrapolation where applicable.
        """
        dt_start = datetime.now()
        tim_start = time.time()
        arr_output = self._extrapolation()
        tim_duration = time.time() - tim_start

        arr_output.meta['extrapolator_start_time'] = dt_start.isoformat()
        arr_output.meta['extrapolator_duration'] = tim_duration

        if self.filepath:
            arr_output.save(self.filepath)
        return arr_output



class AnalyticalModel(object):
    """
    Common class for the development of anylitical models of magnetic fields.
    Use the models to evaluate the accuracy of an extrapolation routine with
    the figures of merit.
    """
    def __init__(self):
        dim = 16
        # 
        header = {'ZNAXIS': 3, 'ZNAXIS1': dim, 'ZNAXIS2': dim, 'ZNAXIS3': dim}
        X, Y, Z = np.zeros([dim, dim, dim]), np.zeros([dim, dim, dim]), np.zeros([dim, dim, dim])
        self.field = Map3D(X, Y, Z, header)
        
        magnetogram = np.zeros([dim, dim])
        magnetogram_header  = {'ZNAXIS': 2, 'ZNAXIS1': dim, 'ZNAXIS2': dim}
        self.magnetogram = sunpy.map.Map((magnetogram, magnetogram_header))

    def generate(self):
        """
        Calculate the vector field and return.
        """
        return self.map

    def to_magnetogram(self):
        """
        Calculate the LoS vector field as a sunpy map and return.
        """
        return self.magnetogram


class Map3D(object):
    """
    A basic data structure for holding a 3D numpy array of floats or 3-float
    vectors and metadata.
    The structure can be saved/loaded (using pickle ATM).
    """
    def __init__(self, data, meta, **kwargs):
        self.data = data
        self.meta = meta
    
    @property
    def is_scalar(self):
        """
        Returns true if data is a volume of scalar values (3D array) or false
        if it is a volume of vector values (4D array).
        """
        return (True if self.data.ndim is 3 else False)
    
    @property
    def xrange(self):
        """Return the X range of the image in arcsec from edge to edge."""
        # Original code from sunpy.map
        #xmin = self.center['x'] - self.shape[1] / 2. * self.scale['x']
        #xmax = self.center['x'] + self.shape[1] / 2. * self.scale['x']
        
        xmin, xmax = -20.0, +20.0
        return [xmin, xmax]

    @property
    def yrange(self):
        """Return the Y range of the image in arcsec from edge to edge."""
        # Original code from sunpy.map
        #ymin = self.center['y'] - self.shape[0] / 2. * self.scale['y']
        #ymax = self.center['y'] + self.shape[0] / 2. * self.scale['y']

        ymin, ymax = -20.0, +20.0
        return [ymin, ymax]

    @property
    def yrange(self):
        """Return the Z range of the image in arcsec from edge to edge."""
        zmin, zmax = 0.0, +40.0
        return [zmin, zmax]
    
# #### I/O routines #### #
    @classmethod
    def load(self, filepath):
        """
        Load a Map3D instance using pickle.
        """
        loaded = pickle.load( open( filepath, "rb" ) )
        return loaded

    def save(self, filepath, filetype='auto', **kwargs):
        """
        Saves the Map3D object to a file.

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
        pickle.dump(self, open( filepath, "wb" ), **kwargs)




from sunpy.util import expand_list


class Map3DCube:
    """
    A basic data structure for holding a list of Map3D objects.
    """
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
        """
        Overriding indexing operation.  If the key results in a single map,
        then a map object is returned.  This allows functions like enumerate to
        work.  Otherwise, a mapcube is returned.
        """
        if isinstance(self.maps[key], Map3D):
            return self.maps[key]
        else:
            return Map3DCube(self.maps[key])

    def __len__(self):
        """
        Return the number of maps in a mapcube.
        """
        return len(self.maps)
    
    
    def all_maps_same_shape(self):
        """
        Tests if all the 3D maps have the same shape.
        """
        return np.all([m.data.shape == self.maps[0].data.shape for m in self.maps])


class Map3DComparer(object):
    """
    | Class for comparrison of vector fields.
    | There are two classification of test:
    | * **Mono**: returns a value for a given vector field. Can be normalized to
    the benchmark field.
    | * **Binary**: requires comparrison between two vector fields.
    | By default:
    | * Benchmark field is the first/original vector field. This is
    used as the baseline for comparrison. This can be changed using the
    benchmark=n kwarg.
    | * Normalise will be set to false.
    | Individual tests can be run and return results for imediate viewing
    (using astropy.table).
    | Likewise, compare_all can be used to run the whole series of tests.
    | Note: all vector fields must be of the same shape.
    """
    def __init__(self, map3D, *args, **kwargs):
        self.maps_list = map3D + expand_list(args)
        self.benchmark = kwargs.get('benchmark', 0) # Defaults to the first vector field in the list
        self.normalise = kwargs.get('normalise', False)
        
        # An empty table for the results:
        N = len(self.maps_list)
        t1, t2, t3, t4, t5, t6, t7 = [None] * N, [None] * N, [None] * N, [None] * N, [None] * N, [None] * N, [None] * N
        self.results = Table([t1, t2, t3, t4, t5, t6, t7], names=('l-infinity norm', 'test 2', 'test 3', 'test 4', 'test 5', 'test 6', 'test 7'), meta={'name': 'Results Table'})
        self.results_normalised = Table([t1, t2, t3, t4, t5, t6, t7], names=('l-infinity norm', 'test 2', 'test 3', 'test 4', 'test 5', 'test 6', 'test 7'), meta={'name': 'Results Table'})

        for m in self.maps:
            if not isinstance(m, Map3D):
                raise ValueError(
                         'CompositeMap expects pre-constructed map3D objects.')


    def L_infin_norm(map_field, benchmark):
        """
        l-infinity norm of the vector field.
        For vector field :math:`\bfx` this would be:
        \left \| \mathbf{x} \right \|_\infty = \sqrt[\infty]{\Sigma_i x_i^\infty} \approx \textup{max}(|x_i|)
        (the malue of the maximum component)
        
        # From: https://rorasa.wordpress.com/2012/05/13/l0-norm-l1-norm-l2-norm-l-infinity-norm/
        """
        
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
            
        if self.normalise:
            self.results_normalised
        else:
            self.results



############
#
#  Examples
#
#######
if __name__ == '__main__':
    aNewMap3D = Map3D()
    
    # Testing a preprocessor
    
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

    aMap2D = sunpy.map.Map('C://git//solarextrapolation//solarextrapolation//data//example_data_(100x100)__01_hmi.fits')
    aPrePro = PreZeros(aMap2D.submap([0,10], [0,10]))
    aPreProData = aPrePro.preprocess()
    # aPreProData = aMap2D.submap([0,10], [0,10])
    
    # Some checks:
    #aPreProData.data # Should be a 2D zeros array.
    #aPreProData.meta
    #aPreProData.meta['preprocessor_routine']
    #aPreProData.meta['preprocessor_start_time']

    
    # Testing an extrapolator
    
    class ExtZeros(Extrapolators):
        def __init__(self, map_magnetogram, **kwargs):
            super(ExtZeros, self).__init__(map_magnetogram, **kwargs)

        def _extrapolation(self):
            # Adding in custom parameters to the meta
            self.meta['extrapolator_routine'] = 'Zeros Extrapolator'

            arr_4d = np.zeros([self.map_boundary_data.data.shape[0], self.map_boundary_data.data.shape[0], self.z, 3])
            return Map3D(( arr_4d, self.meta ))

    aExt = ExtZeros(aPreProData, filepath='C://Users/Alex/solarextrapolation/solarextrapolation/3Dmap.m3d')
    aMap3D = aExt.extrapolate()
    
    # Some checks:
    #aMap3D.data # Should be a 4D zeros array.
    #aMap3D.meta
    #aMap3D.meta['extrapolator_routine']
    #aMap3D.meta['extrapolator_start_time']
    
    # Testing a Map3DCube
    
    aMapCube = Map3DCube(aMap3D, aMap3D)
    aMapCube[0]
    aMapCube[0].data
    aMapCube[0].meta
    aMapCube[1].data
    aMapCube[1].meta

    
    
    


    #class Greens_Potential(Extrapolators):
    #    super(Extrapolators, self).__init__(map_magnetogram)
