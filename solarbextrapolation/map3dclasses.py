# -*- coding: utf-8 -*-

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
from copy import deepcopy
#from sunpy.sun._constants import physical_constants as con

# SunPy imports
import sunpy.map
from sunpy.sun import constants, sun
from sunpy.time import parse_time, is_time
from astropy.table import Table
import astropy.units as u

# Internal imports
from .utilities import *

__all__ = ["Map3D", "Map3DCube", "Map3DComparer"]

class Map3D(object):
    """
    A basic data structure for holding a 3D numpy array of floats or 3-float
    vectors and metadata.
    The structure can be saved/loaded (using pickle ATM).

    Parameters
    ----------

    data : `numpy.array`
        The numpy array containing the numerical data.
    meta : `dictionary`
        The container for additional information about the data in this object.
        Where:
        * x/y/zrange: the max/min spacial positions along the given axis.
        * x/yobsrange: the observational data range, often in arcsec.
        * cdelt1/2/3: the size of each pixel in each axis.
        * unit1/2/3: the spacial units in each axis.
        * naxis1/2/3: the number of pixels in each axis.

    """
    def __init__(self, data, meta, **kwargs):
        self.data = data
        self.meta = meta
        self.xrange = kwargs.get('xrange', [ 0, data.shape[1] ] * u.pixel)
        self.yrange = kwargs.get('yrange', [ 0, data.shape[0] ] * u.pixel)
        self.zrange = kwargs.get('zrange', [ 0, data.shape[2] ] * u.pixel)
        self.xobsrange = kwargs.get('xobsrange', self.xrange)
        self.yobsrange = kwargs.get('yobsrange', self.yrange)

        # Add some general properties to the metadata dictionary
        self.meta['xrange'] = self.xrange
        self.meta['yrange'] = self.yrange
        self.meta['zrange'] = self.zrange
        self.meta['cdelt1'] = ((self.xrange[1] - self.xrange[0]) / self.data.shape[1]).value
        self.meta['cdelt2'] = ((self.yrange[1] - self.yrange[0]) / self.data.shape[0]).value
        self.meta['cdelt3'] = ((self.zrange[1] - self.zrange[0]) / self.data.shape[2]).value
        # Note: should be reversed to fortran array indexing
        self.meta['cunit1'] = self.xrange.unit
        self.meta['cunit2'] = self.yrange.unit
        self.meta['cunit3'] = self.zrange.unit
        self.meta['naxis1'] = self.data.shape[1]
        self.meta['naxis2'] = self.data.shape[0]
        self.meta['naxis3'] = self.data.shape[2]
        if kwargs.get('date_obs', False):
            self.meta['date-obs'] = kwargs.get('date_obs')
        self.meta['rsun_ref'] = kwargs.get('rsun_ref', constants.radius.value)
        if kwargs.get('dsun_obs', False):
            self.meta['dsun_obs'] = kwargs.get('dsun_obs')
        if kwargs.get('bunit', False):
            self.meta['bunit'] = kwargs.get('bunit')

        # For alignment with the boundary data
        self.meta['xobsrange'] = self.xobsrange
        self.meta['yobsrange'] = self.yobsrange


    @property
    def is_scalar(self, **kwargs):
        """
        Returns true if data is a volume of scalar values (3D array) or false
        if it is a volume of vector values (4D array).
        """
        return (True if self.data.ndim is 3 else False)

    @property
    def units(self, **kwargs):
        """
        Image coordinate units along the x, y and z axes (cunit1/2/3).
        """

        # Define a triple, a named tuple object for returning values
        Triple = namedtuple('Triple', 'x y z')

        return Triple(u.Unit(self.meta.get('cunit1', 'pix')),
                      u.Unit(self.meta.get('cunit2', 'pix')),
                      u.Unit(self.meta.get('cunit3', 'pix')))

    @property
    def scale(self, **kwargs):
        """
        Image scale along the x, y and z axes in units/pixel (cdelt1/2/3)
        """
        # Define a triple, a named tuple object for returning values
        from collections import namedtuple
        Triple = namedtuple('Triple', 'x y z')

        '''
        return Triple(u.Unit(self.meta.get('cdelt1', 'arcsec')),
                      u.Unit(self.meta.get('cdelt1', 'arcsec')),
                      u.Unit(self.meta.get('cdelt2', 'arcsec')))
        '''
        return Triple(self.meta.get('cdelt1', 1.) * self.units.x / u.pixel,
                      self.meta.get('cdelt2', 1.) * self.units.y / u.pixel,
                      self.meta.get('cdelt3', 1.) * self.units.z / u.pixel)

    @property
    def rsun_meters(self, **kwargs):
        """Radius of the sun in meters"""
        return u.Quantity(self.meta.get('rsun_ref', constants.radius), 'meter')

    @property
    def date(self, **kwargs):
        """Image observation time"""
        time = parse_time(self.meta.get('date-obs', 'now'))
        if time is None:
            warnings.warn("Missing metadata for observation time. Using current time.",                               Warning)
        return parse_time(time)

    @property
    def dsun(self, **kwargs):
        """
        The observer distance from the Sun.
        """
        dsun = self.meta.get('dsun_obs', None)

        if dsun is None:
            warnings.warn("Missing metadata for Sun-spacecraft separation: assuming Sun-Earth distance",
                                   Warning)
            dsun = sun.sunearth_distance(self.date).to(u.m)

        return u.Quantity(dsun, 'm')

# #### I/O routines #### #
    @classmethod
    def load(self, filepath, **kwargs):
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

    def __getitem__(self, key, **kwargs):
        """
        Overriding indexing operation.  If the key results in a single map,
        then a map object is returned.  This allows functions like enumerate to
        work.  Otherwise, a mapcube is returned.
        """
        if isinstance(self.maps[key], Map3D):
            return self.maps[key]
        else:
            return Map3DCube(self.maps[key])

    def __len__(self, **kwargs):
        """
        Return the number of maps in a mapcube.
        """
        return len(self.maps)


    def all_maps_same_shape(self, **kwargs):
        """
        Tests if all the 3D maps have the same shape.
        """
        return np.all([m.data.shape == self.maps[0].data.shape for m in self.maps])


class Map3DComparer(object):
    """
    | Class for comparrison of vector fields.
    | There are two classification of test:
    | * **Mono**: returns a value for a given vector field. Can be normalized to the benchmark field.
    | * **Binary**: requires comparrison between two vector fields.
    | By default:
    | * Benchmark field is the first/original vector field. This is used as the baseline for comparrison. This can be changed using the ``benchmark=n`` kwarg.
    | * Normalise will be set to false.
    | Individual tests can be run and return results for imediate viewing (using astropy.table).
    | Likewise, compare_all can be used to run the whole series of tests.
    | Note: all vector fields must be of the same shape.

    """
    def __init__(self, map3D, *args, **kwargs):
        # Use all the user parameters
        self.maps_list = map3D + expand_list(args)
        self.benchmark = kwargs.get('benchmark', 0) # Defaults to the first vector field in the list
        self.normalise = kwargs.get('normalise', False)

        # The table to store the test results
        self.results = Table(names=('extrapolator routine', 'extrapolation duration', 'fig of merit 1'), meta={'name': '3D field comparison table'}, dtype=('S24', 'f8', 'f8'))
        t['time (ave)'].unit = u.s

        # An empty table for the results:
        #N = len(self.maps_list)
        #t1, t2, t3, t4, t5, t6, t7 = [None] * N, [None] * N, [None] * N, [None] * N, [None] * N, [None] * N, [None] * N
        #self.results = Table([t1, t2, t3, t4, t5, t6, t7], names=('l-infinity norm', 'test 2', 'test 3', 'test 4', 'test 5', 'test 6', 'test 7'), meta={'name': 'Results Table'})
        #self.results_normalised = Table([t1, t2, t3, t4, t5, t6, t7], names=('l-infinity norm', 'test 2', 'test 3', 'test 4', 'test 5', 'test 6', 'test 7'), meta={'name': 'Results Table'})

        # Ensure that the input maps are all the same type and shape.
        for m in self.maps_list:#self.maps:
            # Check that this is a Map3D object.
            if not isinstance(m, Map3D):
                raise ValueError(
                         'Map3DComparer expects pre-constructed map3D objects.')

            # Compare the shape of this Map3D to the first in the Map3D list.
            if not m.data.shape == self.maps_list[0]:
                raise ValueError(
                         'Map3DComparer expects map3D objects with identical dimensions.')


    def _normalise():
        """
        Return the normalised table.
        """
        # Get the benchmark extrapolation result.
        row_benchmark = self.results[self.benchmark]

        # Create a copy of the table
        tbl_output = deepcopy(self.results)

        for row in tbl_output:
            for val, val_benchmark in zip(row, row_benchmark):
                # If the value is a float then normalise.
                if type(val) == np.float64 or type(val) == np.float32 or type(val) == np.float16:
                    val = val / val_benchmark



    def L_infin_norm(map_field, benchmark, **kwargs):
        """
        l-infinity norm of the vector field.
        For vector field :math:`\bfx` this would be:

        .. math::

           \| \mathbf{x} \| \infty = \sqrt[\infty]{\Sigma_i x_i^\infty} \approx \text{max}(|x_i|)

        (the malue of the maximum component)

        From: https://rorasa.wordpress.com/2012/05/13/l0-norm-l1-norm-l2-norm-l-infinity-norm/
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

    def compare_all(self, **kwargs):
        """
        Compare all of the given vector fields and return the results as an
        astropy.table.
        """
        #num_tests = 1
        #num_maps = len(self.maps)
        #arr_data = np.zeros([num_tests, num_maps])

        # For each given 3D field, run all the tests and add a row to the table.
        for map3D in self.maps:
            # Get the data
            arr_data = map3D.data

            # Store the results from each test for this field.
            lis_results = [ map3D.meta.get('extrapolator_routine', 'Unknown Routine'),
                            map3D.meta.get( 'extrapolator_duration', 0.0 ) ]

            # Run through all the tests and append results to the list.
            lis_results.append(self.L_infin_norm(arr_data))

            # Now add the results to the table.
            self.results.add_row(lis_results)


        if self.normalise:
            self.results_normalised
        else:
            self.results
