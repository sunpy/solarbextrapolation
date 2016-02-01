# -*- coding: utf-8 -*-
"""
===============================================
Defining a Custom Preprocessor and Extrapolator
===============================================

Here you will be creating trivial preprocessor and and exztrqapolatoirs
following the API.

"""

################################################################################
# You start by importing the necessary modules.

# General imports
import sunpy.map as mp
import numpy as np
from mayavi import mlab # Necessary for visulisation

# Module imports
from solarbextrapolation.preprocessors import Preprocessors
from solarbextrapolation.extrapolators import Extrapolators
from solarbextrapolation.map3dclasses import Map3D
from solarbextrapolation.visualisation_functions import visualise

###########################################################################
# Preprocessor
# Defining a trivial preprocessor that returns a zeros map for any given input
# map.
class PreZeros(Preprocessors):
    def __init__(self, map_magnetogram):
        super(PreZeros, self).__init__(map_magnetogram)

    def _preprocessor(self):
        # Adding in custom parameters to the meta
        self.meta['preprocessor_routine'] = 'Zeros Preprocessor'

        # Creating the trivial zeros map of the same shape as the input map
        map_output = mp.Map((np.zeros(self.map_input.data.shape),
                                    self.meta))

        # Outputting the map.
        return map_output

###########################################################################
# Make an input map that we will run the preprocessor on.
# This will be changed to using the sample HMI image.
#aMap2D = mp.Map('C://git//solarextrapolation//solarextrapolation//data//example_data_(100x100)__01_hmi.fits')
from solarbextrapolation.example_data_generator import generate_example_data, dummyDataToMap
import astropy.units as u
aMap2D = arr_Data = dummyDataToMap(generate_example_data([ 20, 20 ],u.Quantity([ -10.0, 10.0 ] * u.arcsec),u.Quantity([ -10.0, 10.0 ] * u.arcsec)), u.Quantity([ -10.0, 10.0 ] * u.arcsec), u.Quantity([ -10.0, 10.0 ] * u.arcsec))

###########################################################################
# Instansiate the preprocessor and process the input map.
aPrePro = PreZeros(aMap2D.submap([0, 10]*u.arcsec, [0, 10]*u.arcsec))
aPreProMap = aPrePro.preprocess()


###########################################################################
# You can plot the preprocessed map using peek.
aPreProMap.peek()

###########################################################################
# You can also access the metadata of the preprocessor like any map:
print "preprocessor_routine: " + str(aPreProMap.meta['preprocessor_routine'])
print "preprocessor_duration: " + str(aPreProMap.meta['preprocessor_duration'])





###########################################################################
# Extrapolator
# Defining a trivial extrapolator that returns a volume of one vectors.
class ExtOnes(Extrapolators):
    def __init__(self, map_magnetogram, **kwargs):
        super(ExtOnes, self).__init__(map_magnetogram, **kwargs)

    def _extrapolation(self):
        # Adding in custom parameters to the meta
        self.meta['extrapolator_routine'] = 'Ones Extrapolator'

        #arr_4d = np.ones([self.map_boundary_data.data.shape[0], self.map_boundary_data.data.shape[0], self.z, 3])
        arr_4d = np.ones(self.shape.tolist() + [3])
        return Map3D(arr_4d, self.meta)

###########################################################################
# Instansiate the preprocessor and extrapolate.
aExt = ExtOnes(aPreProMap, zshape=10)
aMap3D = aExt.extrapolate()

###########################################################################
# You can visulise the field using MayaVi.
fig = visualise(aMap3D,
                boundary=aPreProMap,
                show_boundary_axes=False,
                show_volume_axes=False,
                debug=False)
mlab.show()

"""

# aPreProData = aMap2D.submap([0,10], [0,10])

# Some checks:
#aPreProData.data # Should be a 2D zeros array.
#aPreProData.meta
#aPreProData.meta['preprocessor_routine']
#aPreProData.meta['preprocessor_start_time']

###########################################################################
# Testing an extrapolator


# Define trivial extrapolator
class ExtZeros(Extrapolators):
    def __init__(self, map_magnetogram, **kwargs):
        super(ExtZeros, self).__init__(map_magnetogram, **kwargs)

    def _extrapolation(self):
        # Adding in custom parameters to the meta
        self.meta['extrapolator_routine'] = 'Zeros Extrapolator'

        arr_4d = np.zeros([self.map_boundary_data.data.shape[0],
                           self.map_boundary_data.data.shape[0], self.z, 3])
        return Map3D((arr_4d, self.meta))


aExt = ExtZeros(
    aPreProData,
    filepath='C://Users/Alex/solarextrapolation/solarextrapolation/3Dmap.m3d')
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
"""