# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 18:01:09 2015

@author: alex_

Some basic examples of using the extrapolation API.
"""


###########################################################################
#                          Testing a preprocessor                         #
###########################################################################

# Define a trivial preprocessor
class PreZeros(Preprocessors):
    def __init__(self, map_magnetogram):
        super(PreZeros, self).__init__(map_magnetogram)

    def _preprocessor(self):
        # Adding in custom parameters to the meta
        self.meta['preprocessor_routine'] = 'Zeros Preprocessor'
        
        # Creating the trivial zeros map of the same shape as the input map
        map_output = sunpy.map.Map((np.zeros(self.map_input.data.shape), 
                                        self.meta))
            
        # Outputting the map.
        return map_output

# Generate the boundary
aMap2D = sunpy.map.Map('C://git//solarextrapolation//solarextrapolation//data//example_data_(100x100)__01_hmi.fits')
aPrePro = PreZeros(aMap2D.submap([0,10], [0,10]))
aPreProData = aPrePro.preprocess()
# aPreProData = aMap2D.submap([0,10], [0,10])
    
# Some checks:
#aPreProData.data # Should be a 2D zeros array.
#aPreProData.meta
#aPreProData.meta['preprocessor_routine']
#aPreProData.meta['preprocessor_start_time']

###########################################################################
#                          Testing an extrapolator                        #
###########################################################################
    
# Define trivial extrapolator
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
