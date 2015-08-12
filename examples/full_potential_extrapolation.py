# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 11:26:13 2015

@author: alex_
"""

import numpy as np
import sunpy.map as mp
from sunpy.net import vso
from astropy import units as u

# Module Imports
from solarbextrapolation import *
#from potential_field_extrapolator import *
#from utilities import *
#from example_data_generator import *
#from visualisation_functions import *

if __name__ == '__main__':
    # Download the HMI data from VSO
    
    # create a new VSOClient instance
    client = vso.VSOClient()

    # build our query, this can return one item, or a list of them to DL (matching the filters).
    result_hmi = client.query(
        # The following are filters for what data I want.    
        vso.attrs.Time((2011, 2, 14, 20, 34, 0), (2011, 2, 14, 21, 0, 0)), # Time range.
        vso.attrs.Instrument('HMI'), # Helioseismic and Magnetic Imager.
        vso.attrs.Physobs('LOS_magnetic_field'), # Physical observables
        vso.attrs.Sample(4000 * u.s) # Only take a shot every $var seconds.
        # More observables at http://sdac.virtualsolar.org/cgi/show_details?keyword=PHYSOBS
    )
    
    
    result_aia = client.query(
        # The following are filters for what data I want.    
        vso.attrs.Time((2011, 2, 14, 20, 34, 0), (2011, 2, 14, 21, 0, 0)), # Time range.
        vso.attrs.Instrument('AIA'), # Helioseismic and Magnetic Imager.
        vso.attrs.Physobs('intensity'), # Physical observables
        vso.attrs.Sample(4000 * u.s) # Only take a shot every $var seconds.
        # More observables at http://sdac.virtualsolar.org/cgi/show_details?keyword=PHYSOBS
    )

    # Save the HMI and AIA data to a fits files.
    data_hmi = client.get(result_hmi, methods=('URL-FILE_Rice', 'URL-FILE')).wait()
    data_aia = client.get(result_aia, methods=('URL-FILE_Rice', 'URL-FILE')).wait()
    
    # Cropping into the active region within the HMI map
    str_vol_filepath = data_hmi[0][0:-5] + '_Bxyz.npy'
    xrange = u.Quantity([50,    300] * u.arcsec)
    yrange = u.Quantity([-350, -100] * u.arcsec)
    zrange = u.Quantity([0,     250] * u.arcsec)
    xrangeextended = u.Quantity([ xrange.value[0] - 50, xrange.value[1] + 50 ] * xrange.unit)
    yrangeextended = u.Quantity([ yrange.value[0] - 50, yrange.value[1] + 50 ] * yrange.unit)

    # Open the map and create a cropped version for the extrapolation.
    map_hmi = mp.Map(data_hmi[0])
    map_hmi_cropped = map_hmi.submap(xrange, yrange)
    dimensions = u.Quantity([20, 20] * u.pixel)
    map_hmi_cropped_resampled = map_hmi_cropped.resample(dimensions, method='linear')

    # Open the map and create a cropped version for the visualisation.
    map_boundary = mp.Map(data_hmi[0])
    print 'xrange: ' + str(map_boundary.xrange) + '\nyrange: ' + str(map_boundary.yrange)
    map_boundary = mp.Map(data_aia[0]).rotate()
    print 'xrange: ' + str(map_boundary.xrange) + '\nyrange: ' + str(map_boundary.yrange)
    
    
    map_boundary_cropped = map_boundary.submap(xrangeextended, yrangeextended)
    print 'xrange: ' + str(map_boundary_cropped.xrange) + '\nyrange: ' + str(map_boundary_cropped.yrange)
    
    # Only extrapolate if we don't have a saved version
    if not os.path.isfile(str_vol_filepath):
        aPotExt = PotentialExtrapolator(map_hmi_cropped_resampled, filepath=str_vol_filepath, zshape=20, zrange=zrange)
        aMap3D = aPotExt.extrapolate()
    aMap3D = Map3D.load(str_vol_filepath)    
    
    # Visualise this
    visualise(aMap3D, boundary=map_boundary_cropped, scale=1.0*u.Mm, boundary_unit=1.0*u.arcsec, show_boundary_axes=False, show_volume_axes=True, debug=False)
