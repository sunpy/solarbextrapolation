# -*- coding: utf-8 -*-
"""
=====================
Extrapolation Example
=====================

This is a full example
"""

import numpy as np
import sunpy.map as mp
from astropy import units as u

# Module Imports
from solarbextrapolation.classes import *
from potential_field_extrapolator import *
from utilities import *
from example_data_generator import *
from visualisation_functions import *

# 2015
str_folder = "C://fits//"
str_map_filepath = str_folder + "2015-01-04__19-41-12__02_aia.fits"
str_map_filepath = str_folder + "2015-01-04__19-41-12__01_hmi.fits"
str_vol_filepath = str_folder + "2015-01-04__19-41-12__03_Bxyz.npy"
str_ext_filepath = str_folder + "2015-01-04__19-41-12__04_extrapolator.ext"
lis_cropping = [[-200,200],[-250,150],[0,400], 'data']
xrange = u.Quantity([-200, 200] * u.arcsec)
yrange = u.Quantity([-250, 150] * u.arcsec)
zrange = u.Quantity([0,    400] * u.arcsec)

# 2014
str_map_filepath = str_folder + "2014-01-06__07-28-36__02_aia.fits"
str_map_filepath = str_folder + "2014-01-06__07-28-36__01_hmi.fits"
str_vol_filepath = str_folder + "2014-01-06__07-28-36__03_Bxyz.npy"
str_ext_filepath = str_folder + "2014-01-06__07-28-36__04_extrapolator.ext"
lis_cropping = [[-550,-200],[-245,105],[0,300], 'data']
xrange = u.Quantity([-550, -200] * u.arcsec)
yrange = u.Quantity([-245,  105] * u.arcsec)
zrange = u.Quantity([0,     300] * u.arcsec)

# 2011
str_map_filepath = str_folder + "2011-02-14__20-35-25__02_aia"
str_map_filepath = str_folder + "2011-02-14__20-35-25__01_hmi.fits"
str_vol_filepath = str_folder + "2011-02-14__20-35-25__03_Bxyz.npy"
str_ext_filepath = str_folder + "2011-02-14__20-35-25__04_extrapolator.ext"
lis_cropping = [[50,300],[-350,-100],[0,250], 'data']
xrange = u.Quantity([50,    300] * u.arcsec)
yrange = u.Quantity([-350, -100] * u.arcsec)
zrange = u.Quantity([0,     250] * u.arcsec)

# Open the map and create a cropped version for the visualisation.
aMap2D = mp.Map(str_map_filepath)
#aMap2D_cropped = aMap2D.submap(lis_cropping[0] * u.arcsec, lis_cropping[1] * u.arcsec)
aMap2D_cropped = aMap2D.submap(xrange, yrange)
dimensions = u.Quantity([30, 30] * u.pixel)
aMap2D_cropped_resampled = aMap2D_cropped.resample(dimensions, method='linear')
aMap2D_visulisation = aMap2D.submap([lis_cropping[0][0] - 50, lis_cropping[0][1] + 50] * u.arcsec, [lis_cropping[1][0] - 50, lis_cropping[1][1] + 50] * u.arcsec)

# Only extrapolate if we don't have a saved version
if not os.path.isfile(str_vol_filepath):
    aPotExt = PotentialExtrapolator(aMap2D_cropped_resampled, filepath=str_vol_filepath, zshape=30, zrange=zrange)
    aMap3D = aPotExt.extrapolate()
aMap3D = Map3D.load(str_vol_filepath)

# Visualise this
visualise(aMap3D, boundary=aMap2D_visulisation, scale=1.0*u.Mm, boundary_unit=1.0*u.arcsec, show_boundary_axes=False, show_volume_axes=True, debug=False)
