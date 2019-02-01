# -*- coding: utf-8 -*-
"""
===============================================
Defining and Run a Custom Analytical Model
===============================================

Here you will be creating trivial analytical model following the API.

"""

##############################################################################
# You can start by importing the necessary module components.

# Module imports
from solarbextrapolation.map3dclasses import Map3D
from solarbextrapolation.analyticalmodels import AnalyticalModel
from solarbextrapolation.visualisation_functions import visualise

##############################################################################
# You also need the ability to convert astropyunits, manipulate numpy arrays
# and use MayaVi for visualisation.

# General imports
import astropy.units as u
import numpy as np
from mayavi import mlab


##############################################################################
# You are going to try and define a 3D cuboid grid of 20x22x20 with ranges in
# arcseconds, these parameters can be stored in the following lists and astropy
# quantities.

# Input parameters:
qua_shape   = u.Quantity([   20,   20,  20] * u.pixel)
qua_x_range = u.Quantity([ -80.0,  80 ] * u.Mm)
qua_y_range = u.Quantity([ -80.0,  80 ] * u.Mm)
qua_z_range = u.Quantity([   0.0, 120 ] * u.Mm)

##############################################################################
# From the above parameters you can derive the grid step size and total size in
# each dimension.

"""
# Derived parameters (make SI where applicable)
x_0 = x_range[0].to(u.m).value
Dx = (( x_range[1] - x_range[0] ) / ( tup_shape[0] * 1.0 )).to(u.m).value
x_size = Dx * tup_shape[0]
y_0 = y_range[0].to(u.m).value
Dy = (( y_range[1] - y_range[0] ) / ( tup_shape[1] * 1.0 )).to(u.m).value
y_size = Dy * tup_shape[1]
z_0 = z_range[0].to(u.m).value
Dz = (( z_range[1] - z_range[0] ) / ( tup_shape[2] * 1.0 )).to(u.m).value
z_size = Dy * tup_shape[2]
"""

##############################################################################
# You can define this analytical model as a child of the AnalyticalModel class.
class AnaOnes(AnalyticalModel):
    def __init__(self, **kwargs):
        super(AnaOnes, self).__init__(**kwargs)

    def _generate_field(self, **kwargs):
        # Adding in custom parameters to the metadata
        self.meta['analytical_model_routine'] = 'Ones Model'

        # Generate a trivial field and return (X,Y,Z,Vec)
        arr_4d = np.ones(self.shape.value.tolist() + [3])
        self.field = arr_4d

        # Extract the LoS Magnetogram from this:
        self.magnetogram.data = arr_4d[:,:,0,2]

        # Now return the vector field.
        return Map3D( arr_4d, self.meta )


##############################################################################
# You can instansiate a copy of the new  analytical model.
aAnaMod = AnaOnes(shape=qua_shape, xrange=qua_x_range, yrange=qua_y_range, zrange=qua_z_range)

##############################################################################
# Note: you could use default ranges and grid shape using aAnaMod = AnaOnes().

##############################################################################
# You can now calculate the vector field.
aMap3D = aAnaMod.generate()

##############################################################################
# You can now see the 2D boundary data used for extrapolation.
aMap2D = aAnaMod.to_los_magnetogram()
aMap2D.peek()

##############################################################################
# You also visulise the 3D vector field:
fig = visualise(aMap3D,
                show_boundary_axes=False,
                show_volume_axes=False,
                debug=False)
mlab.show()

# Note: you can add boundary axes using:
"""
fig = visualise(aMap3D,
                show_boundary_axes=False,
                boundary_units=[1.0*u.arcsec, 1.0*u.arcsec],
                show_volume_axes=True,
                debug=False)
"""
