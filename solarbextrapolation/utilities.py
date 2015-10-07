# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 15:32:51 2015

@author: alex_
"""

from astropy import units as u
import numpy as np
import sunpy.map as mp
from copy import deepcopy

__all__ = ["decompose_ang_len", "si_this_map_OLD", "si_this_map"]

def decompose_ang_len(qua_input, **kwargs):
    """
    Function to help decompose quantities that have an equivilence between angles
    and length, such as photospheric observational angles and object sizes.
    The function uses an equivilence to convert into either length or angle
    physical_type units.

    Parameters
    ----------

    qua_input : `astropy.units.quantity.Quantity`
        The quantity you wish decomposed.

    working_units : `astropy.units.quantity.Quantity`, optional
        Unit that will be used for internal working and the returned quantity.
    Ensure that it is of the correct physical type, angle or length.

    equivalencies : astropy equivilence,
        Equivilence used to relate the length and angle units.

    """
    # Parameters
    working_units = kwargs.get('working_units', u.m) * 1.0
    equivalence = kwargs.get('equivalencies', u.dimensionless_angles())

    # Do nothing if the input is dimensionless.
    if qua_input.unit is u.Quantity(1.0).unit:
        return qua_input.decompose()
    else:
        # Components of the quantity
        value = qua_input.value
        length_unit = 0.0 * u.m
        length_exponent = 0.0
        angle_unit = 0.0 * u.radian
        angle_exponent = 0.0

        # If we have at least 1 base, populate from the first base
        if len(qua_input.unit.bases) > 0:
            if qua_input.unit.bases[0].physical_type is u.m.physical_type:
                length_unit     = 1.0 * qua_input.unit.bases[0]
                length_exponent = qua_input.unit.powers[0]

                # convert to SI (meter here)
                length_unit = length_unit.to(u.m)
            elif qua_input.unit.bases[0].physical_type is u.radian.physical_type:
                angle_unit     = 1.0 * qua_input.unit.bases[0]
                angle_exponent = qua_input.unit.powers[0]

                # Convert to SI (radian here)
                angle_unit = angle_unit.to(u.radian)

        # If we have 2 bases, populate from the second base
        if len(qua_input.unit.bases) > 1:
            if qua_input.unit.bases[1].physical_type is u.m.physical_type:
                length_unit     = 1.0 * qua_input.unit.bases[1]
                length_exponent = qua_input.unit.powers[1]

                # Convert to SI (meter here)
                length_unit = length_unit.to(u.m)
            elif qua_input.unit.bases[1].physical_type is u.radian.physical_type:
                angle_unit     = 1.0 * qua_input.unit.bases[1]
                angle_exponent = qua_input.unit.powers[1]

                # Convert to SI (radian here)
                angle_unit = angle_unit.to(u.radian)

        # Convert the incompatible base to the working units using the equivilence
        if working_units.unit.physical_type is u.m.physical_type:
            angle_unit = angle_unit.to(working_units, equivalencies=equivalence)

            # Strip out the units, so the output doesn't have squared lenth units
            #angle_unit = angle_unit.value # Kept in-case it causes bugs
        elif working_units.unit.physical_type is u.radian.physical_type:
            length_unit = length_unit.to(working_units, equivalencies=equivalence)

            # Strip out the units, so the output doesn't have squared length units
            #length_unit = length_unit.value # Kept in-case it causes bugs
        # The quantity to return
        quantity =  value * length_unit ** length_exponent * angle_unit ** angle_exponent
        # Change to the working unit if not dimensionless
        if quantity.unit.physical_type is not (u.m / u.m).decompose().physical_type:
            quantity.to(working_units)
        return quantity.decompose()


def si_this_map_OLD(map):
    """
    Basic function to create a deep copy of a map but with all units in SI.
    """
    # Find out the value units and convert this and data to SI
    units = 1.0 * u.Unit(map.meta['bunit']).to(u.Tesla) * u.Tesla
    data = deepcopy(map.data) * units.value

    # ATM I don't convert the x-axis and y-axis to SI

    # Modify the map header to reflect all these changes
    meta = deepcopy(map.meta)
    meta['bunit']   = units.unit
    meta['datamax'] = data.max()
    meta['datamin'] = data.min()
    #meta['cdelt1'] = 0.504295 # Following modified if we convert x/y-axes
    #meta['cdelt2'] = 0.504295
    #meta['cunit1'] = 'arcsec'
    #meta['cunit2'] = 'arcsec'
    #meta['crpix1'] = data.shape[1] / 2.0 + 0.5, # central x-pixel
    #meta['crpix2'] = data.shape[0] / 2.0 + 0.5, # cnetral y-pixel
    #meta['CRVAL1'] = 0.000000
    #meta['CRVAL2'] = 0.000000

    # Return the modified map
    return mp.Map((data, meta))


def si_this_map(map):
    """
    Basic function to create a deep copy of a map but with all units in SI.
    """
    # Find out the value units and convert this and data to SI
    units = 1.0 * u.Unit(map.meta['bunit']).to(u.Tesla) * u.Tesla
    data = deepcopy(map.data) * units.value

    # Setup the arc to length equivilence
    obs_distance = map.dsun - map.rsun_meters
    radian_length = [ (u.radian, u.meter, lambda x: obs_distance * x, lambda x: x / obs_distance) ]

    # Convert the x-axis and y-axis to SI
    cdelt1 = (float(map.meta['cdelt1']) * u.Unit(map.meta['cunit1'])).to(u.meter, equivalencies=radian_length)
    cdelt2 = (float(map.meta['cdelt2']) * u.Unit(map.meta['cunit2'])).to(u.meter, equivalencies=radian_length)
    crpix1 = (float(map.meta['crpix1']) * u.Unit(map.meta['cunit1'])).to(u.meter, equivalencies=radian_length)
    crpix2 = (float(map.meta['crpix2']) * u.Unit(map.meta['cunit2'])).to(u.meter, equivalencies=radian_length)

    # Modify the map header to reflect all these changes
    meta = deepcopy(map.meta)
    meta['bunit']   = 'Tesla' #units.unit
    meta['datamax'] = data.max()
    meta['datamin'] = data.min()
    # Following modified if we convert x/y-axes
    meta['cdelt1'] = str(cdelt1.value)
    meta['cdelt2'] = str(cdelt2.value)
    meta['cunit1'] = str(cdelt1.unit)
    meta['cunit2'] = str(cdelt2.unit)
    meta['crpix1'] = str(crpix1.value)
    meta['crpix2'] = str(crpix2.value)
    #meta['CRVAL1'] = 0.000000 # Reference data coordinates
    #meta['CRVAL2'] = 0.000000

    # Return the modified map
    return mp.Map((data, meta))




