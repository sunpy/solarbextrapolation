# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 15:32:51 2015

@author: alex_
"""

from astropy import units as u
import numpy as np
import sunpy.map as mp
from copy import deepcopy

def decompose_ang_len(qua_input, **kwargs):
    """
    Function to help decompose quantities that have an equivilence between angles
    and length, such as photospheric observational angles and object sizes.
    The function uses an equivilence to cobnvert into eithe length or angle
    physical_type units.

    Parameters
    ----------
    working_units : astropy unit, a unit that will be used for internal working
    and the returned quantity. Ensure that it is of the correct physical type,
    angle or length.
        A 2d list or ndarray containing the map data
    equivalencies : astropy equivilence, an equivilence used to relate the
    length and angle units.
    
    """
    # Parameters
    working_units = kwargs.get('working_units', u.m)
    equivalence = kwargs.get('equivalencies', u.dimensionless_angles())

    # Components of the quantity
    value = qua_input.value
    length_unit = 0.0 * u.m
    length_exponent = 0.0
    angle_unit = 0.0 * u.radian
    angle_exponent = 0.0
    
    # If we have at least 1 base, populate from the first base
    if len(qua_input.unit.bases) > 0:
        if qua_input.unit.bases[0].physical_type is u.m.physical_type:
            length_unit = 1.0 * qua_input.unit.bases[0]
            length_exponent = qua_input.unit.powers[0]
            
            # Convert to SI (meter here)
            length_unit = length_unit.to(u.m)
        elif qua_input.unit.bases[0].physical_type is u.radian.physical_type:
            angle_unit = 1.0 * qua_input.unit.bases[0]
            angle_exponent = qua_input.unit.powers[0]
            
            # Convert to SI (radian here)
            angle_unit = angle_unit.to(u.radian)
    
    # If we have 2 bases, populate from the second base
    if len(qua_input.unit.bases) > 1:
        if qua_input.unit.bases[1].physical_type is u.m.physical_type:
            length_unit = 1.0 * qua_input.unit.bases[1]
            length_exponent = qua_input.unit.powers[1]
            
            # Convert to SI (meter here)
            length_unit = length_unit.to(u.m)
        elif qua_input.unit.bases[1].physical_type is u.radian.physical_type:
            angle_unit = 1.0 * qua_input.unit.bases[1]
            angle_exponent = qua_input.unit.powers[1]
            
            # Convert to SI (radian here)
            angle_unit = angle_unit.to(u.radian)
        
    # Convert the incompatible base to the working units using the equivilence
    if working_units.physical_type is u.m.physical_type:
        angle_unit = angle_unit.to(working_units, equivalencies=equivalence)
        
        # Strip out the units, so the output doesn't have squared lenth units
        #angle_unit = angle_unit.value
    elif working_units.physical_type is u.radian.physical_type:
        length_unit = length_unit.to(working_units, equivalencies=equivalence)

        # Strip out the units, so the output doesn't have squared length units
        #length_unit = length_unit.value
    # The quantity to return
    quantity =  value * length_unit ** length_exponent * angle_unit ** angle_exponent
    # Change to the working unit if not dimensionless
    if quantity.unit.physical_type is not (u.m / u.m).decompose().physical_type:
        quantity.to(working_units)
    return quantity.decompose()
    

def si_this_map(map):
    # Find out the value units and convert this and data to SI
    units = 1.0 * u.Unit(map.meta['bunit']).to(u.Tesla) * u.Tesla
    print units
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



