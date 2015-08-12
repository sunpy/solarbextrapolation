# -*- coding: utf-8 -*-
"""
=================================
Downloading HMI Data with SunPy
=================================

Before we can extrapolate a magnetic field, we need some boundary data!
"""
import astropy.units as u
import matplotlib.pyplot as plt

import sunpy.map
from sunpy.net import vso


vc = vso.VSOClient()

res = vc.query(vso.attrs.Time('2012/1/1T00:00:00','2012/1/1T00:01:00'),
               vso.attrs.Instrument('HMI'), vso.attrs.Physobs('LOS_magnetic_field'),
               vso.attrs.Sample(1*u.min))

print(res)

files = vc.get(res).wait()

print(files)

hmimap = sunpy.map.Map(files)
hmimap.peek()

plt.show()
