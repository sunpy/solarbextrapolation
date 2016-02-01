

.. _sphx_glr_auto_examples_plot_potential_extrapolation_of_hmi_fits_file.py:


===========================
HMI FITS File Extrapolation
===========================

Example of extrapolating from a HMI fitts file using the potential
extrapolator and visualising.



.. code-block:: pytb

    Traceback (most recent call last):
      File "C:\Users\alex_\Anaconda\lib\site-packages\sphinx_gallery\gen_rst.py", line 467, in execute_script
        exec(code_block, example_globals)
      File "<string>", line 25, in <module>
      File "C:\Users\alex_\Anaconda\lib\site-packages\sunpy-0.6-py2.7.egg\sunpy\map\map_factory.py", line 239, in __call__
        data_header_pairs, already_maps = self._parse_args(*args, **kwargs)
      File "C:\Users\alex_\Anaconda\lib\site-packages\sunpy-0.6-py2.7.egg\sunpy\map\map_factory.py", line 202, in _parse_args
        raise ValueError("File not found or invalid input")
    ValueError: File not found or invalid input





.. code-block:: python

    # General imports
    import numpy as np
    import sunpy.map as mp
    from astropy import units as u
    from mayavi import mlab
    import os

    # Module imports
    from solarbextrapolation.map3dclasses import Map3D
    from solarbextrapolation.extrapolators import PotentialExtrapolator
    from solarbextrapolation.visualisation_functions import visualise

    # Cropping into the active region within the HMI map
    str_vol_filepath = 'C:\\git\\solarbextrapolation\\examples\\2011-02-14__20-35-25__02_Bxyz.npy'
    xrange = u.Quantity([50, 300] * u.arcsec)
    yrange = u.Quantity([-350, -100] * u.arcsec)
    zrange = u.Quantity([0, 250] * u.arcsec)
    xrangeextended = u.Quantity([xrange.value[0] - 50, xrange.value[1] + 50] *
                                xrange.unit)
    yrangeextended = u.Quantity([yrange.value[0] - 50, yrange.value[1] + 50] *
                                yrange.unit)

    # Open the map and create a cropped version for the extrapolation.
    map_hmi = mp.Map(
        'C:\\git\\solarbextrapolation\\examples\\2011-02-14__20-35-25__01_hmi.fits')
    map_hmi_cropped = map_hmi.submap(xrange, yrange)
    dimensions = u.Quantity([100, 100] * u.pixel)
    map_hmi_cropped_resampled = map_hmi_cropped.resample(dimensions,
                                                         method='linear')

    # Open the map and create a cropped version for the visualisation.
    #map_boundary = mp.Map('C:\\git\\solarbextrapolation\\examples\\2011-02-14__20-35-25__02_aia.fits') # For AIA
    map_boundary = mp.Map(
        'C:\\git\\solarbextrapolation\\examples\\2011-02-14__20-35-25__01_hmi.fits'
    )  # For HMI

    map_boundary_cropped = map_boundary.submap(xrangeextended, yrangeextended)

    # Only extrapolate if we don't have a saved version
    if not os.path.isfile(str_vol_filepath):
        aPotExt = PotentialExtrapolator(map_hmi_cropped_resampled,
                                        filepath=str_vol_filepath,
                                        zshape=dimensions[0].value,
                                        zrange=zrange)
        aMap3D = aPotExt.extrapolate()
    aMap3D = Map3D.load(str_vol_filepath)
    print('\nextrapolation duration: ' + str(np.round(aMap3D.meta['extrapolator_duration'], 3)) + ' s\n')

    # Visualise this
    visualise(aMap3D,
              boundary=map_boundary_cropped,
              scale=1.0 * u.Mm,
              boundary_unit=1.0 * u.arcsec,
              show_boundary_axes=False,
              show_volume_axes=True,
              debug=False)
    mlab.show()

**Total running time of the script:**
(0 minutes 0.000 seconds)



.. container:: sphx-glr-download

    **Download Python source code:** :download:`plot_potential_extrapolation_of_hmi_fits_file.py <plot_potential_extrapolation_of_hmi_fits_file.py>`


.. container:: sphx-glr-download

    **Download IPython notebook:** :download:`plot_potential_extrapolation_of_hmi_fits_file.ipynb <plot_potential_extrapolation_of_hmi_fits_file.ipynb>`
