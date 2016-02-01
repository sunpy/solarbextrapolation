

.. _sphx_glr_auto_examples_plot_gaussian_example_data.py:


=========================================
Generating Example Gaussian Boundary Data
=========================================

In this example you will be generating some example data and extrapolate this
using the basic potential extrapolator.

You can start by importing the necessary module components.


.. code-block:: python


    # Module imports
    from solarbextrapolation.example_data_generator import generate_example_data, dummyDataToMap







You also need the ability to convert astropyunits.


.. code-block:: python

    import astropy.units as u







You need to define the parameters of the eare, includsing the x and y ranges
as astropy quantities with angular or distance units and the grid shape.


.. code-block:: python


    # Input parameters:
    arr_grid_shape = [ 20, 22 ]         # [ y-size, x-size ]
    qua_xrange = u.Quantity([ -10.0, 10.0 ] * u.arcsec)
    qua_yrange = u.Quantity([ -11.0, 11.0 ] * u.arcsec)







The generated data will consist of a 2D space with 2 Gaussian spots, one
positive and one negative, on a background of 0.0.
solarbextrapolation.example_data_generator provides many ways to achieve this,
including letting it randomly generate the position, magnitude and size of
each spot/pole.


.. code-block:: python


    # To randomly generate 2 poles simply don't add any pole parameters:
    arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange)
    # Note: each time you run this pole positions/magnitudes will change.







We can now convert this into a a sunpy map object:


.. code-block:: python

    aMap = dummyDataToMap(arr_Data, qua_xrange, qua_yrange)







We can see this map using peek:


.. code-block:: python

    aMap.peek()




.. image:: /auto_examples\images\sphx_glr_plot_gaussian_example_data_001.png
    :align: center




 To manually position poles, simply build lists of parameters for each pole.
 It's often easiest to use percentage units for location/size, wheer we compare
 to the maps region.
arrA0 = [ Position, size, Max Magnitude ]


.. code-block:: python

    arrA0 = [ u.Quantity([ 25, 25 ] * u.percent), 10.0 * u.percent,  0.2 * u.T ]
    arrA1 = [ u.Quantity([ 75, 75 ] * u.percent), 10.0 * u.percent, -0.2 * u.T ]

    # To generate and view:
    arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange, arrA0, arrA1)
    aMap = dummyDataToMap(arr_Data, qua_xrange, qua_yrange)
    aMap.peek()




.. image:: /auto_examples\images\sphx_glr_plot_gaussian_example_data_002.png
    :align: center




But absolute positioning using the map range units is also possible


.. code-block:: python

    arrA2 = [ u.Quantity([ -6,  6 ] * u.arcsec), 2 * u.arcsec, -0.2 * u.T ]
    arrA3 = [ u.Quantity([  6, -7 ] * u.arcsec), 2 * u.arcsec,  0.2 * u.T ]

    # To generate and view:
    arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange, arrA2, arrA3)
    aMap = dummyDataToMap(arr_Data, qua_xrange, qua_yrange)
    aMap.peek()




.. image:: /auto_examples\images\sphx_glr_plot_gaussian_example_data_003.png
    :align: center




You can add as many poles as you want:


.. code-block:: python

    arr_Data = generate_example_data(arr_grid_shape, qua_xrange, qua_yrange, arrA0, arrA1, arrA2, arrA3)
    aMap = dummyDataToMap(arr_Data, qua_xrange, qua_yrange)
    aMap.peek()




.. image:: /auto_examples\images\sphx_glr_plot_gaussian_example_data_004.png
    :align: center




And being a map you can use all the normal SunPy functions, such as saving
the map using aMap.save(filepath).

**Total running time of the script:**
(0 minutes 0.911 seconds)



.. container:: sphx-glr-download

    **Download Python source code:** :download:`plot_gaussian_example_data.py <plot_gaussian_example_data.py>`


.. container:: sphx-glr-download

    **Download IPython notebook:** :download:`plot_gaussian_example_data.ipynb <plot_gaussian_example_data.ipynb>`
