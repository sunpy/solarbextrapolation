

.. _sphx_glr_auto_examples_plot_define_and_run_trivial_preprocessor_and_extrapolator.py:


===============================================
Defining a Custom Preprocessor and Extrapolator
===============================================

Here you will be creating trivial preprocessor and and exztrqapolatoirs
following the API.


You start by importing the necessary modules.


.. code-block:: python


    # General imports
    import sunpy.map as mp
    import numpy as np
    from mayavi import mlab # Necessary for visulisation

    # Module imports
    from solarbextrapolation.preprocessors import Preprocessors
    from solarbextrapolation.extrapolators import Extrapolators
    from solarbextrapolation.map3dclasses import Map3D
    from solarbextrapolation.visualisation_functions import visualise







Preprocessor
Defining a trivial preprocessor that returns a zeros map for any given input
map.


.. code-block:: python

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







 Make an input map that we will run the preprocessor on.
 This will be changed to using the sample HMI image.
aMap2D = mp.Map('C://git//solarextrapolation//solarextrapolation//data//example_data_(100x100)__01_hmi.fits')


.. code-block:: python

    from solarbextrapolation.example_data_generator import generate_example_data, dummyDataToMap
    import astropy.units as u
    aMap2D = arr_Data = dummyDataToMap(generate_example_data([ 20, 20 ],u.Quantity([ -10.0, 10.0 ] * u.arcsec),u.Quantity([ -10.0, 10.0 ] * u.arcsec)), u.Quantity([ -10.0, 10.0 ] * u.arcsec), u.Quantity([ -10.0, 10.0 ] * u.arcsec))







Instansiate the preprocessor and process the input map.


.. code-block:: python

    aPrePro = PreZeros(aMap2D.submap([0, 10]*u.arcsec, [0, 10]*u.arcsec))
    aPreProMap = aPrePro.preprocess()








You can plot the preprocessed map using peek.


.. code-block:: python

    aPreProMap.peek()




.. image:: /auto_examples\images\sphx_glr_plot_define_and_run_trivial_preprocessor_and_extrapolator_001.png
    :align: center




You can also access the metadata of the preprocessor like any map:


.. code-block:: python

    print "preprocessor_routine: " + str(aPreProMap.meta['preprocessor_routine'])
    print "preprocessor_duration: " + str(aPreProMap.meta['preprocessor_duration'])








.. code-block:: pytb

    Traceback (most recent call last):
      File "C:\Users\alex_\Anaconda\lib\site-packages\sphinx_gallery\gen_rst.py", line 467, in execute_script
        exec(code_block, example_globals)
      File "<string>", line 1
        print "preprocessor_routine: " + str(aPreProMap.meta['preprocessor_routine'])
                                     ^
    SyntaxError: invalid syntax




Extrapolator
Defining a trivial extrapolator that returns a volume of one vectors.


.. code-block:: python

    class ExtOnes(Extrapolators):
        def __init__(self, map_magnetogram, **kwargs):
            super(ExtOnes, self).__init__(map_magnetogram, **kwargs)

        def _extrapolation(self):
            # Adding in custom parameters to the meta
            self.meta['extrapolator_routine'] = 'Ones Extrapolator'

            #arr_4d = np.ones([self.map_boundary_data.data.shape[0], self.map_boundary_data.data.shape[0], self.z, 3])
            arr_4d = np.ones(self.shape.tolist() + [3])
            return Map3D(arr_4d, self.meta)







Instansiate the preprocessor and extrapolate.


.. code-block:: python

    aExt = ExtOnes(aPreProMap, zshape=10)
    aMap3D = aExt.extrapolate()







You can visulise the field using MayaVi.


.. code-block:: python

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




.. code-block:: pytb

    Traceback (most recent call last):
      File "C:\Users\alex_\Anaconda\lib\site-packages\sphinx_gallery\gen_rst.py", line 467, in execute_script
        exec(code_block, example_globals)
      File "<string>", line 17
        """

    # aPreProData = aMap2D.submap([0,10], [0,10])

    # Some checks:
    #aPreProData.data # Should be a 2D zeros array.
    #aPreProData.meta
    #aPreProData.meta['preprocessor_routine']
    #aPreProData.meta['preprocessor_start_time']
       

                                             

              
                                               
                 
                                         
                                               ^
    SyntaxError: EOF while scanning triple-quoted string literal




Testing an extrapolator


.. code-block:: python



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


.. code-block:: pytb

    Traceback (most recent call last):
      File "C:\Users\alex_\Anaconda\lib\site-packages\sphinx_gallery\gen_rst.py", line 467, in execute_script
        exec(code_block, example_globals)
      File "<string>", line 36
        """
          ^
    SyntaxError: EOF while scanning triple-quoted string literal




**Total running time of the script:**
(0 minutes 0.261 seconds)



.. container:: sphx-glr-download

    **Download Python source code:** :download:`plot_define_and_run_trivial_preprocessor_and_extrapolator.py <plot_define_and_run_trivial_preprocessor_and_extrapolator.py>`


.. container:: sphx-glr-download

    **Download IPython notebook:** :download:`plot_define_and_run_trivial_preprocessor_and_extrapolator.ipynb <plot_define_and_run_trivial_preprocessor_and_extrapolator.ipynb>`
