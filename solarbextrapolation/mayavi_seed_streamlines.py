# -*- coding: utf-8 -*-
"""
This module contains a custom streamlining class derived from the MayaVi2
streamlining class, modified to accept an array of seed points for visulaisation
using mayavi.

.. warning::
    The documentation for this class cannot be built on Read The Docs, it is possible to build it locally.

You can use this class thus:

Create a new Streamline instance and add it to a pipeline
"""

import numpy as np
from tvtk.api import tvtk
from traits.api import Instance, TraitPrefixList, Trait, Array

import mayavi
from mayavi.modules.streamline import Streamline

__all__ = ['SeedStreamline']

class SeedStreamline(Streamline):
    """
    This class is a modification of the mayavi Streamline class that accepts
    an array of seed points as a input rather than a widget.

    Examples
    --------
    Create a new Streamline instance and add it to a pipeline

    >>> from solarbextrapolation.mayavi_seed_streamlines import SeedStreamline
    >>> import numpy as np
    >>> seeds = [[1, 2, 5], [3, 4, 5]]
    >>> field_lines = SeedStreamline(seed_points = np.array(seeds)) 
    >>> myvectorfield.add_child(field_lines) #doctest: +SKIP
    """

    seed_points = Array(allow_none=False)
    seed = Instance(tvtk.PolyData, args=())
    update_mode = Trait('interactive', TraitPrefixList(['interactive',
                                                         'semi-interactive',
                                                         'non-interactive']),
                         desc='the speed at which the poly data is updated')

    def setup_pipeline(self):
        """Override this method so that it *creates* the tvtk
        pipeline.

        This method is invoked when the object is initialized via
        `__init__`.  Note that at the time this method is called, the
        tvtk data pipeline will *not* yet be setup.  So upstream data
        will not be available.  The idea is that you simply create the
        basic objects and setup those parts of the pipeline not
        dependent on upstream sources and filters.  You should also
        set the `actors` attribute up at this point.
        """
        # Create and setup the default objects.
        self.seed = tvtk.PolyData(points=self.seed_points)
        self.stream_tracer = tvtk.StreamTracer(maximum_propagation=2000,
                                               integration_direction='backward',
                                               compute_vorticity=False,
                                               integrator_type='runge_kutta4',
                                               )
        self.ribbon_filter = tvtk.RibbonFilter()
        self.tube_filter = tvtk.TubeFilter()

        self.actor = mayavi.components.actor.Actor()
        # Setup the actor suitably for this module.
        self.actor.property.line_width = 2.0

    def update_pipeline(self):
        """Override this method so that it *updates* the tvtk pipeline
        when data upstream is known to have changed.

        This method is invoked (automatically) when any of the inputs
        sends a `pipeline_changed` event.
        """
        mm = self.module_manager
        if mm is None:
            return

        src = mm.source
        self.stream_tracer.input = src.outputs[0]
        #self.seed.inputs = [src]

        # Setup the radius/width of the tube/ribbon filters based on
        # given input.
        if self._first:
            b = src.outputs[0].bounds
            l = [(b[1]-b[0]), (b[3]-b[2]), (b[5]-b[4])]
            length = np.sqrt(l[0]*l[0] + l[1]*l[1] + l[2]*l[2])
            self.ribbon_filter.width = length*0.0075
            self.tube_filter.radius = length*0.0075
            self._first = False

        self._streamline_type_changed(self.streamline_type)
        # Set the LUT for the mapper.
        self.actor.set_lut(mm.scalar_lut_manager.lut)

        self.pipeline_changed = True

    def _seed_points_changed(self, old, new):
        self.seed = tvtk.PolyData(points=self.seed_points)

    def _stream_tracer_changed(self, old, new):
        if old is not None:
            old.on_trait_change(self.render, remove=True)
        seed = self.seed
        if seed is not None:
            new.source = seed
        new.on_trait_change(self.render)
        mm = self.module_manager
        if mm is not None:
            new.input = mm.source.outputs[0]

        # A default output so there are no pipeline errors.  The
        # update_pipeline call corrects this if needed.
        self.outputs = [new.output]

        self.update_pipeline()

    def _seed_changed(self, old, new):
        st = self.stream_tracer
        if st is not None:
            st.source = new#.poly_data
        #self._change_components(old, new)
