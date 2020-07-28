=========
Analyzers
=========

All analyses performed by AOSAT are performed by so-called analyzers.
Each analyzer does its own topical analysis - what each one does exactly and how is described
in each ones own description:

* `Pupil fragmentation analyzer <../analyzers/frg_analyzer>`_
* `Residual phase analyzer <../analyzers/phs_analyzer>`_
* `PSF analyzer <../analyzers/psf_analyzer>`_
* `Spatial power spectrum analyzer <../analyzers/sps_analyzer>`_
* `Temporal variance contrast analyzer <../analyzers/sps_analyzer>`_
* `Zernike analyzer <../analyzers/zrn_analyzer>`_

The analyzer class
==================

There is a master class, called ``dmy_analyzer`` which implements all the
methods common to all analyzers.  In addition to these methods, each analyzer will
have a number of properties availabler after the ``finalize`` method has been called,
which contain the analysis results. Again, details are described for each analyzer
under the links given above.

In general, each analyzer implements the following methods:

.. code-block::
  python

  class dmy_analyzer():

      def __init__(self,sd):
          self.sd = sd

      def feed_frame(self,frame,nframes):
          pass

      def finalize(self):
          pass

      def make_report(self):
          report = "This is a dummy with no actual functionality"
          return(report)

      def make_plot(self,fig=None,index='111',plotkwargs={},subplotkwargs={}):
          fig=plt.figure()
          return(fig)

Analyzers are instantiated by passing a valid setup dictionary
`describing the simulation setup <setup>`_.
The calling sequence is generally ``feed_frame(frame, nframes)`` during the
analysis (where nframes is always the total number expected),  ``finalize()``
which produces the various properties that carry the result of the analysis,
and optionally ``make_plot()`` and ``make_report()`` afterwards.

``make_plot`` usually takes two sets of kwarg dictionaries, ``**plotkwargs`` and ``**subplotkwargs``.
``**plotkwargs`` are passed to the actual plotting command, while ``**subplotkwargs``
influence the setup of the (sub)plot produced - tick and axis labels,, sizes, etc.

Writing your own analyzer
=========================

When you require a particular analysis not provided by an existing analyzer, it's
usually a good idea to write your own.  Doing so is easily facilitated by AOSAT.
Generally, derive your new class from ``dmy_analyzer`` and overwrite the methods required.

In general you definitely need the ``feed_frame``, which performs the actual analysis of
each frame as it comes in, and the ``finalize``, which will do some conclusive analysis once
the set is complete and fill the properties defined early on.

So if you were to rite e.g. an analyzer that simply tracks the peak-to-valley of the residual
phase, you could define the new analyzer like so:

.. code-block::
  python

  class ptv_analyzer(dmy_analzer):

  def __init__(self,sd):
    self.sd       = sd
    self.ptv_max  = None
    self.ptv_mean = None
    self.ptv_std  = None
    self.ptv_arr  = None
    self._ffed    = 0

The ``_ffed`` property is used by many of AOSAT's internal analyzers
to track the number of frames fed so far. It is incerased by 1 every time
``feed_frame`` is called.

The feed_frame method could then look like:

.. code-block::
  python

  def feed_frame(self,frame,nframes):
    if self._ffed == 0:
      self.ptv_arr = np.zeros(nframes) # make array to hold PTV values
    self.ptv_arr[self._ffed] = np.max(frame[self.sd['wnz']]) - np.min(frame[self.sd['wnz']])

You will notice the analyzer makes use of a useful property stored in the setup dictionary,
in this case 'wnz', which holds the indices of the pixels where the pupil transmission is
not zero (wnz = where not zero). The setup dictionary holds a ton of these useful things, to find out
what is there you are kindly invited to take a look at the sources, or to set a breakpoint
in your analyzer and do a :command:`dir(self.sd)`. Most things are named appropriately...

Once the analysis is over it is time to finalize:

.. code-block::
  python

  def finalize(self):
    self.ptv_max  = np.max(self.ptv_arr)
    self.ptv_mean = np.mean(self.ptv_arr)
    self.ptv_std  = np.std(self.ptv_arr)

In principle, your analyzer is done now.  The properties available after a run are
the time series itself (in ``ptv_array``), the maximum value (in ``ptv_max``),
the mean, and the standard deviation (guess where...).

If you are kind, you could provide a plot, and a textual report.  The ``make_plot()``
method is expected to return a matplotlib figure instance, but actually you should check
whether a figure already exists and insert a subplot at a desired index if appropriate.
AOSAT's internal analyzers achieve this as follows:

.. code-block::
  python

  def make_plot(self,fig=None,index='111',plotkwargs={},subplotkwargs={}):
    if fig is None:
      fig = plt.figure()

    ##
    ## default appearance
    ##
    if 'xlabel' not in subplotkwargs:
        subplotkwargs['xlabel'] = r'time [ms]'
    if 'ylabel' not in subplotkwargs:
        subplotkwargs['ylabel'] = r'PtV [rad]'
    if 'title' not in subplotkwargs:
        subplotkwargs['title'] = 'Peak-to-valley over time'

    if 'color' not in plotkwargs:
        plotkwargs['color']='blue'

    tv = np.arange(len(self.ptv_arr))/self.sd['cfg']['loopfreq']*1000 # temporal vector in ms
    ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))
    ax.plot(tv,self.ptv_arr,**plotkwargs)

    return(fig)

Subplot positioning via the index argument needs to be taken care of from the outside.

Imports and Supporting CUDA
---------------------------

In the header of your analyzer module you would of course do all the necessary
imports you need.

.. code-block::
  python

  from aosat import aosat_cfg
  from aosat import fftx
  from aosat import frameserver
  from aosat import util

  import scipy
  ...

In order to support GPU operation, it is recommended to do it the same way as the
internal analyzers:

.. code-block::
  python

  from pip._internal.utils.misc import get_installed_distributions
  if any(["cupy" in str(f) for f in get_installed_distributions()]):
      import cupy as np
  else:
      import numpy as np

This way, your statements like :command:`np.mean` will always work, whether there's
CUDA support or not. If you use arrays returned by an external library, chances are it's
a numpy array, so in order to ensure it get's transferred to the GPU if neede, use np.array:

.. code-block::
  python

  new_array = np.array(somelib.function(util.ensure_numpy(orig_array)))

Here, ``orig_array`` and ``new_array`` are on the GPU if CUDA support is active,
and the ``somelib.function()`` receives and returns numpy arrays. In case of
no CUDA support, nothing is changed.

Note that calling ``util.ensure_numpy()`` and ``np.array()`` forces a sync, so
a lot of these will slow down operation.  The goal of using CUDA is of course to
do *all* array operations on the GPU. if you rely heavily on external libraries,
it might be better to not use CUDA support.
