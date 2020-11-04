===========
Frameserver
===========

AOSAT operates on residual phase screens. All input files need to be supplied
as FITS [#f0]_ -type images or cubes.  The two key inputs are

The pupil representing the telescope aperture
  AOSAT operates with transmission numbers, so if your 2D pupil map contains a 1.0
  at a given location, it's transparent, whereas if it contains a 0.0 it's opaque.
  intermediate values are welcome!

The residual phase screens
  A standing term in the AO simulation community, tis is a set of frames that each
  contain a 2D map of the residual wavefront phase at the corresponding time.

Both are 2D spatial representations, their scale in terms of pixels per metre must match,
as must the size of them, and of course their locations in the 2D plane. in short, the full
spatial coordinate system of the inputs must be identical.


Zero-padding
============

Depending on the simulation tool used, the size of the mask may just suffice the actual aperture,
or provide ample opaque space around it.  For the sake of the analyses performed, we
recommend to provide a zero padding such that individual frame sizes are a power of two, and
to be more precise it should be the power of two after the next one that's larger
than the pupil diameter in pixels.  So if the pupil diameter is 720 pixels, the frames
should have a dimension of 2048x2048 pixels (Next power of two: 1024, the one after: 2048).
This will make the FFTs used to transform to and from pupil and focal plane provide decent
sampling in both planes.  Larger sizes will only cost memory and speed.
If your ``pupilmask`` and residual frames (remember, their scales and sizes need to be identical!) do not
provide sufficient zero padding, you can enforce it by using the ``embed_frame = [2048,2048]``
(your mileage may vary) key in the setup file.

Reading a sequence of frames
=============================

AOSAT has an internal frameserver which reads the residual phase frames and serves them
to the analyzers one by one.  Frames should be equidistant timewise, their distance in time
used to derive the ``loopfreq`` key.

The internal frameserver is directed essentially by two keys in the setup file

.. code::

  screen_dir
  screen_fpattern

``screen_dir`` is the directory where to look for the frames (as always relative to the setup file!),
and ``screen_fpattern`` is a wildcard pattern that is applied to identify matching frames inside ``screen_dir``.

In order to provide some flexibility to read phase screens (frames) from different
simulation tools, the frameserver operates according to the following rules:

Matching files will be sorted accroding to an integer in the filename found in an intelligent way.
  This essentially means that your filenames should contain one, and only one number which is an integer,
  and is unique for every file. The order of these integers should reflect the temporal order of the frames,
  but does not have to be consecutive.

Matching files should contain one or more frames
  The FITS files matching the wildcard pattern may contain  2D or 3D data. If a 3D structure is found,
  it is assumed that it represents a temporal series of 2D frames stacked along NAXIS3.

In summary, the following sequence would be a good match for ``screen_fpattern = '*rwf*.fits'``:

.. code::

  mytool_rwf_1.fits
  mytool_rwf_2.fits
  mytool_rwf_03.fits
  mytool_rwf_004.fits
  mytool_rwf_5.fits   // this could be a cube of 100 frames
  mytool_rwf_200.fits // another cube of this time 200 frames
  mytool_rwf_1000.fits
  mytool_rwf_2000.fits

Yes, the files would be sorted in that order by the integer matching, and the frames in the cube would
served one by one, so the full sequence would have 306 frames.



Other frameserver options
=========================

The following options either also influence the frame server directly, or are indispensable for
identifying the frame properties:

loopfreq
  One over the temporal distance of subsequent frames

phaseunit
  The frameserver feeds the enalyzers with frames where the phase is given in unit of radian at
  the analysis wavelength. In order to do so, it must know the unit in which residual phase is given
  in the files. Anything astropy.units understands is fine.

start_skip
  Number of frames to skip before starting to serve to the analyzers. useful if the AO loop
  needed some time to stabilize before results become interesting.

skip_step
  can be set to i in oder to serve only every ith frame. usually discouraged, but can be useful
  for a quick test run ahead of the actual analysis.



Replacing the frameserver - tying into a simulation directly
============================================================

The frameserver is an integral part of AOSAT, but in the end it only serves 2D maps of residual
phases in units of radian to the individual analyzers.

Circumstances are thinkable where you might want to do that yourself, e.g. in order to serve
the frames directly out of a given simulation tool and not write them to disk at all.

In such a case you would first setup a simulation in the usual way, then serve frames in
an appropriate way from whatever tool you are using.

.. code-block::
  python

  import aosat
  from aosat.analyze import frg_analyzer, setup
  aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(my_file)
  sd=setup()
  analyzers=[frg_analyzer(sd)] # there may be others

Then, in your simulation main loop could look something like:

.. code-block::
  python

  def ao_loop(ao_param1, ao_param2, ..., analyzers):
    ##
    ## do the setup
    ##

    for loop_step in range(nsteps):
      ##
      ## do the control, produce residual_frame
      ## it should carry the residual phase
      ## in units of radian at the desired analysis
      ## wavelength
      ##
      for analyzer in analyzers:
        analyzer.feed_frame(residual_frame, nsteps)


    ##
    ## after the loop finished, evaluate the analyzers
    ##
    figurelist = []
    reportlist = []
    for analyzer in analyzers:
      analyzer.finalize()
      ##
      ## the following is a proposal, you could do whatever you want
      ##
      figurelist.append(analyzer.make_plot())
      reportlist.append(analyzer.make_report())


.. rubric:: Footnotes

.. [#f0] Flexible Image Transport System. Most simulation tools do or can easily
  produce this type of output which is simple-structured, human readable and
  flawlessly implemented in just about every astronomical piece of software.
  There is no plan to ever support HDF5.  `Here's why <https://cyrille.rossant.net/moving-away-hdf5/>`_.
