===========
Setup Files
===========

Setup files are the key element to steer AOSAT analyses.  The setup file
contains key value pairs of the format

.. code::

  key = value

One per line.  The content of th efile is turned into an actual configuration plus
an analysis setup like so:

.. code::

  >>> aosat.aosat_cfg.CFG_SETTINGS = aosat_cfg.configure('my_file.setup')
  >>> sd = aosat.analyze.setup()

After this operation, :command:`sd` will contain the setup dictionary that the
individual analyzers need when being set up.

The following key/value pairs are being interpreted by AOSAT:

.. csv-table:: Setup file keys
  :widths: 1, 1, 2, 5
  :header-rows: 1

  Key, type, Meaning, Default
  an_lambda, float, Analysis wavelength [m], 3.0e-6
  aosat_logfile, string, "Name of log file [#f0]_", 'aosat.log'
  aosat_loglevel, string, "ERROR, WARNING, INFO, or DEBUG", 'INFO'
  divide_by_phase_mask, bool, Divide the residual phase by the mask [#f1]_, True
  embed_frame, "[int,int]", Size of support to embed frames in, --
  L0, float, Outer scale [m], 25.0
  loopfreq, float, Frequency of frames [Hz] (not the actual loop...), 1000.0
  ppm, float, Pixels per metre scale, 10.0
  phaseunit, string, "Unit of residual phase, anything in astropy.units", 'micron'
  pupilmask, string, Path to pupil mask file [#f0]_, 'examples/ExampleAnalyze/yao_pupil.fits'
  rm_glob_pist, bool, Remove global piston from residual phases, True
  screen_dir, string, Path to residual screens [#f0]_, 'examples/ExampleClosedLoop'
  screen_fpattern, string, Pattern to match to identify residual screen files, '\*rwf\*.fits'
  skipstep, int, "when i, every ith frame will be used ", 1
  startskip, int, No. of frames to skip at beginning, 0
  totalnumber, int, No. of frames to read (all if not set), --
  ts_basefilename, string, Basename for tearsheet output, 'ts_test'
  ts_title, string, Title for tearsheet, 'Example TS'
  zterms, int, No. of Zernike terms for analysis, 20

Notes on individual keys:

divide_by_phase_mask
  Some simulation tools multiply the pupil mask on to the residual phase to produce
  something resembling zero padding. In case the pupil mask is non-binary and
  contains grey values to represent apodizing or for anti-aliasing purposes, this
  must be divided out before analysis.

embed_frame
  Depending on your simulation tools, your pupil mask and residual frames may be zero padded
  or not. if they are not, it is recommended to embed them in arrays the size of which is the
  next after the next power of two to your frames. I.e. if your frames are e.g. 370 pixels
  across, and the aperture extends right to the edge, you should embed them in an array
  of 1024x1024 pixels.  This is achieved by setting ``embed_frame = [1024,1024]``

L0
  this is needed only if an informative plot of the open loop spectrum is desired by the
  spatial power spectrum analyzer.

loopfreq
  This should actually be one over the temporal separation of two consecutive frames.
  The actual loop frequency is not of interest.




.. rubric:: Footnotes

.. [#f0] relative to setup files
.. [#f1] for simulation tools that brutally multiply the mask onto the phase even though it may contain grey values
