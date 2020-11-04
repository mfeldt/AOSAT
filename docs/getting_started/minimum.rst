===============
A Minimum Setup
===============

The key to AOSAT is the setup file. This file is needed for every analysis,
and it contains the description of your simulation that is needed in order to
perform meaningful analyses.

A detailed description of all setup file parameters and their possible defaults
is given in :doc:`Setup Files<../general_concept/setup>`.

To get started, assume you have run a simulation and it produced phase screens
in :abbr:`FITS (Flexible Image Transport System)` files
in a subdirectory of you current working directory. This subdirectory could e.g.
be called :command:`residual_phasescreens`. In addition, AOSAT needs to know about
the telescope pupil.  It, too, must be stored in a fits file, and the scale
in pixels per metre must naturally be the same as for the residual phase screens.
Let's assume this is named :command:`telescope_pupil.fits` and resides in your current directory.

Then, the first lines of your setup files are already defined:

.. code-block::

  pupilmask     = 'telescope_pupil.fits'
  screen_dir    = 'residual_phasescreens'

The locations are given relative to the setup file!
The setup file is read line by line and turned into a dictionary, it's not executed as python script!
In order to read the phase screens, AOSAT needs to know the filename pattern to
look for:

.. code-block::

  screen_fpattern  = '*residual_phase*.fits'

AOSAT will read all files matching the pattern that can be found in the :command:`screen_dir`
defined above.  The will be sorted in a somewhat intelligent way according to integer
numbers found in the filenames.

Individual files may contain cubes of data along the NAXIS3 of the FITS file, these will
be served sequentially one after the other as individual phase screens. Can be of variable lengths,
and can also be intermixed with flat files.


In addition, AOSAT needs to know about key properties, such as the pixel scale, and the
simulation frequency:

.. code::

  ppm      = 10.0 # pixels pr metre
  loopfreq = 1000.0 # Hz - actually it's the frequency of phase screens on disk

The last two mandatory parameters are the wavelength at which the analysis is desired,
and the unit the residual phases are given in. these can be any unit understood by
astropy.units.

.. code::

  phaseunit       = "micron"               # unit of phase screens
  an_lambda       = 3.0e-6                 # analysis wavelength in metres

The full minimum setup file looks like this then:

.. code::

  pupilmask        = 'telescope_pupil.fits'
  screen_dir       = 'residual_phasescreens'
  screen_fpattern  = '*residual_phase*.fits'
  ppm              = 10.0 # pixels pr metre
  loopfreq         = 1000.0 # Hz - actually it's the frequency of phase screens on disk
  phaseunit        = "micron"               # unit of phase screens
  an_lambda        = 3.0e-6                 # analysis wavelength in metres

create this file in your favourite editor, save it in the appropriate location
(remember, :command:`screen_dir` and :command:`pupilmask` are defined relative to the file's location!),
and use it to do the setup for your simulation:

.. code::

  >>> import aosat
  >>> aosat.aosat_cfg.CFG_SETTINGS = aosat_cfg.configure('my_file.setup')
