============
PSF Analyzer
============

Description
===========

Analyzing the :abbr:`PSF (Point Spread Function)` is done in a straight forward way: The residual phase of stored in each screen is transformed into a PSF by means of the well known

.. code-block:: python

  np.abs(fftForward(aperture*np.exp(1j*phase)))**2

``aperture`` is represented by the telescope mirror defined in the :doc:`setup file <../general_concept/setup>`.
A number of additional analyses is carried out in the resulting PSF, and the input phase.

* the Strehl ratio is measured on the phase :math:`\phi` as :math:`S = e^{-\sigma_\phi^2}`.
* the Strehl ratio is measured on the PSF as :math:`S = I_{peak}/I_{ref, peak}`, where :math:`I` is the PSF's intensity distribution, and :math:`I_{ref}` is the intensity distribution of a reference PSF resulting from a perfectly flat wavefront.
* The tip and tilt excursion of the PSF is measured on the phase :math:`\phi` by means of a least squares fit of a flat wavefront to the phase
* The tip and tilt excursion of the PSF is measured on the PSF by means of fitting a 2D Gaussian to it.



Plot caption
============

When called on its own, or on a figure with sufficient available subplot space, ``frg_anaylzer.makeplot()`` will produce a figure like so:

.. image:: psf.png
  :width: 50%

The caption for the figure could be:

*Resulting time-averaged PSF in units of peak intensity.  Additionally shown is the Strehl ratio derived from the peak intensity, denoted as "SR(PSF)", and derived from the wavefront quality, denoted as "SR(WF)". The tip-tilt statistics are shown in the lower part, also derived from the PSF directly as well as from the wavefronts.*  :math:`Q_{90}` *means that 90% of TT excursions are smaller than the quoted value.*

Resulting properties
====================

``psf_analyzer`` exposes the following properties after ``psf_analyzer.finalize()`` has been called:

.. csv-table:: psf_analyzer porperties
  :widths: 1, 3, 5
  :header-rows: 1

  Property, type, Explanation
  **psf**, 2D float array, image of the last psf analyzed
  **psf_mean**, 2D float array, image of time-averaged PSF
  **sr_wf**, float, Strehl ratio determined from individual wavefronts
  **ttx**, 1D float array of length n_frames, vector of individual tip excursions from PSF (in milli arcsec)
  **tty**, 1D float array of length n_frames, vector of individual tilt excursions from PSF (in milli arcsec)
  **ttilt**, 1D float array of length n_frames, vector of individual total excursions from WF (in milli arcsec)
  **ttjit**, float, standard deviation of **ttilt**
  **ttq90**, float, 90% quantile of **ttilt**
  **ttjit_psf**, float, standard deviation of sqrt(**ttx**\*\*2 + **tty**\*\*2)
  **ttq90_psf**, float, 90% quantile of sqrt(**ttx**\*\*2 + **tty**\*\*2)
