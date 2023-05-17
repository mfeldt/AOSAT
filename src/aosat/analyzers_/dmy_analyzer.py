import os


import sys

if sys.version_info >= (3, 8):
    from importlib.metadata import distributions as get_installed_distributions
else:
    import importlib_metadata.distributions as get_installed_distributions

if any(["cupy" in f.metadata["Name"] for f in get_installed_distributions()]):
    import cupy as np
else:
    import numpy as np

import scipy
    
from astropy import units
from astropy.io import fits as pyfits
import pandas as pd
import copy

from aosat import aosat_cfg
from aosat import fftx
from aosat import frameserver
from aosat import util
from scipy import ndimage
from poppy import zernike
from importlib import reload

import pdb
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.image as mpimg

from astropy.modeling import models, fitting


import logging
logger = logging.getLogger(__name__)

class dmy_analyzer():
    """Dummy analyzer.

    Analyzers do the actual analysis of residual wavefronts. The concept of AOSAT is
    to have one analyzer for each analysis task.  Analyzers are generally constructed,
    then they are fed frames one by one, followed by a 'finalize' stage.  After the
    analysis is completed, the 'make_report' and 'make_plot' methods supply the relevant
    information.

    This is a dummy class with no actually functionality for documentation purposes.
    It might be used to derive subclasses, but these can also be written independently.

    All analyzers should provide the methods of this one, and all methods
    should accept the parameters of the ones implemented here.  Accepting
    additional optional parameters is fine, but not supplying those should not imply
    fatal behaviour.



    Parameters
    ----------
    sd : setup dictionary
        created by setup()

    Attributes
    ----------

    sd
      the setup dictionary used for construction

    Additional attributes may be available depending on the actual analysis performed.

    """

    def __init__(self,sd):
        """Construct the analyzer

        Parameters
        ----------
        sd : setup dictionary
            created by setup()

        Returns
        -------
        nothing

        Examples
        -------

        >>> sd = setup()
        >>> dan = dmy_analyzer(sd)


        """
        self.sd = sd

    def feed_frame(self,frame,nframes):
        """Accept a single frame for analysis



        Parameters
        ----------
        frame : 2D residual wavefront
            residual wavefront of a single timestep of the AO simulation.
            These residual WF frames are fed one by one, calculations should
            take into account this circumstance.
        nframes: int
            total number of frames to be fed


        Returns
        -------
        nothing

        Examples
        -------

        >>> sd = setup()
        >>> dan = dmy_analyzer(sd)
        >>> rf  = np.zeros((1024,1024))
        >>> dan.feed_frame(rf,1)


        """
        pass


    def finalize(self):
        """Finish the analysis

        To be called after the last frame has been fed to execute
        final calculations

        Parameters
        ----------
        none

        Returns
        -------
        nothing

        Examples
        -------

        >>> sd=setup()
        >>> dan = dmy_analyzer(sd)
        >>> rf = np.zeros((1024,1024))
        >>> dan.feed_frame(rf,1)
        >>> dan.finalize()
        """
        pass

    def make_report(self):
        """Generate a report about the result.

        Parameters
        ----------
        None

        Returns
        -------
        report
            String containing the report

        Examples
        -------

        >>> sd=setup()
        >>> dan = dmy_analyzer(sd)
        >>> rf = np.zeros((1024,1024))
        >>> dan.feed_frame(rf,1)
        >>> dan.finalize()
        >>> dan.make_report()
        'This is a dummy with no actual functionality'
        """
        report = "This is a dummy with no actual functionality"
        return(report)

    def make_plot(self,fig=None,index='111',plotkwargs={},subplotkwargs={}):
        """Generate a plot about the result.

        Parameters
        ----------
        fig : matplolib figure
            If a figure pre-exists to which the plot should be attached,
            this should be provided.
        index : matplotlib subplot index code.
            If a pre-existing figure is supplied, this must be the index of the subplots
            to be attached to the figure. (the default is '111').
        plotkwargs : dictionary
            Keyword arguments supplied to the actual plot command (the default is {}).
        subplotkwargs : type
            Keyword arguments supplied to the add_subplot command (the default is {}).

        Returns
        -------
        matplotlib figure instance
            Either the pre-existing one with an attached subplot, or newly created

        Examples
        -------

        >>> sd=setup()
        >>> dan = dmy_analyzer(sd)
        >>> rf = np.zeros((1024,1024))
        >>> dan.feed_frame(rf,1)
        >>> dan.finalize()
        >>> rplot = dan.make_plot()


        """
        fig=plt.figure()
        return(fig)
