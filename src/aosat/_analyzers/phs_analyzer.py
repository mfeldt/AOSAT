

import os

from pip._internal.utils.misc import get_installed_distributions
if any(["cupy" in str(f) for f in get_installed_distributions()]):
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


class phs_analyzer():
    """Analyze phase error

    Parameters
    ----------
    sd : dictionary
        Analysis setup dictionary

    Examples
    -------
    >>> sd={}
    >>> an=phs_analyzer

    Attributes
    ----------
    Available after calling finalize()!
    rms : float
        mean RMS of all wavefronts in nm.
    lastphase : 2D array
        Last residual phasescreen (in nm)

    """
    def __init__(self,sd):
        self.sd        = sd
        self.rmst      = None
        self.rms       = 0.0
        self.lastphase = None
        self._ffed = 0

    def feed_frame(self,frame,nframes):
        if self._ffed == 0:
            self.rmst = np.zeros(nframes)
        if self._ffed == (nframes-2):
            self.lastphase = frame/2/np.pi*self.sd['cfg']['an_lambda']*1e9 # in nm
        self.rmst[self._ffed] = np.std(frame[self.sd['wnz']]/2/np.pi*self.sd['cfg']['an_lambda']*1e9) # in nm
        self._ffed += 1

    def finalize(self):
        self.rms=np.mean(self.rmst)

    def make_report(self):
        report =  "##\n##\n"
        report += "## reporting analyzer: %s\n##\n##\n" % self.__class__.__name__
        report += "## Mean RMS of wavefront:\n\n"
        report += "std_wf    = %s # nm\n\n\n" % self.rms
        return(report)

    def make_plot(self,fig=None,index=111,plotkwargs={},subplotkwargs={}):
        if fig is None:
            fig = plt.figure()

        ## plot region
        sdim  = self.lastphase.shape[0]
        x, y  = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]
        r     = (x**2+y**2)**0.5
        plts  = int(min([max((r*self.sd['tel_mirror']).flatten())*1.1,sdim/2])/2)*2
        sreg = slice(int(sdim/2-plts),int(sdim/2+plts))

        ## what to show
        mphase = np.min(self.lastphase[self.sd['wnz']])
        sphase = copy.copy(self.lastphase)
        #sphase[self.sd['wnz']] -=  mphase
        pp     = sphase[self.sd['wnz']].flatten()
        pmin  = np.percentile(pp,5)
        pmax  = np.percentile(pp,95)
        sphase[self.sd['tel_mirror']==0] = pmin
        sshow = sphase[sreg,sreg]


        ##
        ## default appearance
        ##
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = r'$\delta$x [m]'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = r'$\delta$y [m]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = 'Last Residual Phase\n '


        if 'cmap' not in plotkwargs:
            plotkwargs['cmap'] = 'nipy_spectral'
        if 'extent' not in plotkwargs:
            plotkwargs['extent']=[-(1.0*plts)/self.sd['ppm'],(1.0*plts)/self.sd['ppm'],-(1.0*plts)/self.sd['ppm'],(1.0*plts)/self.sd['ppm']]
        if 'origin' not in plotkwargs:
            plotkwargs['origin'] = 'lower'
        if 'vmin' not in plotkwargs:
            plotkwargs['vmin'] = pmin
        if 'vmax' not in plotkwargs:
            plotkwargs['vmax'] = pmax


        ##
        ##  create (only) subplot
        ##
        logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargs))
        logger.debug("Plot keyword args:\n"+aosat_cfg.repString(plotkwargs))
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))#
        im=ax.imshow(util.ensure_numpy(sphase[sreg,sreg]), **plotkwargs)
        ax.text(0.5,0.5,"Mean RMS WF:\n %.1f nm"%self.rms,color='white',transform=ax.transAxes,size=5,ha='center',va='center')

        dividerpi = make_axes_locatable(ax)
        caxpi = dividerpi.append_axes("right", size="10%", pad=0.05)
        caxpi.tick_params(axis='both', which='major', labelsize=6)
        caxpi.tick_params(axis='both', which='minor', labelsize=5)
        t = util.ensure_numpy(pmin+(np.arange(11))/10.0*(pmax-pmin))
        cbarpi = plt.colorbar(util.ensure_numpy(im), cax=caxpi, ticks=t,label=r'Phase [nm]')

        return(fig)
