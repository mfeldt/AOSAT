
import os


from pip._internal.utils.misc import get_installed_distributions

if any(["cupy" in str(f) for f in get_installed_distributions()]):
    import cupy as np
else:
    import numpy as np
#import numpy as np
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

class zrn_analyzer():
    """Short summary.

    Parameters
    ----------
    sd : type
        Description of parameter `sd`.

    Attributes
    ----------
    ffed : type
        Description of attribute `ffed`.
    modes : 1D float NDarray of length zterms
        Time averaged mean amplitude of each mode
    dmodes : 1D float NDarray of length zterms
        Standard deviation of each mode amplitude
    modest : 2D float NDarray of shape (zterms, nframes)
        Individual modal expansion for each frame
    _ffed : type
        Description of attribute `_ffed`.
    self,frame : type
        Description of attribute `self,frame`.
    sd

    """

    def __init__(self,sd):
        self.sd      = sd
        self.ffed    = 0
        self.modes   = np.zeros(sd['cfg']['zterms'])
        self.dmodes  = np.zeros(sd['cfg']['zterms'])
        self.modest  = None
        self._ffed   = 0

    def feed_frame(self,frame,nframes):
        if self._ffed == 0:
            self.modest = np.zeros((self.sd['cfg']['zterms'],nframes))

        self.modest[:,self._ffed] = util.basis_expand(frame/2/np.pi*self.sd['cfg']['an_lambda']*1e9,self.sd['zernike_basis'],self.sd['tel_mirror']) #in nm

        self._ffed +=1

    def finalize(self):
        self.modes  = np.mean(self.modest,axis=1)
        self.dmodes = np.std(self.modest,axis=1)


    def make_plot(self,fig=None,index=111,plotkwargs={},subplotkwargs={}):

        from numpy import arange
        if fig is None:
            fig = plt.figure()

        ##
        ## default appearance
        ##
        ylim = util.ensure_numpy(max(np.maximum(np.abs(self.modes+self.dmodes),np.abs(self.modes-self.dmodes))[1:])*1.1).item()
        if not 'color' in plotkwargs:
            plotkwargs['color'] = 'blue'
        if 'xlim' not in subplotkwargs:
            subplotkwargs['xlim'] =(1,self.sd['cfg']['zterms']-1)
        if 'ylim' not in subplotkwargs:
            subplotkwargs['ylim'] = (-1*ylim,ylim)
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = 'Term #'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = r'Amplitude and Standard Deviation [nm]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = r'Zernike Expansion'

        if 'color' not in plotkwargs:
            plotkwargs['color']='blue'
        ##
        ##  create (only) subplot
        ##
        logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargs))
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))
        ax.fill_between(arange(len(self.modes)),util.ensure_numpy(self.modes-self.dmodes),util.ensure_numpy(self.modes+self.dmodes),alpha=0.5,**plotkwargs)
        ax.plot(arange(len(self.modes)),util.ensure_numpy(self.modes),**plotkwargs)
        ax.text(0.75,0.95,r'Amplitude',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        ax.text(0.75,0.9,r'Std. Dev.',transform=ax.transAxes,size=6,ha='left',alpha=0.5,color=plotkwargs['color'])

        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        return(fig)

    def make_report(self):
        from numpy import array2string
        report =  "##\n##\n"
        report += "## reporting analyzer: %s\n##\n##\n" % self.__class__.__name__
        report += "## Basis expansion values (mean) [nm]:\n\n"
        report += "modes_mean    = %s\n" % (array2string(util.ensure_numpy(self.modes),separator=',',precision=2))
        report += "## Basis expansion values (standard deviation) [nm]:\n\n"
        report += "modes_std     = %s\n" % (array2string(util.ensure_numpy(self.dmodes),separator=',',precision=2))
        return(report)
