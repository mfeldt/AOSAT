

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


class frg_analyzer():
    """Analyze the island effect, i.e. the impact of pupil
    fragmentation.

    Parameters
    ----------
    sd : setup dictionary

    Examples
    -------

    >>> sd={}
    >>> ### never do the below, it's just to get the example work!
    >>> sd['fragmask'] = np.arange(6)
    >>> a1 = frg_analyzer(sd)

    Attributes
    ----------
    available after calling finalize()

    ffed : int
        Number of frames fed into analyzer.
    piston : numpy array (float)
        Array holding the mean piston value for each pupil
        fragment (in nm).
    dpiston : numpy array (float)
        Array holding the standard deviation of the piston
        value for each pupil fragment (in nm).
    ttx : numpy array (float)
        Array holding mean tip deviation for each fragment
        (in mas)
    dttx : numpy array (float)
        Array holding the standard deviation of tip
        deviation for each pupil fragment (in mas)
    tty : numpy array (float)
        Array holding mean tilt deviation for each fragment
        (in mas)
    dtty : numpy array (float)
        Array holding the standard deviation of tilt
        deviation for each pupil fragment (in mas)
    pistframe : numpy array
        Frame holding the worst piston pattern across the
        top (1-tile) (see below) part of the simulated
        sequence
    tile : float
        Fractional tile above which the anylzer looks
        for the worst piston frame in the sequence.


    """
    def __init__(self,sd):
        num_fragments = np.max(sd['fragmask']).item()
        self.sd      = sd
        self.ffed    = 0
        self.piston  = np.zeros(num_fragments)
        self.dpiston = np.zeros(num_fragments)
        self.pistont = None
        self.ttx     = np.zeros(num_fragments)
        self.dttx    = np.zeros(num_fragments)
        self.ttxt    = None
        self.tty     = np.zeros(num_fragments)
        self.dtty    = np.zeros(num_fragments)
        self.ttyt    = None
        self.wfrag   = []
        self.x       = None
        self.y       = None
        self.A       = []
        self.pistframe = None
        self.tile      = 0.8
        self.worstid   = 0

        for i in range(num_fragments):
            self.wfrag.append(np.where(sd['fragmask'] == (i+1)))


    def feed_frame(self,frame,nframes):
        num_fragments = np.max(self.sd['fragmask']).item()

        if self.ffed == 0:
            logger.debug("Initializing variables...")
            sdim = frame.shape[0]
            self.pistont = np.zeros((nframes,num_fragments))
            self.ttxt = np.zeros((nframes,num_fragments))
            self.ttyt = np.zeros((nframes,num_fragments))
            self.x, self.y = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]/self.sd['ppm']
            for i in range(num_fragments):
                self.A.append(np.c_[self.x[self.wfrag[i]], self.y[self.wfrag[i]], self.x[self.wfrag[i]]*0.0+1.0])

            logger.debug("Done!")

        for i in range(num_fragments):
            ## tilt from WF
            C,_,_,_ = np.linalg.lstsq(self.A[i], frame[self.wfrag[i]])
            self.ttxt[self.ffed,i]    = C[0]
            self.ttyt[self.ffed,i]    = C[1]
            self.pistont[self.ffed,i] = C[2]
        self.ffed+=1

    def finalize(self,tile=0.8):
        for i in range(np.max(self.sd['fragmask']).item()):
            self.piston[i]  = self.pistont[:,i].mean()/2/np.pi*self.sd['cfg']['an_lambda']*1e9 # nano metres
            self.dpiston[i] = self.pistont[:,i].std()/2/np.pi*self.sd['cfg']['an_lambda']*1e9 # nano metres
            self.ttx[i]     = self.ttxt[i].mean()/2/np.pi*self.sd['cfg']['an_lambda']/np.pi*180*3600*1000 # milli arc sec
            self.dttx[i]    = self.ttxt[i].std()/2/np.pi*self.sd['cfg']['an_lambda']/np.pi*180*3600*1000 # milli arc sec
            self.tty[i]     = self.ttyt[i].mean()/2/np.pi*self.sd['cfg']['an_lambda']/np.pi*180*3600*1000 # milli arc sec
            self.dtty[i]    = self.ttyt[i].std()/2/np.pi*self.sd['cfg']['an_lambda']/np.pi*180*3600*1000 # milli arc sec

        self.tile = tile
        npist = len(self.pistont[:,0])
        lpa = int(npist*(1-tile))

        pistmin = np.amin(self.pistont,axis=1)
        pistmax = np.amax(self.pistont,axis=1)
        wmaxpist = np.argmax(pistmax[-lpa:]-pistmin[-lpa:])

        pistframe=self.sd['tel_mirror']*0
        for i in range(np.max(self.sd['fragmask']).item()):
            xc = np.mean(self.x[np.where(self.sd['fragmask']==i)])
            yc = np.mean(self.y[np.where(self.sd['fragmask']==i)])
            tphase = (self.ttxt[-lpa:,i][wmaxpist]*(self.x-xc) + self.ttyt[-lpa:,i][wmaxpist]*(self.y-yc) + self.pistont[-lpa:,i][wmaxpist])/2/np.pi*self.sd['cfg']['an_lambda']*1e9 # frame in nm
            pistframe[self.wfrag[i]] = tphase[self.wfrag[i]]
        self.pistframe = pistframe

    def make_plot(self,fig=None,index=111,plotkwargs={},subplotkwargs={}):
        if fig is None:
            fig = plt.figure()


        ## plot region
        sdim  = self.pistframe.shape[0]
        x, y  = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]
        r     = (x**2+y**2)**0.5
        plts  = int(min([max((util.ensure_numpy(r*self.sd['tel_mirror'])).flatten())*1.1,sdim/2])/2)*2
        sreg = slice(int(sdim/2-plts),int(sdim/2+plts))

        ## what to show
        mphase = np.mean(self.pistframe[self.sd['wnz']])
        sphase = copy.copy(self.pistframe)
        sphase[self.sd['wnz']] -=  mphase
        pp     = sphase[self.sd['wnz']].flatten()
        pmin  = np.percentile(pp,5).item()
        pmax  = np.percentile(pp,95).item()
        sshow = sphase[sreg,sreg]
        sshow[np.where(self.sd['tel_mirror'][sreg,sreg] == 0)] = pmin
        sshow = util.ensure_numpy(sshow)

        ##
        ## default appearance
        ##
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = r'$\delta$x [m]'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = r'$\delta$y [m]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = 'Strongest Residual Fragment Piston\n '


        if 'cmap' not in plotkwargs:
            plotkwargs['cmap'] = 'nipy_spectral'
        if 'extent' not in plotkwargs:
            plotkwargs['extent']=[(-1.0*plts)/self.sd['ppm'],(1.0*plts)/self.sd['ppm'],-(1.0*plts)/self.sd['ppm'],(1.0*plts)/self.sd['ppm']]
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
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))#

        im=ax.imshow(sshow, **plotkwargs)
        x, y =  (np.mgrid[0:sdim,0:sdim]-sreg.start)/(sreg.stop-sreg.start)#np.mgrid[0:sdim,0:sdim]/sdim
        for i in range(np.max(self.sd['fragmask']).item()):
            #pdb.set_trace()
            xc = np.mean(x[self.wfrag[i]])
            yc = np.mean(y[self.wfrag[i]])
            ax.text(yc,xc,'%.2f' % self.dpiston[i],color='white',transform=ax.transAxes,size=5,ha='center')

        ax.text(0.5,0.5,r'$\sigma_{segpist}$[nm]' '\n' '(full seq.)',color='white',transform=ax.transAxes,size=5,ha='center',va='center')
        ax.text(0.5,1.05,r'(in last %d percent of sequence)' % int((1-self.tile)*100),color='black',transform=ax.transAxes,size=9,ha='center',va='center')

        dividerpi = make_axes_locatable(ax)
        caxpi = dividerpi.append_axes("right", size="10%", pad=0.05)
        caxpi.tick_params(axis='both', which='major', labelsize=6)
        caxpi.tick_params(axis='both', which='minor', labelsize=5)
        t = util.ensure_numpy(pmin+(np.arange(11))/10.0*(pmax-pmin))
        cbarpi = plt.colorbar(im, cax=caxpi, ticks=t,label=r'Phase [nm]')

        return(fig)



    def make_report(self):

        report =  "##\n##\n"
        report += "## reporting analyzer: %s\n##\n##\n" % self.__class__.__name__
        report += "## Fragment IDs:\n"
        report += "frag_id = %s\n##\n" % np.arange(np.max(self.sd['fragmask']))
        report += "## Average piston of fragments [nm]:\n"
        report += "pist_mean = %s\n##\n" % self.piston
        report += "## Standard deviation of piston of fragments [nm]:\n"
        report += "pist_std = %s\n##\n" % self.dpiston
        report += "## Average tip of fragments [mas]:\n"
        report += "ttx_mean = %s\n##\n" % self.ttx
        report += "## Standard deviation of tip of fragments [mas]:\n"
        report += "ttx_std = %s\n##\n" % self.dttx
        report += "## Average tilt of fragments [mas]:\n"
        report += "tty_mean = %s\n##\n" % self.tty
        report += "## Standard deviation of tilt of fragments [mas]:\n"
        report += "tty_std = %s\n##\n\n\n" % self.dtty
        return(report)
