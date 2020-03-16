

import os
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

class psf_analyzer():
    def __init__(self,sd):
        self.sd = sd
        self.ffed     = 0
        self.psf      = None
        self.psf_mean = None
        self.strehl   = 0.0
        self.sr_wf    = None
        self.ttx      = None
        self.tty      = None
        self.ttx_wf   = None
        self.ttx_wf   = None
        self.A        = None
        self.M        = None
        self.lmf      = fitting.LevMarLSQFitter()
        self.x        = None
        self.y        = None
        xx,yy = np.indices([32,32])
        self.sdim     = 0
        self.xx       = xx
        self.yy       = yy

    def feed_frame(self,frame,nframes):
        in_field  = self.sd['tel_mirror']*np.exp(1j*frame)
        frame_dev = fftx.FFTprepare(in_field)
        fftframe  = fftx.FFTshift(fftx.FFTforward(self.sd['fft_plan'],self.sd['fft_out'],frame_dev))
        psf       = np.abs(fftframe)**2
        if self.ffed == 0:
            logger.debug("Initializing variables...")
            self.psf       = psf*0
            self.sdim      = psf.shape[0]
            self.sr_wf     = np.zeros(nframes)
            self.ttx       = np.zeros(nframes)
            self.tty       = np.zeros(nframes)
            self.ttilt     = np.zeros(nframes)
            self.ttjit     = 0.0
            self.ttq90     = 0.0
            self.ttjit_psf = 0.0
            self.ttq90_psf = 0.0
            self.tty_wf    = np.zeros(nframes)
            self.x, self.y = np.mgrid[-self.sdim/2:self.sdim/2,-self.sdim/2:self.sdim/2]/self.sd['ppm']
            self.A         = np.c_[self.x[self.sd['wnz']], self.y[self.sd['wnz']], self.x[self.sd['wnz']]*0.0+1.0]
            self.M         = models.Gaussian2D(np.max(psf), 16, 16, 3, 3)
            logger.debug("Done!")

        ## Strehl from WF
        self.sr_wf[self.ffed] = np.exp(-1*np.std(frame[self.sd['wnz']])**2)
        logger.debug("Found SR: %s"%self.sr_wf[self.ffed])

        ## mean PSF (and variance)
        self.psf += psf/nframes

        logger.debug("Solving tilt...")
        ## tilt from WF
        C,_,_,_ = scipy.linalg.lstsq(self.A, frame[self.sd['wnz']])
        logger.debug("Constructing fitted WF")
        tphase = (C[0]*self.x + C[1]*self.y + C[2])/2/np.pi*self.sd['cfg']['an_lambda'] # frame in m
        self.ttilt[self.ffed] = (np.max(tphase[self.sd['wnz']]) - np.min(tphase[self.sd['wnz']]))/(self.sd['pupildiam']/self.sd['ppm'])/np.pi*180.0*3600.0 # arc sec
        logger.debug("Found tilt of %s arc sec", self.ttilt[self.ffed])
        logger.debug("Fitting position...")
        ## tip and tilt from PSF fitted position
        result_lmf = self.lmf(self.M, self.xx, self.yy, psf[int(self.sdim/2-16):int(self.sdim/2+16),int(self.sdim/2-16):int(self.sdim/2+16)])
        self.ttx[self.ffed] = (result_lmf.x_mean - 16.0)*self.sd['aspp']
        self.tty[self.ffed] = (result_lmf.y_mean - 16.0)*self.sd['aspp']
        logger.debug("Found tt excursion of %.3f, %.3f mas from PSF!"%(self.ttx[self.ffed]*1000,self.tty[self.ffed]*1000))
        self.ffed += 1


    def make_report(self):
        report =  "##\n##\n"
        report += "## reporting analyzer: %s\n##\n##\n" % self.__class__.__name__
        report += "## Strehl from average's PSF peak intensity:\n"
        report += "s_psf = %s\n##\n" % self.strehl
        report += "## Strehl from averaged wavefront Strehls:\n"
        report += "s_wf = %s\n##\n" % self.sr_wf.mean()
        report += "##\n##All tt values in mas\n##\n"
        report += "## TT jitter from average's PSF peak position:\n"
        report += "ttjit_psf = %s\n##\n" % self.ttjit_psf
        report += "## TT jitter from wavefront tilts:\n"
        report += "ttjit_wf = %s\n##\n" % self.ttjit
        report += "## TT Q90 excursion from average's PSF peak position:\n"
        report += "ttjit_psf = %s\n##\n" % self.ttq90_psf
        report += "## TT Q90 excursion from wavefront tilts:\n"
        report += "ttjit_wf = %s\n##\n" % self.ttq90
        return(report)

    def make_plot(self,fig=None,index=111,plotkwargs={},subplotkwargs={}):

        if fig is None:
            fig = plt.figure()

        ##
        ## default appearance
        ##
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = r'$\delta$RA[mas]'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = r'$\delta$DEC[mas]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = 'PSF'

        # size of extraction area
        plts = int(min([np.round(self.sd['crad']/self.sd['aspp']/2.0)*2,self.sdim/2]))
        sd2 = int(self.sdim/2)

        if 'cmap' not in plotkwargs:
            plotkwargs['cmap'] = 'nipy_spectral'
        if 'norm' not in plotkwargs:
            plotkwargs['norm'] = LogNorm(vmin=1e-7,vmax=0.5)
        if 'extent' not in plotkwargs:
            plotkwargs['extent'] = [plts*self.sd['aspp']*1000,-plts*self.sd['aspp']*1000,-plts*self.sd['aspp']*1000,plts*self.sd['aspp']*1000]
        if 'origin' not in plotkwargs:
            plotkwargs['origin'] = 'lower'
        ##
        ##  create (only) subplot
        ##
        logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargs))
        logger.debug("Plot keyword args:\n"+aosat_cfg.repString(plotkwargs))
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))

        #   pdb.set_trace()
        im = ax.imshow(self.psf[sd2-plts:sd2+plts,sd2-plts:sd2+plts]/np.max(self.psf), **plotkwargs)

        ## useful info
        ax.text(0.7,0.95,'SR (PSF) = %.3f' % self.strehl,transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.7,0.9,'SR (WF) = %.3f' % self.sr_wf.mean(),transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.05,0.1,r'$\sigma_{tt}$(PSF) = %.3f mas' % self.ttjit_psf,transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.05,0.05,r'$\sigma_{tt}$(WF) = %.3f mas' % self.ttjit,transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.65,0.1,r'$Q90_{tt}$(PSF) = %.3f mas' % self.ttq90_psf,transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.65,0.05,r'$Q90_{tt}$(WF) = %.3f mas' % self.ttq90,transform=ax.transAxes,size=6,ha='left',color='white')

        ## color bar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="10%", pad=0.05)
        cax.tick_params(axis='both', which='major', labelsize=6)
        cax.tick_params(axis='both', which='minor', labelsize=5)
        t = [1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]
        cbar = plt.colorbar(im, cax=cax, ticks=t,label='Intensity')

        return(fig)


    def finalize(self):
        frame_dev = fftx.FFTprepare(self.sd['tel_mirror']*np.exp(1j*0.0))
        fftframe  = fftx.FFTshift(fftx.FFTforward(self.sd['fft_plan'],self.sd['fft_out'],frame_dev))
        psf_ref   = np.abs(fftframe)**2
        logger.debug("Max of reference PSF: %s" % np.max(psf_ref))
        logger.debug("Max of measured PSF:  %s" % np.max(self.psf))
        self.strehl   = np.max(self.psf)/np.max(psf_ref)

        ttmas = self.ttilt*1000.0

        self.ttjit = ((ttmas**2).sum()/len(ttmas))**0.5 # tip-tilt jitter in mas
        self.ttq90 = np.sort(ttmas)[int(len(ttmas)*0.9)]

        ttmas_psf      = ((self.ttx**2+self.tty**2)**0.5)*1000.0
        self.ttjit_psf = ((ttmas_psf**2).sum()/len(ttmas_psf))**0.5
        self.ttq90_psf = np.sort(ttmas_psf)[int(len(ttmas_psf)*0.9)]
