

import os
import numpy as np
import scipy
from astropy import units
from astropy.io import fits as pyfits
import pandas as pd
import copy
import pdb

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
from statsmodels.tsa.stattools import acf


import logging
logger = logging.getLogger(__name__)


class tvc_analyzer():
    """Short summary.

    Parameters
    ----------
    sd : setup dictionary
        must be fully set up!

    ctype : string
        Can be either 'raw' or 'icor'. When 'raw', the analysis is
        executed on the raw, non-coronagraphic PSF. When 'icor', an ideal coronagraph
        is simulated according to  (the default is 'raw').

    Examples
    -------

    >>> sd={}
    >>> a1 = tvc_analyzer(sd)

    Attributes
    ----------
    available after 'finalize'

    variance : 2D array
        variance of the PSF
    contrast : 2D array
        5 sigma contrast of the PSF
    rvec : array
        Vector of radii where contrast values are given
    cvecmean : array
        Mean contrast at locatstepoions in rvec. Same length as rvec.
    cvecmin : type
        Minimum contrast at locations in rvec. Same length as rvec.
    cvecmax : type
        Maximum contrast at locations in rvec. Same length as rvec.
    ctype : string
        Flag signalling which analysis is to be performed.
    sd    : set-up dictionary

    """
    def __init__(self,sd,ctype='nocor'):
        self.ctype     = ctype
        self.sd        = sd
        self.variance  = (0,None,None)
        self.variance2 = (0,None,None)
        self.contrast  = None
        self.rcontrast = None
        self.mean      = None
        self.rvec      = None
        self.cvecmean  = None
        self.cvecmin   = None
        self.cvecmax   = None
        self._ffed     = 0
        self.ntracks   = self.sd['cfg']['ntracks']
        self.ptrak     = None
        self.ppos      = np.zeros(self.ntracks,dtype=np.int32)
        self.corrlen   = 0.0

    def feed_frame(self,frame,nframes):
        in_field  = self.sd['tel_mirror']*np.exp(1j*frame)
        this_strehl   = np.exp(-1*np.std(frame[self.sd['wnz']])**2) ## remember: Frames are in radians!
        logger.debug("Found Strehl ratio: %s" % this_strehl)
        frame_dev = fftx.FFTprepare(in_field- (this_strehl**0.5)*self.sd['tel_mirror']*int(self.ctype=='icor'))
        fftframe  = fftx.FFTshift(fftx.FFTforward(self.sd['fft_plan'],self.sd['fft_out'],frame_dev))
        psf       = np.abs(fftframe)**2
        if self.ctype =='icor':
            frame_dev2 = fftx.FFTprepare(in_field)
            fftframe2  = fftx.FFTshift(fftx.FFTforward(self.sd['fft_plan'],self.sd['fft_out'],frame_dev2))
            psf2       = np.abs(fftframe2)**2
        if self.variance == (0,None,None):
            self.variance=(0,psf*0,psf*0)
            self.variance2=(0,psf*0,psf*0)
        self.variance = util.rolling_variance(self.variance,psf)
        if self.ctype=='icor':
            self.variance2 = util.rolling_variance(self.variance2,psf2)
        if self._ffed == 0:
            posvec = (np.random.random(self.ntracks)*len(self.sd['wnz'][0])).astype(int)
            self.ppos = (self.sd['wnz'][0][posvec],self.sd['wnz'][1][posvec])
            self.ptrak=np.zeros((nframes,self.ntracks))
        self.ptrak[self._ffed] = frame[self.ppos] # track phase at 6 random locations
        self._ffed +=1
    def make_plot(self,fig=None,index=111,plotkwargs={},subplotkwargs={}):

        if fig is None:
            fig = plt.figure()

        ##
        ## default appearance
        ##
        if not 'color' in plotkwargs:
            plotkwargs['color'] = 'red'
        if 'ylim' not in subplotkwargs:
            subplotkwargs['ylim'] = (1e-8,1)
        if 'xlim' not in subplotkwargs:
            subplotkwargs['xlim'] =(0,self.sd['crad']*1.2)
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = 'Sep. [arc sec.]'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = 'Contrast [Peak]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = r'Simple Single-Pixel %s Contrast' % self.ctype
        if 'yscale' not in subplotkwargs:
            subplotkwargs['yscale'] = 'log'
        #if 'aspect' not in subplotkwargs:
        #    subplotkwargs['aspect'] = 1.0
        #subplotkwargs['adjustable'] = 'box'
        ##
        ##  create (only) subplot
        ##
        logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargs))
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))
        ax.fill_between(self.rvec,self.cvecmin,self.cvecmax,alpha=0.15,**plotkwargs)
        ax.plot(self.rvec,self.cvecmean,**plotkwargs)
        if self.corrlen > 0:
            num_frames_per_hour = 3600.0/(self.corrlen/self.sd['loopfreq'])
            plotkwargsP = copy.copy(plotkwargs)
            plotkwargsP['linestyle'] = 'dotted'
            ax.plot(self.rvec,self.cvecmean/(num_frames_per_hour**0.5),**plotkwargsP)

        ax.text(0.5,0.9,'No photon noise,',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        ax.text(0.5,0.85,'due to PSF variation only!',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        ax.text(0.5,0.95,r'5 $\sigma$ TVC',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        plotkwargs2 = copy.copy(plotkwargs)
        plotkwargs2['color']='black'
        ax.fill_between(self.rvec,self.rcvecmin,self.rcvecmax,alpha=0.15,**plotkwargs2)
        ax.plot(self.rvec,self.rcvecmean,**plotkwargs2)
        if self.corrlen > 0:
            plotkwargsP = copy.copy(plotkwargs2)
            plotkwargsP['linestyle'] = 'dotted'
            ax.plot(self.rvec,self.rcvecmean/(num_frames_per_hour**0.5),**plotkwargsP)
        ax.text(0.5,0.8,'Raw (PSF profile)',transform=ax.transAxes,size=6,ha='left',color=plotkwargs2['color'])
        if self.corrlen>0:
            ax.text(0.5,0.75,'Dotted lines: est. best 1hr ADI contrast',transform=ax.transAxes,size=6,ha='left',color=plotkwargs2['color'])
        #ax.set_aspect(subplotkwargs['aspect'],'box')

        return(fig)

    def make_report(self):

        rep_vec_r     = np.arange(50)*self.sd['aspp']

        rep_vec_cmean = np.interp(rep_vec_r,self.rvec,self.cvecmean)
        rep_vec_cmin  = np.interp(rep_vec_r,self.rvec,self.cvecmin)
        rep_vec_cmax  = np.interp(rep_vec_r,self.rvec,self.cvecmax)

        rep_vec_rcmean = np.interp(rep_vec_r,self.rvec,self.rcvecmean)
        rep_vec_rcmin  = np.interp(rep_vec_r,self.rvec,self.rcvecmin)
        rep_vec_rcmax  = np.interp(rep_vec_r,self.rvec,self.rcvecmax)

        report =  "##\n##\n"
        report += "## reporting analyzer: %s\n##\n##\n" % self.__class__.__name__
        report += "## contrast type: %s\n" % self.ctype
        report += "## Interpolated contrast values:\n\n"
        report += "tvc_%s_sepvec    = %s\n" % (self.ctype,np.array2string(rep_vec_r,separator=',',precision=4))
        report += "tvc_%s_ctrstmean = %s\n" % (self.ctype,np.array2string(rep_vec_cmean,separator=','))
        report += "tvc_%s_ctrstmin  = %s\n" % (self.ctype,np.array2string(rep_vec_cmin,separator=','))
        report += "tvc_%s_ctrstmax  = %s\n\n\n" % (self.ctype,np.array2string(rep_vec_cmax,separator=','))
        report += "raw_%s_sepvec    = %s\n" % (self.ctype,np.array2string(rep_vec_r,separator=',',precision=4))
        report += "raw_%s_ctrstmean = %s\n" % (self.ctype,np.array2string(rep_vec_rcmean,separator=','))
        report += "raw_%s_ctrstmin  = %s\n" % (self.ctype,np.array2string(rep_vec_rcmin,separator=','))
        report += "raw_%s_ctrstmax  = %s\n##\n" % (self.ctype,np.array2string(rep_vec_rcmax,separator=','))
        report += "## estimated correlation time [ms] \n##\n"
        report += "estim_corr_time   = %s \n\n\n" % (1000*self.corrlen/self.sd['loopfreq'])
        return(report)


    def finalize(self):

        variance      = self.variance[2]/self.variance[0]
        mean          = self.variance[1]#*int(self.ctype == 'nocor') + self.variance2[1]*int(self.ctype=='icor')
        sigma         = variance**0.5
        max_no_cor    = np.max(self.variance[1]*int(self.ctype == 'nocor') + self.variance2[1]*int(self.ctype=='icor'))
        self.mean     = mean

        self.contrast = 5*sigma/max_no_cor
        self.rcontrast = mean/max_no_cor
        ra,pixels,rord = util.rad_order(self.contrast)

        self.rvec  = ra[pixels][rord]*self.sd['aspp']
        c          = pd.Series(self.contrast[pixels][rord])
        r          = pd.Series(self.rcontrast[pixels][rord])

        self.cvecmean = c.rolling(window=50).mean()
        self.cvecmin  = c.rolling(window=50).min()
        self.cvecmax  = c.rolling(window=50).max()
        self.rcvecmean = r.rolling(window=50).mean()
        self.rcvecmin  = r.rolling(window=50).min()
        self.rcvecmax  = r.rolling(window=50).max()

        max_corr_len=0.0
        for i in range(self.ntracks):
            acff,conf    = acf(self.ptrak[:,i],fft=True,alpha=0.05,nlags=min([500,len(self.ptrak[:,i])-1]))
            where_uncorr = np.where(acff < (conf[:,1] - acff))[0]
            if len(where_uncorr) > 0:
                max_corr_len = np.max([max_corr_len,np.min(where_uncorr)])
        self.corrlen = max_corr_len*1.0
