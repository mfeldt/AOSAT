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



##
##
## Simple PSF Analyzer
##
##


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
        Mean contrast at locations in rvec. Same length as rvec.
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
        if 'aspect' not in subplotkwargs:
            subplotkwargs['aspect'] = 1.0
        subplotkwargs['adjustable'] = 'box'

        ##
        ##  create (only) subplot
        ##
        logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargs))
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))
        ax.fill_between(self.rvec,self.cvecmin,self.cvecmax,alpha=0.15,**plotkwargs)
        ax.plot(self.rvec,self.cvecmean,**plotkwargs)
        ax.text(0.5,0.9,'No photon noise,',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        ax.text(0.5,0.85,'due to PSF variation only!',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        ax.text(0.5,0.95,r'5 $\sigma$ TVC',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        plotkwargs2 = copy.copy(plotkwargs)
        plotkwargs2['color']='black'
        ax.fill_between(self.rvec,self.rcvecmin,self.rcvecmax,alpha=0.15,**plotkwargs2)
        ax.plot(self.rvec,self.rcvecmean,**plotkwargs2)
        ax.text(0.5,0.8,'Raw (PSF profile)',transform=ax.transAxes,size=6,ha='left',color=plotkwargs2['color'])
        ax.set_aspect(subplotkwargs['aspect'],'box')
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
        report += "raw_%s_ctrstmax  = %s\n\n\n" % (self.ctype,np.array2string(rep_vec_rcmax,separator=','))

        return(report)

    def finalize(self):

        variance      = self.variance[2]/self.variance[0]
        mean          = self.variance[1]*int(self.ctype == 'nocor') + self.variance2[1]*int(self.ctype=='icor')
        sigma         = variance**0.5
        self.mean     = mean
        self.contrast = 5*sigma/np.max(self.variance[1]*int(self.ctype == 'nocor') + self.variance2[1]*int(self.ctype=='icor'))
        self.rcontrast = mean/np.max(self.variance[1]*int(self.ctype == 'nocor') + self.variance2[1]*int(self.ctype=='icor'))
        r,pixels,rord = util.rad_order(self.contrast)

        self.rvec  = r[pixels][rord]*self.sd['aspp']
        c          = pd.Series(self.contrast[pixels][rord])
        r          = pd.Series(self.rcontrast[pixels][rord])

        self.cvecmean = c.rolling(window=50).mean()
        self.cvecmin  = c.rolling(window=50).min()
        self.cvecmax  = c.rolling(window=50).max()
        self.rcvecmean = r.rolling(window=50).mean()
        self.rcvecmin  = r.rolling(window=50).min()
        self.rcvecmax  = r.rolling(window=50).max()

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
        plts  = int(np.min([max((r*self.sd['tel_mirror']).flatten())*1.1,sdim/2])/2)*2
        sreg = slice(int(sdim/2-plts),int(sdim/2+plts))

        ## what to show
        mphase = np.min(self.lastphase[self.sd['wnz']])
        sphase = copy.copy(self.lastphase)
        sphase[self.sd['wnz']] -=  mphase
        pp     = sphase[self.sd['wnz']].flatten()
        pmin  = np.percentile(pp,5)
        pmax  = np.percentile(pp,95)
        sshow = sphase[sreg,sreg]


        ##
        ## default appearance
        ##
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = r'$\delta$x[m]'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = r'$\delta$y[m]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = 'Last Residual Phase\n '


        if 'cmap' not in plotkwargs:
            plotkwargs['cmap'] = 'nipy_spectral'
        if 'extent' not in plotkwargs:
            plotkwargs['extent']=[(1.0*plts)/self.sd['ppm'],-(1.0*plts)/self.sd['ppm'],-(1.0*plts)/self.sd['ppm'],(1.0*plts)/self.sd['ppm']]
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
        im=ax.imshow(sshow*self.sd['tel_mirror'][sreg,sreg], **plotkwargs)
        ax.text(0.5,0.5,"Mean RMS WF:\n %.1f nm"%self.rms,color='white',transform=ax.transAxes,size=5,ha='center',va='center')

        dividerpi = make_axes_locatable(ax)
        caxpi = dividerpi.append_axes("right", size="10%", pad=0.05)
        caxpi.tick_params(axis='both', which='major', labelsize=6)
        caxpi.tick_params(axis='both', which='minor', labelsize=5)
        t = pmin+(np.arange(11))/10.0*(pmax-pmin)
        cbarpi = plt.colorbar(im, cax=caxpi, ticks=t,label=r'Phase [nm]')

        return(fig)


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
        self.sd      = sd
        self.ffed    = 0
        self.piston  = np.zeros(np.max(sd['fragmask']))
        self.dpiston = np.zeros(np.max(sd['fragmask']))
        self.pistont = None
        self.ttx     = np.zeros(np.max(sd['fragmask']))
        self.dttx    = np.zeros(np.max(sd['fragmask']))
        self.ttxt    = None
        self.tty     = np.zeros(np.max(sd['fragmask']))
        self.dtty    = np.zeros(np.max(sd['fragmask']))
        self.ttyt    = None
        self.wfrag   = []
        self.x       = None
        self.y       = None
        self.A       = []
        self.pistframe = None
        self.tile      = 0.8
        self.worstid   = 0

        for i in range(np.max(sd['fragmask'])):
            self.wfrag.append(np.where(sd['fragmask'] == (i+1)))


    def feed_frame(self,frame,nframes):
        if self.ffed == 0:
            logger.debug("Initializing variables...")
            sdim = frame.shape[0]
            self.pistont = np.zeros((nframes,np.max(self.sd['fragmask'])))
            self.ttxt = np.zeros((nframes,np.max(self.sd['fragmask'])))
            self.ttyt = np.zeros((nframes,np.max(self.sd['fragmask'])))
            self.x, self.y = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]/self.sd['ppm']
            for i in range(np.max(self.sd['fragmask'])):
                self.A.append(np.c_[self.x[self.wfrag[i]], self.y[self.wfrag[i]], self.x[self.wfrag[i]]*0.0+1.0])

            logger.debug("Done!")

        for i in range(np.max(self.sd['fragmask'])):
            ## tilt from WF
            C,_,_,_ = scipy.linalg.lstsq(self.A[i], frame[self.wfrag[i]])
            self.ttxt[self.ffed,i]    = C[0]
            self.ttyt[self.ffed,i]    = C[1]
            self.pistont[self.ffed,i] = C[2]
        self.ffed+=1

    def finalize(self,tile=0.8):
        for i in range(np.max(self.sd['fragmask'])):
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
        for i in range(np.max(self.sd['fragmask'])):
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
        plts  = int(np.min([max((r*self.sd['tel_mirror']).flatten())*1.1,sdim/2])/2)*2
        sreg = slice(int(sdim/2-plts),int(sdim/2+plts))

        ## what to show
        mphase = np.mean(self.pistframe[self.sd['wnz']])
        sphase = copy.copy(self.pistframe)
        sphase[self.sd['wnz']] -=  mphase
        pp     = sphase[self.sd['wnz']].flatten()
        pmin  = np.percentile(pp,5)
        pmax  = np.percentile(pp,95)
        sshow = sphase[sreg,sreg]
        sshow[np.where(self.sd['tel_mirror'][sreg,sreg] == 0)] = pmin


        ##
        ## default appearance
        ##
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = r'$\delta$x[m]'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = r'$\delta$y[m]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = 'Strongest Residual Fragment Piston\n '


        if 'cmap' not in plotkwargs:
            plotkwargs['cmap'] = 'nipy_spectral'
        if 'extent' not in plotkwargs:
            plotkwargs['extent']=[(1.0*plts)/self.sd['ppm'],-(1.0*plts)/self.sd['ppm'],-(1.0*plts)/self.sd['ppm'],(1.0*plts)/self.sd['ppm']]
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
        for i in range(np.max(self.sd['fragmask'])):
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
        t = pmin+(np.arange(11))/10.0*(pmax-pmin)
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






class zrn_analyzer():
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

        if fig is None:
            fig = plt.figure()

        ##
        ## default appearance
        ##
        if not 'color' in plotkwargs:
            plotkwargs['color'] = 'blue'
        if 'xlim' not in subplotkwargs:
            subplotkwargs['xlim'] =(1,self.sd['cfg']['zterms'])
        if 'ylim' not in subplotkwargs:
            subplotkwargs['ylim'] =(0,np.max(np.maximum(self.modes**2,self.dmodes**2)[1:])*1.1)
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = 'Term #'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = r'Amplitude$^2$ and Variance [nm$^2$]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = r'Zernike Expansion'

        if 'color' not in plotkwargs:
            plotkwargs['color']='blue'
        ##
        ##  create (only) subplot
        ##
        logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargs))
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))
        ax.fill_between(np.arange(len(self.modes)),np.repeat(0.0,len(self.modes)),self.dmodes**2,alpha=0.5,**plotkwargs)
        ax.plot(np.arange(len(self.modes)),self.modes**2,**plotkwargs)
        ax.text(0.75,0.95,r'Amplitude$^2$',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        ax.text(0.75,0.9,r'Variance',transform=ax.transAxes,size=6,ha='left',alpha=0.5,color=plotkwargs['color'])

        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        return(fig)

    def make_report(self):

        report =  "##\n##\n"
        report += "## reporting analyzer: %s\n##\n##\n" % self.__class__.__name__
        report += "## Basis expansion values (mean) [nm]:\n\n"
        report += "modes_mean    = %s\n" % (np.array2string(self.modes,separator=',',precision=2))
        report += "## Basis expansion values (standard deviation) [nm]:\n\n"
        report += "modes_std     = %s\n" % (np.array2string(self.dmodes,separator=',',precision=2))
        return(report)






def getFontProperties():
    """Helper function.

    Parameters
    ----------
    none

    Returns
    -------
    font
        font for matplotlib

    """
    font0 = FontProperties()
    font = font0.copy()
    font.set_family('monospace')
    font.set_size(5)
    return(font)




def setup():
    """Set up an analysis run

    Parameters
    ----------
    none.

    Relevant settings are taken from the config settings.
    In particular the following parameters are sued:


    CFG_SETTINGS['pupilmask']
        Path to FITS file containing the telescope aperture,
        relative to the path to the config file itself

    CFG_SETTINGS['ppm']
        float giving the pixels per metre in the aperture plane

    CFG_SETTINGS['an_lambda']
        analysis wavelength in metres

    CFG_SETTINGS['zterms']
        int giving the number of zernike terms to analyze


    Returns
    -------
    setup dictionary
        Dictionary describing the analysis to be performed
    Examples
    -------

    >>> sd=setup()

    """
    logger.info("Setting up analysis...")
    setup_dict = {}

    setup_dict['tel_mirror'] = pyfits.getdata(os.path.join(aosat_cfg.CFG_SETTINGS['setup_path'],aosat_cfg.CFG_SETTINGS['pupilmask']))
    ext_x  = np.max(np.where(setup_dict['tel_mirror'].sum(axis=1) != 0))-np.min(np.where(setup_dict['tel_mirror'].sum(axis=1) != 0))+1
    ext_y  = np.max(np.where(setup_dict['tel_mirror'].sum(axis=0) != 0))-np.min(np.where(setup_dict['tel_mirror'].sum(axis=0) != 0))+1
    setup_dict['pupildiam']  = max([ext_x,ext_y])

    logger.debug("Found pupil diameter of %s pixels" % setup_dict['pupildiam'])

    setup_dict['sdim']          = setup_dict['tel_mirror'].shape[0] # always needs to be squared...
    setup_dict['divpm']         = setup_dict['tel_mirror'] * 1.0
    setup_dict['divpm'][np.where(setup_dict['tel_mirror'] == 0.0)] = 1.0 # pupil mask for division if necessary
    setup_dict['wnz']           = np.where(setup_dict['tel_mirror'] != 0.0)
    setup_dict['fragmask']      = ndimage.label(setup_dict['tel_mirror']>1e-6)[0]

    logger.debug("Found %s independent pupil fragments." % np.max(setup_dict['fragmask']))
    logger.info("Setting up FFTs...")
    setup_dict['fft_plan']      = fftx.FFTmakeplan(setup_dict['tel_mirror'])
    setup_dict['fft_out']       = fftx.FFTmakeout(setup_dict['tel_mirror'])
    setup_dict['cfg']           = aosat_cfg.CFG_SETTINGS
    setup_dict['ppm']           = aosat_cfg.CFG_SETTINGS['ppm'] # pix per metre
    setup_dict['aspp']          = aosat_cfg.CFG_SETTINGS['ppm']*3600/(np.pi/180.0)/setup_dict['sdim'] *aosat_cfg.CFG_SETTINGS['an_lambda'] # arc sec per pix
    setup_dict['crad']          = aosat_cfg.CFG_SETTINGS['an_lambda']/(5.0/setup_dict['ppm'])/np.pi*180*3600 # estimate of control radius assuming
                                                                                                             # a sampling of 5 pixels per actuator, in arc sec
    setup_dict['dl']            = aosat_cfg.CFG_SETTINGS['an_lambda']/(setup_dict['pupildiam']/setup_dict['ppm'])*3600*180/np.pi # diffraction limit


    x, y     = np.mgrid[-setup_dict['sdim']/2:setup_dict['sdim']/2,-setup_dict['sdim']/2:setup_dict['sdim']/2]/setup_dict['ppm']
    setup_dict['wtrk']          = np.where((setup_dict['tel_mirror'] != 0 ) * ((x-x.astype(int)) ==0) * ((y-y.astype(int))==0))

    logger.info("Making Modal basis with %s terms ..." % aosat_cfg.CFG_SETTINGS['zterms'])
    #setup_dict['zernike_basis'] = zernike.zernike_basis(nterms=aosat_cfg.CFG_SETTINGS['zterms'],npix=int(setup_dict['pupildiam']))
    mreg = slice (int(setup_dict['sdim']/2-setup_dict['pupildiam']/2),int(setup_dict['sdim']/2+setup_dict['pupildiam']/2))
    setup_dict['zernike_basis'] = zernike.arbitrary_basis(setup_dict['tel_mirror'][mreg,mreg],nterms=aosat_cfg.CFG_SETTINGS['zterms'],outside=0.0)


    return(setup_dict)

def sizeTearsheetLabels(f):
    ca = f.gca()
    ca.tick_params(axis='both', which='major', labelsize=6)
    ca.tick_params(axis='both', which='minor', labelsize=5)
    ca.yaxis.label.set_size(8)
    ca.xaxis.label.set_size(8)

def makeTearsheetFigure(analyzers):
    """Helper function generating the tearsheet figure

    Parameters
    ----------
    analyzers : type
        A list of analyzers to be executed and then
        called one by one to produce a plot on a "tearsheet".


    Returns
    -------
    none

    """
    num_plots  = len(analyzers)
    num_pages  = int(num_plots/6)+1
    logger.debug("Preparing %s plots on %s pages!" % (num_plots,num_pages))
    plots_done = 1

    filename = aosat_cfg.CFG_SETTINGS["ts_basefilename"]+".pdf"
    with PdfPages(filename) as pdf:

        f =plt.figure(figsize=(8.27,11.69))

        font = getFontProperties()
        logo = mpimg.imread(os.path.join(os.path.dirname(os.path.abspath(__file__)),'img','aosat_logo.png'))

        alignment = {'horizontalalignment': 'left', 'verticalalignment': 'center'}

        ##
        ## fundamental data
        ##
        c_axarr = f.add_subplot('321')#axarr[0,0]
        sizeTearsheetLabels(f)
        c_axarr.axis('off')
        c_axarr.text(0,0.7,"CONFIG DATA:\n"+aosat_cfg.CFG_SETTINGS['ts_title']+"\n\n"+aosat_cfg.repString(aosat_cfg.CFG_SETTINGS),fontproperties=font,transform=c_axarr.transAxes,**alignment)
        c_axarr.imshow(logo,alpha=0.2,zorder=1)
        ##
        ## plot result of all analyzers
        ##
        while plots_done <= num_plots:
            if plots_done%6 == 0:
                ##
                ## new page
                ##
                f.suptitle(r'AOSAT Simulation Tear Sheet - %s @ %s$\mu$m' % (aosat_cfg.CFG_SETTINGS['ts_title'],aosat_cfg.CFG_SETTINGS['an_lambda']*1e6))
                plt.tight_layout()
                pdf.savefig(f)
                plt.close()
                plt.clf()
                f = plt.figure(figsize=(8.27,11.69))

            ppage = plots_done % 6
            index = int('32'+str(int(ppage)+1))
            logger.debug("Plots on page:     %s"% ppage)
            logger.debug("Plots done so far: %s"% plots_done)
            logger.debug("Current index:     %s" % index)
            skwa={}
            f = analyzers[plots_done-1].make_plot(fig=f,index=index,subplotkwargs=skwa)
            sizeTearsheetLabels(f)
            plots_done +=1
        f.suptitle(r'AOSAT Simulation Tear Sheet - %s @ %s$\mu$m' % (aosat_cfg.CFG_SETTINGS['ts_title'],aosat_cfg.CFG_SETTINGS['an_lambda']*1e6),y=1.0)
        plt.tight_layout()
        pdf.savefig(f)
        plt.close()
        ##

def makeTearsheetReport(analyzers):
    """Helper function to generate a textual report
    of the tearsheet


    Parameters
    ----------
    analyzers : list of analyzers
        Executed and then called one by one to produce the
        report.


    Returns
    -------
    nothing



    """

    filename = aosat_cfg.CFG_SETTINGS["ts_basefilename"]+".txt"

    report = ""

    report += "## CONFIG DATA:\n## "+aosat_cfg.CFG_SETTINGS["ts_title"]+"\n##\n"+aosat_cfg.repString(aosat_cfg.CFG_SETTINGS,prefix='## ')
    for analyzer in analyzers:
        report += analyzer.make_report()
    with open(filename, "w") as text_file:
        text_file.write(report)








def run(analyzers):
    """Run the analysis using a set of Analyzers.

    All frames of residual wavefronts are served to each
    analyzer one by one. When finished, the 'finalize()'
    method of each analyzer is called.

    Thus, all analyzers should contain the relevant informative
    attributes upon completion, and be ready to have the
    'make_report()' and 'make_plot()' methods called.


    Parameters
    ----------
    analyzers : list of analyzers


    Returns
    -------
    nothing


    """


    fs = frameserver.frameServer()
    logger.info("Analyzing...")
    for (frame,tfnum,totnum) in fs:
        if tfnum == 1:
            #rd = initResult(sd,totnum)
            util.printProgressBar(logger,1, totnum, prefix = 'Analyzing: ', suffix = 'Complete', length = 50)
            barsteps = (np.arange(20)/20.0*totnum).astype(int)
        for analyzer in analyzers:
            analyzer.feed_frame(frame,totnum)
        if tfnum in barsteps:
            util.printProgressBar(logger,tfnum, totnum, prefix = 'Analyzing: ', suffix = 'Complete', length = 50)
    util.printProgressBar(logger,totnum, totnum, prefix = 'Analyzing: ', suffix = 'Complete', length = 50)
    logger.info("Preparing reports...")
    for analyzer in analyzers:
        analyzer.finalize()

def tearsheet(config_file):
    """Generate a 'tearsheet', i.e. a collection
    of informative plots for a simulation

    Parameters
    ----------
    config_file : string
        (Full path to) configuration file for the analysis.
        In addition to the standard descriptive parameters, tearsheet
        evaluates the following:

        ts_basefilename  - base filename for the plot and report file.
                           If it contains path information, it's relative
                           to the config file!

        ts_title         - title of the tearsheet written on top.

        A copy of the config file's content is printed in the first
        plot space on the tearsheet. So in principle you can add more
        parameters than are being evaluated!

    Returns
    -------
    nothing.

    Examples
    -------

    >>> tearsheet('examples/example_analyze_closed_loop.setup')

    """

    ##
    ## prepare config
    ##
    aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(config_file)
    logger.debug("\n"+aosat_cfg.repString(aosat_cfg.CFG_SETTINGS))
    reload(fftx)

    ##
    ## set up analysis
    ##
    sd=setup()
    logger.debug("\n"+aosat_cfg.repString(sd))

    ##
    ## add all available analyzers, then run
    ##
    analyzers=[psf_analyzer(sd), frg_analyzer(sd), phs_analyzer(sd), zrn_analyzer(sd),tvc_analyzer(sd,ctype='icor'), tvc_analyzer(sd)]
    run(analyzers)

    ##
    ## plot
    ##
    makeTearsheetFigure(analyzers)
    makeTearsheetReport(analyzers)



if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
