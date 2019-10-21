import os
import numpy as np
import scipy
from astropy import units
from astropy.io import fits as pyfits
import pandas as pd

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

    This is a dummy class with no actually functionality for ocumentation purposes.
    It might be used to derive subclasses, but these can also be written independently.

    All analyzers should provide the methods of this one, and all methods
    should accept the aparemters of the ones implemented here,  Accepting
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
            total number of rames to be fed


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
        'This is a dummy with no actuctual functionality'
        """
        report = "This is a dummy with no actuctual functionality"
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
        self.psf += psf

        logger.debug("Solving tilt...")
        ## tilt from WF
        C,_,_,_ = scipy.linalg.lstsq(self.A, frame[self.sd['wnz']])#/2/np.pi*self.sd['cfg']['an_lambda']) # frame in m
        tphase = C[0]*self.x + C[1]*self.y + C[2]
        self.ttilt[ffed] = np.max(tphase[self.sd['wnz']]) - np.min(tphase[self.sd['wnz']])

        logger.debug("Fitting position...")
        ## tip and tilt from PSF fitted position
        result_lmf = self.lmf(self.M, self.xx, self.yy, psf[int(self.sdim/2-16):int(self.sdim/2+16),int(self.sdim/2-16):int(self.sdim/2+16)])
        self.ttx[self.ffed] = (result_lmf.x_mean - 16.0)*self.sd['aspp']
        self.tty[self.ffed] = (result_lmf.y_mean - 16.0)*self.sd['aspp']
        logger.debug("Found tt excursion of %.3f, %.3f from PSF!"%(self.ttx[self.ffed],self.tty[self.ffed]))

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


    def make_plot(self,fig=None,index=111,plotkwargs={},subplotkwargs={}):

        if fig is None:
            fig = plt.figure()

        ##
        ## default appearance
        ##
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = '$\delta$RA[mas]'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = '$\delta$DEC[mas]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = 'PSF'

        # size of extraction area
        plts = int(min([np.round(sd['crad']/1000.0*1.2/sd['aspp']/2.0)*2,self.sdim/2]))

        if 'map' not in plotkwargs:
            plotkwargs['map'] = 'nipy_spectral'
        if 'norm' not in plotkwargs:
            plotkwargs['norm'] = LogNorm(vmin=1e-7,vmax=0.5)
        if 'extent' not in plotkwargs:
            plotkwargs['extent'] = [plts*sd['aspp']*1000,-plts*['aspp']*1000,-plts*sd['aspp']*1000,plts*sd['aspp']*1000]
        if 'origin' not in plotkwargs:
            plotkwargs['origin'] = 'lower'
        ##
        ##  create (only) subplot
        ##
        logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargs))
        logger.debug("Plot keyword args:\n"+aosat_cfg.repString(plotkwargs))
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))

        im = ax.imshow(psf[self.sdim/2-plts:self.sdim/2+plts,self.sdim/2-plts:self.sdim/2+plts]/np.max(psf),
                           **plotkwargs)

        ## useful info
        ax.text(0.7,0.9,'SR (PSF) = %.3f' % self.strehl,transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.7,0.85,'SR (WF) = %.3f' % self.sr_wf,transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.7,0.8,r'$\sigma_{tt}$(PSF) = %.3f mas' % self.ttjit_psf,transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.7,0.75,r'$\sigma_{tt}$(WF) = %.3f mas' % self.ttjit,transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.7,0.8,r'$\Q90_{tt}$(PSF) = %.3f mas' % self.ttq90_psf,transform=ax.transAxes,size=6,ha='left',color='white')
        ax.text(0.7,0.75,r'$\Q90_{tt}$(WF) = %.3f mas' % self.ttq90,transform=ax.transAxes,size=6,ha='left',color='white')

        ## color bar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="10%", pad=0.05)
        cax.tick_params(axis='both', which='major', labelsize=6)
        cax.tick_params(axis='both', which='minor', labelsize=5)
        t = [1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]
        cbar = plt.colorbar(im, cax=cax, ticks=t,label='Intensity')




    def finalize(self):
        frame_dev = fftx.FFTprepare(self.sd['tel_mirror']*np.exp(1j*0.0))
        fftframe  = fftx.FFTshift(fftx.FFTforward(self.sd['fft_plan'],self.sd['fft_out'],frame_dev))
        psf_ref   = np.abs(tfftframe)**2
        logger.debug("Max of reference PSF: %s" % np.max(psf_ref))
        logger.debug("Max of measured PSF:  %s" % np.max(self.variance[1]))
        self.psf_mean = self.psf/self.ffed
        self.strehl   = np.max(self.psf_mean)/np.max(psf_ref)

        ttmas = self.ttilt/(37.0)/np.pi*180.0*3600.0*1000.0

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
    Examples should be written in doctest format, and
    should illustrate how to use the function/class.
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
    def __init__(self,sd,ctype='raw'):
        self.ctype     = ctype
        self.sd        = sd
        self.variance  = (0,None,None)
        self.variance2 = (0,None,None)
        self.contrast  = None
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
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = 'Sep. [arc sec.]'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = 'Contrast [Peak]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = r'Simple Single-Pixel 5$\sigma$ %s Contrast' % self.ctype
        if 'yscale' not in subplotkwargs:
            subplotkwargs['yscale'] = 'log'

        ##
        ##  create (only) subplot
        ##
        logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargs))
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))
        ax.fill_between(self.rvec,self.cvecmin,self.cvecmax,alpha=0.15,**plotkwargs)
        ax.plot(self.rvec,self.cvecmean,**plotkwargs)
        ax.text(0.5,0.95,'No photon noise,',transform=ax.transAxes,size=6,ha='left')
        ax.text(0.5,0.9,'due to PSF variation only!',transform=ax.transAxes,size=6,ha='left')
        return(fig)

    def make_report(self):

        rep_vec_r     = np.arange(30)*self.sd['aspp']

        rep_vec_cmean = np.interp(rep_vec_r,self.rvec,self.cvecmean)
        rep_vec_cmin  = np.interp(rep_vec_r,self.rvec,self.cvecmin)
        rep_vec_cmax  = np.interp(rep_vec_r,self.rvec,self.cvecmax)

        report =  "##\n##\n"
        report += "## reporting analyzer: %s\n##\n##\n" % self.__class__.__name__
        report += "## contrast type: %s\n" % self.ctype
        report += "## Interpolated contrast values:\n\n"
        report += "tvc_%s_sepvec    = %s\n" % (self.ctype,np.array2string(rep_vec_r,separator=',',precision=4))
        report += "tvc_%s_ctrstmean = %s\n" % (self.ctype,np.array2string(rep_vec_cmean,separator=','))
        report += "tvc_%s_ctrstmin  = %s\n" % (self.ctype,np.array2string(rep_vec_cmin,separator=','))
        report += "tvc_%s_ctrstmax  = %s\n\n\n" % (self.ctype,np.array2string(rep_vec_cmax,separator=','))
        return(report)

    def finalize(self):

        variance      = self.variance[2]/self.variance[0]
        sigma         = variance**0.5

        self.contrast = 5*sigma/np.max(self.variance[1]*int(self.ctype == 'raw') + self.variance2[1]*int(self.ctype=='icor'))

        r,pixels,rord = util.rad_order(self.contrast)

        self.rvec  = r[pixels][rord]*self.sd['aspp']
        c          = pd.Series(self.contrast[pixels][rord])
        self.cvecmean = c.rolling(window=50).mean()
        self.cvecmin  = c.rolling(window=50).min()
        self.cvecmax  = c.rolling(window=50).max()





def getFontProperties():
    """Helper function.

    Parameters
    ----------
    none

    Returns
    -------
    font
        fonr for matplotlib

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
        float giving th epixels per metre in the aperture plane

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
    ext_x  = np.max(np.where(setup_dict['tel_mirror'].sum(axis=1) != 0))-np.min(np.where(setup_dict['tel_mirror'].sum(axis=1) != 0))
    ext_y  = np.max(np.where(setup_dict['tel_mirror'].sum(axis=0) != 0))-np.min(np.where(setup_dict['tel_mirror'].sum(axis=0) != 0))
    setup_dict['pupildiam']  = max([ext_x,ext_y])

    logger.debug("Found pupil diameter of %s pixels" % setup_dict['pupildiam'])

    setup_dict['sdim']          = setup_dict['tel_mirror'].shape[0] # always needs to be squared...
    setup_dict['divpm']         = setup_dict['tel_mirror'] * 1.0
    setup_dict['divpm'][np.where(setup_dict['tel_mirror'] == 0.0)] = 1.0 # pupil mask for division if necessary
    setup_dict['wnz']           = np.where(setup_dict['tel_mirror'] != 0.0)[0]
    setup_dict['fragmask']      = ndimage.label(setup_dict['tel_mirror']>1e-6)[0]

    logger.debug("Found %s independent pupil fragments." % np.max(setup_dict['fragmask']))
    logger.info("Setting up FFTs...")
    setup_dict['fft_plan']      = fftx.FFTmakeplan(setup_dict['tel_mirror'])
    setup_dict['fft_out']       = fftx.FFTmakeout(setup_dict['tel_mirror'])
    setup_dict['cfg']           = aosat_cfg.CFG_SETTINGS
    setup_dict['ppm']           = aosat_cfg.CFG_SETTINGS['ppm'] # pix per metre
    setup_dict['aspp']          = aosat_cfg.CFG_SETTINGS['ppm']*3600/(np.pi/180.0)/setup_dict['sdim'] *aosat_cfg.CFG_SETTINGS['an_lambda'] # arc sec per pix
    setup_dict['dl']            = aosat_cfg.CFG_SETTINGS['an_lambda']/(setup_dict['pupildiam']/setup_dict['ppm'])*3600*180/np.pi # diffraction limit


    x, y     = np.mgrid[-setup_dict['sdim']/2:setup_dict['sdim']/2,-setup_dict['sdim']/2:setup_dict['sdim']/2]/setup_dict['ppm']
    setup_dict['wtrk']          = np.where((setup_dict['tel_mirror'] != 0 ) * ((x-x.astype(int)) ==0) * ((y-y.astype(int))==0))

    logger.info("Making Zernike basis with %s terms ..." % aosat_cfg.CFG_SETTINGS['zterms'])
    setup_dict['zernike_basis'] = zernike.zernike_basis(nterms=aosat_cfg.CFG_SETTINGS['zterms'],npix=int(setup_dict['pupildiam']))



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
    method of each analyzer is calledself.

    Thus, all analyzers should contain the relvant informative
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
        (Full path to) configuration file for the analysisself.
        In addition to the standard descriptive parameters, tearsheet
        evaluates the following:

        ts_basefilename  - base filename for the plot and report file.
                           If it contains path information, it's relative
                           to the config file!

        ts_title         - title of the tearsheet written on topself.

        A copy of the config file's content is printed in the first
        plot space on the tearsheet. So in pronciple you can add more
        parameters than are being evaluated!

    Returns
    -------
    nothing.

    Examples
    -------

    >>> tearsheet('examples/example_analyze.setup')

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
    analyzers=[psf_analyzer(sd), tvc_analyzer(sd,ctype='icor'), tvc_analyzer(sd)]
    run(analyzers)

    ##
    ## plot
    ##
    makeTearsheetFigure(analyzers)
    makeTearsheetReport(analyzers)



if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
