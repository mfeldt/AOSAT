import os
import numpy as np
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
    Examples should be written in doctest format, and
    should illustrate how to use the function/class.
    >>> sd={}
    >>> a1 = tvc_analyzer(sd)

    Attributes
    ----------
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

    def feed_frame(self,frame):
        in_field  = self.sd['tel_mirror']*np.exp(1j*frame)
        this_strehl   = np.exp(-1*np.std(frame[self.sd['wnz']])**2) ## remember: Frames are in radians!
        frame_dev = fftx.FFTprepare(in_field- (this_strehl**0.5)*self.sd['tel_mirror']*(self.ctype=='icor'))
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

        ax = fig.add_subplot(index,**subplotkwargs)
        ax.fill_between(self.rvec,self.cvecmin,self.cvecmax,alpha=0.15,**plotkwargs)
        ax.plot(self.rvec,self.cvecmean,**plotkwargs)
        ax.text(0.5,0.95,'No photon noise,',transform=ax.transAxes,size=6,ha='left')
        ax.text(0.5,0.9,'due to PSF variation only!',transform=ax.transAxes,size=6,ha='left')

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
        self.contrast = 5*sigma/np.max(self.variance[1]*int(self.ctype=='raw')+self.variance2[1]*int(self.ctype=='icor'))

        r,pixels,rord = util.rad_order(self.contrast)

        self.rvec  = r[pixels][rord]*self.sd['aspp']
        c          = pd.Series(self.contrast[pixels][rord])
        self.cvecmean = c.rolling(window=50).mean()
        self.cvecmin  = c.rolling(window=50).min()
        self.cvecmax  = c.rolling(window=50).max()

        self.make_report()
        self.make_plot()




def getFontProperties():
    font0 = FontProperties()
    font = font0.copy()
    font.set_family('monospace')
    font.set_size(5)
    return(font)


def assignMaskFromStringOrValue(sd,keyword,defval=0.0):
    """Assigns a setup dictionary keyword with a vale
    (if defval is a number) or a mask read from a FITS file
    (if defval is a string)

    Parameters
    ----------
    sd : setup dictionary

    keyword : the keyword in the dictionary to define/overwrite
        Description of parameter `keyword`.
    defval : string or number
        If it is a string will be interpreted as the name of FITS file
        that contains the mask that will be assigned as value to the key.
        If it is an umber, it will be assigned directly (the default is 0.0).

    Returns
    -------
    type
        swetup dictionary with key/value-pair set.
    Examples
    -------

    >>> sd={}
    >>> sd = assignMaskFromStringOrValue(sd,'Pupil')
    >>> sd['Pupil']
    0.0

    """
    if keyword in aosat_cfg.CFG_SETTINGS:
        if type(keyword) == str:
            logger.info("Reading mask %s from %s!" % (keyword,os.path.join(aosat_cfg.CFG_SETTINGS['setup_path'],aosat_cfg.CFG_SETTINGS[keyword])))
            sd[keyword] = pyfits.getdata(os.path.join(aosat_cfg.CFG_SETTINGS['setup_path'],aosat_cfg.CFG_SETTINGS[keyword]))
        else:
            sd[keyword] = aosaot_cfg.CFG_SETTINGS[keyword]*1.0
    else:
        sd[keyword] = sd['tel_mirror']*0.0+defval
    return(sd)

def setupC2(sd):
    """Generate a Lyot coronagraph with 2l/D

    Parameters
    ----------
    sd : setup dictionary


    Returns
    -------
    setup dictionary


    Examples
    -------
    Examples should be written in doctest format, and
    should illustrate how to use the function/class.
    >>> sd = {'tel_mirror':1.0}
    >>> sd = setupC2(sd)
    >>> sd['c2_pp1_phase']
    0.0


    """

    sd = assignMaskFromStringOrValue(sd,"c2_pp1_ampl",defval=1.0)
    sd = assignMaskFromStringOrValue(sd,"c2_pp1_phase",defval=0.0)
    sd = assignMaskFromStringOrValue(sd,"c2_fs_ampl",defval=0.0)
    if sd['c2_fs_ampl'].sum()==0:
        # make a 2l/D field stop
        x, y     = np.mgrid[-sd['sdim']/2:sd['sdim']/2,-sd['sdim']/2:sd['sdim']/2]
        r        = np.sqrt(x**2+y**2)
        sd['c2_fs_ampl'] = 1- (r*sd['aspp'] < 2*sd['dl'])*1.0 # coronagraph mask
    sd = assignMaskFromStringOrValue(sd,"c2_fs_phase",defval=0.0)
    sd = assignMaskFromStringOrValue(sd,"c2_ls_ampl",defval=0.0)
    if sd['c2_ls_ampl'].sum()==0:
        sd['c2_ls_ampl'] = sd['tel_mirror'] * 1.0 # Lyot stop equal to aperture
    sd = assignMaskFromStringOrValue(sd,"c2_ls_phase",defval=0.0)

    return(sd)


def setup():
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
    setup_dict['dl']            = aosat_cfg.CFG_SETTINGS['an_lambda']/setup_dict['pupildiam']*3600*180/np.pi # diffraction limit

    setup_dict = setupC2(setup_dict)

    x, y     = np.mgrid[-setup_dict['sdim']/2:setup_dict['sdim']/2,-setup_dict['sdim']/2:setup_dict['sdim']/2]/setup_dict['ppm']
    setup_dict['wtrk']          = np.where((setup_dict['tel_mirror'] != 0 ) * ((x-x.astype(int)) ==0) * ((y-y.astype(int))==0))

    logger.info("Making Zernike basis with %s terms ..." % aosat_cfg.CFG_SETTINGS['zterms'])
    setup_dict['zernike_basis'] = zernike.zernike_basis(nterms=aosat_cfg.CFG_SETTINGS['zterms'],npix=int(setup_dict['pupildiam']))



    return(setup_dict)

def initResult(sd,n):
    result_dict = {}
    result_dict['psf_avg'] = sd['tel_mirror']*0.0
    result_dict['psf_var'] = sd['tel_mirror']*0.0
    result_dict['psf_corop_avg'] = sd['tel_mirror']*0.0
    result_dict['psf_corop_var'] = sd['tel_mirror']*0.0
    result_dict['psf_coro2_avg'] = sd['tel_mirror']*0.0
    result_dict['psf_coro2_var'] = sd['tel_mirror']*0.0
    result_dict['ampl_track']    = np.zeros((len(sd['wtrk'][0]),n,3),dtype=np.complex128)
    result_dict['tt_track']      = np.zeros((n,2))
    result_dict['zernike_track']  = np.zeros((n,sd['zernike_basis'].shape[0]))
    result_dict['piston_track']   = np.zeros((n,np.max(sd['fragmask'])))
    result_dict['power_spectrum'] = np.zeros(int(sd['sdim']/2),dtype=np.float64)

    return(result_dict)

def makeTearsheetFigure(analyzers):
    num_plots  = len(analyzers)
    num_pages  = int(num_plots/6)+1
    logger.debug("Preparing %s plots on %s pages!" % (num_plots,num_pages))
    plots_done = 1

    filename = aosat_cfg.CFG_SETTINGS["ts_basefilename"]+".pdf"
    with PdfPages(filename) as pdf:

        f, axarr = plt.subplots(3,2, figsize=(8.27,11.69))
        for ax in axarr.flatten():
            ax.tick_params(axis='both', which='major', labelsize=6)
            ax.tick_params(axis='both', which='minor', labelsize=5)
            ax.yaxis.label.set_size(8)
            ax.xaxis.label.set_size(8)
        font = getFontProperties()
        logo = mpimg.imread(os.path.join(os.path.dirname(os.path.abspath(__file__)),'aosat_logo.png'))

        alignment = {'horizontalalignment': 'left', 'verticalalignment': 'center'}
        f.suptitle(r'AOSAT Simulation Tear Sheet - %s @ %s$\mu$m' % (aosat_cfg.CFG_SETTINGS['ts_title'],aosat_cfg.CFG_SETTINGS['an_lambda']*1e6))
        ##
        ## fundamental data
        ##
        c_axarr = axarr[0,0]
        c_axarr.axis('off')
        c_axarr.text(0,0.7,"CONFIG DATA:\n"+aosat_cfg.CFG_SETTINGS['ts_title']+"\n\n"+aosat_cfg.repString(aosat_cfg.CFG_SETTINGS),fontproperties=font,transform=axarr[0,0].transAxes,**alignment)
        c_axarr.imshow(logo,alpha=0.2,zorder=1)
        ##
        ## plot result of all analyzers
        ##
        while plots_done <= num_plots:
            if plots_done%6 == 0:
                ##
                ## new page
                ##
                plt.tight_layout()
                pdf.savefig(f)
                plt.close()
                plt.clf()
                f, axarr = plt.subplots(3,2, figsize=(8.27,11.69))
                for ax in axarr.flatten():
                    ax.tick_params(axis='both', which='major', labelsize=6)
                    ax.tick_params(axis='both', which='minor', labelsize=5)
                    ax.yaxis.label.set_size(8)
                    ax.xaxis.label.set_size(8)

            ppage = plots_done % 6
            index = '32'+str(int(ppage)+1)
            logger.debug("Plots on page:     %s"% ppage)
            logger.debug("Plots done so far: %s"% plots_done)
            logger.debug("Current index:     %s" % index)

            dummy = analyzers[plots_done-1].make_plot(fig=f,index=index)
            plots_done +=1
        pdf.savefig(f)
        plt.close()
        ##

def makeTearsheetReport(analyzers):
    filename = aosat_cfg.CFG_SETTINGS["ts_basefilename"]+".txt"

    report = ""

    report += "## CONFIG DATA:\n## "+aosat_cfg.CFG_SETTINGS["ts_title"]+"\n##\n"+aosat_cfg.repString(aosat_cfg.CFG_SETTINGS,prefix='## ')
    for analyzer in analyzers:
        report += analyzer.make_report()
    with open(filename, "w") as text_file:
        text_file.write(report)








def run(analyzers):


    fs = frameserver.frameServer()
    logger.info("Analyzing...")
    for (frame,tfnum,totnum) in fs:
        if tfnum == 1:
            #rd = initResult(sd,totnum)
            util.printProgressBar(logger,1, totnum, prefix = 'Analyzing: ', suffix = 'Complete', length = 50)
            barsteps = (np.arange(20)/20.0*totnum).astype(int)
        for analyzer in analyzers:
            analyzer.feed_frame(frame)
        if tfnum in barsteps:
            util.printProgressBar(logger,tfnum, totnum, prefix = 'Analyzing: ', suffix = 'Complete', length = 50)
    util.printProgressBar(logger,totnum, totnum, prefix = 'Analyzing: ', suffix = 'Complete', length = 50)
    logger.info("Preparing reports...")
    for analyzer in analyzers:
        analyzer.finalize()

def tearsheet(setup_file):

    ##
    ## prepare config
    ##
    aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(setup_file)
    reload(fftx)

    ##
    ## set up analysis
    ##
    sd=setup()

    ##
    ## add all available analyzers, tzhen run
    ##
    analyzers=[tvc_analyzer(sd), tvc_analyzer(sd,ctype='icor')]
    run(analyzers)

    ##
    ## plot
    ##
    makeTearsheetFigure(analyzers)
    makeTearsheetReport(analyzers)
