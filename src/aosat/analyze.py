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


##
##
## Simple PSF Analyzer
##
##

from aosat._analyzers.frg_analyzer import frg_analyzer
from aosat._analyzers.zrn_analyzer import zrn_analyzer
from aosat._analyzers.dmy_analyzer import dmy_analyzer
from aosat._analyzers.phs_analyzer import phs_analyzer
from aosat._analyzers.tvc_analyzer import tvc_analyzer
from aosat._analyzers.psf_analyzer import psf_analyzer

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
    In particular the following parameters are used:

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
    setup_dict['tel_mirror'] = np.nan_to_num(np.array(1.0*pyfits.getdata(os.path.join(aosat_cfg.CFG_SETTINGS['setup_path'],aosat_cfg.CFG_SETTINGS['pupilmask']))))

    if 'embed_frame' in aosat_cfg.CFG_SETTINGS:
        des_x_size = aosat_cfg.CFG_SETTINGS["embed_frame"][0]
        des_y_size = aosat_cfg.CFG_SETTINGS["embed_frame"][1]
        xs = setup_dict['tel_mirror'].shape[1]
        ys = setup_dict['tel_mirror'].shape[0]

        embd = ((int((des_x_size-xs)/2),
                int(des_x_size-xs-int((des_x_size-xs)/2))),
                (int((des_y_size-ys)/2),
                int(des_y_size-ys-int((des_y_size-ys)/2))))
        setup_dict['tel_mirror'] = np.pad(setup_dict['tel_mirror'],embd,'constant')
        logger.debug("Embedded mirror in %s,%s array" % (aosat_cfg.CFG_SETTINGS['embed_frame'][0],aosat_cfg.CFG_SETTINGS['embed_frame'][1]))

    ext_x  = np.atleast_1d(np.max(np.where(setup_dict['tel_mirror'].sum(axis=1) != 0)[0])-np.min(np.where(setup_dict['tel_mirror'].sum(axis=1) != 0)[0])+1).item()
    ext_y  = np.atleast_1d(np.max(np.where(setup_dict['tel_mirror'].sum(axis=0) != 0)[0])-np.min(np.where(setup_dict['tel_mirror'].sum(axis=0) != 0)[0])+1).item()
    setup_dict['pupildiam']  = max([ext_x,ext_y])

    logger.debug("Found pupil diameter of %s pixels" % setup_dict['pupildiam'])

    setup_dict['sdim']          = setup_dict['tel_mirror'].shape[0] # always needs to be squared...
    setup_dict['divpm']         = setup_dict['tel_mirror'] * 1.0
    setup_dict['divpm'][np.where(setup_dict['tel_mirror'] == 0.0)] = 1.0 # pupil mask for division if necessary
    setup_dict['wnz']           = np.where(setup_dict['tel_mirror'] != 0.0)
    setup_dict['fragmask']      = np.array(ndimage.label(util.ensure_numpy(setup_dict['tel_mirror'])>1e-6)[0])

    logger.debug("Found %s independent pupil fragments." % np.max(setup_dict['fragmask']))
    logger.info("Setting up FFTs...")
    setup_dict['fft_plan']      = fftx.FFTmakeplan(setup_dict['tel_mirror'])
    setup_dict['fft_out']       = fftx.FFTmakeout(setup_dict['tel_mirror'])
    setup_dict['cfg']           = aosat_cfg.CFG_SETTINGS
    setup_dict['ppm']           = aosat_cfg.CFG_SETTINGS['ppm'] # pix per metre
    setup_dict['aspp']          = aosat_cfg.CFG_SETTINGS['ppm']*3600/(np.pi/180.0)/setup_dict['sdim'] *aosat_cfg.CFG_SETTINGS['an_lambda'] # arc sec per pix
    setup_dict['crad']          = aosat_cfg.CFG_SETTINGS['an_lambda']/(5.0/setup_dict['ppm'])/np.pi*180*3600 # estimate of control radius assuming                                                                                                     # a sampling of 5 pixels per actuator, in arc sec
    setup_dict['dl']            = aosat_cfg.CFG_SETTINGS['an_lambda']/(setup_dict['pupildiam']/setup_dict['ppm'])*3600*180/np.pi # diffraction limit
    setup_dict['loopfreq']      = aosat_cfg.CFG_SETTINGS['loopfreq'] # loop frequency

    x, y     = np.mgrid[-setup_dict['sdim']/2:setup_dict['sdim']/2,-setup_dict['sdim']/2:setup_dict['sdim']/2]/setup_dict['ppm']
    setup_dict['wtrk']          = np.where((setup_dict['tel_mirror'] != 0 ) * ((x-x.astype(int)) ==0) * ((y-y.astype(int))==0))

    logger.info("Making Modal basis with %s terms ..." % aosat_cfg.CFG_SETTINGS['zterms'])
    #setup_dict['zernike_basis'] = zernike.zernike_basis(nterms=aosat_cfg.CFG_SETTINGS['zterms'],npix=int(setup_dict['pupildiam']))
    mreg = slice (int(setup_dict['sdim']/2-setup_dict['pupildiam']/2),int(setup_dict['sdim']/2+setup_dict['pupildiam']/2))
    setup_dict['zernike_basis'] = np.array(zernike.arbitrary_basis(util.ensure_numpy(setup_dict['tel_mirror'][mreg,mreg]),nterms=aosat_cfg.CFG_SETTINGS['zterms'],outside=0.0))
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

        If set to None, no new configuration will be preformed.

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
    if config_file is not None:
        aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(config_file)
    aosat_cfg.configureLogging(aosat_cfg.CFG_SETTINGS)
    logger.debug("\n"+aosat_cfg.repString(aosat_cfg.CFG_SETTINGS))
    logger.debug("\n"+aosat_cfg.repString(aosat_cfg.LOG_SETTINGS))

    reload(fftx)

    ##
    ## set up analysis
    ##
    sd=setup()
    logger.debug("\n"+aosat_cfg.repString(sd))

    ##
    ## add all available analyzers, then run
    ##
    analyzers=[zrn_analyzer(sd)]#psf_analyzer(sd), frg_analyzer(sd), phs_analyzer(sd), zrn_analyzer(sd),tvc_analyzer(sd,ctype='icor'), tvc_analyzer(sd)]
    run(analyzers)

    ##
    ## plot
    ##
    makeTearsheetFigure(analyzers)
    makeTearsheetReport(analyzers)



if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
