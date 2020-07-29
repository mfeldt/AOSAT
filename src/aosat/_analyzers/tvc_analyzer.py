

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
from aosat import analyze
from scipy import ndimage
from poppy import zernike
from importlib import reload

import matplotlib.pyplot as plt
import matplotlib as mp
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
        self.max_no_cor = 0.0

    def feed_frame(self,frame,nframes):
        in_field  = self.sd['tel_mirror']*np.exp(1j*frame)
        #this_strehl   = np.exp(-1*np.std(frame[self.sd['wnz']])**2) ## remember: Frames are in radians!
        #logger.debug("Found Strehl ratio: %s" % this_strehl)
        frame_dev = fftx.FFTprepare(in_field- self.sd['tel_mirror']*int(self.ctype=='icor'))
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




    def make_plot(self,fig=None,index=121,plotkwargs={},subplotkwargs={},plotkwargsC={},subplotkwargsC={}):

        if fig is None:
            fig =plt.figure(figsize=(8.27,11.69/3.0))
            no_fig=True
        else:
            no_fig=False

        pidx  = index - (index//10)*10
        nrows = index // 100
        ncols = (index - nrows*100)//10
        #import pdb; pdb.set_trace()
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
        ##  create first subplot
        ##
        logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargs))
        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))

        rep_vec_r     = util.ensure_numpy(np.arange(1000)/1000*max(self.rvec))

        from numpy import interp, array2string

        rep_vec_cmean = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.cvecmean))
        rep_vec_cmin  = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.cvecmin))
        rep_vec_cmax  = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.cvecmax))

        rep_vec_rcmean = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.rcvecmean))
        rep_vec_rcmin  = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.rcvecmin))
        rep_vec_rcmax  = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.rcvecmax))



        ax.fill_between(rep_vec_r,rep_vec_cmin,rep_vec_cmax,alpha=0.15,**plotkwargs)
        ax.plot(util.ensure_numpy(self.rvec),util.ensure_numpy(self.cvecmean),**plotkwargs)
        ax.plot(rep_vec_r,rep_vec_cmean,**plotkwargs)

        if self.corrlen > 0:
            num_frames_per_hour = 3600.0/(self.corrlen/self.sd['loopfreq'])
            plotkwargsP = copy.copy(plotkwargs)
            plotkwargsP['color'] = 'green'
            ax.plot(rep_vec_r,rep_vec_cmean/(num_frames_per_hour**0.5),**plotkwargsP)

        ax.text(0.5,0.9,'No photon noise,',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        ax.text(0.5,0.85,'due to PSF variation only!',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        ax.text(0.5,0.95,r'5 $\sigma$ TVC',transform=ax.transAxes,size=6,ha='left',color=plotkwargs['color'])
        plotkwargs2 = copy.copy(plotkwargs)
        plotkwargs2['color']='black'
        ax.fill_between(rep_vec_r,rep_vec_rcmin,rep_vec_rcmax,alpha=0.15,**plotkwargs2)
        ax.plot(rep_vec_r,rep_vec_rcmean,**plotkwargs2)

        # if self.corrlen > 0:
        #
        #     plotkwargsP = copy.copy(plotkwargs2)
        #     plotkwargsP['color'] = 'green'
        #     ax.plot(util.ensure_numpy(self.rvec),util.ensure_numpy(self.rcvecmean)/(num_frames_per_hour**0.5),**plotkwargsP)
        ax.text(0.5,0.8,'Raw (PSF profile)',transform=ax.transAxes,size=6,ha='left',color=plotkwargs2['color'])
        # if self.corrlen>0:
        ax.text(0.5,0.75,'Estd. best 1hr ADI contrast',transform=ax.transAxes,size=6,ha='left',color=plotkwargsP['color'])
        #ax.set_aspect(subplotkwargs['aspect'],'box')

        if self.ctype=='icor':
            ##
            ## add coronagraphic psf
            ##
            ##
            ## default appearance
            ##

            if pidx > (nrows*ncols)-1 and self.ctype=='icor':
                logger.error("Not enough free plot positions on page - not adding TVC PSF plot!")
                return(fig)

            if 'xlabel' not in subplotkwargsC:
                subplotkwargs['xlabel'] = r'$\delta$RA [mas]'
            if 'ylabel' not in subplotkwargsC:
                subplotkwargs['ylabel'] = r'$\delta$DEC [mas]'

            # size of extraction area
            sdim = self.mean.shape[0]
            plts = int(min([np.around(self.sd['crad']/self.sd['aspp']/2.0)*2,sdim/2]))
            sd2 = int(sdim/2)
            if 'cmap' not in plotkwargsC:
                plotkwargsC['cmap'] = 'nipy_spectral'
            if 'norm' not in plotkwargsC:
                plotkwargsC['norm'] = LogNorm(vmin=1e-7,vmax=0.5)
            if 'extent' not in plotkwargsC:
                plotkwargsC['extent'] = [plts*self.sd['aspp']*1000,-plts*self.sd['aspp']*1000,-plts*self.sd['aspp']*1000,plts*self.sd['aspp']*1000]
            if 'origin' not in plotkwargsC:
                plotkwargsC['origin'] = 'lower'
            if 'title' not in subplotkwargsC:
                subplotkwargsC['title'] = 'Coronagraphic PSF'
            ##
            ##  create (only) subplot
            ##
            logger.debug("Subplot keyword args:\n"+aosat_cfg.repString(subplotkwargsC))
            logger.debug("Plot keyword args:\n"+aosat_cfg.repString(plotkwargsC))
            nindex = nrows*100+ncols*10+pidx+1
            ax2 = fig.add_subplot(nindex,**subplotkwargsC,label=str(index*2))
            im = ax2.imshow(util.ensure_numpy(self.mean[sd2-plts:sd2+plts,sd2-plts:sd2+plts]/np.max(self.max_no_cor)), **plotkwargsC)
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes("right", size="10%", pad=0.05)
            cax.tick_params(axis='both', which='major', labelsize=6)
            cax.tick_params(axis='both', which='minor', labelsize=5)
            t = [1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]
            cbar = plt.colorbar(im, cax=cax, ticks=t,label='Intensity')

        if no_fig:
            analyze.sizeTearsheetLabels(fig)
            plt.tight_layout()
        return(fig)

    def make_report(self):

        rep_vec_r     = util.ensure_numpy(np.arange(50)*self.sd['aspp'])

        from numpy import interp, array2string

        rep_vec_cmean = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.cvecmean))
        rep_vec_cmin  = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.cvecmin))
        rep_vec_cmax  = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.cvecmax))

        rep_vec_rcmean = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.rcvecmean))
        rep_vec_rcmin  = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.rcvecmin))
        rep_vec_rcmax  = interp(rep_vec_r,util.ensure_numpy(self.rvec),util.ensure_numpy(self.rcvecmax))

        report =  "##\n##\n"
        report += "## reporting analyzer: %s\n##\n##\n" % self.__class__.__name__
        report += "## contrast type: %s\n" % self.ctype
        report += "## Interpolated contrast values:\n\n"
        report += "tvc_%s_sepvec    = %s\n" % (self.ctype,array2string(rep_vec_r,separator=',',precision=4))
        report += "tvc_%s_ctrstmean = %s\n" % (self.ctype,array2string(rep_vec_cmean,separator=','))
        report += "tvc_%s_ctrstmin  = %s\n" % (self.ctype,array2string(rep_vec_cmin,separator=','))
        report += "tvc_%s_ctrstmax  = %s\n\n\n" % (self.ctype,array2string(rep_vec_cmax,separator=','))
        report += "raw_%s_sepvec    = %s\n" % (self.ctype,array2string(rep_vec_r,separator=',',precision=4))
        report += "raw_%s_ctrstmean = %s\n" % (self.ctype,array2string(rep_vec_rcmean,separator=','))
        report += "raw_%s_ctrstmin  = %s\n" % (self.ctype,array2string(rep_vec_rcmin,separator=','))
        report += "raw_%s_ctrstmax  = %s\n##\n" % (self.ctype,array2string(rep_vec_rcmax,separator=','))
        report += "## estimated correlation time [ms] \n##\n"
        report += "estim_corr_time   = %s \n\n\n" % (1000*self.corrlen/self.sd['loopfreq'])
        return(report)


    def finalize(self):

        variance      = self.variance[2]/self.variance[0]
        mean          = self.variance[1]#*int(self.ctype == 'nocor') + self.variance2[1]*int(self.ctype=='icor')
        sigma         = variance**0.5
        self.max_no_cor    = np.max(self.variance[1]*int(self.ctype == 'nocor') + self.variance2[1]*int(self.ctype=='icor'))
        self.mean     = mean

        self.contrast = 5*sigma/self.max_no_cor
        self.rcontrast = mean/self.max_no_cor
        ra,pixels,rord = util.rad_order(self.contrast)

        self.rvec  = ra[pixels][rord]*self.sd['aspp']
        c          = pd.Series(util.ensure_numpy(self.contrast[pixels][rord]))
        r          = pd.Series(util.ensure_numpy(self.rcontrast[pixels][rord]))

        self.cvecmean = c.rolling(window=50).mean()
        self.cvecmin  = c.rolling(window=50).min()
        self.cvecmax  = c.rolling(window=50).max()
        self.rcvecmean = r.rolling(window=50).mean()
        self.rcvecmin  = r.rolling(window=50).min()
        self.rcvecmax  = r.rolling(window=50).max()

        max_corr_len=0.0
        for i in range(self.ntracks):
            acff,conf    = acf(util.ensure_numpy(self.ptrak[:,i]),fft=True,alpha=0.05,nlags=min([500,len(self.ptrak[:,i])-1]))
            where_uncorr = np.where(np.array(acff) < (np.array(conf[:,1]) - np.array(acff)))[0]
            if len(where_uncorr) > 0:
                max_corr_len = max([max_corr_len,np.min(where_uncorr)])
        self.corrlen = max_corr_len.item()*1.0
