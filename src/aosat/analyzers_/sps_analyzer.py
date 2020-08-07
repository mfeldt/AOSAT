from aosat.analyzers_.dmy_analyzer import dmy_analyzer

from aosat import aosat_cfg
from aosat import fftx
from aosat import frameserver
from aosat import util

from astropy import units
from pip._internal.utils.misc import get_installed_distributions
if any(["cupy" in str(f) for f in get_installed_distributions()]):
    import cupy as np
else:
    import numpy as np

import matplotlib.pyplot as plt

import logging
logger = logging.getLogger(__name__)


def powerSpectrum(ps,mask=None,radialIndex=None):
    nX             = ps.shape[0]
    nWaveNumber    = int(nX / 2)
    x, y           = np.mgrid[-nX/2:nX/2,-nX/2:nX/2]
    waveNumberDist = np.sqrt(x**2+y**2)
    if mask is None:
        fftArray       = np.fft.fftshift(np.fft.fft2(ps,norm='ortho'))
    else:
        fftArray       = np.fft.fftshift(np.fft.fft2(ps*mask,norm='ortho'))
    fftArray       = (fftArray * np.conj(fftArray)).astype(np.float)

    outArray       = np.arange(nWaveNumber)*1.0
    if radialIndex is None:
        radialIndex = []
        for i in range(nWaveNumber):
            radialIndex.append(np.where((waveNumberDist >= (i -0.5)) * (waveNumberDist <= (i+0.5))))
    for i in range(nWaveNumber):
        outArray[i] = fftArray[radialIndex[i]].mean()
    return(outArray,radialIndex)


class sps_analyzer(dmy_analyzer):

    def __init__(self,sd):
        self.sd           = sd
        self.mask         = util.apodize_mask(self.sd['tel_mirror'] != 0)          # apodized mask
        self.pcf          = 1./frameserver.getUnitCnv(units.Unit('nm'),self.sd['cfg']['an_lambda'])   # phase conversion factor to nm
        xsize=self.sd['tel_mirror'].shape[0]
        self.pix_scale    = self.sd['ppm']                                         # spatial frequency
        self.f_sampling   = self.sd['ppm']#1./self.pix_scale                                      # frequency sampling
        self.delta_f      = self.f_sampling/xsize                                  # e.g. 10/4096 is delta_f
        f_spatial         = np.arange(xsize/2) + 1.
        self.f_spatial    = f_spatial * self.delta_f                               # spatial frequency vector
        self.radial_index = None
        self.psd          = np.zeros(int(xsize/2))
        self.ps_psd       = None


    def feed_frame(self,frame,nframes):
        this_psd,ri      = powerSpectrum(frame*self.pcf,mask=self.mask,radialIndex=self.radial_index)
        self.radialIndex = ri
        self.psd += this_psd/nframes



    def finalize(self):
        self.ps_psd    = self.psd * 1./self.pix_scale


    def make_plot(self,fig=None,index=111,plotkwargs={},subplotkwargs={}):

        if fig is None:
            fig = plt.figure()

        #if 'ylim' not in subplotkwargs:
        #    subplotkwargs['ylim'] = (1e-8,1)
        #if 'xlim' not in subplotkwargs:
        #    subplotkwargs['xlim'] =(0,self.sd['crad']*1.2)
        if 'xlabel' not in subplotkwargs:
            subplotkwargs['xlabel'] = 'Spatial frequency. [1/m]'
        if 'ylabel' not in subplotkwargs:
            subplotkwargs['ylabel'] = r'Power Spectral Density [nm$^2$/m$^{-1}$]'
        if 'title' not in subplotkwargs:
            subplotkwargs['title'] = r'Spatial PSD'
        if 'yscale' not in subplotkwargs:
            subplotkwargs['yscale'] = 'log'

        ol_psd = self.f_spatial**(-11.0/3.0)
        ol_psd = ol_psd /ol_psd[-1]*self.ps_psd[-1] # poor man's fit...

        if 'L0' in self.sd:
            logger.debug("Plotting Karman spectrum!")
            vk_psd    = (self.f_spatial**2+(1/self.sd['L0'])**2)**(-11.0/6.0)
            vk_psd    = vk_psd /vk_psd[-1]*self.ps_psd[-1] # poor man's fit...

        ax = fig.add_subplot(index,**subplotkwargs,label=str(index*2))
        ax.loglog(util.ensure_numpy(self.f_spatial[1:]),util.ensure_numpy(self.ps_psd[1:]),color='k',**plotkwargs)
        ax.loglog(util.ensure_numpy(self.f_spatial[1:]),util.ensure_numpy(ol_psd[1:]),color='r',**plotkwargs)
        if 'L0'in self.sd:
            ax.loglog(util.ensure_numpy(self.f_spatial[1:]),util.ensure_numpy(vk_psd[1:]),color='b',**plotkwargs)

        ax.text(0.5,0.8,r'Kolmogorov, $k^{-11/3}$',color='r',size=6,transform=ax.transAxes)
        if 'L0' in self.sd:
            ax.text(0.5,0.75,r'von Kármán, $(k^2+(1/L_0)^2)^{-11/6}$,',color='b',size=6,transform=ax.transAxes)
            ax.text(0.5,0.7,r'$L_0=%s$m'%self.sd['L0'],color='b',size=6,transform=ax.transAxes)
        return(fig)



    def make_report(self):

        from numpy import array2string

        report =  "##\n##\n"
        report += "## reporting analyzer: %s\n##\n##\n" % self.__class__.__name__
        report += "## Spatial frequencies sampled:\n"
        report += "f_spatial = %s\n##\n" % array2string(util.ensure_numpy(self.f_spatial),separator=',',precision=4)
        report += "## PSD in nm^2/m^-1:\n"
        report += "f_spatial = %s\n##\n" % array2string(util.ensure_numpy(self.ps_psd),separator=',',precision=4)
        return(report)
