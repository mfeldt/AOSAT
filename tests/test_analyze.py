import pytest, math
import aosat
from aosat import util
from aosat import analyze
from aosat import aosat_cfg
import os

from pip._internal.utils.misc import get_installed_distributions
if any(["cupy" in str(f) for f in get_installed_distributions()]):
    import cupy as np
else:
    import numpy as np

from astropy import units
from scipy import ndimage

from astropy.io import fits as pyfits
import numpy.polynomial.polynomial as poly

import pdb

def test_psf_analyze_strehl():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    ## phase screen in rad at 1 micron
    ps = np.array(pyfits.getdata(os.path.join(bdir,'ExampleClosedLoop','metis_370P_35L_rwf800.fits'))*2*np.pi)
    tm = np.array(pyfits.getdata(os.path.join(bdir,'ExampleAnalyze','yao_pupil.fits'))*1.0)
    wp = np.where(tm != 0)

    ## Strehl at 1 micron
    sr   = np.exp(-1*(ps)[wp].std()**2).item()
    psf  = np.abs(np.fft.fft2(tm*np.exp(1j*ps)))**2
    psfr = np.abs(np.fft.fft2(tm*np.exp(1j*ps*0)))**2
    srp  = (np.max(psf)/np.max(psfr)).item()

    sd = analyze.setup()

    ## make analyzer and runs
    an = analyze.psf_analyzer(sd)
    an.feed_frame(ps,2)
    an.feed_frame(ps,2)
    an.finalize()
    assert np.abs(an.sr_wf.mean() - sr) < 1e-5
    assert np.abs(an.strehl - srp) < 1e-5


def polyfit2d(x, y, f, deg):
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = poly.polyvander2d(util.ensure_numpy(x), util.ensure_numpy(y), util.ensure_numpy(deg))
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(np.array(vander), np.array(f),rcond=-1)[0]
    return c.reshape(util.ensure_numpy(deg+1))


def test_psf_analyze_ttwf():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    tm = np.array(pyfits.getdata(os.path.join(bdir,'ExampleAnalyze','yao_pupil.fits')))
    wp = np.where(tm != 0)
    sdim = tm.shape[0]


    sd = analyze.setup()


    x,y = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]/sd['ppm']

    ps = x/3600.0/180.0*np.pi/1000 # tilted by 1mas
    ps *= 1e6*2*np.pi         # to micron to rad

    ## independent 2D plane fit
    c = polyfit2d(x,y,ps,[1,1])

    tip  = c[0,1] # in rad per metre
    tilt = c[1,0] # in rad per metre

    tilt *= 1e-6/2/np.pi/np.pi*180*3600.0 # in arc seconds
    tip  *= 1e-6/2/np.pi/np.pi*180*3600.0 # in arc seconds


    texc = (tip**2+tilt**2)**0.5

    ## make analyzer and runs
    an = analyze.psf_analyzer(sd)
    an.feed_frame(ps,1)
    an.finalize()

    assert np.abs(an.ttjit - texc*1000.0) < 2e-3


def test_psf_analyze_ttpsf():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    tm = np.array(pyfits.getdata(os.path.join(bdir,'ExampleAnalyze','yao_pupil.fits')))
    wp = np.where(tm != 0)
    sdim = tm.shape[0]


    sd = analyze.setup()


    x,y = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]/sd['ppm']

    ps = x/3600.0/180.0*np.pi/1000.0 # tilted by 1mas
    ps *= 1e6*2*np.pi         # to micron to rad


    ## make analyzer and runs
    an = analyze.psf_analyzer(sd)
    an.feed_frame(ps,1)
    an.finalize()


    assert np.abs(an.ttjit_psf - 1.0) < 2e-3

def test_tvc():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    tm = np.array(pyfits.getdata(os.path.join(bdir,'ExampleAnalyze','yao_pupil.fits'))*1.0)
    wp = np.where(tm != 0)
    sdim = tm.shape[0]

    sd = analyze.setup()

    ##
    ## make phasescreens
    ##
    x,y = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]/sd['ppm']
    ps = x/3600.0/180.0*np.pi/1000 # tilted by 1mas
    ps *= 1e6*2*np.pi         # to micron to rad
    ps2 = ps/2 # titlted by 0.5mas

    ##
    ## run analyzer
    ##
    an = analyze.tvc_analyzer(sd)
    an.feed_frame(ps,2)
    an.feed_frame(ps2,2)
    an.finalize()

    ##
    ## make PSFs and compute TVC
    ##
    psf  = np.abs(np.fft.fftshift(np.fft.fft2(tm*np.exp(1j*ps))))**2
    psf2 = np.abs(np.fft.fftshift(np.fft.fft2(tm*np.exp(1j*ps2))))**2
    s5   = 5*np.abs(psf-psf2)/2 # 5 sigma distance of 2 values
    s5c  = s5 / np.max((0.5*(psf+psf2))) # contrast to peak
    assert np.abs(s5c - an.contrast).sum() < 1e-12


def test_analyze_npupils():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')
    aosat.aosat_cfg.CFG_SETTINGS['setup_path'] = bdir
    aosat.aosat_cfg.CFG_SETTINGS['pupilmask'] = 'frag_pupil.fits'
    sd = analyze.setup()
    assert np.max(sd['fragmask']) == 6


#==========================================================#
def test_frg_worst():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    tm = pyfits.getdata(os.path.join(bdir,'frag_pupil.fits'))
    fr = np.array(ndimage.label(tm>1e-6)[0])

    p1 = np.array(tm *0.0)
    p2 = np.array(tm*0.0)
    p2[np.where(fr == 2)]=1.0

    ##
    ## make analyzer
    ##
    aosat.aosat_cfg.CFG_SETTINGS['setup_path'] = bdir
    aosat.aosat_cfg.CFG_SETTINGS['pupilmask'] = 'frag_pupil.fits'
    sd = analyze.setup()
    an = analyze.frg_analyzer(sd)

    ## feed frames
    for i in range(45):
        an.feed_frame(p1,100)
        an.feed_frame(p2*(1+(np.random.randn()*0.001)),100)
    an.feed_frame(p2*2,100) # highest frame index is 90
    for i in range(9):
        an.feed_frame(p1,100)

    ev = 2.0/2/np.pi*1e-6*1e9 # piston of 2 rad at 1mum in nm
    ##
    ## finalize and check
    ##
    an.finalize()
    assert np.abs(np.max(an.pistframe) - ev).item() < 1e-4


#==========================================================#
def test_frg_tilt():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    tm = pyfits.getdata(os.path.join(bdir,'frag_pupil.fits'))
    fr = np.array(ndimage.label(tm>1e-6)[0])


    ##
    ## make analyzer
    ##
    aosat.aosat_cfg.CFG_SETTINGS['setup_path'] = bdir
    aosat.aosat_cfg.CFG_SETTINGS['pupilmask'] = 'frag_pupil.fits'
    sd = analyze.setup()
    an = analyze.frg_analyzer(sd)


    ##cd test_psf_analyze_strehl
    ## make frames
    ##
    sdim=tm.shape[0]
    x,y = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]/sd['ppm']

    p1 = x/3600.0/180.0*np.pi/1000.0 # tilted by 1mas in x
    p1 *= 1e6*2*np.pi         # to micron to rad


    p2 = p1*1.0
    p2[np.where(fr == 2)]+=1.0

    ## feed frames
    for i in range(5):
        an.feed_frame(p1,100)
        an.feed_frame(p2*(1+(np.random.randn()*0.001)),100)


    ##
    ## finalize and check
    ## tolerance 5 micro arc sec
    ##
    an.finalize()
    from numpy import repeat,abs
    assert abs(util.ensure_numpy(an.ttx) - repeat(1.0,6)).sum().item() < 5e-3
    assert abs(util.ensure_numpy(an.tty) - repeat(0.0,6)).sum().item() < 5e-3


def test_frg_all():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    tm = pyfits.getdata(os.path.join(bdir,'frag_pupil.fits'))
    fr = np.array(ndimage.label(tm>1e-6)[0])


    ##
    ## make analyzer
    ##
    aosat.aosat_cfg.CFG_SETTINGS['setup_path'] = bdir
    aosat.aosat_cfg.CFG_SETTINGS['pupilmask'] = 'frag_pupil.fits'
    sd = analyze.setup()
    an = analyze.frg_analyzer(sd)


    ##cd test_psf_analyze_strehl
    ## make frames
    ##
    sdim=tm.shape[0]
    x,y = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]/sd['ppm']

    p1 = x/3600.0/180.0*np.pi/1000.0 # tilted by 1mas in x
    p1 *= 1e6*2*np.pi         # to micron to rad


    p2=x*0.0
    p2[np.where(sd['fragmask'] == 1)] = 1.0
    p2[np.where(sd['fragmask'] == 2)] = 2.0
    p2[np.where(sd['fragmask'] == 3)] = 3.0
    p2[np.where(sd['fragmask'] == 4)] = 1.0 + p1[np.where(sd['fragmask'] == 4)]
    p2[np.where(sd['fragmask'] == 5)] = 2.0 +2*p1[np.where(sd['fragmask'] == 5)]

    for i in range(10):
        an.feed_frame(p2,10)
    an.finalize()

    mean_pist = np.zeros(6)
    for i in range(1,7,1):
        mean_pist[i-1] = p2[np.where(sd['fragmask']==i)].mean()
        assert np.abs(an.pistont[0,i-1]-mean_pist[i-1])<1e-5

    assert np.all(an.pistframe-p2/2/np.pi*an.sd['cfg']['an_lambda']*1e9) < 1e-3



def test_phs():


    sd = analyze.setup()
    an = analyze.phs_analyzer(sd)


    ##
    ## make frames
    ##
    sdim=sd['tel_mirror'].shape[0]

    p1 = sd['tel_mirror']*0+1.0
    p1[::2] = 0.0 # every other row 0, so SDEV=0.5

    ## feed frames
    nframes=10
    for i in range(nframes):
        an.feed_frame(p1,nframes)


    ##
    ## finalize and check
    ##
    an.finalize()
    assert np.abs(an.rms - 0.5/2/np.pi*sd['cfg']['an_lambda']*1e9).item()<5e-5


def test_zrn_modevec():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')
    aosat.aosat_cfg.CFG_SETTINGS['setup_path'] = bdir
    aosat.aosat_cfg.CFG_SETTINGS['pupilmask'] = 'ExampleAnalyze/yao_pupil.fits'
    sd = analyze.setup()
    an = analyze.zrn_analyzer(sd)

    mv = np.array([0.1,0.2,0.1,0.2,0.5,0.4,0.3,0.2,0.1,0.5,0.1,0.2,0.1,0.2,0.5,0.4,0.3,0.2,0.1,0.5])
    #mv=np.repeat(0,20)
    #mv[0]=1
    assert len(mv) == sd['zernike_basis'].shape[0]
    wsize    = sd['tel_mirror'].shape
    bsize    = sd['zernike_basis'][0].shape
    breg     = slice(int(wsize[0]/2-bsize[0]/2),int(wsize[0]/2+bsize[0]/2))
    aperture = sd['tel_mirror'][breg,breg]
    wgood    = np.where((aperture != 0.0) & np.isfinite(sd['zernike_basis'][1]))
    ngood    = (wgood[0]).size
    p1       = sd['tel_mirror']*0.0
    for i in range(len(mv)):
        p1[breg,breg] += mv[i]*sd['zernike_basis'][i]
    ##
    ## feed frames
    ##
    nframes=10
    for i in np.arange(nframes):
        an.feed_frame(p1*i,nframes)


    ##
    ## finalize and check
    ##
    an.finalize()
    #final set of modes in nm
    fmv  = mv*np.mean(np.arange(10))/2/np.pi*1e3
    fdmv = mv*np.std(np.arange(10))/2/np.pi*1e3
    assert not any(np.abs(fmv-an.modes)>1e-3)
    assert not any(np.abs(fdmv-an.dmodes)>1e-3)
