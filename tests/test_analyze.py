import pytest, math
import aosat
from aosat import util
from aosat import analyze
from aosat import aosat_cfg
import os
import numpy as np
from astropy import units
from astropy.io import fits as pyfits
import numpy.polynomial.polynomial as poly

import pdb

def test_psf_analyze_strehl():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    ## phase screen in rad at 1 micron
    ps = pyfits.getdata(os.path.join(bdir,'ExampleClosedLoop','metis_370P_35L_rwf800.fits'))*2*np.pi
    tm = pyfits.getdata(os.path.join(bdir,'ExampleAnalyze','yao_pupil.fits'))
    wp = np.where(tm != 0)

    ## Strehl at 1 micron
    sr   = np.exp(-1*(ps)[wp].std()**2)
    psf  = np.abs(np.fft.fft2(tm*np.exp(1j*ps)))**2
    psfr = np.abs(np.fft.fft2(tm*np.exp(1j*ps*0)))**2
    srp  = np.max(psf)/np.max(psfr)

    sd = analyze.setup()

    ## make analyzer and runs
    an = analyze.psf_analyzer(sd)
    an.feed_frame(ps,2)
    an.feed_frame(ps,2)
    an.finalize()

    assert an.sr_wf.mean() == sr
    assert an.strehl == srp


def polyfit2d(x, y, f, deg):
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = poly.polyvander2d(x, y, deg)
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f,rcond=-1)[0]
    return c.reshape(deg+1)


def test_psf_analyze_ttwf():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    tm = pyfits.getdata(os.path.join(bdir,'ExampleAnalyze','yao_pupil.fits'))
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

    assert np.abs(an.ttjit - texc*1000.0) < 1e-4


def test_psf_analyze_ttpsf():
    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    tm = pyfits.getdata(os.path.join(bdir,'ExampleAnalyze','yao_pupil.fits'))
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

    tm = pyfits.getdata(os.path.join(bdir,'ExampleAnalyze','yao_pupil.fits'))
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
