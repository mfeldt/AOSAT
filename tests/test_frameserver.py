#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import aosat
from aosat import frameserver
from aosat import aosat_cfg
from aosat import analyze
import os
import numpy as np
from astropy import units
from astropy.io import fits as pyfits

import pdb

__author__ = "Markus Feldt"
__copyright__ = "Markus Feldt"
__license__ = "mit"


def test_fs():
    cfg_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','src','aosat','examples','example_frameserver.setup')
    aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(cfg_file)
    fs = frameserver.frameServer()
    for i in range(10):
        assert next(fs)[0].mean() == i+1

def test_fs2():
    cfg_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','src','aosat','examples','example_frameserver.setup')
    aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(cfg_file)
    ##
    ## skip first two frames, and serve every second
    ##
    ## we thus expect frames 3,5,7, and 9 (starting to count with 1)
    ##
    aosat_cfg.CFG_SETTINGS['startskip'] = 2
    aosat_cfg.CFG_SETTINGS['skipstep'] = 2

    exp_result = [3,5,7,9]
    rc = 0
    fs = frameserver.frameServer()

    for f in fs:
        assert f[0].mean() == exp_result[rc]
        rc +=1

def test_uc():
    tunit     = units.Unit('m')
    an_lambda = 3.0e-6 # always in metres!
    pcf = frameserver.getUnitCnv(tunit,an_lambda)
    assert pcf == 2094395.1023931953

def test_uc2():
    tunit     = units.Unit('deg')
    an_lambda = 3.0e-6 # always in metres!
    pcf = frameserver.getUnitCnv(tunit,an_lambda)
    assert pcf == 0.017453292519943295

def test_uc3():
    cfg_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','src','aosat','examples','example_analyze_closed_loop.setup')
    aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(cfg_file)
    fs = frameserver.frameServer()
    frame = next(fs)[0] ## should be in rad at 3 mum


    bdir = os.path.join(os.path.dirname(os.path.abspath(aosat.__file__)),'examples')

    ## phase screen in rad at 1 micron
    ps = pyfits.getdata(os.path.join(bdir,'ExampleClosedLoop','metis_370P_35L_rwf800.fits'))/3*2*np.pi
    tm = pyfits.getdata(os.path.join(bdir,'ExampleAnalyze','yao_pupil.fits'))
    wp = np.where(tm != 0)

    sigma_ref = ps[wp].std()
    sigma_fs  = frame[wp].std()
    assert sigma_ref == pytest.approx(sigma_fs.item(),1e-4)
