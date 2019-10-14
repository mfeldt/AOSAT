import pytest, math
from aosat import util
from aosat import aosat_cfg
import os
import numpy as np
from astropy import units

def test_rv():
    a=np.random.randn(100)
    b=(0.,0.,0.)
    for i in range(100):
        b = util.rolling_variance(b,a[i])
    assert b[0] == 100
    assert math.isclose(b[1],a.mean())
    assert math.isclose(b[2]/b[0], a.var())

def test_zernike():
    zb = util.zernike_basis(10,512,None)
    s = np.array([zb[i]*i/1e9 for i in range(10)]).sum(axis=0)
    tm = np.isfinite(s)*1.0
    v=util.basis_expand(s,zb,tm)
    assert not any(np.abs(v-np.arange(10)/1e9)>2e-11)
