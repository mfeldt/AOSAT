import pytest, math
from aosat import util
from aosat import aosat_cfg
import os

# from pip._internal.utils.misc import get_installed_distributions
# if any(["cupy" in str(f) for f in get_installed_distributions()]):
#     import cupy as np
# else:
#     import numpy as np
import numpy as np

from astropy import units
import pdb

def test_rv():
    a=np.random.randn(100)
    b=(0.,0.,0.)
    for i in range(100):
        b = util.rolling_variance(b,a[i])
    assert b[0] == 100
    assert math.isclose(b[1],a.mean())
    assert math.isclose(b[2]/b[0], a.var())

def test_zernike():
    x,y = np.mgrid[-256:256,-256:256]
    r=(x**2+y**2)**0.5
    tm = (r<=255)*1.0
    zb = util.zernike_basis(10,512,tm)
    s = np.array([zb[i]*i/1e1 for i in range(10)]).sum(axis=0)
    v=util.basis_expand(s,zb,tm)
    assert not any(np.abs(v-np.arange(10)/1e1)>3e-3)
