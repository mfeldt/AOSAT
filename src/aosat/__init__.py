# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = 'unknown'
finally:
    del get_distribution, DistributionNotFound

__all__ = ["aosat_cfg","fftx","frameserver","analyze","analyzers_"]

#import aosat_cfg

import logging, logging.config




from . import aosat_cfg

repdict = aosat_cfg.configure(None)
aosat_cfg.configureLogging(repdict)
logger = logging.getLogger(__name__)
logger.debug('Completed configuring logger()!')

from . import fftx
from . import frameserver
from . import analyze
