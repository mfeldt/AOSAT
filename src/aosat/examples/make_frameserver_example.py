import numpy as np
from astropy.io import fits as pyfits
import os

ecube = np.zeros((10,1024,1024))
for i in range(10):
    ecube[i] += i+1

pyfits.writeto(os.path.join('ExampleFrameserver','exmpl_frame1.fits'),ecube[0],overwrite=True)
pyfits.writeto(os.path.join('ExampleFrameserver','exmpl_frame2.fits'),ecube[1:8],overwrite=True)
pyfits.writeto(os.path.join('ExampleFrameserver','exmpl_frame9.fits'),ecube[8:],overwrite=True)

config = """
screen_fpattern = "exmpl*.fits"          # filename pattern
screen_dir      = "ExampleFrameserver"   # directory for phasescreens, relative to *this* file!
phaseunit       = "radian"               # unit of phase screens
skipstep        = 1                      # step between consecutive frames, 1=use all
startskip       = 0                      # frames to skip before evaluation starts, 0= start from beginning
totalnumber     = 10                     # number of frames to analyze in total
an_lambda       = 1.0                    # analysis wavelength in micrometres
rm_glob_pist    = False                  # don't remove average piston for this test
"""

with open("example_frameserver.setup", "w") as text_file:
    text_file.write(config)
