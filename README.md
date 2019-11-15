<img src="/src/aosat/img/aosat_logo.png" width="350px"/>


Adaptive Optics Simulation Analysis Tool


Description
===========

Extensive analysis of sets of residual wavefront phasescreens as output by many adaptive optics simulation tools


Installing while still in development:
======================================

1. Make a new empty directory, e.g. like this:
```
mkdir AOSAT
cd AOSAT
```

2. Clone the repository
```
git clone https://github.com/mfeldt/AOSAT.git
```

3. Create and activate a virtual environment (if you don't know/have virtualenv, use pip to install it)
```
virtualenv -p python3.6 venv
source venv/bin/activate
```
(This assumes using bash, there's also venv/bin/activate.csh and a few others)

4. Change to the repository and install:
```
cd AOSAT
python setup.py install
```

That's it, python should install the package and all required dependencies!

Verifying the installation:
===========================

While still verifying that the installation is fine, you can do the following:

1. run the test suite
```
python setup.py test
```

2. Try the individual files
```
cd src/aosat
python fftx.py
python aosat_cfg.py
python util.py
python analyze.py
```
Ideally everything should terminate without failures. Beware it may take a while.


New Tearsheet:
==============
The last of the above tests already produces a tearsheet from the examples (there are
only 2 analyzers as of yet).

This is made from the examples provided using a config/setup file that looks as follows:

```
loopfreq        = 100.0
ppm             = 16.0
screen_fpattern = "yao_rwf*.fits"        # filename pattern
screen_dir      = "ExampleAnalyze"       # directory for phasescreens, realtive to *this* file!
phaseunit       = "micron"               # unit of phase screens, must be a valid astropy.units unit
skipstep        = 1                      # step between consecutive frames, 1=use all
startskip       = 0                      # frames to skip before evaluation starts, 0= start from beginning
an_lambda       = 1.0e-6                 # analysis wavelength in metres
pupilmask       = 'ExampleAnalyze/yao_pupil.fits' # relative to this file!
aosat_fft_forcenumpy = True
ts_basefilename = 'ts_test'              # if path contained, it's relative to *this* file!
ts_title        = 'Example TS'
```

This is essentially modeled after the old "report" files from simmetis_tearsheet.
The primary information to make a tearsheet should be clear from this example.
The screen_fpattern should allow to find all phasesecreen files, each must carry a unique number
in the filename, they will be sorted numerically whatever the number format.
Files can contain multiple screens along NAXIS3.

Everything else in the setup file should of course fit the format of these files!

To run a tearsheet of your own, prepare a file with content similar to the above (but matching
your circumstances), and run in a python shell:
```
import aosat
aosat.analyze.tearsheet('path/to/mysetupfile')
```

Editing / Changing Code:
========================
There is a unit test suite that can be run from the setup.py file:

```
python setup.py test
```
In addition to the tests already provided by running the individual modules like in no. 2 under installation,
this performs a number of tests ensuring that the functionality, logic, and unit conversions of the code work correctly.

So whenever you change something, make sure you don't break the tests! If you have to, because something was seriously wrong,
think thrice and talk to someone about it!


For the adventurous:
====================

You can try to install a working openCL or CUDA environment on your machine, accompanied
by either pyopencl (with pybin11) or pycuda, plus reikna.

Then you can create a configuration file with the following contents
```
aosat_fft_forcenumpy = False
```
and read it with 
```
aosat.aosat_cfg.configure('/path/to/myfile')
```
If you're lucky, you will  then see some interesting log messages about some GPU features
being selected for doing ffts. Try
```
aosat.fftx.fftSpeedTest()
```
to see if you actually gain.
If yes, maybe enjoy another tearsheet?


Note
====

This project has been set up using PyScaffold 3.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
