<img src="/src/aosat/img/aosat_logo.png" width="350px"/>


Adaptive Optics Simulation Analysis Tool


Description
===========

Extensive analysis of sets of residual wavefront phasescreens as output by many adaptive optics simulation tools


Installing while still in development:
====================

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
(This assumes using bash, theres also venv/bin/activate.csh and a few others)

4. Change to the repository and install:
```
cd AOSAT
python setup.py install
```

That's it, python should install th epackage and all required dependencies!

New Tearsheet
---------------
* mehr Doku folgt ab morgen


Note
====

This project has been set up using PyScaffold 3.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
