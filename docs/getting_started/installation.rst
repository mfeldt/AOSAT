============
Installation
============


Installing while still in development
=====================================

1. Make a new empty directory, e.g. like this:

.. code::

  $ mkdir AOSAT
  $ cd AOSAT


2. Clone the repository

.. code::

  $ git clone https://github.com/mfeldt/AOSAT.git


.. _virtualenv:

3. Create and activate a virtual environment (if you don't know/have virtualenv, use pip to install it)

.. code::

  $ virtualenv -p python3.6 venv
  $ source venv/bin/activate

(This assumes using bash, there's also venv/bin/activate.csh and a few others)

4. Change to the repository and install:

.. code-block::

 $ cd AOSAT
 $ python setup.py install


That's it, python should install the package and all required dependencies!

Verifying the installation
==========================

While still verifying that the installation is fine, you can do the following:

1. run the test suite

.. code-block::

  $ python setup.py test


2. Try the individual files

.. code-block::

  $ cd src/aosat
  $ python fftx.py
  $ python aosat_cfg.py
  $ python util.py
  $ python analyze.py

Ideally everything should terminate without failures. Beware it may take a while.

Preparing for GPU computation
=============================

AOSAT supports working with CUDA on NVIDIA cards via CUPY.
By default, cupy is not required, and all computations are executed on the CPU.

Since CUDA can on modern (and expensive...) cards gain afactor of about 5 in
computation speed, you might wish to install nvidia drivers, CUDA, and cupy.

Currently, the latest available cupy version is cupy-cuda102. However,
due to a bug in the 10.x versions of CUDA, `SVDs are not possible on large arrays
<https://github.com/cupy/cupy/issues/2351>`_.,and thus these versions don't work
with AOSAT. Until version cupy-cuda110 becomes available (and CUDA 11 has been
confirmed to not suffer from the issue), you need to install CUDA9.2 and
cupy-cuda92.

A rough outline of the installation process is thus:

1. Verify you have a CUDA-capable NVIDIA card
2. Make sure you have the appropriate graphics driver installed and active
3. Install `CUDA 9.2 <https://developer.nvidia.com/cuda-92-download-archive>`_.
4. Install cupy_cuda92

With the first three, you're on your own, the last one should be done by
entering a command prompt in the active :ref:`virtualenv created in step 3<virtualenv>` of the
original installation and issue

.. code-block::

  $ pip install cupy-cuda92


Once the above steps are completed successfully, AOSAT will silently perform
all array operations on the GPU.  It is recommended to run all tests described
in the _verification. step above again.  You should notice a hefty speed increase.

To switch off using cupy-cuda, uninstall it from the virtualenv.


CUDA and OpenCL via reikna
==========================

Support is in principle implemented, but discouraged.
There is no gain in computation time on most machines due to frequent synchronizations.
Support will likely be scrapped in forthcoming versions.
