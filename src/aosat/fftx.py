# -*- coding: utf-8 -*-
"""
This selects the best (i.e. fastest) FFT routine to be used in AOSAT.

The module will look first for CUDA, then for OpenCL, and fall back
to numpy FFTs if neither is available.
"""

import logging

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

##
##
## Check for GPU support and select FFT function
##
##

import numpy as np
from pip._internal.utils.misc import get_installed_distributions
from packaging import version

from aosat import aosat_cfg

avl_pkgs = [f.project_name for f  in get_installed_distributions()]
avl_vrss = [f.version for f  in get_installed_distributions()]

cuda_vrs = avl_vrss[avl_pkgs.index('pycuda')] if 'pycuda' in avl_pkgs else "0.0"
opcl_vrs = avl_vrss[avl_pkgs.index('pyopencl')] if 'pyopencl' in avl_pkgs else "0.0"
FORCE_NUMPY = aosat_cfg.CFG_SETTINGS['aosat_fft_forcenumpy']

if (FORCE_NUMPY or (version.parse(opcl_vrs) < version.parse("2018.2.4") and version.parse(cuda_vrs) < version.parse("0.94"))):
    ##
    ## use numpy fftx
    ##
    log.debug("Neither CUDA nor OpenCL found or not desired!")
    log.info("Using numpy FFTs!")
    def FFTmakeplan(samplearr):
        """Make a plan for FFTs

        Parameters
        ----------
        arrshape : numpy array
            the shape of the arrays to be transformed

        Returns
        -------
        plan
            FFT plan to be executed upon calls.
            In the current case where numpy FFTs
            are used, plan is None!

        Examples
        -------

        >>> plan = FFTmakeplan(np.zeros((1024,1024)).shape)
        >>> plan == None
        True

        """
        return(None)

    def FFTmakeout(inarr):
        """Make an output for FFTs

        Parameters
        ----------
        inarr : numpy array
            Then input array

        Returns
        -------
        plan
            FFT plan to be exectued upon calls.
            In the current case where numpy FFTs
            are used, plan is None!

        Examples
        -------

        >>> out = FFTmakeout(np.zeros((1024,1024)))
        >>> out == None
        True

        """
        return(None)


    def FFTprepare(inarr):
        """Prepare an array for input int FFTforward or
        FFTinverse

        Parameters
        ----------
        inarr : numpy array
            The array to be transforemed

        Returns
        -------
        numpy array
            In the current case of numpy FFTs, the array itself
        Examples
        -------

        >>> arr  = np.zeros((1024,1024))
        >>> iarr = FFTprepare(arr)
        >>> np.array_equal(arr,iarr)
        True

        """
        return(inarr)
    def FFTforward(plan,outarr,inarr):
        """Perform forward FFT.

        Parameters
        ----------
        plan : fft plan
            In this case of numpy FFTs, plan should be None.
        outarr : device-suitable output array
            The device-suitable output. Should be None.
        inarr : numpy array
            The array to be transformed. Should be square.

        Returns
        -------
        numpy array
            The transformed array

        Examples
        -------

        >>> arr = np.zeros((1024,1024))
        >>> oarr = FFTforward(None,None,arr)
        >>> oarr.dtype
        dtype('complex128')
        >>> oarr.shape
        (1024, 1024)
        """
        return(np.fft.fft2(inarr))
    def FFTinverse(plan,outarr, inarr):
        """Perform inverse FFT.

        Parameters
        ----------
        plan : fft plan
            In this case of numpy FFTs, plan should be None.
        outarr : device-suitable output array
            The device-suitable output. Should be None.
        inarr : numpy array
            The array to be transformed. Should be square.

        Returns
        -------
        numpy array
            The transformed array

        Examples
        -------

        >>> arr = np.random.rand(1024,1024)
        >>> farr = FFTforward(None,None,arr)
        >>> oarr = FFTinverse(None,None,farr)
        >>> np.abs((np.real(oarr)-arr)).sum() < 2e-10
        True
        >>> oarr.dtype
        dtype('complex128')
        >>> oarr.shape
        (1024, 1024)
        """
        return(np.fft.ifft2(inarr))
    def FFTshift(inarr):
        """Shift array transformed by FFT such that
        the zero component is centered

        Parameters
        ----------
        inarr : numpy array
            Output of an FFT transform
        Returns
        -------
        numpy array
            Centered array.
        Examples
        -------
        >>> arr = np.random.rand(1024,1024)
        >>> farr = np.abs(FFTshift(FFTforward(None,None,arr)))
        >>> np.unravel_index(np.argmax(farr, axis=None), farr.shape)
        (512, 512)

        """
        return(np.fft.fftshift(inarr))
else:
    ##
    ## use one of the CLUDA FFTs
    ##
    if version.parse(cuda_vrs) >= version.parse("0.94"):
        log.debug("CUDA version %s found" % cuda_vrs)
        log.info("Using CUDA FFTs!")
        from reikna.cluda import dtypes, cuda_api
        FFTapi = cuda_api()
    else:
        log.debug("OpenCL version %s found" % opcl_vrs)
        log.info("Using OpenCL FFTs!")
        from reikna.cluda import dtypes, ocl_api
        FFTapi = ocl_api()


    import reikna.fft as cludaFFT

    FFTthr = FFTapi.Thread.create()

    def FFTmakeplan(samplearr):
        """Make a plan for FFTs

        Parameters
        ----------
        samplearr : numpy array
            an array of the same type and size as the ones
            you're going to transform a lot later on

        Returns
        -------
        plan
            FFT plan to be exectued upon calls.

        Examples
        -------

        >>> plan = FFTmakeplan(np.zeros((1024,1024),dtype=np.float64))
        >>> plan
        <reikna.core.computation.ComputationCallable object at ...

        """
        fft = cludaFFT.FFT(samplearr.astype(np.complex128))
        return(fft.compile(FFTthr))

    def FFTmakeout(inarr):
        """Make a suitable output object for FFTs

        Parameters
        ----------
        inarr : numpy array, preferrably np.complex64
            Sample of input array, same size and type which you
            are going to transform a lot later on.

        Returns
        -------
        FFTout object
            The output object which must be passed to FFTplan
            later on. Will not be used explicitely

        Examples
        -------
        >>> inarr = np.zeros((1024,1024))
        >>> outobj = FFTmakeout(inarr)
        >>> type(outobj)
        <class 'reikna.cluda.ocl.Array'>

        """
        return(FFTthr.array(inarr.shape, dtype=np.complex128))

    def FFTprepare(inarr):
        """Prepare an input array to suit machine/GPU architecture

        Parameters
        ----------
        inarr : numpy array
            The input array

        Returns
        -------
        inarr_dev
                device dependent input array representation

        Examples
        -------
        >>> arr = np.random.rand(1024,1024)
        >>> arr_dev = FFTprepare(arr)
        >>> type(arr_dev)
        <class 'reikna.cluda.ocl.Array'>

        """
        return(FFTthr.to_device(inarr.astype(np.complex128)))
    def FFTforward(plan,outarr,inarr_dev):
        """Perform forward FFT

        Parameters
        ----------
        plan : fft plan
            The plan produced by FFTmakeplan.
        outarr : reikna.cluda.ocl.Array
            Output array, produced by FFTmakeout
        inarr_dev : reikna.cluda.ocl.Array
            Device-suitable input array representation,
            produced by FFTprepare

        Returns
        -------
        numpy array
            FFT transfor of input array

        Examples
        -------
        >>> arr = np.random.rand(1024,1024)
        >>> plan = FFTmakeplan(arr)
        >>> arr_dev = FFTprepare(arr)
        >>> arr_out = FFTmakeout(arr)
        >>> arr_fft = FFTforward(plan,arr_out,arr_dev)
        >>> ref_fft = np.fft.fft2(arr)
        >>> print(np.testing.assert_almost_equal(arr_fft,ref_fft))
        None

        You can also time the execution to see if it's really faster than the standard numpy FFT:

        >>> import time
        >>> start=time.time()
        >>> for i in range(100):
        ...    arr_fft = FFTforward(plan,arr_out,arr_dev)
        ...
        >>> end=time.time()
        >>> print("100 CLUDA FFTs: ",end-start)
        100 CLUDA FFTs: ...

        >>> start=time.time()
        >>> for i in range(100):
        ...   ref_fft = np.fft.fft2(arr)
        ...
        >>> end=time.time()
        >>> print("100 numpy FFTs: ",end-start)
        100 numpy FFTs: ...

        If that's not the case, you can force the use of numpy's FFTs by putting the line
        force_np_fft = True
        in any configuration or report file!

        """
        plan(outarr, inarr_dev, inverse=0)
        return(outarr.get())
    def FFTinverse(plan,outarr,inarr_dev):
        """Perform inverse FFT

        Parameters
        ----------
        plan : fft plan
            The plan produced by FFTmakeplan.
        outarr : reikna.cluda.ocl.Array
            Output array, produced by FFTmakeout
        inarr_dev : reikna.cluda.ocl.Array
            Device-suitable input array representation,
            produced by FFTprepare

        Returns
        -------
        numpy array
            FFT inverse transform of input array

        Examples
        -------
        >>> arr = np.random.rand(1024,1024)
        >>> plan = FFTmakeplan(arr)
        >>> arr_dev = FFTprepare(arr)
        >>> arr_out = FFTmakeout(arr)
        >>> arr_fft = FFTinverse(plan,arr_out,arr_dev)
        >>> ref_fft = np.fft.ifft2(arr)
        >>> print(np.testing.assert_almost_equal(arr_fft,ref_fft))
        None

        """
        plan(outarr, inarr_dev, inverse=1)
        return(outarr.get())

    def FFTshift(inarr):
        """Shift array transformed by FFT such that
        the zero component is centered

        Parameters
        ----------
        inarr : numpy array
            Output of an FFT transform
        Returns
        -------
        numpy array
            Centered array.
        Examples
        -------
        >>> arr = np.random.rand(1024,1024)
        >>> plan = FFTmakeplan(arr)
        >>> arr_dev = FFTprepare(arr)
        >>> arr_out = FFTmakeout(arr)
        >>> farr    = FFTshift(FFTforward(plan,arr_out,arr_dev))
        >>> np.unravel_index(np.argmax(farr, axis=None), farr.shape)
        (512, 512)

        """
        return(np.fft.fftshift(inarr))

def fftSpeedTest(max_res=13):
    """Test of speed of currently selected FFT vs. numpy.

    This is intended to verify that the selectd OpenCL / CUDA
    FFTs are indeed faster than numpy.  For older GPUs, that
    may not always be the case.

    When called, fftSpeedtest will test a series of array sizes
    and print the ratio of execution times against standardnumpy
    implementations. If you see a lot of number smaller than one,
    in particular in the line for the array size of interest,
    make sure that the golbal aosat_cfg.CFG_SETTINGS contains the
    entry CFG_SETTINGS['aosat_fft_forcenumpy'] = True .
    In aosat, this can be achieved by calling

    aosat.aosat_cfg.configure(setupfile)

    where the setupfile contains the line

    aosat_fft_forcenumpy = True

    The setup (or "report") fiel is also an argument of many
    helper and convenience functions.

    After an execution of configure, issue a

    reload(aosat.fftx)

    This will redefine the fftx.FFT* functions!


    Parameters
    ----------


    Returns
    -------
    nothing


    Examples
    -------
    >>> fftSpeedTest(max_res=4)
    array size =     4 x     4 : gpu speedup = ...


    """
    from time import time

    dtype       = np.complex128
    resolutions = range(2,max_res)
    Nloops      = 20
    rtol        = 1e-7
    atol        = 0
    for n in resolutions:
        shape, axes  = (2**n,2**n), (0,1)
        data         = np.random.rand(*shape).astype(dtype)
        plan         = FFTmakeplan(data)
        rtime, ntime = 0., 0.
        for nloop in range(Nloops):

            data = np.random.rand(*shape).astype(dtype)

            # forward
            t0 = time()
            data_dev = FFTprepare(data)

            fwd = FFTforward(plan,data_dev,data_dev)
            rtime += time() - t0
            t0 = time()
            fwd_ref = np.fft.fft2(data, axes=axes).astype(dtype)
            ntime += time() - t0
            actualf = np.real(fwd * np.conj(fwd))
            desiredf =  np.real(fwd_ref * np.conj(fwd_ref))

            # inverse
            t0 = time()
            data_dev = FFTprepare(data)
            inv = FFTinverse(plan,data_dev,data_dev)
            rtime += time() - t0
            t0 = time()
            inv_ref = np.fft.ifft2(data, axes=axes).astype(dtype)
            ntime += time() - t0
            actuali = np.real(inv * np.conj(inv))
            desiredi =  np.real(inv_ref * np.conj(inv_ref))

            np.testing.assert_allclose(desiredf, actualf, rtol, atol)
            np.testing.assert_allclose(desiredi, actuali, rtol, atol)

        print ('array size = %5d x %5d : gpu speedup = %g' % (2**n, 2**n, ntime / rtime))

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
