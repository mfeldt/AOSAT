#from pip._internal.utils.misc import get_installed_distributions
# if any(["cupy" in str(f) for f in get_installed_distributions()]):
#     import cupy as np
# else:
#     import numpy as np
import numpy as np

from poppy import zernike
import logging
logger = logging.getLogger(__name__)

from astropy.io import fits as pyfits
import os
dn = dn = os.path.dirname(os.path.realpath(__file__))
pupil_file = os.path.join(dn, 'examples/ExampleAnalyze/yao_pupil.fits')

# tel_mirror =  pyfits.getdata('examples/ExampleAnalyze/yao_pupil.fits')
tel_mirror = pyfits.getdata(pupil_file)

def ensure_numpy(a):
    if 'cupy' in str(type(a)):
        return(np.asnumpy(a))
    else:
        return(a)


def printProgressBar (logger,iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '#'):
    """
    Call in a loop to create terminal progress bar
    @params:
        logger      - Required  : logger instance to send output to
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """

    xxxdot = 3
    if decimals >= 1:
      xxxdot = 4 + decimals

    percent = ("{0:" + str(xxxdot) +"." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    logger.info("\r%s |%s| %s%% %s" % (prefix, bar, percent, suffix))
    # Print New Line on Complete
    if iteration == total:
        logger.info("%s Done!" % prefix)

def rad_order(frame,outer_radius=np.inf):
    """Compute the radial order of pixels from the centre.

    Parameters
    ----------
    frame : 2D array, square(!)
        Input frame of the shape for which the 2D order should be calculated
    outer_radius : float
        Radius (in pixels) out to which the order should be calculated.
        Default is to calculate fora all pixels

    Returns
    -------
    tuple (r,pixels,rord)
        r - array pof the same shape as frame where each pixel contains
            the radial distance to the centre

        pixels - indices of the pixels closer to the centre than
                 outer_radius.

        rord - Containing radial order of the pixels. Can be used
               to index e.g. the original array as follows:
               frame[pixels][rord]

    Examples
    -------

    >>> a = np.zeros((512,512))
    >>> o = rad_order(a)


    """

    sdim=frame.shape[0] ##

    x, y = np.mgrid[-sdim/2:sdim/2,-sdim/2:sdim/2]
    r    = np.sqrt(x**2+y**2)
    pixels = np.where(r<outer_radius)
    rord    = np.argsort(r[pixels])
    return((r,pixels,rord))

def rolling_variance(old,newValue):
    """Rolling variance computation according to  Welford's algorithm.
    https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_Online_algorithm

    Parameters
    ----------
    old : tuple containing the triplet
        count, mean, M2. Initialize with (0,0,0)
        on the first call
    newValue : new element to be added to the sample

    Returns
    -------
    type
        tuple containing the triplet
        count, mean, M2. To retrieve variance and
        sample variance compute M2/count and M2/(count-1),
        respectively

    Examples
    -------

    >>> a=np.arange(10)
    >>> b=(0.,0.,0.)
    >>> for i in range(10):
    ...     b=rolling_variance(b,a[i])
    >>> b
    (10.0, 4.5, 82.5)


    """
    (count, mean, M2) = old
    count += 1
    delta = newValue - mean
    mean += delta / count
    delta2 = newValue - mean
    M2 += delta * delta2

    return (count, mean.item(), M2.item())

def zernike_basis(nterms,npix,tel_mirror):
    """Generate a Zernike basis

    Note that even though called "zernike", it actually calls
    poppy's "arbitrary_basis", which uses a Gram-Schmidt
    orthonormalization to adapt to the actual pupil!

    For this reason, the aperture mask (tel_mirror) must be
    supplied!

    Parameters
    ----------
    nterms : int
        Number of basis functions to generate
    npix : int
        Number of pixels in each screen (will be npix x npix)
    tel_mirror : 2D array
        Aperture mask representing the telescope mirror.
        If larger than npix x npix the central area will be used.


    Returns
    -------
    3D array of shape (nterms,npix,npix)
        Containing the Zernike basis screens

    Examples
    -------
    Examples should be written in doctest format, and
    should illustrate how to use the function/class.
    >>> zb = zernike_basis(10,512,tel_mirror)
    >>> zb.shape
    (10, 512, 512)


    """
    if tel_mirror is not None:
        tms = tel_mirror.shape[0]
        tm = tel_mirror[int(tms/2-npix/2):int(tms/2+npix/2),int(tms/2-npix/2):int(tms/2+npix/2)]
    else:
        logger.error("Aperture (tel_mirror) must be supplied!")
    #zbase = zernike.zernike_basis(nterms=nterms,npix=npix,aperture=tm)
    zbase = np.array(zernike.arbitrary_basis(ensure_numpy(tm),nterms=nterms,outside=0.0))
    return(zbase)


def basis_expand(wf, basis, tel_mirror):
    """Axpand a wavefront into a basis set, operating
    on an aperture

    Parameters
    ----------
    wf : ndarray, 2D
        input wavefront
    basis : ndarray, 3D
        stack of basis functions to expand the wavefront in.
        Must be of dimension (m,n,n), when wf
        contains m elements. If n is smaller than the
        dimension of wf, wf will be cropped. In such a case,
        make sure that the outer regions of wf contain zero
        padding only!
    tel_mirror : ndarray, 2D
        Mask for region where fit is performed.
        Regions where mask == 0 will be ignored,
        values should be 1 or 0.

    Returns
    -------
    type
        Coefficients of the basis expansion
    Examples
    -------

    >>> zb = zernike_basis(10,512,tel_mirror)
    >>> s = np.array([zb[i]*i/10.0 for i in range(10)]).sum(axis=0)
    >>> tm = np.isfinite(s)*1.0
    >>> basis_expand(s,zb,tm)
    array([...

    """
    wsize    = wf.shape
    bsize    = basis[0].shape
    opd      = wf[int(wsize[0]/2-bsize[0]/2):int(wsize[0]/2+bsize[0]/2),int(wsize[0]/2-bsize[1]/2):int(wsize[0]/2+bsize[1]/2)]
    aperture = tel_mirror[int(wsize[0]/2-bsize[0]/2):int(wsize[0]/2+bsize[0]/2),int(wsize[0]/2-bsize[0]/2):int(wsize[0]/2+bsize[0]/2)]
    wgood    = np.where((aperture != 0.0) & np.isfinite(basis[1]))
    ngood    = (wgood[0]).size
    coeffs = [(opd * b)[wgood].sum() / ngood for b in basis]
    return np.array(coeffs)



if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
