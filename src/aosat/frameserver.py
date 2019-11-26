

import re,os,glob, errno
import numpy as np
from astropy import units
from astropy.io import fits as pyfits
from aosat import aosat_cfg

import logging
logger = logging.getLogger(__name__)

def sortFiles(file_list):
    """Sort filenames according to inter number in every name

    If filenames contain more than one integer, the last one
    will be used for sorting!
    The list will be sorted according to these numbers in
    ascending order

    Parameters
    ----------
    file_list : list of strings (filenames)
        A list with all the filenames to read residual
        phase screens from.

    Returns
    -------
    type
        Sorted list.

    Examples
    -------

    >>> list = ['res_phase0125.fits','res_phase1835.fits','res_phase12.fits']
    >>> sortFiles(list)
    array(['res_phase12.fits', 'res_phase0125.fits', 'res_phase1835.fits'],
    ...

    """
    numlist = [int(re.findall(r'\d+', fname)[-1]) for fname in file_list]
    order = np.argsort(numlist)
    fnames = np.array(file_list)
    return(fnames[order])


def getUnitCnv(sunit,an_lambda):
    """Helper function for unit conversion
    """
    angle_length = [(units.meter, units.rad, lambda x: x/(an_lambda)*2*np.pi , lambda x: x / 2/np.pi*units.micron)]
    pcf = (1.0*sunit).to(units.rad,equivalencies = angle_length)
    logger.debug("Phase conversion factor to radians is %s." % pcf)
    return(pcf.value)

def frameServer():
    """Serve phase screens, one after the other

    The frame server is configured by the configuration stored in
    aosat_cfg.CFG_SETTINGS.
    In particular by the keywords

    CFG_SETTINGS['startskip']
        Number of screens to skip before commencing analysis. good
        for avoiding loop closure troubles
    CFG_SETTINGS['skipstep']
        Number of frames to skip over between analyzed frames.
        Can speed up execution during trial runs, or avoid
        high-correlation frames when the simulation frequency is much
        higher than the AO loop frequency
    CFG_SETTINGS['totalnumber']
        Total number of frames to be analyzed, should be equal
        or lower to the number of available frames. If the key
        does not exist, will be all frames available
    CFG_SETTINGS['an_lambda']
        The wavelength at which the analysis will be carried out.
        Used to convert frames to unit rad before serving them.
        I repeat: ALL FRAMES SERVED WILL BE IN UNITS OF RADIANS!
    CFG_SETTINGS['phaseunit']
        Unit of the residual phases stored in the fits files.
        Must be an astropy.units.Unit() readable string for an angle
        or a length unit. Phase screens are always converted to radians
        vi an_lambda!
    CFG_SETTINGS['setup_path']
        Path to the setup file used last for configuration. Defaults
        to current working directory if never invoked.
    CFG_SETTINGS['screen_dir']
        Path to the directory where the screens are in,
        RELATIVE to the above!
    CFG_SETTINGS['screen_fpattern']
        Wildcard pattern that matches all files where phase screens are
        in, but nothing else. Must not contain (parts of) paths!


    Returns
    -------
    type
        (numpy array, int, int)

    Examples
    -------
    A simple usage on a provided example which contains ten
    phase screens, each one an array consisting of 1024x1024
    pixels filled with it's own sequential number would be:

    >>> import os
    >>> from aosat import aosat_cfg
    >>> cfg_file = os.path.join(os.path.dirname(
    ... os.path.realpath(__file__)),'examples','example_frameserver.setup')
    >>> aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(cfg_file)
    >>> fs = frameServer()
    >>> [f[0].mean() for f in fs]
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]


    A slightly more complex example involving skipping
    the first two frames and serving only every second would
    be

    >>> import os
    >>> from aosat import aosat_cfg
    >>> cfg_file = os.path.join(os.path.dirname(
    ... os.path.realpath(__file__)),'examples','example_frameserver.setup')
    >>> aosat_cfg.CFG_SETTINGS = aosat_cfg.configure(cfg_file)
    >>> aosat_cfg.CFG_SETTINGS['startskip'] = 2
    >>> aosat_cfg.CFG_SETTINGS['skipstep'] = 2
    >>> fs = frameServer()
    >>> [f[0].mean() for f in fs]
    [3.0, 5.0, 7.0, 9.0]

    """

    skipstep    = aosat_cfg.CFG_SETTINGS['skipstep']
    startskip   = aosat_cfg.CFG_SETTINGS['startskip']
    an_lambda   = aosat_cfg.CFG_SETTINGS['an_lambda']
    setup_path  = aosat_cfg.CFG_SETTINGS['setup_path']
    ##
    ## phase conversion factor
    ##
    logger.info("Phase screen unit is %s." % aosat_cfg.CFG_SETTINGS['phaseunit'])
    pcf = getUnitCnv(units.Unit(aosat_cfg.CFG_SETTINGS['phaseunit']),an_lambda)

    ##
    ## file list
    ##
    file_list = sortFiles(glob.glob(os.path.join(aosat_cfg.CFG_SETTINGS['setup_path'],
                                                  aosat_cfg.CFG_SETTINGS['screen_dir'],
                                                  aosat_cfg.CFG_SETTINGS['screen_fpattern'])))
    if len(file_list) == 0:
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), os.path.join(aosat_cfg.CFG_SETTINGS['setup_path'],
                                                          aosat_cfg.CFG_SETTINGS['screen_dir'],
                                                          aosat_cfg.CFG_SETTINGS['screen_fpattern']))
    ##
    ## division mask
    ##
    tel_mirror_divpm =  pyfits.getdata(os.path.join(aosat_cfg.CFG_SETTINGS['setup_path'],aosat_cfg.CFG_SETTINGS['pupilmask']))
    tel_mirror_divpm[tel_mirror_divpm == 0.0] = 1.0
    if 'divide_phase_by_mask' in aosat_cfg.CFG_SETTINGS:
        if aosat_cfg.CFG_SETTINGS == True:
            pass
        else:
            tel_mirror_divpm = tel_mirror_divpm *0.+1.0
    else:
        tel_mirror_divpm = tel_mirror_divpm *0.+1.0

    frame_num = 0
    total_num = 0
    embd = None
    # if "shift_frame" in aosat_cfg.CFG_SETTINGS:
    #     shft = aosat_cfg.CFG_SETTINGS["shift_frame"]
    # else:
    #     shft = None
    start_frame_index  = np.zeros(len(file_list)).astype(int)
    num_frames_in_file = np.zeros(len(file_list)).astype(int)
    logger.info("Scanning file headers...")
    for i in range(len(file_list)):
        hdr = pyfits.getheader(file_list[i])
        start_frame_index[i] = total_num

        if "embed_frame" in aosat_cfg.CFG_SETTINGS and embd is None:
            des_x_size = aosat_cfg.CFG_SETTINGS["embed_frame"][0]
            des_y_size = aosat_cfg.CFG_SETTINGS["embed_frame"][1]

            embd = ((int((des_x_size-hdr['NAXIS1'])/2),
                    int(des_x_size-hdr['NAXIS1']-int((des_x_size-hdr['NAXIS1'])/2))),
                    (int((des_y_size-hdr['NAXIS2'])/2),
                    int(des_y_size-hdr['NAXIS2']-int((des_y_size-hdr['NAXIS2'])/2))))


        if hdr['NAXIS'] == 2:
            total_num += 1
            num_frames_in_file[i] = 1
        else:
            if hdr['NAXIS'] == 3:
                total_num += hdr['NAXIS3']
                num_frames_in_file[i] = hdr['NAXIS3']
            else:
                raise ValueError("Wrong format of FITS file, should have 2 or 3 axes, with 3rd axis storing multiple phase screens if exists!")

    if 'totalnumber' in aosat_cfg.CFG_SETTINGS:
        totalnumber = aosat_cfg.CFG_SETTINGS['totalnumber']
        if totalnumber + startskip > total_num:
            logger.warning("Desired total number of %s is larger than available number of frames %s, adjusting!" % (totalnumber,total_num))
            totalnumber=total_num-startskip
    else:
        totalnumber = total_num-startskip
    logger.info("Found %s residual phase screens, serving screen %s to %s in steps of %s, total number will be %s!" % (total_num, startskip+1, startskip+totalnumber*skipstep,skipstep, totalnumber ))

    ##
    ## read and serve the frames
    ##
    this_frame_num = startskip
    next_step = 0
    frames_served = 0
    for i in range(len(file_list)):
        perc_done = frames_served / totalnumber
        if perc_done >= next_step:
            logger.debug("Serving from file %s, %s%% done! %s" % (os.path.basename(file_list[i]),perc_done*100,   '#'*int(30*perc_done)))
            next_step += 0.0333
        if this_frame_num >= start_frame_index[i] and this_frame_num < (start_frame_index[i] + num_frames_in_file[i]) and frames_served < totalnumber:
            logger.debug("Reading file %s, start index is %s, number of frames in file is %s!" % (file_list[i],start_frame_index[i],num_frames_in_file[i]))
            data = pyfits.getdata(file_list[i])
            this_index = this_frame_num - start_frame_index[i]
            logger.debug("Overall frame %s has index %s!" % (this_frame_num,this_index))
            while this_index < (num_frames_in_file[i]) and frames_served < totalnumber:
                logger.debug("Serving frame %s as %s!" % (this_index,frames_served))
                if num_frames_in_file[i] > 1:
                    frame = data[this_index] * pcf / tel_mirror_divpm
                else:
                    frame = data*pcf / tel_mirror_divpm
                if embd is not None:
                    logger.debug("Embedding frame in %s,%s array!"%embd)
                    frame = np.pad(frame,embd,'constant')

                this_frame_num += skipstep
                this_index = this_frame_num - start_frame_index[i]
                frames_served += 1
                yield(np.nan_to_num(frame),frames_served,totalnumber)

    return




if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
