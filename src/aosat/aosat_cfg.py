"""
These are a number of default parameters to be used of nothing can be found elsewhere

"""

import collections
import os


fpars = collections.OrderedDict([('an_lambda',     3e-6),  # analysis wavelength in metres
                                 ('ppm',           10.0),  # pixel per metre in pupil plane
                                 ('pup_diam',      37.0),  # pupil diameter in metres
                                 ('L0',            25.0),  # outer scale in metres
                                 ('crad',          700.0),  # AO control radius [mas]
                                 ('zterms',        20),    # number of Zernike terms for decomposition
                                 ('aosat_logfile',   'aosat.log'),  # log file
                                 ('aosat_loglevel',  'DEBUG'),      # log level
                                 ('aosat_fft_forcenumpy',  True),   # force usage of numpy FFTs
                                 ('setup_path','.'),                # path to setup file
                                 ('pupilmask','examples/ExampleAnalyze/yao_pupil.fits'),
        ])


LOG_SETTINGS = {
    'version': 1,
    'root': {
        'level': 'NOTSET',
        'handlers': ['console', 'file'],
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'detailed',
            'stream': 'ext://sys.stdout',
        },
        'file': {
            'class': 'logging.handlers.RotatingFileHandler',
            'formatter': 'detailed',
            'filename': 'logs/MyProject.log',
            'mode': 'a',
            'maxBytes': 10485760,
            'backupCount': 5,
        },
    },
    'formatters': {
        'detailed': {
            'format': '%(asctime)s %(name)s %(module)-17s %(funcName)s line:%(lineno)-4d ' \
            '%(levelname)-8s %(message)s',
        },
        'normal': {
            'format': '%(asctime)s %(name)s ' \
            '%(levelname)-8s %(message)s',
        },
    },
    'disable_existing_loggers':False
}


def configure(setupfile):
    """Set up a configuration 'repdict' from default or coniguration file.

    The general format of configuration files is

    name = value

    pairs. If a # appears in a line, the remainder is ignored.
    Also [ and ] are ignored. This is historical as we have some
    auto-generated report files from simulation software that
    contains array constructors.

    The default parameters are always present and need to be overwritten
    if you want to change them. So make sure that your configuration file
    always takes care of

    an_lambda     = 3e-6,  # analysis wavelength in metres
    ppm           = 10.0,  # pixel per metre in pupil plane
    pup_diam      = 37.0,  # pupil diameter in metres
    L0            = 25.0,  # outer scale in metres
    crad          = 700.0,  # AO control radius [mas]
    zterms        = 200,     # number of Zernike terms for decomposition
    aosat_logfile  = 'aosat.log', # log file
    aosat_loglevel = 'DEBUG'      # log level

    If you don't like them!

    There can be parameters that are never used. A report can be generated
    by calling 'repstring(repdict)'

    Parameters
    ----------
    setupfile : <class 'str'> or <class 'NoneType'>
        Name of coniguration file to parse. If None,
        the default configuration repdict is returned
    Returns
    -------
    collections.OrderedDict
        Ordered dictionary of the AOSAT configuration


    Examples
    -------

    The most simple example is a call like


    >>> configure(None)
    OrderedDict([('...


    When reading a configuration file, these files should have the format of
    a series of simple assignments, one per line, e.g:

    parfile                  = 'metis_370P_35L.par'
    sim_name                 = 'amp_v_gain_pyr74_mag=02.00_seeing=0.43'
    seeing                   =   0.43
    starmag                  =   2.00
    optical_throughput       = 0.2243
    aosat_logfile            = None
    aosat_loglevel           = 'INFO'

    When reading this file from the examples diretory, the resulting repdict will change
    to contain the new loglevel 'INFO':

    >>> import os
    >>> repdict = configure(os.path.join(os.path.dirname(os.path.abspath(__file__)),'examples','example_report.py'))
    >>> repdict
    OrderedDict([('...
    >>> repdict['aosat_loglevel']
    'INFO'

    """
    repdict=collections.OrderedDict(fpars)
    if setupfile is not None:
        path_to_setup = os.path.dirname(os.path.abspath(setupfile))
        repdict['setup_path'] = path_to_setup
        lines = open(setupfile,'r').readlines()
        for line in lines:
            newline=line.split('#')[0]
            if len(newline)>0:
                if (not 'zeros' in newline) and ('=' in newline):
                    newline = "repdict['"+newline.split('=')[0].replace(" ","")+"']="+"=".join(newline.split('=')[1:])#
                    exec(newline)
    return(repdict)


CFG_SETTINGS = configure(None)

def repString(repdict,prefix = ''):
    """returns a summary string of the repdict

    Parameters
    ----------
    repdict : collections.OrderedDict
        AOSAT configuration ordered dictionary
    Returns
    -------
    type
        Representative string for use e.g. in plots

    Examples
    -------
    >>> import os
    >>> repdict = configure(os.path.join(os.path.dirname(os.path.abspath(__file__)),'examples','example_report.py'))
    >>> print(repString(repdict))
               an_lambda : 3e-06...


    """
    rep_string=''
    maxklen    = max(map(len, repdict.keys()))
    for key,val in repdict.items():
            newline = prefix + '{:>{width}}'.format(key,width=maxklen)+' : '+ ("%s" % val)+"\n"
            if len(newline)>50:
                newline = newline[:50]+"...\n"
            rep_string += newline
    return(rep_string)



def configureLogging(repdict):
    """Reconfigures the log settings.

    Need to call
    logging.config.dictConfig(aosat_cfg.LOG_SETTINGS)
    after.

    Not intended for explicit use.


    Parameters
    ----------
    repdict : collections.OrderedDict
        AOSAT configuration ordered dictionary

    Returns
    -------
    nothing

    """

    if repdict['aosat_logfile'] is None:
        LOG_SETTINGS['handlers'].pop('file',None) # remove the file handlers
        LOG_SETTINGS['root']['handlers'] = ['console']
        if repdict['aosat_loglevel'] == 'DEBUG':
            LOG_SETTINGS['handlers']['console']['formatter']='detailed'
        else:
            LOG_SETTINGS['handlers']['console']['formatter']='normal'
    else:
        LOG_SETTINGS['handlers']['file'] = {
            'class': 'logging.handlers.RotatingFileHandler',
            'mode': 'a',
            'maxBytes': 10485760,
            'backupCount': 5,
        }
        LOG_SETTINGS['handlers']['file']['filename'] = repdict['aosat_logfile']
        if repdict['aosat_loglevel'] == 'DEBUG':
            LOG_SETTINGS['handlers']['file']['formatter']='detailed'
        else:
            LOG_SETTINGS['handlers']['file']['formatter']='normal'
    LOG_SETTINGS['root']['level']=repdict['aosat_loglevel']



if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
