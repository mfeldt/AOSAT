import numpy as np

#----------------------------------------------------------
#Tue Jan  8 21:54:21 2019
parfile                  = 'metis_370P_35L.par'
sim_name                 = 'ExampleClosedLoop'
seeing                   =   0.43
starmag                  =   2.00
optical_throughput       = 0.2243
photons_pi_star          = 59006168.0 # per subap in SHS, total in pyr
photons_pi_sky           = 6554.3 # per subap in SHS, total in pyr
regpar                   = 0.1000
gain                     = 0.4000
calamp                   =    0.2500
threshold                =  1.000
loopfreq                 = 1000.0
zenithdist               =  0
pixscale                 = 0.7261 # input value, irrelevant for pyramid
fssize                   = 1.8000
resultfile               = '/data/feldt/yao_sims/NewBase//amp_v_gain_pyr74_mag=02.00_seeing=0.43_gain=0.40_calamp=000.25_regpar=0.100_threshold=1_loopfreq=1000_zenithdist=00/metis_370P_35L.res'
nact                     = '3982.0'
nsub                     = '4160.0'
lambdae                  = np.zeros([ 4])
strehl                   = np.zeros([ 4, 1]) # nlambda, ntarget
fwhm                     = np.zeros([ 4, 1]) # nlambda, ntarget
e50                      = np.zeros([ 4, 1]) # nlambda, ntarget
xpos                     = np.zeros([ 1]) # ntarget
ypos                     = np.zeros([ 1]) #  ntarget
lambdae[ 0]             =  10.00
xpos[ 0]                =   0.00
ypos[ 0]                =   0.00
strehl[ 0, 0]          =  0.998
fwhm[ 0, 0]            = 55.965
e50[ 0, 0]             = 95.729
lambdae[ 1]             =   3.70
xpos[ 0]                =   0.00
ypos[ 0]                =   0.00
strehl[ 1, 0]          =  0.984
fwhm[ 1, 0]            = 20.493
e50[ 1, 0]             = 35.420
lambdae[ 2]             =   3.00
xpos[ 0]                =   0.00
ypos[ 0]                =   0.00
strehl[ 2, 0]          =  0.976
fwhm[ 2, 0]            = 16.528
e50[ 2, 0]             = 28.719
lambdae[ 3]             =   2.20
xpos[ 0]                =   0.00
ypos[ 0]                =   0.00
strehl[ 3, 0]          =  0.957
fwhm[ 3, 0]            = 12.120
e50[ 3, 0]             = 22.169
