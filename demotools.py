################################################################
# imports

import numpy as np
from scipy import interpolate

print('Generating interpolators...')

################################################################
# bookkeeping for Shankar full-aperture counts (230 GHz)

input_table = './tables/source_counts_nu=230GHz_Shankar.txt'

# read in the table
logtheta_r, logsigma_thresh, logN = np.loadtxt(input_table,unpack=True)

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N = 10.0**logN

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N = N.reshape((len(theta_r),len(sigma_thresh)))

# create an interpolation function
interp_SC_Sh_230 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N),kx=1,ky=1,s=0.0)
def source_counts_Sh_230(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC_Sh_230.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for KH full-aperture counts (230 GHz)

input_table = './tables/source_counts_nu=230GHz_KH.txt'

# read in the table
logtheta_r, logsigma_thresh, logN = np.loadtxt(input_table,unpack=True)

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N = 10.0**logN

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N = N.reshape((len(theta_r),len(sigma_thresh)))

# create an interpolation function
interp_SC_KH_230 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N),kx=1,ky=1,s=0.0)
def source_counts_KH_230(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC_KH_230.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for RV full-aperture counts (230 GHz)

input_table = './tables/source_counts_nu=230GHz_RV.txt'

# read in the table
logtheta_r, logsigma_thresh, logN = np.loadtxt(input_table,unpack=True)

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N = 10.0**logN

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N = N.reshape((len(theta_r),len(sigma_thresh)))

# create an interpolation function
interp_SC_RV_230 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N),kx=1,ky=1,s=0.0)
def source_counts_RV_230(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC_RV_230.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for Shankar baseline-projected counts (230 GHz)

input_table = './tables/source_counts_nu=230GHz_Shankar_baseline_projection.txt'

# read in the table
logtheta_r, logsigma_thresh, logN0, logN1, logN2 = np.loadtxt(input_table,usecols=(0,1,2,3,4),unpack=True)

# remove nans
logN0[~np.isfinite(logN0)] = -20.0
logN1[~np.isfinite(logN1)] = -20.0
logN2[~np.isfinite(logN2)] = -20.0

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N0 = 10.0**logN0
N1 = 10.0**logN1
N2 = 10.0**logN2

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N0 = N0.reshape((len(theta_r),len(sigma_thresh)))
N1 = N1.reshape((len(theta_r),len(sigma_thresh)))
N2 = N2.reshape((len(theta_r),len(sigma_thresh)))

# create interpolation function
interp_SC0_SH_230 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N0),kx=1,ky=1,s=0.0)
def source_counts_baseline_proj_Sh_230(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC0_SH_230.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for KH baseline-projected counts (230 GHz)

input_table = './tables/source_counts_nu=230GHz_KH_baseline_projection.txt'

# read in the table
logtheta_r, logsigma_thresh, logN0, logN1, logN2 = np.loadtxt(input_table,usecols=(0,1,2,3,4),unpack=True)

# remove nans
logN0[~np.isfinite(logN0)] = -20.0
logN1[~np.isfinite(logN1)] = -20.0
logN2[~np.isfinite(logN2)] = -20.0

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N0 = 10.0**logN0
N1 = 10.0**logN1
N2 = 10.0**logN2

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N0 = N0.reshape((len(theta_r),len(sigma_thresh)))
N1 = N1.reshape((len(theta_r),len(sigma_thresh)))
N2 = N2.reshape((len(theta_r),len(sigma_thresh)))

# create interpolation function
interp_SC0_KH_230 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N0),kx=1,ky=1,s=0.0)
def source_counts_baseline_proj_KH_230(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC0_KH_230.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for Shankar full-aperture counts (86 GHz)

input_table = './tables/source_counts_nu=86GHz_Shankar.txt'

# read in the table
logtheta_r, logsigma_thresh, logN = np.loadtxt(input_table,unpack=True)

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N = 10.0**logN

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N = N.reshape((len(theta_r),len(sigma_thresh)))

# create an interpolation function
interp_SC_Sh_86 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N),kx=1,ky=1,s=0.0)
def source_counts_Sh_86(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC_Sh_86.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for KH full-aperture counts (86 GHz)

input_table = './tables/source_counts_nu=86GHz_KH.txt'

# read in the table
logtheta_r, logsigma_thresh, logN = np.loadtxt(input_table,unpack=True)

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N = 10.0**logN

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N = N.reshape((len(theta_r),len(sigma_thresh)))

# create an interpolation function
interp_SC_KH_86 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N),kx=1,ky=1,s=0.0)
def source_counts_KH_86(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC_KH_86.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for RV full-aperture counts (86 GHz)

input_table = './tables/source_counts_nu=86GHz_RV.txt'

# read in the table
logtheta_r, logsigma_thresh, logN = np.loadtxt(input_table,unpack=True)

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N = 10.0**logN

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N = N.reshape((len(theta_r),len(sigma_thresh)))

# create an interpolation function
interp_SC_RV_86 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N),kx=1,ky=1,s=0.0)
def source_counts_RV_86(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC_RV_86.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for KH full-aperture counts (345 GHz)

input_table = './tables/source_counts_nu=345GHz_KH.txt'

# read in the table
logtheta_r, logsigma_thresh, logN = np.loadtxt(input_table,unpack=True)

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N = 10.0**logN

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N = N.reshape((len(theta_r),len(sigma_thresh)))

# create an interpolation function
interp_SC_KH_345 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N),kx=1,ky=1,s=0.0)
def source_counts_KH_345(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC_KH_345.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for RV full-aperture counts (345 GHz)

input_table = './tables/source_counts_nu=345GHz_RV.txt'

# read in the table
logtheta_r, logsigma_thresh, logN = np.loadtxt(input_table,unpack=True)

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N = 10.0**logN

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N = N.reshape((len(theta_r),len(sigma_thresh)))

# create an interpolation function
interp_SC_RV_345 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N),kx=1,ky=1,s=0.0)
def source_counts_RV_345(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC_RV_345.ev(log_theta_r,log_sigma_thresh)

################################################################
# bookkeeping for Shankar full-aperture counts (345 GHz)

input_table = './tables/source_counts_nu=345GHz_Shankar.txt'

# read in the table
logtheta_r, logsigma_thresh, logN = np.loadtxt(input_table,unpack=True)

# convert to linear scale
theta_r = 10.0**logtheta_r
sigma_thresh = 10.0**logsigma_thresh
N = 10.0**logN

# collapse to original dimensions
theta_r = np.unique(theta_r)
sigma_thresh = np.unique(sigma_thresh)
N = N.reshape((len(theta_r),len(sigma_thresh)))

# create an interpolation function
interp_SC_Sh_345 = interpolate.RectBivariateSpline(np.log10(theta_r),np.log10(sigma_thresh),np.log10(N),kx=1,ky=1,s=0.0)
def source_counts_Sh_345(log_theta_r,log_sigma_thresh):
    return 10.0**interp_SC_Sh_345.ev(log_theta_r,log_sigma_thresh)

################################################################
# define master function

print('Generated!')

def source_counts(theta_r,sigma_thresh,freq='230',counts_type='full',BHMF_type='KH'):
    """ A function to compute SMBH shadow counts, based on
        computations from Pesce et al. (2021)
       
        Inputs:
            theta_r: angular resolution, in uas
            sigma_thresh: sensitivity, in Jy
            freq: either '86', '230', or '345'
            counts_type: either 'baseline' or 'full'
            BHMF_type: either 'Shankar', 'KH', or 'RV'
        
        Output: source counts at the specified resolution
                and sensitivity
    """

    # take logs
    log_theta_r = np.log10(theta_r)
    log_sigma_thresh = np.log10(sigma_thresh)

    # do some basic bounds checking
    if np.any((log_theta_r > 2.0) | (log_theta_r < -2.0)):
        print('WARNING: angular resolution is outside of computed range')
    if np.any((log_sigma_thresh > 1.0) | (log_sigma_thresh < -9.0)):
        print('WARNING: sensitivity is outside of computed range')

    # select the appropriate interpolation function
    if (freq == '230'):
        if (counts_type == 'baseline'):
            if (BHMF_type == 'Shankar'):
                out = source_counts_baseline_proj_Sh_230(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'KH'):
                out = source_counts_baseline_proj_KH_230(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'RV'):
                out = source_counts_baseline_proj_RV_230(log_theta_r,log_sigma_thresh)
        if (counts_type == 'full'):
            if (BHMF_type == 'Shankar'):
                out = source_counts_Sh_230(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'KH'):
                out = source_counts_KH_230(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'RV'):
                out = source_counts_RV_230(log_theta_r,log_sigma_thresh)
    elif (freq == '86'):
        if (counts_type == 'baseline'):
            if (BHMF_type == 'Shankar'):
                out = source_counts_baseline_proj_Sh_86(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'KH'):
                out = source_counts_baseline_proj_KH_86(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'RV'):
                out = source_counts_baseline_proj_RV_86(log_theta_r,log_sigma_thresh)
        if (counts_type == 'full'):
            if (BHMF_type == 'Shankar'):
                out = source_counts_Sh_86(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'KH'):
                out = source_counts_KH_86(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'RV'):
                out = source_counts_RV_86(log_theta_r,log_sigma_thresh)
    elif (freq == '345'):
        if (counts_type == 'baseline'):
            if (BHMF_type == 'Shankar'):
                out = source_counts_baseline_proj_Sh_345(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'KH'):
                out = source_counts_baseline_proj_KH_345(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'RV'):
                out = source_counts_baseline_proj_RV_345(log_theta_r,log_sigma_thresh)
        if (counts_type == 'full'):
            if (BHMF_type == 'Shankar'):
                out = source_counts_Sh_345(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'KH'):
                out = source_counts_KH_345(log_theta_r,log_sigma_thresh)
            if (BHMF_type == 'RV'):
                out = source_counts_RV_345(log_theta_r,log_sigma_thresh)
            
    return out

