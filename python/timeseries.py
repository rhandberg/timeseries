#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 02 13:05:17 2015

@author: f041080
"""

from __future__ import division, with_statement, print_function
import os
import numpy as np
import scipy.stats as st
from bottleneck import nanmedian, nanmin, nanmax
import scipy.interpolate as npint
from copy import deepcopy as dc
import subprocess
import tempfile
from scipy.integrate import trapz


def gapfill(time, flux):

	maxgap = 1/24

	dt = nanmedian(np.diff(time))

    # Do linear interpolation of the points:
    indx2 = np.isfinite(flux)

	df = np.diff(time[indx2])
	indx = np.where( (df > 1.5*dt) & (df < maxgap) )[0]

    f = npint.interp1d(time[indx2], flux[indx2], kind='linear')

    flux_gapfill = dc(flux)
    flux_gapfill[indx] = f(time[indx])
    return flux_gapfill


#------------------------------------------------------------------------------
def gapfill(time, flux, maxgap=1.0):
    # Create list of points to be interpolated:
    indx = [ i for i,f in enumerate(flux[1:-1], start=1) if not np.isfinite(f) and np.isfinite(flux[i+1]) and np.isfinite(flux[i-1])]
    #indx = [ i for i,f in enumerate(flux) if not np.isfinite(f) and np.isfinite(flux[i+1]) and np.isfinite(flux[i-1])]

    indx = []
    for k in xrange(len(time)):
        if not np.isfinite(flux[k]):


    # Do linear interpolation of the points:
    indx2 = np.isfinite(flux)
    f = npint.interp1d(time[indx2], flux[indx2], kind='linear')

    flux_gapfill = dc(flux)
    flux_gapfill[indx] = f(time[indx])
    return flux_gapfill

#------------------------------------------------------------------------------
def window(time, flux, sigma, width=50, ofac=1, use_weights=True):

    nyquist = 0.5e6/(86400*np.median(np.diff(time)))
    dnu = 1e6/(86400.0*( nanmax(time)-nanmin(time) ))

    # Commstring to send to WINDOW program:
    if sigma is None or np.allclose(sigma, 1):
        commstring_window = "%.16e\n%.16e\n%.16e\n"%(nyquist/2, width, dnu/ofac)
    else:
        commstring_window = "%d\n%.16e\n%.16e\n%.16e\n"%(int(use_weights), nyquist/2, width, dnu/ofac)

    # Calculate highly oversampled windowfunction and perform integral:
    # Get temp file for input:
    fd_in, input_timeseries = tempfile.mkstemp(text=True)
    fd_out, output_window = tempfile.mkstemp(text=True)
    try:
        # Write file with input timeseries:
        with open(input_timeseries, 'w') as fid:
            for k in xrange(len(time)):
                fid.write("%.16e  %.16e  %.16e\n"%(time[k], flux[k], sigma[k]))

        # Oversample the windowfunction by a factor 10:
        cmd = r'window -tday -quiet "%s" "%s"'%(input_timeseries, output_window)
        output = subprocess.Popen(cmd, cwd='.', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        ret = output.communicate(commstring_window)
        if not ret[1] == '':
            raise Exception('WINDOW failed with error: %s' % ret[1])

        # Load oversampled windowfunction:
        window_nu, window_P = np.loadtxt(output_window, usecols=(0,1), unpack=True)
    except:
        raise
    finally:
        # Delete input file again:
        os.close(fd_in)
        os.remove(input_timeseries)
        os.close(fd_out)
        os.remove(output_window)

    return window_nu, window_P

#------------------------------------------------------------------------------
def wps(file_tmseries, use_weights=True):
    """Weighted power spectrum"""
    # If needed, load timeseries and find nyquist and mean spacing:
    Time, Flux, FluxErr = np.loadtxt(file_tmseries, usecols=(0,1,2), comments='#', unpack=True)
    difft = np.diff(Time*86400.0)
    nyquist = 0.5e6/nanmedian(difft)
    dnu = 1e6/(86400.0*( nanmax(Time)-nanmin(Time) ))

    if nyquist < 300:
        width = 50
    else:
        width = 300

    invalid_weights = np.allclose(FluxErr, 1)

    #print(nyquist, dnu, width)

    # Commstring to send to WINDOW program:
    if invalid_weights:
        commstring_window = "%.16e\n%.16e\n%.16e\n"%(nyquist/2, width, dnu/10)
    else:
        commstring_window = "%d\n%.16e\n%.16e\n%.16e\n"%(int(use_weights), nyquist/2, width, dnu/10)

    # Calculate highly oversampled windowfunction and perform integral:
    # Get temp file for input:
    fd_in, file_window_high = tempfile.mkstemp(text=True)
    try:
        # Oversample the windowfunction by a factor 10:
        cmd = r'window -tday -quiet "%s" "%s"'%(file_tmseries, file_window_high)
        output = subprocess.Popen(cmd, cwd='.', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        ret = output.communicate(commstring_window)
        if not ret[1] == '':
            raise Exception('WINDOW failed with error: %s' % ret[1])

        # Load oversampled windowfunction:
        window_nu, window_P = np.loadtxt(file_window_high, usecols=(0,1), unpack=True)
    except:
        raise
    finally:
        # Delete input file again:
        os.close(fd_in)
        os.remove(file_window_high)

    # Load oversampled windowfunction:
    Pint = trapz(window_P, window_nu)

    if invalid_weights:
        commstring_spec = "%.16e\n%.16e\n%.16e\n%.16e\n"%(0, nyquist, Pint, Pint)
    else:
        commstring_spec = "%d\n%.16e\n%.16e\n%.16e\n%.16e\n"%(int(use_weights), 0, nyquist, Pint, Pint)

    # Calculate powerspectrum:
    # Get temp file for output:
    fd_in, file_spec_out = tempfile.mkstemp(text=True)
    try:
        # Call spec-program:
        # Run at fundamental sampling
        cmd = r'spec -tday -pd -quiet "%s" "%s"'%(file_tmseries, file_spec_out)
        output = subprocess.Popen(cmd, cwd='.', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        ret = output.communicate(commstring_spec)
        if not ret[1] == '':
            raise Exception('SPEC failed with error: %s' % ret[1])

        # Load the results:
        nu, P = np.loadtxt(file_spec_out, usecols=(0,1), unpack=True)
    except:
        raise
    finally:
        os.close(fd_in)
        os.remove(file_spec_out)

    # Covert to the same scaling as the unweighted spectrum (rms-scaling):
    P /= 2.0

    return nu, P

#------------------------------------------------------------------------------
def spec(time, flux, sigma, nuspan=None, ofac=1, use_weights=True):

    if nuspan is None:
        numin = 0
        numax = 0.5e6/(86400*np.median(np.diff(time)))
    else:
        numin = nuspan[0]
        numax = nuspan[1]

    dnu = 1e6/(86400*(np.max(time) - np.min(time)))
    dnu /= ofac

    if sigma is None or np.allclose(sigma, 1):
        commstring = "%.16e\n%.16e\n%.16e\n"%(numin, numax, dnu)
    else:
        commstring = "%d\n%.16e\n%.16e\n%.16e\n"%(int(use_weights), numin, numax, dnu)

    fd_in, input_timeseries = tempfile.mkstemp(text=True)
    fd_out, output_powerspectrum = tempfile.mkstemp(text=True)
    try:
        # Write file with input timeseries:
        with open(input_timeseries, 'w') as fid:
            for k in xrange(len(time)):
                fid.write("%.16e  %.16e  %.16e\n"%(time[k], flux[k], sigma[k]))

        # Call spec-program:
        cmd = r'spec -tday -p -quiet "%s" "%s"'%(input_timeseries, output_powerspectrum)
        output = subprocess.Popen(cmd, cwd='.', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        ret = output.communicate(commstring)
        if not ret[1] == '':
            raise Exception('SPEC failed with error: %s' % ret[1])

        # Load spectrum:
        nu, P = np.loadtxt(output_powerspectrum, unpack=True, usecols=(0,1))

    except:
        raise

    finally:
        os.close(fd_in)
        os.remove(input_timeseries)
        os.close(fd_out)
        os.remove(output_powerspectrum)

    return nu, P

#------------------------------------------------------------------------------
def clean(time, flux, sigma, npeaks=1, nuspan=None, use_weights=True):

    if nuspan is None:
        numin = 0
        numax = 0.5e6/(86400*np.median(np.diff(time)))
    else:
        numin = nuspan[0]
        numax = nuspan[1]

    if sigma is None or np.allclose(sigma, 1):
        commstring = "%.16e\n%.16e\n%d\n"%(numin, numax, npeaks)
    else:
        commstring = "%d\n%.16e\n%.16e\n%d\n"%(int(use_weights), numin, numax, npeaks)

    fd_in, input_timeseries = tempfile.mkstemp(text=True)
    fd_out, output_timeseries = tempfile.mkstemp(text=True)
    try:
        # Write file with input timeseries:
        with open(input_timeseries, 'w') as fid:
            for k in xrange(len(time)):
                fid.write("%.16e  %.16e  %.16e\n"%(time[k], flux[k], sigma[k]))

        # Call clean-program:
        cmd = r'clean -tday "%s" "%s"'%(input_timeseries, output_timeseries)
        output = subprocess.Popen(cmd, cwd='.', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        ret = output.communicate(commstring)
        if not ret[1] == '':
            raise Exception('CLEAN failed with error: %s' % ret[1])

        # Load final timeseries:
        time, flux, sigma = np.loadtxt(output_timeseries, unpack=True)
    except:
        raise

    finally:
        os.close(fd_in)
        os.remove(input_timeseries)
        os.close(fd_out)
        os.remove(output_timeseries)

    return time, flux, sigma

#------------------------------------------------------------------------------
def cdpp_powerspectrum(nu, P):

    # Load SG filter:
    nu_gs, P_gs = np.loadtxt('tranfunc.txt', usecols=(1,2), unpack=True)
    P_gs = kf.smooth(P_gs, 21)
    pint = interp1d(nu_gs, P_gs, kind='linear', bounds_error=False, fill_value=0)

    # Caluclate CDPP:
    Tav = 6.5*60*60
    dnu = nu[1] - nu[0]
    intp = P * np.sinc(1e-6*nu*Tav)**2 * pint(nu)
    cdpp = np.sqrt(dnu*np.sum(intp))

    return cdpp

#------------------------------------------------------------------------------
def cleanpeaks(time, flux, sigma, alpha=0.10, nmax=30, use_weights=True):

    # Calculate nyquist frequency for this timeseries:
    nyquist = 0.5e6/(86400*np.median(np.diff(time)))
    dof = 2

    # Loop through the harmonics to try to remove:
    for harmonic in (1, 2, 3, 4, 5): #0.5,
        # Set frequencies to CLEAN between:
        nucen = harmonic * 47.23394501
        numin = nucen - 1.0
        numax = nucen + 1.0

        # If outside of the spectrum, skip it:
        if (nucen > nyquist): break

        # Width of the interval we are cleaning:
        dspan = numax-numin

        # Start cleaning peaks in this interval:
        n = 1
        while (n <= nmax):
            # Calculate power spectrum in small segment around interval
            # and estimate the background level and set the threshold
            # that peaks should be above to be significant:
            nu_test, P_test = spec(time, flux, sigma, nuspan=(numin-dspan, numax+dspan), ofac=1, use_weights=use_weights)
            background = 1.4238281250000002 * np.median(P_test)

            # Calculate power spectrum in interval to see if there are any
            # significant peaks, and if not stop:
            nu_test, P_test = spec(time, flux, sigma, nuspan=(numin, numax), ofac=10, use_weights=use_weights)

            #imax = np.argmax(P_test)
            Chi2 = (P_test.max()/background)*dof
            p_val = st.chi2.sf(Chi2, dof)
            if not p_val < alpha: break

            # Clean the highest peak away:
            time, flux, sigma = clean(time, flux, sigma, npeaks=1, nuspan=(numin, numax), use_weights=use_weights)
            #time, flux, sigma = clean(time, flux, sigma, npeaks=1, nuspan=(nu_test[imax-10], nu_test[imax+10]), use_weights=use_weights)

            n += 1

    # Remove the clean.log that was created:
    try:
        os.remove('clean.log')
    except:
        pass

    return time, flux, sigma
