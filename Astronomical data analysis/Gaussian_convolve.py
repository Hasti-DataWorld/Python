import numpy as np 
from matplotlib.pylab import *
from astropy import units as u
from astropy.io import fits 
from astropy.io import ascii
from astropy.table import Table, column
import os
import math
from astropy.convolution import convolve
from astropy.modeling import functional_models
from astropy.convolution import Gaussian1DKernel
from astropy import constants as consts
from scipy.interpolate import interp1d
from PyAstronomy import pyasl

"""
This code reads the spectrum in velocity space
and
Convolve the spectral feature with a Gaussian kernel
and
Plots the outcome

by HNateghi.
"""



# read velocity profile
synth_data = Table.read('.dat', format='ascii')
vlos = (synth_data['profile_vlos']).data
flux = (synth_data['profile_flux']).data

# make sure they have same pix size, width
vlos_inter = np.linspace(vlos[0], vlos[len(vlos)-1], 20000, endpoint=True)
flux_inter = np.interp(vlos_inter, vlos, flux)
disp = np.round(vlos_inter[300] - vlos_inter[299], 6)

x = np.arange(-30, 30, disp)
sig = 20  # this is in velocity

# make a gaussian emissin, normalized
gauss_ker = 1 / sig * np.sqrt(2 * np.pi) * np.exp(-x ** 2 / (sig**2) * 2.)
#plt.plot(x, gauss_ker)
result = convolve(flux_inter, gauss_ker, boundary='extend',normalize_kernel=True)
plt.plot(vlos_inter,flux_inter, 'b-', label="generated profile")
plt.plot(vlos_inter,result, 'm-', label="convolved profile, Gauss_sig ="+ str(sig)+"  kms")
plt.legend(loc=0)
plt.show()

#--- make a gaussian absorption, normalized

def emit2abs(emitspec, continuum=-1):
    """
    Takes an emission spectrum and "reverses" it, that is it:
    Makes a default flat continuum of 1 and subtracts exp(line)
    Inputs:
    emitspec (1D arraylike): the emission spectrum
    continuum (array or -1) : if -1, default of continuum of 1,
    otherwise must be same length as emitspec
    """
    if type(continuum) == int:
        if continuum == -1:
            continuum = np.ones(len(emitspec))
    abs_spec = continuum * np.exp(-1.0 * emitspec)
    return abs_spec
abs_spec = emit2abs(y, continuum=-1)

result_abs = convolve(flux_inter, abs_spec, boundary='extend',normalize_kernel=True)

plt.plot(vlos_inter,flux_inter, 'b-')
plt.plot(vlos_inter,result, 'm-')
plt.plot(vlos_inter,result_abs, 'k-')
plt.show()


#-----------------------------------------
# using a different function from PyAstronomy libarary to compare the outcome
# this function reads the wavelength, not working in velocity space

synth_data = Table.read('.dat', format='ascii')
vlos = (synth_data['profile_vlos']).data
flux = (synth_data['profile_flux']).data

zgal = 0.3815113
ion = 1215.67 
cval = consts.c.to('km/s').value

zlam = ion*(1.0+zgal)
wave = (vlos*zlam/cval)+zlam

wave_inter = np.linspace(wave[0], wave[len(wave)-1], 20000, endpoint=True)
flux_inter = np.interp(wave_inter, wave, flux)

sig = 1   # this is in wavelength
r = pyasl.broadGaussFast(wave_inter, flux_inter, sig, edgeHandling="firstlast")

plt.plot(wave_inter, flux_inter, 'k:', label="generated profile")
plt.plot(wave_inter, r, 'b-', label="convolved")
plt.legend(loc=4)
plt.show()


