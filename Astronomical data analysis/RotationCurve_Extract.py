from astropy.io import fits 
import numpy as np 
from matplotlib.pylab import *
import pandas as pd
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker
from scipy.integrate import quad
from astropy import units as u
from astropy.io import ascii
from astropy.table import Table, column
from astropy import constants as consts
import os

# use these functions to extract rotation curve of a galaxy
# by HNateghi-2021 


### creat an empty table to append rotation curves:    CHANGE THE NAME OF ION !!!!!!!
def table(ion_name):
    tablev = Table({ 'slit_position': [], 'center_velocity': [], 'err_velocity': []}, names=('slit_position', 'center_velocity', 'err_velocity'))
    tablev.write(ion_name+'_rotv.dat', format='ascii', overwrite=False)

    tablew = Table({ 'slit_position': [], 'center_wave': [], 'err_wave': [], 'z': [] ,'err_z': []} ,names=('slit_position', 'center_wave', 'err_wave', 'z', 'err_z'))
    tablew.write(ion_name+'_rotw.dat', format='ascii', overwrite=False)
    return

########## plot the spec for initial guess
def spec_plot(inputspec, z, ion_wave):
    plt.close()
    with fits.open(inputspec) as hdul:
            hdr = hdul[0].header 
            sci = hdul[0].data     # science spectrum
            disper = hdr['CD1_1'] 
            startwave = hdr['CRVAL1']
            pixnum = hdr['NAXIS1'] 
            wave_range = np.arange(0,pixnum)
            wave_delt = wave_range*disper
            wave = startwave + wave_delt
            # since the spape of line*_emiss.fits spectra are (4, 1, 2042):
            flux = sci[1,0,:] # we take the second extension for non-weight flux, change to error correspond to emission. if we read whole spec sci[2,7,:]
            err = sci[3,0,:]
    
    emiss_line = ion_wave*(1+z)

   
    fig, ax = plt.subplots()
    ax.step(wave[0:flux.shape[0]], flux, 'k')
    ax.axvline(x=emiss_line, color='m')
    ax.set_xlabel('Angstrom')
    ax.set_ylabel("Flux")

    show()
    return

    ######## calculate the center of emission line in each aperture spectrum to get the rotation curve

def emiss_fit(inputspec, shift, z, ion_wave, ion_name, init_guess):     # init_guess = [h, c, sig]
    with fits.open(inputspec) as hdul:
            hdr = hdul[0].header 
            sci = hdul[0].data     # science spectrum
            disper = hdr['CD1_1'] 
            startwave = hdr['CRVAL1']
            pixnum = hdr['NAXIS1'] 
            wave_range = np.arange(0,pixnum)
            wave_delt = wave_range*disper
            wavelength = startwave + wave_delt
            # since the spape of line*_emiss.fits spectra are (4, 1, 2042):
            flux = sci[1,0,:] # we take the second extension for non-weight flux, change to error correspond to emission. if we read whole spec sci[2,7,:]
            err = sci[3,0,:]
            wave = wavelength[0:flux.shape[0]]

    emiss_line = ion_wave*(1+z)
    cont_emiss = [[emiss_line - 70, emiss_line - 30],[emiss_line + 30, emiss_line + 70]]

    indcont_emiss = ((wave > cont_emiss[0][0]) & (wave < cont_emiss[0][1])) |((wave > cont_emiss[1][0]) & (wave < cont_emiss[1][1]))
    indrange_emiss = (wave > cont_emiss[0][0]) & (wave < cont_emiss[1][1])
    flux_norm = np.zeros((flux.shape[0]))
    err_norm = np.zeros((err.shape[0]))

    for i in range(flux.shape[0]):
        linecoeff = np.polyfit(wave[indcont_emiss], flux[indcont_emiss],2)
        flux_norm = flux[indrange_emiss]/np.polyval(linecoeff, wave[indrange_emiss])
        err_norm  = err[indrange_emiss]/np.polyval(linecoeff, wave[indrange_emiss])

    wave_emiss = wave[indrange_emiss]
    conti = np.polyval(linecoeff, wave[indrange_emiss])

    def gauss(x, h,c,sig):
        gauss_func = h * np.exp( - (x - c)**2.0 / (2.0 * sig**2.0) )
        return gauss_func

    popt, pcov = curve_fit(gauss, wave_emiss, flux_norm-1, p0=init_guess , sigma= err_norm-1, absolute_sigma=True)

    amplitude, emiss_center, sigma = popt[0], popt[1], popt[2]

    sderr = np.sqrt(np.diag(pcov))
    ERRamplitude, ERRemiss_center, ERRsigma = sderr[0], sderr[1], sderr[2]  
        
    cval = consts.c.to('km/s').value
    vel = (emiss_center-(ion_wave*(1+z)))/(ion_wave*(1+z))*cval
    vel_err= ERRemiss_center/(ion_wave*(1+z))*cval

    z = (popt[1] - ion_wave)/ion_wave
    zerr = 1/ion_wave * ERRemiss_center

    table = Table({ 'slit_position': [shift], 'center_velocity': [vel], 'err_velocity': [vel_err]} ,names=('slit_position', 'center_velocity', 'err_velocity'))
    with open(ion_name+'_rotv.dat', mode='a') as f:
    #    f.seek(0, os.SEEK_END)  # Some platforms don't automatically seek to end when files opened in append mode
        table.write(f, format='ascii.no_header') 

    table = Table({ 'slit_position': [shift], 'center_wave': [emiss_center], 'err_wave': [ERRemiss_center], 'z': [z] ,'err_z': [zerr]} ,names=('slit_position', 'center_wave', 'err_wave', 'z', 'err_z'))
    with open(ion_name+'_rotw.dat', mode='a') as f:
    #    f.seek(0, os.SEEK_END)  # Some platforms don't automatically seek to end when files opened in append mode
        table.write(f, format='ascii.no_header')     


    #x_fit  = np.linspace(wave_emiss[0],wave_emiss[-1], 1000)
    plt.step( wave_emiss, flux_norm, c='k', where='mid' )
    plt.step( wave_emiss, err_norm, c='k', where='mid', alpha=0.3 )
    plt.plot(wave_emiss, (gauss(x_fit, amplitude, emiss_center, sigma)+1), 'r-')
    show()

    return










