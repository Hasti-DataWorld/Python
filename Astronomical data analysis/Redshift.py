from astropy.io import fits 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker
from scipy.integrate import quad
from astropy import units as u
from astropy.io import ascii
from astropy.table import Table, column
from astropy import constants as consts
import os

# this function calculates the redshift of the object using emission or absorption line
# by HNateghi-2021

def redshift(specTable, roughZ, ion_wave,ion_name):  

    data = Table.read(specTable)
    # science spectrum
    wave = data['wave'] 
    flux = data['flux']
    err = data['err'] 

    emiss_line = ion_wave*(1+roughZ)

    fig, ax = plt.subplots()
    ax.step(wave, flux, 'k')
    ax.axvline(x=emiss_line, color='m')
    ax.set_xlabel('Angstrom')
    ax.set_ylabel("Flux")

    plt.show()

    guess_h = input ("Enter init_guess_h:")
    guess_c = input ("Enter init_guess_c:")
    guess_sig = input ("Enter init_guess_sig:")
    

    in_h = float(guess_h)
    in_c = float(guess_c)
    in_sig = float(guess_sig)


    cont_emiss = [[emiss_line - 90, emiss_line - 40],[emiss_line + 40, emiss_line + 90]]

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

    popt, pcov = curve_fit(gauss, wave_emiss, flux_norm-1, p0=[in_h, in_c, in_sig] , sigma= err_norm, absolute_sigma=True)
    # note that we shouldn't subtract 1 from error, subtraction doesn't matter in error propogation
    amplitude, emiss_center, sigma = popt[0], popt[1], popt[2]

    sderr = np.sqrt(np.diag(pcov))
    ERRamplitude, ERRemiss_center, ERRsigma = sderr[0], sderr[1], sderr[2]  

    z = (popt[1] - ion_wave)/ion_wave
    zerr = 1/ion_wave * ERRemiss_center

    #cval = consts.c.to('km/s').value
    #vel = (emiss_center-(ion_wave*(1+z)))/(ion_wave*(1+z))*cval
    #vel_err= ERRemiss_center/(ion_wave*(1+z))*cval

    print('z = '+ion_name+str(z))
    print('zerr = '+ion_name+str(zerr))

    plt.step( wave_emiss, flux_norm, c='k', where='mid' )
    plt.step( wave_emiss, err_norm, c='k', where='mid', alpha=0.3 )
    plt.plot(wave_emiss, (gauss(wave_emiss, amplitude, emiss_center, sigma)+1), 'r-')
    #plt.plot(wave_emiss, (gauss(wave_emiss, amplitude, emiss_center, sigma)+1+err_norm), 'c')
    #plt.plot(wave_emiss, (gauss(wave_emiss, amplitude, emiss_center, sigma)+1-err_norm), 'g')
    plt.show()

    return

