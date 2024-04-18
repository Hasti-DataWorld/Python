"""
These functions help to grab the chunk of spectra that you need with the
absorption or emission line and then normalize the continnum to
proceed with the analysis

by HNateghi-2022
"""


from astropy.io import fits 
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker
from scipy.integrate import quad
from astropy import units as u
from astropy.io import ascii
from astropy.table import Table, column
import os

plt.close()

def close_figure(event): # close the plot window with 'esc' key
    if event.key == 'escape':
        plt.close(event.canvas.figure)

def norm_ion(inputspec, z, ion_name, ion, pol_order):  # normalize the chunk of spec including ion of interest

    wave, flux, err = np.genfromtxt(inputspec, skip_header=1, usecols=(0,1,2), unpack=True)

    ion_abs = ion * (1+z)

    fig, ax1 = plt.subplots()

    ax1.step(wave, flux, 'b')
    ax1.axvline(ion_abs, c='r')
    ax1.set_ylim(np.min(flux),5e-14)
    plt.gcf().canvas.mpl_connect('key_press_event', close_figure)
    plt.show()


    min_boundleft = input ("Enter the minimum of the left regoin:")
    max_boundleft = input ("Enter the maximum of the left regoin:")
    min_boundright = input ("Enter the minimum of the right regoin:")
    max_boundright = input ("Enter the maximum of the right regoin:")
    
    Lmin = float(min_boundleft)
    Lmax = float(max_boundleft)
    Rmin = float(min_boundright)
    Rmax = float(max_boundright)
    
      
    left_range = input ("Enter the min wave of the spec chunk:")
    right_range = input ("Enter the max wave of the spec chunk:")
    range_min = float(left_range)
    range_max = float(right_range)
    range_ion = [range_min, range_max] 
    indrange_ion = (wave > range_min) & (wave < range_max)
    
    cont_ion = [[Lmin, Lmax],[Rmin, Rmax]]
    
    indcont_ion = ((wave > cont_ion[0][0]) & (wave < cont_ion[0][1])) |((wave > cont_ion[1][0]) & (wave < cont_ion[1][1]))

    linecoeff = np.polyfit(wave[indcont_ion], flux[indcont_ion], pol_order)
    flux_norm = flux[indrange_ion]/np.polyval(linecoeff, wave[indrange_ion])
    err_norm  = err[indrange_ion]/np.polyval(linecoeff, wave[indrange_ion])


    dataout = Table([wave[indrange_ion], flux_norm, err_norm], names=['lambda', 'flux_norm', 'err_norm'])
    ascii.write(dataout, ion_name+'_z'+ str(z) +'_norm_'+inputspec, overwrite=True)

    fig, ax2 = plt.subplots(figsize=(7,6)) 
    ax2.step(wave[indrange_ion], flux_norm, 'k')
    ax2.step(wave[indrange_ion], err_norm, 'g', alpha=0.4)
    ax2.axhline(1, c='m')
    ax2.set_ylim(-2,4)
    ax2.set_xlabel('wave')
    ax2.set_ylabel("Norm_flux")
    plt.show()

    return

def norm_ion_unbinned(inputspec, z, ion_name, ion, pol_order):  
    spec = ascii.read(inputspec, delimiter=',')
    wave = spec['col1']
    flux = spec['col2']
    err = spec['col3']

    ion_abs = ion * (1+z)

    fig, ax1 = plt.subplots()

    ax1.step(wave, flux, 'b')
    ax1.axvline(ion_abs, c='r')
    plt.gcf().canvas.mpl_connect('key_press_event', close_figure)
    plt.show()

    min_boundleft = input ("Enter the minimum of the left regoin:")
    max_boundleft = input ("Enter the maximum of the left regoin:")
    min_boundright = input ("Enter the minimum of the right regoin:")
    max_boundright = input ("Enter the maximum of the right regoin:")
    
    Lmin = float(min_boundleft)
    Lmax = float(max_boundleft)
    Rmin = float(min_boundright)
    Rmax = float(max_boundright)
    
    cont_ion = [[Lmin, Lmax],[Rmin, Rmax]]
    indcont_ion = ((wave > cont_ion[0][0]) & (wave < cont_ion[0][1])) |((wave > cont_ion[1][0]) & (wave < cont_ion[1][1]))
    indrange_ion = (wave > cont_ion[0][0]) & (wave < cont_ion[1][1])

    linecoeff = np.polyfit(wave[indcont_ion], flux[indcont_ion], pol_order)
    flux_norm = flux[indrange_ion]/np.polyval(linecoeff, wave[indrange_ion])
    err_norm  = err[indrange_ion]/np.polyval(linecoeff, wave[indrange_ion])

    dataout = Table([wave[indrange_ion], flux_norm, err_norm], names=['lambda', 'flux_norm', 'err_norm'])
    ascii.write(dataout, 'test_unbinned.dat', overwrite=True)

    fig, ax2 = plt.subplots(figsize=(7,6)) 
    ax2.step(wave[indrange_ion], flux_norm, 'k')
    ax2.axhline(1, c='m')
    ax2.set_xlabel('wave')
    ax2.set_ylabel("Norm_flux")
    plt.show()

    return


