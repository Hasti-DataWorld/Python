"""
This code reads the full spectrum, finds the metal and hydrogen absorption lines 
and create a plot of the spectrum chunks

can make the plots in wavelength and velocity space

by HNateghi-2023
"""



from astropy.io import fits 
import numpy as np 
from matplotlib.pylab import *
import pandas as pd
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker
from scipy.integrate import quad
from astropy import units as u
from astropy import constants
from astropy.io import ascii
from astropy.table import Table, column
from astropy.table import vstack
from matplotlib import rc
import os
import matplotlib.cm as cm
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from astropy import constants as consts

plt.rcParams['xtick.labelsize'] = 5
plt.rcParams['ytick.labelsize'] = 5
cval = consts.c.to('km/s').value

# read the full spectrum
wave, flux, error = np.genfromtxt('./.dat',usecols=(0,1,2),skip_header=1, unpack=True)

zabs = 0.0587548

ions = np.array([1215.67,1025.74,972.54,949.74,937.80,930.75,926.23,923.15,920.96,919.35,
                    918.13,917.18,1334.54,1036.34 ,1548.20,1550.78,1083.99,989.80,1238.82,
                    1242.80,1302.17,1039.23,988.77,1031.93,1037.62,989.87,1190.42,1193.29,
                    1260.42,1304.37,1526.71,1206.50,1393.76,1402.77, 3934.78,
                    2796.35,2803.53, 977.02, 834.47, 833.33, 832.76, 832.93, 702.33, 1012.49, 1062.664])

ions_name = ["HI","HI","HI","HI","HI","HI","HI","HI","HI","HI","HI","HI","CII","CII",
            "CIV","CIV","NII","NIII","NV","NV","OI","OI","OI","OVI","OVI","SiII",
            "SiII","SiII","SiII","SiII","SiII","SiIII","SiIV","SiIV","CaII","MgII","MgII",
            "CIII", "OII", "OII", "OII", "OIII", "OIII", "SIII", "SIV"]



detection_wave = []
detection_name = []

for i in range(len(ions)):
    ion_obs = ions[i] * (1+zabs) 

    if ion_obs > np.min(wave) and ion_obs < np.max(wave):
        detection_wave.append(ions[i])
        detection_name.append(ions_name[i])


detection_wave_arr = np.array(detection_wave)

square = np.sqrt(len(detection_wave_arr))
check_square = len(detection_wave_arr)%square
non_square = 1

if check_square == 0:
    nrows = int(np.sqrt(len(detection_wave_arr)))
    ncols = int(np.sqrt(len(detection_wave_arr)))
else:
    nrows = int(np.sqrt(len(detection_wave_arr)))+non_square
    ncols = int(np.sqrt(len(detection_wave_arr)))

plt.close()
fig, ax = plt.subplots(nrows, ncols, sharex=True, sharey=False, figsize = (9,8))

def velspace(wave,flux,ion,zabs):
    zlam = ion*(1.0+zabs)

    vel = np.empty(len(wave))
    for i in range(len(wave)):
        velcalc = cval*(wave[i]-zlam)/zlam
        vel[i] = velcalc

    vind500, = np.where((vel>=-500.0) & (vel<=500.0))

    velarr = vel[vind500[0]:vind500[-1]]
    fluxarr = flux[vind500[0]:vind500[-1]]

    if velarr[0] > -500:
        velarr = np.insert(velarr,0,-500.0)
        fluxarr = np.insert(fluxarr,0,1.0)

    if velarr[-1] < 500:
        velarr = np.append(velarr,500.0)
        fluxarr = np.append(fluxarr,1.0)

    return(velarr,fluxarr)


for i in range(nrows):
    for j in range(ncols):
        if(j + ncols*i) == len(detection_wave_arr):
            break 
    
        if i==0:
            
            vel, v_flux = velspace(wave,flux,detection_wave_arr[j],zabs)
            vel, v_err = velspace(wave,error,detection_wave_arr[j],zabs)
            ax[i,j].step(vel, v_flux, c='k', linewidth=0.65)
            ax[i,j].step(vel[v_err!=1], v_err[v_err!=1], c='k', linewidth=0.65, alpha=0.3)
            ax[i,j].axvline(x=0, c='g', ls='--', lw=0.6)
            ax[i,j].set_title(detection_name[j]+str(detection_wave_arr[j]), fontsize=6.8, c='r')
            ax[i,j].tick_params(axis='both', which='major', length=2.5)
            #ax[i,j].set_ylim( np.min(flux_inrange), np.max(flux_inrange))
        else:
            vel, v_flux = velspace(wave,flux,detection_wave_arr[j+ncols*i],zabs)
            vel, v_err = velspace(wave,error,detection_wave_arr[j+ncols*i],zabs)
            ax[i,j].step(vel, v_flux, c='k', linewidth=0.65)
            ax[i,j].step(vel[v_err!=1], v_err[v_err!=1], c='k', linewidth=0.65, alpha=0.3)
            ax[i,j].axvline(x=0, c='g', ls='--', lw=0.6)
            ax[i,j].set_title(detection_name[j+ncols*i]+str(detection_wave_arr[j+ncols*i]), fontsize=6.8, c='r')
            ax[i,j].tick_params(axis='both', which='major', length=2.5)


y = int(np.sqrt(len(detection_wave_arr))/2)
ax[nrows-1,y].set_xlabel(r'Velocity (km/s)', fontsize=9)
ax[y,0].set_ylabel(r'Flux', fontsize=9)

fig.subplots_adjust(wspace=0.2, hspace=0.3)
fig.savefig('./.png', overwrite=True, bbox_inches='tight',dpi=300)
       









