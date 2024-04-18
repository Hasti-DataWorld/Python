from astropy.io import fits
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.table import Column
from astropy.cosmology import FlatLambdaCDM


# 1 -calculate the look back time in universe
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
cosmo.lookback_time(4.5)
ar = np.array([1,0.2,3])
cosmo.lookback_time(ar)


# 2- create a cumulative distribution function out of an array of data
def mkcdf(x):
    # This array does not have to be sorted, this function sorts it

    # New arrays for sorted data and CDF value
    newx = [None]*(len(x)+1)
    newy = [None]*(len(x)+1)

    # Sort x array and place it into newx
    newx = np.sort(x)

    # Put a 0 at the start of the array
    newx = np.insert(newx,0,0)

    # Figure out the values of the cdf
    dy = 1.0/len(x)

    newy[0] = 0

    i = 1
    while i <= len(x):
        newy[i] = newy[i-1] + dy
        i += 1

    return newx,newy




# 3- canlculate the impact parameter using the angular-sep on a 2d spectrum
def impact_param(z,ang_sep):
    
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=69.6, Om0=0.286)
    dc = cosmo.kpc_proper_per_arcmin(z)
    D = dc/60.     #reads the angle in arcsec and convert them to arcmin
    D_kpc = D * ang_sep
    return D_kpc

