import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import constants as consts
from astropy.wcs import WCS
from PyAstronomy import pyasl
from astropy.table import Table
from astropy.table import Column
from astropy.io import ascii

"""
This code snippet applies the helio correction on Keck spectra.

I find the middle exposure from the list of exposures and Keck nightlog. 
I have to read the start time to calculated the correction using middle exp frame. 
Since the time of observation matters for the velocity and position of the Earth 
relative to the sun and we take the middle one, however it is not that significant.

"""

# read the spectrum
spec = fits.open('//')
head = spec[0].header

racenter = head['RA']
decenter = head['DEC']
expstime = head['UT']      # not duration, just start time
expdate = head['DATE-OBS']
coord = " ".join([racenter,decenter])
ra , dec = pyasl.coordsSexaToDeg(coord)
expstart = "T".join([expdate,expstime])


keck = EarthLocation.from_geodetic(lat=19.8283*u.deg, lon=-155.4783*u.deg,
                                   height=4160*u.m)

sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
barycorr = sc.radial_velocity_correction(obstime=Time(expstart),
                                         location=keck)
bcorr = barycorr.to(u.km/u.s)
cval = consts.c.to('km/s').value



with fits.open(inputspec) as hdul:  # reading nonweight finalspec
        hdr = hdul[0].header   
        disper = hdr['CDELT1'] 
        startwave = hdr['CRVAL1']
        pixnum = hdr['NAXIS1'] 
        wave_range = np.arange(0,pixnum)
        wave_delt = wave_range*disper
        wave_vac = startwave + wave_delt
        flux = hdul[0].data

with fits.open(inputspec) as hdulerr: # reading error finalspec
        error = hdulerr[0].data

barywave = wave_vac * (1.0 + (bcorr.value/cval))

dataout = Table([barywave, flux, error], names=['wave', 'flux', 'err'])
dataout.write('heliCorr.fits')
