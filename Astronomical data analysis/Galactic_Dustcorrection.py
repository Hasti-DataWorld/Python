import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy import coordinates as co
from astropy import units as u
import os
import healpy as hp

panstarrs_dust_coefficients = { band : value for ( band, value ) in zip(
    'grizy', [3.62671, 2.71095, 2.07001, 1.59439, 1.33419] ) }
# can replace the constants for different sky surveys


planck_dust_map_filename = '/Users/hnateghi/Desktop/Folders/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
if not os.path.exists( planck_dust_map_filename ):
    print('WARNING: Cannot find Planck dust map with file name:')
    print(' '*8, planck_dust_map_filename)
    print()
    print('This can be downloaded via the Planck Explanatory Supplement wiki:')
    print(' '*8, 'http://wiki.cosmos.esa.int/planckpla')
    print('(Look under Mission Products > CMB and astrophysical component maps,')
    print('and be sure to download the higher resolution Nside=2048 version.')
    print('Note that this is a ~1.6 Gb download!)')
    print()
else :
    print('reading Planck dust map from file:')
    print(' '*8, planck_dust_map_filename)
    planckDustMap = hp.read_map( planck_dust_map_filename, field=2 )
    hp.mollview( planckDustMap, min=0., max=1., fig=1,
                 title='Planck dust map: %s'
                 % planck_dust_map_filename.split('/')[-1] )
print() ; print()


def healpix_coords_from_radec( ra, dec ):
    pos = co.SkyCoord( ra*u.deg, dec*u.deg ).galactic
    theta, phi = np.pi/2. - pos.b.rad, pos.l.rad
    return theta, phi

def Planck_EBV( ra, dec ):
    theta, phi = healpix_coords_from_radec( ra, dec )
    return hp.get_interp_val( planckDustMap, theta, phi )


def dust_correction( extBV, RV=None, band=None ):
    if RV is None and band is None :
        print( 'you must specify either a band, or a value for Rv.' )
        return 0.
    elif band in panstarrs_dust_coefficients.keys():
        RV = panstarrs_dust_coefficients[ band ]
    elif RV is None :
        print( 'could not determine Rv based on specified band:', band )
    return RV * extBV 



data = Table.read('magigroup_psout.fit', format='fits') 

# the coordinates of the targets
ra, dec = data['raStack'], data['decStack']
Ebv = Planck_EBV( ra, dec )



# reading PS1 r-Kron fluxes from Database
for band in 'r':
# turn the correction in mag to a multiplicative scalar for fluxes
    dust_scaling = 10**( 0.4*dust_correction( Ebv, band=band)) 

    gooddata = data[ band+'KronFlux'] > -999.0   

    dust_scaling = np.where( gooddata, dust_scaling, 1 )

    data[band+'KronFlux'] *= ( dust_scaling *10**6 )

for band in 'r':
    data[ band+'KronMag' ] -= dust_correction( Ebv, band=band )
    
# reading PS1 5 bands aperture fluxes from Database
for band in 'grizy':
    dust_scaling = 10**( 0.4*dust_correction( Ebv, band=band)) 

    gooddata2 = data[ band+'c6flxR4'] > -999.0    

    dust_scaling2 = np.where( gooddata2, dust_scaling, 1 )

    
    data[band+'c6flxR4'] *= ( dust_scaling2 *10**6 )
    data[ band+'c6flxR4Err' ] *= ( dust_scaling2 *10**6 )


# saving into fits files
data.write('PSout_corrected.fits')












    
 

