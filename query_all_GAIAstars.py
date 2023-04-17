import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord
import galsim
from astropy.table import Table, Row
import pdb
from astroquery.vizier import Vizier
import os

import instrument as inst
import photometry as phot
import galaxy_params as gp
import psf as psf_module

def get_transmission(band, path=None):
    if path is None:
        path = './'
    if band == 'u':
        return np.genfromtxt(os.path.join(path,'data/instrument/bandpass/u_2023.csv'), delimiter=',')[:, 2][1:]
    elif band == 'b':
        return np.genfromtxt(os.path.join(path,'data/instrument/bandpass/b_2023.csv'), delimiter=',')[:, 2][1:]
    elif band == 'g':
        return np.genfromtxt(os.path.join(path,'data/instrument/bandpass/g_2023.csv'), delimiter=',')[:, 2][1:]
    elif band == 'r':
        return np.genfromtxt(os.path.join(path,'data/instrument/bandpass/r_2023.csv'), delimiter=',')[:, 2][1:]
    elif band == 'nir':
        return np.genfromtxt(os.path.join(path,'data/instrument/bandpass/nir_2023.csv'), delimiter=',')[:, 2][1:]
    elif band == 'lum':
        return np.genfromtxt(os.path.join(path,'data/instrument/bandpass/lum_2023.csv'), delimiter=',')[:, 2][1:]
    else:
        raise ValueError("Invalid band.")

class SciImage:
    '''
    Stores image attributes for star positions
    '''

    def __init__(self):
        self.sci_img_size_x_pix = None
        self.sci_img_size_y_pix = None
        self.sci_img_size_x_arcsec = None
        self.sci_img_size_y_arcsec = None
        self.cluster_x_pix = None
        self.cluster_y_pix = None
        self.height = None
        self.width = None
        self.image = None

    def make_sci_image(self, camera, bandpass, cluster_coord):
        '''
        Get WCS & detector image coordinates for Vizier query
        Issue many reassuring print statements
        '''
        sci_img_size_x_pix = camera.npix_H
        sci_img_size_y_pix = camera.npix_V
        sci_img_size_x_arcsec = sci_img_size_x_pix * bandpass.plate_scale
        sci_img_size_y_arcsec = sci_img_size_y_pix * bandpass.plate_scale

        print("Scicam X: {}, Y: {}".format(sci_img_size_x_pix,
                                        sci_img_size_y_pix))

        cluster_ra_deg = cluster_coord.ra.to(u.deg)
        cluster_dec_deg = cluster_coord.dec.to(u.deg)
        cluster_x_pix = (sci_img_size_x_pix.value/2)
        cluster_y_pix = (sci_img_size_y_pix.value/2)

        print("Cluster ra and dec: {:.2f} and {:.2f}".format(
            cluster_coord.ra, cluster_coord.dec))
        print("Cluster x and y: {:.2f} and {:.2f}".format(
            cluster_x_pix, cluster_y_pix))

        height = ((sci_img_size_y_arcsec.value) * u.arcsec).to(u.deg)
        width = ((sci_img_size_x_arcsec.value) * u.arcsec).to(u.deg)
        print(f'sci image height is {height} width is {width}')
        self.height = height; self.width = width

        sci_img_max_ra = (cluster_ra_deg + (sci_img_size_x_arcsec/2)).to(u.deg)
        sci_img_min_ra = (cluster_ra_deg - (sci_img_size_x_arcsec/2)).to(u.deg)
        sci_img_max_dec = (cluster_dec_deg +
                        (sci_img_size_y_arcsec/2)).to(u.deg)
        sci_img_min_dec = (cluster_dec_deg -
                        (sci_img_size_y_arcsec/2)).to(u.deg)

        print("Sci image min RA: {:.2f}, max RA: {:.2f}".format(
        sci_img_min_ra, sci_img_max_ra))
        print("Sci image min DEC: {:.2f}, max DEC: {:.2f}".format(
        sci_img_min_dec, sci_img_max_dec))

        return


def query_gaia(cluster_coord, sci_im):
    '''
    Query GAIA catalog at position of cluster, keep only well-detected
    stars that lie within detector, also calculate magnitudes
    '''
    '''
    # Columns of the result by default, have to add PQSO and PGal after the fact
    columns=['RA_ICRS', 'DE_ICRS', 'Source', 'e_RA_ICRS', 'e_DE_ICRS', 
                             'Plx', 'e_Plx', 'PM', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'RUWE', 
                             'FG', 'e_FG', 'Gmag', 'FBP', 'e_FBP', 'BPmag', 'FRP', 'e_FRP', 'RPmag', 
                             'BP-RP', 'RV', 'e_RV', 'Vbroad', 'GRVSmag', 'QSO', 'Gal', 'NSS', 'XPcont', 
                             'XPsamp', 'RVS', 'EpochPh', 'EpochRV', 'MCMCGSP', 'MCMCMSC', 'And', 
                             'Teff', 'logg', '__Fe_H_', 'Dist', 'A0', 'HIP', 'PS1', 'SDSS13', 
                             'SKYM2', 'TYC2', 'URAT1', 'AllWISE', 'APASS9', 'GSC23', 'RAVE5', 
                             '_2MASS', 'RAVE6', 'RAJ2000', 'DEJ2000']
    # Kept old query call for historical reasons
    result = Vizier.query_region(coordinates=cluster_coord,
                             height=sci_im.height,
                             width=sci_im.width,
                             catalog='I/355/gaiadr3')
    '''

    # Load in Catalog using Vizier
    catalogs = Vizier(catalog='I/355/gaiadr3')

    # ROW_LIMIT defaults to 50, set to -1 to get all rows
    catalogs.ROW_LIMIT = -1

    # add PQSO and PGal columns to query
    catalogs.columns += ['PQSO', 'PGal']

    # Query catalog
    result = catalogs.query_region(coordinates=cluster_coord,
                               width=sci_im.width,
                               height=sci_im.height)
    

    result = result[0].filled()
    ZP_G = 25.7915509947
    ZP_BP = 25.3861560855

    df_stars = pd.DataFrame()
    df_stars['RA_ICRS'] = result['RA_ICRS']
    df_stars['DE_ICRS'] = result['DE_ICRS']
    df_stars['FG'] = result['FG']
    df_stars['FBP'] = result['FBP']
    df_stars['PGal'] = result['PGal']
    df_stars['PQSO'] = result['PQSO']
    df_stars['BPmag'] = (-2.5 * np.log10(result['FBP'])) + ZP_BP
    df_stars['Gmag'] = (-2.5 * np.log10(result['FG'])) + ZP_G
    df_stars.dropna()

    # Remove stars that have essentially any probability of being Quasars or Galaxies
    notgalorquasar = (df_stars.PGal < 0.01) & (df_stars.PQSO < 0.01)
    df_stars = df_stars[notgalorquasar]

    # Yay, have stars now
    # Clean it up.
    df_stars = df_stars[df_stars['RA_ICRS'] >=
                    (cluster_coord.ra.value-(sci_im.width.value/2))]
    df_stars = df_stars[df_stars['RA_ICRS'] <=
                    (cluster_coord.ra.value+(sci_im.width.value/2))]
    df_stars = df_stars[df_stars['DE_ICRS'] >= (
                    cluster_coord.dec.value-(sci_im.height.value/2))]
    df_stars = df_stars[df_stars['DE_ICRS'] <= (
                    cluster_coord.dec.value+(sci_im.height.value/2))]
    df_stars = df_stars[df_stars['Gmag'] >= 0]
    df_stars = df_stars[df_stars['BPmag'] >= 0]
    df_stars = df_stars.reset_index(drop=True)

    print("Adding {} stars".format(len(df_stars)))

    return df_stars


def plot_gaia(result):
    '''
    fig,ax = plt.subplots(1,1,figsize=(13,9))

    ax.hist(result['BPmag'],bins=100,label='BPmag',log=True,color='blue',range=[10,25],histtype='step',lw=2)
    ax.hist(result['RPmag'],bins=100,label='RPmag',log=True,color='red',range=[7,22],histtype='step',lw=2)
    ax.hist(result['Gmag'],bins=100,label='Gmag',log=True,color='darkgreen',range=[7,22],histtype='step',lw=2)
    ax.legend()
    ax.set_xlabel('Gaia Magnitude')
    ax.set_ylabel('N')
    ax.set_title('All GAIA stars for APRA clusters')
    fig.tight_layout()

    fig.savefig('Gaia_passband_lumfunc_Apra.png')
    '''
    return

# Has been updated to do all bandpasses
def get_BIT_fluxes(df_stars, camera, bandpass, exp_time, telescope):
    '''
    Turn thing into other thing
    '''
    # Store fluxes in BIT filters
    bit_flux_b = []
    bit_flux_lum = []
    bit_flux_u = []
    bit_flux_nir = []
    bit_flux_g = []
    bit_flux_r = []

    # Store transmission
    b_transmission = get_transmission('b', path=None)
    lum_transmission = get_transmission('lum', path=None)
    u_transmission = get_transmission('u', path=None)
    nir_transmission = get_transmission('nir', path=None)
    g_transmission = get_transmission('g', path=None)
    r_transmission = get_transmission('r', path=None)

    # start with b
    for idx in range(len(df_stars)):
        if idx % 50 == 0:
            print(idx)
        gmag = df_stars['BPmag'][idx]
        mean_fnu_gmag = phot.abmag_to_mean_fnu(abmag=gmag)
        mean_flambda = phot.mean_flambda_from_mean_fnu(mean_fnu=mean_fnu_gmag,
                                            bandpass_transmission=b_transmission,
                                            bandpass_wavelengths=bandpass.wavelengths)
        crate_electrons_pix = phot.crate_from_mean_flambda(mean_flambda=mean_flambda,
                                illum_area=telescope.illum_area.value,
                                bandpass_transmission=b_transmission,
                                bandpass_wavelengths=bandpass.wavelengths)
        #crate_adu_pix = crate_electrons_pix / camera.gain.value
        flux_electrons_b = crate_electrons_pix * exp_time.value

        if flux_electrons_b > 4*10**6:
            flux_electrons_b = 4*10**6

        # Save flux & positions
        bit_flux_b.append(flux_electrons_b)


    # Next, do lum
    for idx in range(len(df_stars)):
        if idx % 50 == 0:
            print(idx)
        gmag = df_stars['BPmag'][idx]
        mean_fnu_gmag = phot.abmag_to_mean_fnu(abmag=gmag)
        mean_flambda = phot.mean_flambda_from_mean_fnu(mean_fnu=mean_fnu_gmag,
                                            bandpass_transmission=lum_transmission,
                                            bandpass_wavelengths=bandpass.wavelengths)
        crate_electrons_pix = phot.crate_from_mean_flambda(mean_flambda=mean_flambda,
                                illum_area=telescope.illum_area.value,
                                bandpass_transmission=lum_transmission,
                                bandpass_wavelengths=bandpass.wavelengths)
        #crate_adu_pix = crate_electrons_pix / camera.gain.value
        flux_electrons_lum = crate_electrons_pix * exp_time.value

        if flux_electrons_lum > 4*10**6:
            flux_electrons_lum = 4*10**6

        # Save flux & positions
        bit_flux_lum.append(flux_electrons_lum)
    
    # Next, do u
    for idx in range(len(df_stars)):
        if idx % 50 == 0:
            print(idx)
        gmag = df_stars['BPmag'][idx]
        mean_fnu_gmag = phot.abmag_to_mean_fnu(abmag=gmag)
        mean_flambda = phot.mean_flambda_from_mean_fnu(mean_fnu=mean_fnu_gmag,
                                            bandpass_transmission=u_transmission,
                                            bandpass_wavelengths=bandpass.wavelengths)
        crate_electrons_pix = phot.crate_from_mean_flambda(mean_flambda=mean_flambda,
                                illum_area=telescope.illum_area.value,
                                bandpass_transmission=u_transmission,
                                bandpass_wavelengths=bandpass.wavelengths)
        #crate_adu_pix = crate_electrons_pix / camera.gain.value
        flux_electrons_u = crate_electrons_pix * exp_time.value

        if flux_electrons_u > 4*10**6:
            flux_electrons_u = 4*10**6

        # Save flux & positions
        bit_flux_u.append(flux_electrons_u)

    # Next, do nir
    for idx in range(len(df_stars)):
        if idx % 50 == 0:
            print(idx)
        gmag = df_stars['BPmag'][idx]
        mean_fnu_gmag = phot.abmag_to_mean_fnu(abmag=gmag)
        mean_flambda = phot.mean_flambda_from_mean_fnu(mean_fnu=mean_fnu_gmag,
                                            bandpass_transmission=nir_transmission,
                                            bandpass_wavelengths=bandpass.wavelengths)
        crate_electrons_pix = phot.crate_from_mean_flambda(mean_flambda=mean_flambda,
                                illum_area=telescope.illum_area.value,
                                bandpass_transmission=nir_transmission,
                                bandpass_wavelengths=bandpass.wavelengths)
        #crate_adu_pix = crate_electrons_pix / camera.gain.value
        flux_electrons_nir = crate_electrons_pix * exp_time.value

        if flux_electrons_nir > 4*10**6:
            flux_electrons_nir = 4*10**6

        # Save flux & positions
        bit_flux_nir.append(flux_electrons_nir)

    # Next, do g
    for idx in range(len(df_stars)):
        if idx % 50 == 0:
            print(idx)
        gmag = df_stars['BPmag'][idx]
        mean_fnu_gmag = phot.abmag_to_mean_fnu(abmag=gmag)
        mean_flambda = phot.mean_flambda_from_mean_fnu(mean_fnu=mean_fnu_gmag,
                                            bandpass_transmission=g_transmission,
                                            bandpass_wavelengths=bandpass.wavelengths)
        crate_electrons_pix = phot.crate_from_mean_flambda(mean_flambda=mean_flambda,
                                illum_area=telescope.illum_area.value,
                                bandpass_transmission=g_transmission,
                                bandpass_wavelengths=bandpass.wavelengths)
        #crate_adu_pix = crate_electrons_pix / camera.gain.value
        flux_electrons_g = crate_electrons_pix * exp_time.value

        if flux_electrons_g > 4*10**6:
            flux_electrons_g = 4*10**6

        # Save flux & positions
        bit_flux_g.append(flux_electrons_g)

    # Finally, do r
    for idx in range(len(df_stars)):
        if idx % 50 == 0:
            print(idx)
        gmag = df_stars['BPmag'][idx]
        mean_fnu_gmag = phot.abmag_to_mean_fnu(abmag=gmag)
        mean_flambda = phot.mean_flambda_from_mean_fnu(mean_fnu=mean_fnu_gmag,
                                            bandpass_transmission=r_transmission,
                                            bandpass_wavelengths=bandpass.wavelengths)
        crate_electrons_pix = phot.crate_from_mean_flambda(mean_flambda=mean_flambda,
                                illum_area=telescope.illum_area.value,
                                bandpass_transmission=r_transmission,
                                bandpass_wavelengths=bandpass.wavelengths)
        #crate_adu_pix = crate_electrons_pix / camera.gain.value
        flux_electrons_r = crate_electrons_pix * exp_time.value

        if flux_electrons_r > 4*10**6:
            flux_electrons_r = 4*10**6

        # Save flux & positions
        bit_flux_r.append(flux_electrons_r)


    # Add columns to df stars
    df_stars['bitflux_electrons_b'] = bit_flux_b
    df_stars['bitflux_electrons_lum'] = bit_flux_lum
    df_stars['bitflux_electrons_u'] = bit_flux_u
    df_stars['bitflux_electrons_nir'] = bit_flux_nir
    df_stars['bitflux_electrons_g'] = bit_flux_g
    df_stars['bitflux_electrons_r'] = bit_flux_r

    return df_stars


def main():
    '''
    Obtain GAIA star catalogs for every cluster in the full APRA target list
    Save catalog to file, and make at least one luminosity function plot

    As a bonus, take the current ASCII APRA file and make it a nice table
    with the number of GAIA stars as an additional column
    '''
    n_exp = 1
    exp_time = 1*u.s
    Vizier.ROW_LIMIT = -1

    # allowed filters are u, b, g, r, i and shape (2021) or lum (2019)
    data_path = '/home/wslgeorgios/superbit_photometry'

    apra_name = '/home/wslgeorgios/bitfluxes/APRA_target_clusters.csv'
    apra = Table.read(apra_name,format='ascii',delimiter=',')
    output_table = Table()

    telescope = inst.Telescope('superbit')
    camera = inst.Camera('imx455')
    bandpass = inst.Bandpass('superbit')
    bandpass.wavelengths = camera.wavelengths

    plate_scale = ((206265 * camera.pixel_size.to(u.micron).value)
                        / (1000 * telescope.focal_length.to(u.mm))).value * (u.arcsec/u.pixel)
    bandpass.plate_scale = plate_scale

    ##
    ## Loop over objects in APRA catalog, get star cat
    ##

    #primary_hdu = fits.PrimaryHDU(n)
    names = []; ra = []; dec = []; nstars = []

    for cluster in apra:

        cluster_coord = SkyCoord(cluster['col2'], cluster['col3'],
                        unit=(u.hourangle, u.deg), frame='icrs')
        # Has science image attributes
        sci_image = SciImage()
        sci_image.make_sci_image(camera, bandpass, cluster_coord)

        # This is a pandas data frame
        vizier_result = query_gaia(cluster_coord, sci_image)

        # Convert to flux? Also has all the info from Vizier result
        # I think this is a data frame
        cluster_fluxes = get_BIT_fluxes(vizier_result, camera, bandpass, exp_time, telescope)

        string = 'GAIAstars_{}.csv'.format(cluster['col1'].replace(' ','_'))

        names.append(cluster['col1'])
        ra.append(cluster['col2'])
        dec.append(cluster['col3'])
        nstars.append(len(cluster_fluxes))

        cluster_fluxes.to_csv(string, index_label='NUM')

    output_tab = Table(data=[names, ra, dec, nstars], names=['names', 'ra', 'dec', 'nstars'])
    output_tab.write('apra_vizier_starinfo.fits', format='fits', overwrite=True)



    return 0

if __name__ == '__main__':

    print('executing main')
    rc = main()

    if rc !=0:
        import pdb
        #raise Exception
        pdb.set_trace()
