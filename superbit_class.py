import numpy as np
from astropy.io import fits
from astropy.table import Table
import os, re
from astropy.constants import h, c
from astropy import units as u

###############
#### This code defines a  SuperBIT class holding some useful
#### attributes and a function to do a flux conversion.
#### TO DO: Add in astropy units!
##############

class SuperBIT:

    def __init__(self, gain=None, exptime=None):
        '''
        Class to hold some basic attributes of the SuperBIT
        camera and do flux conversions
        '''
        self.exptime = exptime
        self.gain = gain
        self.mirror_diam = 0.5

        self.band_names = None
        self.mean_qe = None
        self.pivot_lam = None
        self.filter_edges = None
        self.mag_colnames = None
        self.flux_colnames = None
 
        if gain == None:
            self.gain = 0.343
            
        if exptime == None:
            self.exptime = 300.0

        self._set_defaults()
        
        self._make_magname_dict()

        return
    
    
    def _set_defaults(self):
        '''
        v. sloppy with band vs. filter, but they mean the same thing
        '''

        mean_qe = {'u': 0.39063022388059704,
                       'b': 0.8085098086124402,
                       'g': 0.6640071428571429, 
                       'r': 0.5715020134228189,
                       'lum': 0.7241146017699115,
                       'shape': 0.49407140468227423
                       }
            
        pivots = {'u': 395.,
                       'b': 476., 
                       'g': 597., 
                       'r': 640., 
                       'lum': 522., 
                       'shape': 650.
                       }

        filter_edges = {'u': (300., 435.),
                            'b': (365., 575.),
                            'g': (515., 705.),
                            'r': (570., 720.),
                            'lum': (370., 710.),
                            'shape': (530., 830.)
                            }
            
        band_names = ['u', 'b', 'g', 'r', 'lum', 'shape']

        self.band_names = band_names
        self.pivot_lam = pivots
        self.mean_qe = mean_qe
        self.filter_edges = filter_edges
        
        return
     
                       
    def _make_magname_dict(self):
        '''
        Utility script to access magnitude column names
        '''
        
        mag_colnames = dict([('u', 'mag_u'), ('blue', 'mag_b'),
                                 ('g', 'mag_g'), ('b', 'mag_b'),
                                 ('r', 'mag_r'), ('nir', 'mag_nir'),
                                 ('lum', 'mag_lum'),
                                 ('shape', 'shapemag')]
                                )
        
        flux_colnames = dict([('u', 'crates_u'), ('blue', 'crates_b'),
                                  ('g', 'crates_g'), ('b', 'crates_b'),
                                  ('r', 'crates_r'), ('nir', 'crates_nir'),
                                  ('lum', 'crates_lum'), ('shape', 'crates_shape')]
                                 )

        self.mag_colnames = mag_colnames
        self.flux_colnames = flux_colnames

        return

    def get_mean_qe(self, path, qe_name):
        '''
        calculate quantum efficiency in some band.
        qe = Table.read('/Users/j.mccleary/Research/SuperBIT/\
        superbit_photometry/data/instrument/camera/imx455.csv')
        '''

        qe = Table.read(os.path.join(path, qe_name))

        with open('mean_gain_imx455.txt', 'w') as f:
            f.write('# bandpass mean_gain\n')
            for band in self.filter_edges.keys():
                    wg = np.where((qe['wavelength'] > bandict[band][0]) \
                                          & (qe['wavelength'] < bandict[band][1]))
                    mean_gain = np.mean(qe['qe'][wg])
                    print(f'mean gain in band {band} is {mean_gain}')
                    f.write(f'{band:s}\t{mean_gain:.5f} \n')

        f.close()

        return 


    def _photons_to_magab(self, flux_density, pivot_lam):
        '''
        Lightly adapted from bittools photometry routines.
        Converts photon flux density from photons/nm/s/m^2 to AB magnitude at
        a given wavelength

        Args:
            photons_per_second (float): value of photon flux density in
            photons/nm/s/m^2 central_wavelength (float): value of photon
            wavelength in nanometers

        Returns:
            float: value of flux density in AB magnitude
        '''

        g_nu = flux_density * u.photon / u.m ** 2 / u.s / u.nm
        f_nu = g_nu.to(
            u.ABmag, equivalencies=u.spectral_density(pivot_lam * u.nm)
            )

        return f_nu.value  # mag AB

   
    def flux_to_abmag(self, tab, band, colname='flux_auto'):
        '''
        Take Source Extractor FLUX_AUTO, turn it into AB mags!

          - grab ADUs from SEXtractor
          - convert to electrons using gain
          - convert from electrons to photons using mean QE
          - take photons to flam via pivot wavelength 
          - use astropy to convert flam to fnu
          - use ABMAG=−2.5×log(fnu)−48.60 
        '''
        
        if band == 'blue': band = 'b'
        
        assert band in self.band_names, f'band {band} not in {self.band_names}'
        
        band_width = self.filter_edges[band][1] - self.filter_edges[band][0]
        pivot_lam = self.pivot_lam[band]
        this_qe = self.mean_qe[band]
        area = np.pi * (self.mirror_diam/2.0)**2 

        fluxes = tab[colname].data

        crates = fluxes * self.gain / self.exptime

        photons = crates / this_qe / area 

        flux_density = photons/(band_width)
        
        ab_mag = self._photons_to_magab(flux_density, pivot_lam)

        return ab_mag

    
