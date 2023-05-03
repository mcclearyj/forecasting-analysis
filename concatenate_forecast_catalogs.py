from astropy.table import Table, vstack
import numpy as np
import os
import glob
from esutil import htm
from superbit_class import SuperBIT
from argparse import ArgumentParser
import ipdb

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('basedir', default=None,
                            help = 'Path to simulations')
    parser.add_argument('-cosmosdir', default=None,
                            help = 'Path to COSMOS catalogs')
    
    return parser.parse_args()


class COSMOSCats:

    def __init__(self, cosmosdir, shape_cat_name=None, f2023_cat_name=None):
        '''
        A bit unwieldy, but hold COSMOS catalogs in memory for matching against
        joined and lensing-sample galaxy catalogs. Required because the 
        shape band is not in the COSMOS2023 catalog
        '''
        
        self.cosmosdir = cosmosdir
        self.shape_cat_name = shape_cat_name
        self.f2023_cat_name = f2023_cat_name
        self.shape_cat = None
        self.f2023_cat = None
        
        self._load_catalogs()

        
    def _load_catalogs(self):
        '''
        Load in COSMOS catalogs with both 2023 bands and old shape band
        '''
        cosmosdir = self.cosmosdir
        shape_cat_name = self.shape_cat_name
        f2023_cat_name = self.f2023_cat_name

        if shape_cat_name is None:
            shape_cat_name = 'cosmos2015_cam2021_filt2021_crates_fullshape.fits'
            print(f'Using default shape band COSMOS catalog {shape_cat_name}')
            
        if f2023_cat_name is None:
            f2023_cat_name = 'cosmos15_superbit2023_phot_shapes.csv' 
            print(f'Using default 2023 bands COSMOS catalog {f2023_cat_name}')
            
        
        shape_cat_path = os.path.join(cosmosdir, shape_cat_name)
        f2023_cat_path = os.path.join(cosmosdir, f2023_cat_name)
 
        if not os.path.isdir(cosmosdir):
            raise OSError(f'No such directory: {cosmosdir}')
        
        if not os.path.exists(shape_cat_path):
            raise OSError(f'No shape band catalog found at \n{shape_cat_path}')
        else:
            self.shape_cat = Table.read(shape_cat_path)
        
        if not os.path.exists(f2023_cat_path):
            raise OSError(f'No 2023 band catalog found at \n{f2023_cat_path}')
        else:
            self.f2023_cat = Table.read(f2023_cat_path)

            
    def choose_cat(self, bandpass):
        '''
        Do the matching thing, I guess? So one needs band, and one needs catalogs
        This is something that has to be done at the individual cluster realization level, 
        much like the z, etc. histograms
        '''
        
        if bandpass == 'shape':
            return self.shape_cat
        else:
            return self.f2023_cat

    
class ConcatGalCats:

    def __init__(self, basedir, cluster_name, forecast_name, bandpass, vb=False):
        '''
        Concatenate galaxy catalogs, including making a 
        COSMOS catalog-based luminosity function
        
        Note: self.cosmos is a COSMOSCats object
        Note: self.superbit is a SuperBIT object
        '''
        self.basedir = basedir
        self.cluster_name = cluster_name
        self.forecast_name = forecast_name
        self.bandpass = bandpass
        self.cosmos = None
        self.superbit = None
        self.truth_tables = None
        self.vb = vb
        

        
    def match_coords(self, cat1, cat2):
        '''
        Utility function to match cat1 to cat 2 using celestial coordinates
        '''

        # Either 'ra/dec' or 'ALPHAWIN_J2000/DELTAWIN_J2000'!

        if 'ra' in cat1.colnames:
            cat1_ra = cat1['ra']
            cat1_dec =  cat1['dec']
        elif 'ALPHAWIN_J2000' in cat1.colnames:
            cat1_ra = cat1['ALPHAWIN_J2000']
            cat1_dec =  cat1['DELTAWIN_J2000']
        else:
            raise KeyError('non-standard RA/Dec column in cat1')
        
        if 'ra' in cat2.colnames:
            cat2_ra = cat2['ra']
            cat2_dec =  cat2['dec']
        elif 'ALPHAWIN_J2000' in cat2.colnames:
            cat2_ra = cat2['ALPHAWIN_J2000']
            cat2_dec =  cat2['DELTAWIN_J2000']
        else:
            raise KeyError('non-standard RA/Dec column in cat2')

        cat1_matcher = htm.Matcher(16, ra=cat1_ra, dec=cat1_dec)
        
        cat2_ind, cat1_ind, dist = cat1_matcher.match(ra=cat2_ra,
                                                      dec=cat2_dec, 
                                                      maxmatch=1, 
                                                      radius=0.5/3600.
                                                      )
        if self.vb == True:
            print(f'{len(dist)}/{len(cat1)} gals matched to truth')
                      
        return cat1[cat1_ind], cat2[cat2_ind]
                                                          
    
    def _load_truth_tables(self, joined_names):
        '''
        A bit of a kludge to ensure that the truth table
        matches the joined table too! This returns a list of Tables
        '''
        
        forecast_name = self.forecast_name
        all_truth_tables = []
        
        for joined in joined_names:
            
            path = os.path.dirname(joined)
            truth_file = os.path.join(path, f'{forecast_name}_truth.fits')
            this_truth_tab = Table.read(truth_file)
            all_truth_tables.append(this_truth_tab)
            
        return all_truth_tables

    
    def concat_gal_tabs(self, tab_names):
        '''
        Open each table in tab_name, extract the ra/dec, redshifts, SEXtractor ID,
        MAG_AUTO and SNR_WIN, and create a Table() object with the results
        '''
        sb = self.superbit
        mag_col = sb.mag_colnames[self.bandpass]
        flux_col = sb.flux_colnames[self.bandpass]

        # Don't actually need cosmos anymore?
        cosmos = self.cosmos
                
        redshifts = []; ra = []; dec = []; sid = []; mag = []
        snr = []; flux = []; flux_err = []; cosmos_id = []

        for i, tab_name in enumerate(tab_names):

            z = Table.read(tab_name)
            
            # First, (re-)match to truth in RA/Dec
            this_truth = self.truth_tables[i]
            wg  = this_truth['redshift'] > 0
            this_truth = this_truth[wg]

            truth_match, z_match = self.match_coords(this_truth, z)

            # Now get COSMOS ID. Could get   
            # matching_cosmos_entries == self.cosmos[truth_match['cosmos_index']]
            cosmos_id.extend(truth_match['cosmos_index'].data)
            
            # Then get easy bits
            redshifts.extend(truth_match['redshift'].data)
            ra.extend(z_match['ra'].data)
            dec.extend(z_match['dec'].data)
            sid.extend(z_match['id'].data)
            mag.extend(z_match['MAG_AUTO'].data)
            flux.extend(z_match['FLUX_AUTO'].data)
            flux_err.extend(z_match['FLUXERR_AUTO'].data)
            snr.extend(z_match['SNR_WIN'].data)


        zcol = Table.Column(redshifts, name='redshift')
        ra_col = Table.Column(ra, name='ra')
        dec_col = Table.Column(dec, name='dec')
        id_col = Table.Column(sid, name='sex_id')
        mag_col = Table.Column(mag, name='mag_auto')
        flux_col = Table.Column(flux, name='flux_auto')
        fluxerr_col = Table.Column(flux_err, name='fluxerr_auto')
        snr_col = Table.Column(snr, name='snr_win')
        cosmos_id_col = Table.Column(cosmos_id, name='COSMOS_id')

        ztab = Table()
        ztab.meta['n_tables'] = len(tab_names)
        ztab.add_columns([id_col, ra_col, dec_col, zcol, mag_col,  \
                              flux_col, fluxerr_col, snr_col, cosmos_id_col])

        return ztab
      
    
    def make_a_tab(self, sbit):
        '''
        Find all realization galaxy catalogs matching cluster_name and forecast_name, 
        then call get_a_tab() to do the catalog concatenating. Save results to file.
    
        This process will be repeated for both the "joined" (pre-lensing selection) galaxy 
        catalogs and the "annular" (post-lensing selection) catalogs. 
        
        Note: sb is a SuperBIT object defined for flux conversions
        '''
        
        cluster_name = self.cluster_name
        forecast_name = self.forecast_name
        bandpass = self.bandpass
        basedir = self.basedir

        if forecast_name == 'forecast_blue':
            out_name = 'forecast_b'
        else:
            out_name = forecast_name

        # Name the outputs something reasonable
    
        joined_gals_name = f'{cluster_name}_{out_name}_gals_joined_master_cat.fits'
        select_gals_name = f'{cluster_name}_{out_name}_annular_gals_master_cat.fits'

        # Find all matching catalogs
    
        all_joined_names = glob.glob(os.path.join(basedir,forecast_name,\
                                        cluster_name, 'r*/*gals_joined_catalog.fits'))
        all_select_names = glob.glob(os.path.join(basedir,forecast_name, \
                                        cluster_name,'r*/*annular.fits'))
        all_joined_names.sort(); all_select_names.sort()

        assert len(all_joined_names)==len(all_select_names), "oh no, N_annulars =/= N_joined"
                                                
        # Get corresponding truth files
        
        all_truth_tables = self._load_truth_tables(all_joined_names)
        self.truth_tables = all_truth_tables

        # Make and save the galaxy tables
        
        if (len(all_joined_names) != 0):
            
            joined_gals_master = self.concat_gal_tabs(all_joined_names)
            select_gals_master = self.concat_gal_tabs(all_select_names)

            ab_mag_joined = sbit.flux_to_abmag(joined_gals_master, bandpass)
            ab_mag_select = sbit.flux_to_abmag(select_gals_master, bandpass)

            joined_gals_master.add_column(ab_mag_joined, name='ab_mag')
            select_gals_master.add_column(ab_mag_select, name='ab_mag')

            
            joined_gals_master.write(joined_gals_name, format='fits', overwrite=True)        
            select_gals_master.write(select_gals_name, format='fits', overwrite=True)

            print(f'\nDone\n')

        else:
            
            print(f'No joined galaxy catalogs found for {cluster_name} in {forecast_name}')
            ipdb.set_trace()
        return

    
def main(args):
    
    basedir = args.basedir
    cosmosdir = args.cosmosdir

    # CCV: cosmosdir = /users/jmcclear/data/superbit/superbit-metacal/GalSim/data
    # Disco: cosmosdir = /work/mccleary_group/superbit/mock-data-forecasting/mock_catalogs/ 
    
    #bands = ['u', 'b']
    #bandnames =  ['u', 'blue']

    bands =  ['lum', 'shape']
    bandnames = ['lum', 'shape']

    redshifts = ['0.059', '0.3', '0.45']
    masses = ['m4.1e14']

    # Create SuperBIT instance to do flux_auto to ABmag conversion
    sbit = SuperBIT()

    # Create COSMOSCats instance to store COSMOS catalogs -- not needed?
    
    # cosmos_cats = COSMOSCats(cosmosdir)
    
    # There HAS to be a better way to do this
    
    for mass in masses:
        for z in redshifts:
            
            cluster_name = f'cl_{mass}_z{z}'
            print(f'\nWorking on {cluster_name}\n')
            
            for i, band in enumerate(bands):
                
                bandname = bandnames[i]
                forecast_name = f'forecast_{bandname}'
                
                print(f'\nWorking on {cluster_name}: {forecast_name}\n')

                forecast_cats = ConcatGalCats(basedir=basedir,
                                                  cluster_name=cluster_name,
                                                  forecast_name=forecast_name,
                                                  bandpass=band
                                                  )
                
                #forecast_cats.cosmos = cosmos_cats.choose_cat(bandpass=bandpass)
                forecast_cats.superbit = sbit
                forecast_cats.make_a_tab(sbit)

    return 0
    
if __name__ == '__main__':

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('concatenate_forecast_catalogs.py has completed succesfully')
    else:
        print(f'concatenate_forecast_catalogs.py has failed w/ rc={rc}')

                
        
            
