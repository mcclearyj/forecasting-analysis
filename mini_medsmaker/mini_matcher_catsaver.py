from astropy.table import Table, vstack
import numpy as np
import os
import glob
from esutil import htm

import ipdb

class MiniMatcher():

    def __init__(self, basedir, run_name, joined_name=None, annular_name=None, vb=True):
        '''
        Read in a mini-coadd catalog, match to annular and to joined
        galaxy catalogs in realization-level directory, save the output,
        and also save the number of matches to the annular and joined
        catalogs for each each mini-coadd. Then write those all to file
        '''
        self.basedir = basedir
        self.run_name = run_name
        self.joined_name = joined_name
        self.annular_name = annular_name
        self.sexcat_name = None
        self.sex_cat = None
        self.joined_cat = None
        self.annular_cat = None
        self.num_matches = None
        self.vb = True
        self.min_snr = None

        self._load_catalogs()
        self._create_nmatch_dict()


    def _load_catalogs(self):
        '''
        Utility functions to load in annular & joined catalogs.
        '''
        basedir = self.basedir
        run_name = self.run_name

        if self.joined_name == None:
            joined_name = run_name + '_gals_joined_catalog.fits'
            joined_cat_path = os.path.join(basedir, joined_name)
            self.joined_name = joined_name

        if self.annular_name == None:
            annular_name = run_name + '_annular.fits'
            annular_cat_path = os.path.join(basedir, annular_name)
            self.annular_name = annular_name

        if os.path.exists(joined_cat_path) == False:
            raise OSError(f'No joined catalog {joined_cat_path} found')
        else:
            joined_cat = Table.read(joined_cat_path)
            self.joined_cat = joined_cat

        if os.path.exists(annular_cat_path) == False:
            raise OSError(f'No annular catalog {annular_cat_path} found')
        else:
            annular_cat = Table.read(annular_cat_path)
            self.annular_cat = annular_cat

        if self.sexcat_name is not None:
            sexcat_path = os.path.join(basedir, self.sexcat_name)
            self.sex_cat = Table.read(sexcat_path)


    def _create_nmatch_dict(self):
        '''
        Utility function to store the number of galaxies that matched
        annular & joined galaxy catalogs
        '''

        nmatch_dict = {'n_exp': [],
                           'len_sexcat': [],
                           'joined_match': [],
                           'annular_match': []
                           }

        self.num_matches = nmatch_dict


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

        assert len(cat1[cat1_ind]) == len(cat2[cat2_ind])

        return cat1[cat1_ind], cat2[cat2_ind]


    def load_sex_cat(self, sexcat_name, n_exp, min_snr=None):
        '''
        Probably superfluous, but gotta make sure it exists then load it

        Rework this slightly: assume the same sexcat name but now
        n_exp is used to load the mini_coadds folder.

        '''
        assert type(sexcat_name) is str, 'required argument sexcat_name must be a string'
        assert type(n_exp) in [int, float], 'required argument: need to supply a number of exposures'

        if min_snr == None:
            min_snr = self.min_snr
            print(f'Using min_snr = {min_snr} for SExtractor cat')

        if os.path.exists(sexcat_name) == True:
            sexcat_path = sexcat_name
        else:
            sexcat_path = os.path.join(self.basedir, sexcat_name)

        try:
           sex_cat = Table.read(sexcat_path)

           if min_snr is not None:
               wg = sex_cat['SNR_WIN'] > 5
               sex_cat = sex_cat[wg]

           self.sex_cat = sex_cat
           self.sexcat_path = sexcat_path
           self.num_matches['len_sexcat'].append(len(sex_cat))
           self.num_matches['n_exp'].append(n_exp)

        except OSError:
            raise f'No sextractor catalog {sexcat_name} found'

        return


    def match_to_analysis_cats(self):
        '''
        Match analysis objects to mini-coadd catalogs, save number of matched objects
        '''

        sex_cat = self.sex_cat
        joined_cat = self.joined_cat
        annular_cat = self.annular_cat
        sexcat_path = self.sexcat_path

        if self.sex_cat is None:
            raise AttributeError('No sextractor catalog loaded')

        # First, match joined
        matched_sex_cat1, matched_joined_cat = self.match_coords(joined_cat, sex_cat)
        self.num_matches['joined_match'].append(len(matched_joined_cat))

        # Also save the joined catalog
        mini_joined_cat_name = sexcat_path.replace('mock_coadd_cat', 'gals_joined_match')
        self.joined_cat[matched_joined_cat].write(mini_joined_cat_name)

        # Now match annular
        matched_sex_cat1, matched_annular_cat = self.match_coords(annular_cat, sex_cat)
        self.num_matches['annular_match'].append(len(matched_annular_cat))

        # Also save the annular catalog
        mini_annular_cat_name = sexcat_path.replace('mock_coadd_cat', 'annular_match')
        self.annular_cat[matched_annular_cat].write(mini_joined_cat_name)

        return


    def save_num_matches(self, outname=None):
        '''
        Save that num_matches dict to file after
        we are done looping
        '''
        basedir = self.basedir

        if outname == None:
            outname = 'analysis_cat_matches.csv'

        assert len(self.num_matches['n_exp']) == len(self.num_matches['annular_match']), \
          'output dict has uneven column lengths, oh dear'

        dd = Table(self.num_matches)

        outfile_path = os.path.join(basedir, outname)
        dd.write(outfile_path, format='ascii.csv')

        return


    def run(self, sexcat_name):
        '''
        Run the matching for this sextractor catalog
        '''

        # Load sextractor catalog
        self.load_sex_cat(sexcat_name)

        # Match to joined catalog
        self.match_to_analysis_cats()
