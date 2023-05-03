import os
from argparse import ArgumentParser
import subprocess
import json
import ipdb
import pandas as pd
import glob
from astropy.table import Table, vstack
from superbit_class import SuperBIT
from scipy import stats
import numpy as np

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('-mock_dir', type=str,
                        help='Directory containing mock data')
    parser.add_argument('-band_name', type=str,
                        help='Filter [u, g, b, lum, shape]')
    parser.add_argument('-outfile', type=str,
                        help='Name of output MEDS file')
    parser.add_argument('-outdir', type=str, default=None,
                        help='Output directory for MEDS file')
    parser.add_argument('-run_name', action='store', type=str, default=None,
                        help='Name of mock simulation run')
    parser.add_argument('--overwrite', action='store_true', default=False,
                        help='Set to overwrite files')

    return parser.parse_args()

# note that exposure_lists is a list of lists, or basically a list of "exposure_list"

def main(args):
    mock_dir = args.mock_dir
    outfile = args.outfile
    outdir = args.outdir
    run_name = args.run_name
    band_name = args.band_name

    # This also loads the joined and annular catalogs -- note basedir is assumed
    # to be a realization-level thing where annular and joined catalogs live
    # run_name should be forecast_blue or forecast_uv

    exp_list = list(np.arange(3, 37, 3))
    exp_list.extend([1,2,4,5])
    exp_list.sort()

    f = open(os.path.join(mock_dir, f'minimocks_{run_name}_depths.csv'), 'a', encoding="utf-8")
    f.write('# run_name, n_exp, n_obj, median_abmag, mean_abmag, std_abmag, sn10_abmag, argmax_abmag\n')


    for n_exp in exp_list:

        print(f'working on minimocks n_exp={n_exp}')

        minimocks_dir = f'r*/mini_coadds/nexp_{n_exp}'
        outdir = os.path.join(mock_dir, minimocks_dir)
        joined_files = glob.glob(os.path.join(outdir,
                                    f'{run_name}_gals_joined_match.ldac')
                                    )
        print(f'Found: {joined_files}')

        concatenated_tab = Table()

        # Loop over the files and append their data to the concatenated table
        for joined in joined_files:
            table = Table.read(joined)
            table.remove_column('VIGNET')
            concatenated_tab = vstack([concatenated_tab, table])

        # Now convert flux to AB mag
        sbit = SuperBIT()
        ab_mags = sbit.flux_to_abmag(concatenated_tab, band_name, colname='FLUX_AUTO')
        concatenated_tab.add_column(ab_mags, name='ab_mag')

        # Write out the concatenated table to a new file
        outfile_name = f'{run_name}_minimocks_nexp_{n_exp}_match2joinedgals.fits'
        concatenated_tab.write(os.path.join(mock_dir, outfile_name),
                                    format='fits', overwrite=True
                                    )
        print(f'Saved concatenated table to {outfile_name}')


        ##
        ## Get some summary statistics
        ##

        concatenated_tab = concatenated_tab[~np.isnan(concatenated_tab['ab_mag'])]
        len_tab = len(concatenated_tab)

        ab_mag = concatenated_tab['ab_mag']

        sn = concatenated_tab['FLUX_AUTO']/concatenated_tab['FLUXERR_AUTO']

        wg_sn10 = (sn > 9.8) & (sn < 10.2)

        sn10_abmag = np.median(ab_mag[wg_sn10])

        med_ab = np.median(ab_mag)
        mean_ab = np.mean(ab_mag)
        std_ab = np.std(ab_mag)

        depth_b = (stats.mode(ab_mag)).mode[0]

        state = f'{run_name}, {n_exp}, {len_tab}, {med_ab:.2f}, {mean_ab:.2f}, {std_ab:.2f}, {sn10_abmag:.2f}, {depth_b:.2f}\n'
        print(state)
        f.write(state)

    f.close()

    return

if __name__ == '__main__':
    args = parse_args()
    rc = main(args)

    if rc !=0:
        print(f'mini_medsmaker failed w/ return code {rc}!')
