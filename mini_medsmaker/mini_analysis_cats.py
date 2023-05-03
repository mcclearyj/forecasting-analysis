import os
from argparse import ArgumentParser
import subprocess
import json
import ipdb
import pandas as pd
from mini_matcher import MiniMatcher
from mini_matcher_catsaver import MiniMatcher

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('--mock_dir', type=str,
                        help='Directory containing mock data')
    parser.add_argument('--outfile', type=str,
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
    min_snr = args.min_snr
    vb = args.vb

    # This also loads the joined and annular catalogs -- note basedir is assumed
    # to be a realization-level thing where annular and joined catalogs live
    # run_name should be forecast_blue or forecast_uv

    mini_matcher = MiniMatcher(basedir=mock_dir, run_name=run_name)
    exp_list = list(np.arange(3, 37, 3))
    exp_list.extend([1,2,4,5])

    for n_exp in exp_list:

        minimocks_dir = f'mini_coadds/nexp_{n_exp}'
        outdir = os.path.join(mock_dir, minimocks_dir)
        sexcat_file = os.path.join(outdir, f'{run_name}_mock_coadd_cat.ldac')

        mini_matcher.load_sex_cat(sexcat_file, n_exp)

        # Now, do the matching against joined and annular catalogs
        mini_matcher.match_to_analysis_cats()

    return 

if __name__ == '__main__':
    args = parse_args()
    rc = main(args)

    if rc !=0:
        print(f'mini_medsmaker failed w/ return code {rc}!')
