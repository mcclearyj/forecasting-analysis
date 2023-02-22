import os,sys
from pathlib import Path
import glob
import esutil as eu
import meds
from argparse import ArgumentParser
import superbit_lensing.utils as utils
from superbit_lensing.medsmaker.superbit import medsmaker_mocks as medsmaker
import json

import ipdb

# Kept mock_dir so that it would know where to look for the files.
# Initially replaced mock_dir with exposure_list but I think they should be separate instead
# because the exposure_list object is not a path, so it will look in mock_dir for the exposure_list
def parse_args():

    parser = ArgumentParser()

    parser.add_argument('--mock_dir', type=str,
                        help='Directory containing mock data')
    parser.add_argument('--exposure_list', type=list,
                        help='List containing which exposures to make coadds from')
    parser.add_argument('--outfile', type=str,
                        help='Name of output MEDS file')
    parser.add_argument('-outdir', type=str, default=None,
                        help='Output directory for MEDS file')
    parser.add_argument('-fname_base', action='store', type=str, default=None,
                        help='Basename of mock image files')
    parser.add_argument('-run_name', action='store', type=str, default=None,
                        help='Name of mock simulation run')
    parser.add_argument('--meds_coadd', action='store_true', default=False,
                        help='Set to keep coadd cutout in MEDS file')
    parser.add_argument('--overwrite', action='store_true', default=False,
                        help='Set to overwrite files')
    parser.add_argument('--source_select', action='store_true', default=False,
                        help='Set to select sources during MEDS creation')
    parser.add_argument('--vb', action='store_true', default=False,
                        help='Verbosity')

    return parser.parse_args()

# Decided to make two different files, this file is for getting a coadd from your set list of exposures
# Created new object called "exposure_list"
# This argument will be called in the other file that will loop through the 12 lists

def main(args):
    mock_dir = args.mock_dir
    exposure_list = args.exposure_list
    outfile = args.outfile
    outdir = args.outdir
    run_name = args.run_name
    use_coadd = args.meds_coadd
    overwrite = args.overwrite
    source_selection = args.source_select
    vb = args.vb

# This might be unnecessary but I'm not sure, probably will keep it
    if args.outdir is None:
        outdir = mock_dir

# Log for medsmaker
    logfile = 'mini_mocks.log'
    logdir = outdir
    log = utils.setup_logger(logfile, logdir=logdir)
    logprint = utils.LogPrint(log, vb=vb)

# I do not know what this does!
    if args.fname_base is None:
        fname_base = run_name
    else:
        fname_base = args.fname_base

# Instead of using glob.glob to search for the files,
# the science object will just be the exposure list
    exposure_list_string = "".join(exposure_list)
    exposure_list = exposure_list_string.split()
    logprint(f'Science frames: {exposure_list}')
# Shows the path of the outdir
    outfile = os.path.join(outdir, outfile)

    logprint('Setting up configuration...')
    bm = medsmaker.BITMeasurement(
        image_files=exposure_list, data_dir=mock_dir, run_name=run_name,
        log=log, vb=vb
        )

    bm.set_working_dir(path=outdir)
    mask_dir = os.path.join(mock_dir,'mask_files')
    weight_dir = os.path.join(mock_dir,'weight_files')
    bm.set_mask(
        mask_name='forecast_mask.fits', mask_dir=mask_dir
        )
    bm.set_weight(
        weight_name='forecast_weight.fits', weight_dir=weight_dir
        )

    # Combine images, make a catalog.
    logprint('Making coadd & its catalog...')
    bm.make_coadd_catalog(source_selection=source_selection)

    logprint(f'Writing to {outfile}')

    logprint('Done!')

    return 0

if __name__ == '__main__':
    args = parse_args()
    rc = main(args)

    if rc !=0:
        print(f'mini_mocks failed w/ return code {rc}!')
