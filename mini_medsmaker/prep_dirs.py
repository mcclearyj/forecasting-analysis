import os
from argparse import ArgumentParser

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('--exposures_per_list', type=int, default=3,  
                        help='The number of exposures per list (default: 3)')
    parser.add_argument('--mock_dir', type=str, default='.',  
                        help='The directory in which to create the mock directories (default: current directory)')
    return parser.parse_args()


def main(args):
    exposures_per_list = args.exposures_per_list

    make_mini_coadds_dirs(args)

    make_nexp_dirs(args)

    print("Directories made!")


def make_mini_coadds_dirs(args):
    os.makedirs(os.path.join(args.mock_dir, 'mini_coadds'), exist_ok=True)


def make_nexp_dirs(args):
    run_dir = os.path.join(args.mock_dir, 'mini_coadds')
    os.makedirs(run_dir, exist_ok=True)
    for n in range(3, 37, args.exposures_per_list):
        dir_name = f"nexp_{str(n)}"
        os.makedirs(os.path.join(run_dir, dir_name), exist_ok=True)
    # Makes directories for lower nexp runs to go to
    os.makedirs(os.path.join(run_dir, 'nexp_1'), exist_ok=True)
    os.makedirs(os.path.join(run_dir, 'nexp_2'), exist_ok=True)
    os.makedirs(os.path.join(run_dir, 'nexp_4'), exist_ok=True)
    os.makedirs(os.path.join(run_dir, 'nexp_5'), exist_ok=True)


if __name__ == '__main__':
    args = parse_args()
    rc = main(args)
