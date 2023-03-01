import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, vstack, hstack
import sys
import os
import glob
from argparse import ArgumentParser
import pdb


parser = ArgumentParser()

parser.add_argument('-shear_tables',type=str,default=None,
                    help = 'shear tables to read in: xxx_shear_profile_cat.fits')
parser.add_argument('-truth_tables',type=str,default=None,
                    help = 'truth tables to read in: xxx_truth.fits.fits')
parser.add_argument('-outname', type=str, default=None,
                    help = 'output table name')

def make_alpha_tab(alpha_arr, sig_alpha_arr, weight_arr):

    alpha_arr = np.array(alpha_arr)
    sig_alpha_arr = np.array(sig_alpha_arr)
    weight_arr = np.array(weight_arr)

    # Normalize weight array if needed
    w_sum = np.sum(weight_arr)
    if (w_sum != 1):
        norm_weight = weight_arr/w_sum
    else:
        norm_weight = weight_arr

    # Get alpha stats
    mean_alpha = np.mean(alpha_arr)
    std_alpha = np.std(alpha_arr)/np.sqrt(len(alpha_arr))
    weighted_mean_alpha = np.sum(alpha_arr*norm_weight)/np.sum(norm_weight)

    # Make table
    alpha_tab = Table()
    alpha_tab.add_columns([alpha_arr, sig_alpha_arr, norm_weight],\
                          names = ['alpha', 'sig_alpha', 'normed_weight'])

    # Add mean and standard deviation of alpha to metadata of table
    alpha_tab.meta['mean_alpha'] = mean_alpha
    alpha_tab.meta['std_alpha'] = std_alpha
    alpha_tab.meta['wmean_alpha'] = weighted_mean_alpha

    print(f'\nmean alpha is {mean_alpha:.4f} +/- {std_alpha:.4f}')
    print(f'weighted mean alpha = {weighted_mean_alpha}')

    print(f'\n\n\n a = {list(alpha_arr)}')
    print(f'\n\n sig = {list(sig_alpha_arr)}\n')

    return alpha_tab 

def main():
    
    args = parser.parse_args()
    shear_tables = args.shear_tables
    truth_tables = args.truth_tables
    outname = args.outname 

    if shear_tables is None:
        shear_tables = 'r*/*shear_profile_cat.fits'    
    if truth_tables is None:
        truth_tables = 'r*/*truth.fits'
    if outname is None:
        outname = 'alpha_table.fits'

    glob_tables = glob.glob(shear_tables)
    glob_tables.sort()
    print(f'reading tables {glob_tables}')

    glob_truths = glob.glob(truth_tables)
    glob_truths.sort()
    print(f'reading tables {glob_truths}')


    alpha_arr = []
    sig_alpha_arr = []
    weight_arr = []

    for i, tabn in enumerate(glob_tables):
        #truth_name = tabn.replace('shear_profile_cat', 'truth')
        truth_name = glob_truths[i]
        truth_cat = Table.read(truth_name, format='fits')
        nstars = len((truth_cat['obj_class'] == 'star').nonzero()[0])
        weight_arr.append(nstars)

        shear_tab = Table.read(tabn, format='fits')

        try:
            this_alpha = shear_tab.meta['ALPHA']
            #this_sig_alpha = shear_tab.meta['SIG_ALPHA']
            this_sig_alpha = shear_tab.meta['sig_alpha']
            alpha_arr.append(this_alpha)
            sig_alpha_arr.append(this_sig_alpha)
            print(f'{tabn} has shear bias {this_alpha:.4f} +/- {this_sig_alpha:.4f} + {nstars} stars')

        except KeyError as e:
            raise e

    alpha_tab = make_alpha_tab(alpha_arr, sig_alpha_arr, weight_arr)
    alpha_tab.write(outname, format='fits', overwrite=True)

    return 0

if __name__ == "__main__":

    rc = main()

    if rc !=0:
        raise Exception
