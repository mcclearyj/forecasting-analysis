import os
from argparse import ArgumentParser
import subprocess
import json
import ipdb
import pandas as pd
from mini_matcher import MiniMatcher


# Current goals done: 
#   - Generated the list of exposures for the mini_mocks.py file to take (number 1)
#   - Create directory just for this mini_medsmaker run inside each r# directory (number 2)
#       - Named "run_name_mini_coadds" and will be made in each r# directory
#   - Created directories that the coadds will end up going in (number 3)
#       - These go inside the "run_name_mini_coadds" directory
#   - Call the mini_mocks.py file within this one and feed it certain arguments (number 4)
#   - Loop it so that it runs for each list out of 12 lists
#   - Decided that it is best to submit in Sbatch, takes too long even with 18 CPUs
#   - Added lower nexp lists to the list (1, 2, 4, 5)
#       - Edited prep_dirs.py to make nexp_{1,2,4,5}
# TO-DO Goals:
#   - 
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
    parser.add_argument('-min_snr', action='store', type=float, default=5,
                            help='Set an SNR_WIN threshold for sextractor catalog')
    parser.add_argument('--vb', action='store_true', default=False,
                        help='Verbosity')
    parser.add_argument('--exposures_per_list', type=int, default=3,  
                        help='The number of exposures per list (default: 3)')

    return parser.parse_args()

# note that exposure_lists is a list of lists, or basically a list of "exposure_list"

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
    exposures_per_list = args.exposures_per_list

    #1. Generates the list of the exposures
    def generate_exposure_lists(run_name, exposures_per_list=3):
        exposure_lists = []
        for i in range(exposures_per_list, 37, exposures_per_list):
            exposures = []
            for j in range(1, i + 1):
                exposure = mock_dir + "/" + run_name + '_' + str(j).zfill(3) + '.fits'
                exposures.append(exposure)
            exposure_lists.append(exposures)

        return exposure_lists

    exposure_lists = generate_exposure_lists(run_name, exposures_per_list)
    # Makes the lower exposure lists to add to exposure_lists
    new_list_1 = [mock_dir + "/" + run_name + '_' + str(1).zfill(3) + '.fits']
    new_list_2 = [mock_dir + "/" + run_name + '_' + str(1).zfill(3) + '.fits', mock_dir + "/" + run_name + '_' + str(2).zfill(3) + '.fits']
    new_list_4 = [mock_dir + "/" + run_name + '_' + str(1).zfill(3) + '.fits', mock_dir + "/" + run_name + '_' + str(2).zfill(3) + '.fits',
                mock_dir + "/" + run_name + '_' + str(3).zfill(3) + '.fits', mock_dir + "/" + run_name + '_' + str(4).zfill(3) + '.fits']
    new_list_5 = [mock_dir + "/" + run_name + '_' + str(1).zfill(3) + '.fits', mock_dir + "/" + run_name + '_' + str(2).zfill(3) + '.fits',
                 mock_dir + "/" + run_name + '_' + str(3).zfill(3) + '.fits', mock_dir + "/" + run_name + '_' + str(4).zfill(3) + '.fits',
              mock_dir + "/" + run_name + '_' + str(5).zfill(3) + '.fits']
    
    exposure_lists.append(new_list_1)
    exposure_lists.append(new_list_2)
    exposure_lists.append(new_list_4)
    exposure_lists.append(new_list_5)

    # Create an instance of the MiniMatcher object, loading joined & annular cats 
    mini_matcher = MiniMatcher(basedir=mock_dir, run_name=run_name)
    mini_matcher.min_snr = min_snr
    
    #4. Calls the mini_mocks.py file with arguments
    for exposure_list in exposure_lists:
        exposures_per_list = len(exposure_list)
        args.exposure_list = exposure_list
        # Name for each output file will be run_name_coadd_3.fits, etc.
        outfile = str(run_name + "_coadd_" + str(exposures_per_list) + ".fits")
        args.outfile = outfile
        # Outdir parses into each nexp_# directory for each coadd
        outdir = str(mock_dir + "/" + "mini_coadds/" + "nexp_" + str(exposures_per_list))
        args.outdir = outdir
        # Had to convert the exposure list into a string to pass it thorugh subprocess
        # Gets joined again in mini_mocks.py on lines 76 & 77
        exposure_list_string = " ".join(exposure_list)
        subprocess.run(["python", "/work/mccleary_group/vassilakis.g/superbit-metacal/superbit_lensing/medsmaker/scripts/mini_mocks.py",
         "--outfile=" + outfile, "-outdir=" + outdir, "--mock_dir=" + mock_dir, "--exposure_list=" + exposure_list_string,
         "-run_name=" + run_name,
         *filter(lambda x: x is str or x is bytes or x is os.PathLike, vars(args).values())])
         
        # Alrighty, it's time to match. First, load in SExtractor catalog
        real_outfile = os.path.join(outdir, f'{run_name}_mock_coadd_cat.ldac')
        mini_matcher.load_sex_cat(real_outfile, exposures_per_list)
        
        # Now, do the matching against joined and annular catalogs. Be chatty, print to screen
        mini_matcher.match_to_analysis_cats()

    # 5. If all has gone OK, save mini_matcher output to a table. 
    # os.basepath(outdir) should be r{0...29} 
    matched_obj_name = os.path.basename(mock_dir)+ '_mini_coadd_matches.csv'
    mini_matcher.save_num_matches(outname=matched_obj_name)

    #Calculating mini_codd stats
    run_name = args.run_name
    parent_dir = os.path.dirname(mock_dir)
    
    calculate_statistics(run_name, parent_dir)

    outcsv = os.path.join(mock_dir, run_name + '_averages.csv')
    rename_headers(outcsv)        

    return

def calculate_statistics(run_name, parent_dir):
    # Initialize an empty dataframe to store the values from each run
    data = pd.DataFrame(columns=['n_exp', 'len_sexcat', 'joined_match', 'annular_match'])
    
    # Loop through each folder
    for i in range(30):
        # Construct the file path
        file_path = os.path.join(parent_dir, 'r{}'.format(i), 'r{}_mini_coadd_matches.csv'.format(i))
        
        # Read the data from the file
        run_data = pd.read_csv(file_path)
        
        # Append the values to the dataframe
        data = data.append(run_data, ignore_index=True)
        
    # Group the data by the n_exp column
    grouped_data = data.groupby('n_exp')
    
    # Calculate the mean and standard deviation for each group
    mean = grouped_data.mean()
    std = grouped_data.std()
    
    # Add the n_exp column to the mean and std dataframes
    mean['n_exp'] = mean.index
    std['n_exp'] = std.index
    
    # Reset the index of the dataframes so that the n_exp column is the first column
    mean = mean.reset_index(drop=True)
    std = std.reset_index(drop=True)
    
    # Create a new dataframe to store the results
    results = pd.concat([mean, std], axis=1, keys=['mean', 'std'])
    
    # Save the results to a CSV file
    outcsv = os.path.join(parent_dir, run_name + '_averages.csv')
    results.to_csv(outcsv, index=False)

# This function renames the headers of the output csv file
#   - Renamed because above function outputs two header rows, function below joins them

def rename_headers(outcsv):
    df = pd.read_csv(outcsv)
    df.drop(0, axis=0, inplace=True) # Delete first row of the dataframe
    df.rename(columns={df.columns[0]: 'mean_len_sexcat', df.columns[1]: 'mean_joined_match', 
                       df.columns[2]: 'mean_annular_match', df.columns[3]: 'n_exp', 
                       df.columns[4]: 'std_len_sexcat', df.columns[5]: 'std_joined_match', 
                       df.columns[6]: 'std_annular_match', df.columns[7]: 'n_exp'}, inplace=True)
    df.to_csv(outcsv, index=False)

    return 0


if __name__ == '__main__':
    args = parse_args()
    rc = main(args)

    if rc !=0:
        print(f'mini_medsmaker failed w/ return code {rc}!')
