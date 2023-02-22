import pandas as pd
import matplotlib
matplotlib.use('PDF')  # change the backend to PDF
import matplotlib.pyplot as plt
import argparse
import yaml
import ipdb
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str, help='Path to the YAML config file')
    args = parser.parse_args()
    with open(args.config_file, 'r') as f:
        yaml_data = yaml.safe_load(f)
    return yaml_data



def set_rc_params():
    '''
    Set figure parameters
    This should be a config one day
    '''
    plt.rcParams.update({'figure.facecolor':'w'})
    plt.rcParams.update({'axes.linewidth': 1.3})
    plt.rcParams.update({'xtick.labelsize': 16})
    plt.rcParams.update({'ytick.labelsize': 16})
    plt.rcParams.update({'xtick.major.size': 8})
    plt.rcParams.update({'xtick.major.width': 1.3})
    plt.rcParams.update({'xtick.minor.visible': True})
    plt.rcParams.update({'xtick.minor.width': 1.})
    plt.rcParams.update({'xtick.minor.size': 6})
    plt.rcParams.update({'xtick.direction': 'out'})
    plt.rcParams.update({'ytick.major.width': 1.3})
    plt.rcParams.update({'ytick.major.size': 8})
    plt.rcParams.update({'ytick.minor.visible': True})
    plt.rcParams.update({'ytick.minor.width': 1.})
    plt.rcParams.update({'ytick.minor.size':6})
    plt.rcParams.update({'ytick.direction':'out'})
    plt.rcParams.update({'font.family': 'serif'})

    return



def main(args):

    set_rc_params()
    
   # Create the joined plot and axis objects
    fig1, ax_joined = plt.subplots(figsize=(8,4))

    # Call the plot_joined_data function for each CSV file to plot
    for i in range(1, 7):
        if args.get(f'plot_csv{i}', False):
            plot_joined_data(fig1, ax_joined, file_path=args[f'csvpath{i}'], label=args[f'label{i}'], color=args[f'color{i}'],
            joined_xlabel=args['joined_xlabel'], joined_ylabel=args['joined_ylabel'], joined_title=args['joined_title'])

    # Adjust the plot margins to ensure all axes labels are visible
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.9)

    # Create path for pdf export
    pdf_path_j = os.path.join(args['base_dir'], args['pdf_name_j'])

    # Save the joined plot
    fig1.savefig(pdf_path_j)

    print(f"Joined plot printed to {pdf_path_j}")

    # Create the annular plot and axis objects
    fig2, ax_annular = plt.subplots(figsize=(8,4))

    # Call the plot_annular_data function for each CSV file to plot
    for i in range(1, 7):
        if args.get(f'plot_csv{i}', False):
            plot_annular_data(fig2, ax_annular, file_path=args[f'csvpath{i}'], label=args[f'label{i}'], color=args[f'color{i}'],
            annular_xlabel=args['annular_xlabel'], annular_ylabel=args['annular_ylabel'], annular_title=args['annular_title'])

    # Adjust the plot margins to ensure all axes labels are visible
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.9)

    # Create path for pdf export
    pdf_path_a = os.path.join(args['base_dir'], args['pdf_name_a'])

    # Save the joined plot
    fig2.savefig(pdf_path_a)

    print(f"Annular plot printed to {pdf_path_a}")



# Define a function to load and plot data from a CSV file
def plot_joined_data(fig1, ax_joined, file_path, label, color, joined_xlabel, joined_ylabel, joined_title):
    # Load the data from CSV file
    data = pd.read_csv(file_path)

    # Extract the relevant columns
    n_exp = data['n_exp']
    mean_joined_density = data['mean_joined_match']/350
    std_joined_density = data['std_joined_match']/350
    

    # Plot mean_joined_match vs n_exp with line connecting points and error bars
    ax_joined.errorbar(n_exp, mean_joined_density, yerr=std_joined_density, fmt='o', capsize=4, markersize=5, color=color, label=label)
    ax_joined.plot(n_exp, mean_joined_density, '-o', markersize=5, color=color)

    # Set the x and y axis labels and title
    ax_joined.set_xlabel(joined_xlabel)
    ax_joined.set_ylabel(joined_ylabel)
    ax_joined.set_title(joined_title)
    ax_joined.tick_params(axis='both', which='major', labelsize=12)

    # Set the x-axis limits to 0 and 36, with ticks every 3 units
    ax_joined.set_xlim(0, 37)
    ax_joined.set_xticks(range(0, 37, 3))

    # Add a legend to the plot
    ax_joined.legend()

    return fig1, ax_joined



# Define a function to load and plot data from a CSV file
def plot_annular_data(fig2, ax_annular, file_path, label, color, annular_xlabel, annular_ylabel, annular_title):
    # Load the data from CSV file
    data = pd.read_csv(file_path)

    # Extract the relevant columns
    n_exp = data['n_exp']
    mean_annular_density = data['mean_annular_match']/350
    std_annular_density = data['std_annular_match']/350
    
    # Plot mean_annular_match vs n_exp with line connecting points and error bars
    ax_annular.errorbar(n_exp, mean_annular_density, yerr=std_annular_density, fmt='o', capsize=4, markersize=5, color=color, label=label)
    ax_annular.plot(n_exp, mean_annular_density, '-o', markersize=5, color=color)

    # Set the x and y axis labels and title
    ax_annular.set_xlabel(annular_xlabel)
    ax_annular.set_ylabel(annular_ylabel)
    ax_annular.set_title(annular_title)
    ax_annular.tick_params(axis='both', which='major', labelsize=12)

    # Set the x-axis limits to 0 and 36, with ticks every 3 units
    ax_annular.set_xlim(0, 37)
    ax_annular.set_xticks(range(0, 37, 3))

    # Add a legend to the plot
    ax_annular.legend()

    return fig2, ax_annular



if __name__ == '__main__':
    yaml_data = parse_args()

    main(yaml_data)