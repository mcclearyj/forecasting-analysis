import matplotlib.pyplot as plt
import numpy as np
import astropy.table
from astropy.table import Table, hstack, vstack
from astropy.io import fits
import os, re
import ipdb
from argparse import ArgumentParser
from superbit_class import SuperBIT


from matplotlib import ticker, rc
plt.style.use('default')
rc('font',**{'family':'serif'})
rc('text', usetex=True)


import seaborn as sns
sns.set_theme(font="serif", style="darkgrid", font_scale=1.2)

###
###

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('path', default=None,
                            help = 'Path to histogram outputs')
    parser.add_argument('-masses', default=None,
                            help = 'Cluster masses')
    parser.add_argument('-redshifts', default=None,
                            help = 'Cluster_redshifts [default: 0.059, 0.3, 0.45]')
    parser.add_argument('-bands', default=None,
                            help = 'Bandpasses for redshift histograms [default: "blue", "lum", "shape"]')
    parser.add_argument('-distype', default="kde",
                            help='Probability distribution type: "hist" or "kde" [default: "kde"]')
    parser.add_argument('--overwrite', action='store_true', default=False,
                            help='Overwrite existing joined/annular master galaxy catalogs [default: False]')    
    return parser.parse_args()


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
    plt.rcParams.update({'xtick.major.visible': True})
    plt.rcParams.update({'xtick.minor.visible': True})
    plt.rcParams.update({'xtick.minor.width': 1.})
    plt.rcParams.update({'xtick.minor.size': 6})
    plt.rcParams.update({'xtick.direction': 'out'})
    plt.rcParams.update({'xtick.major.visible': True})
    plt.rcParams.update({'ytick.major.width': 1.3})
    plt.rcParams.update({'ytick.major.size': 8})
    plt.rcParams.update({'ytick.minor.visible': True})
    plt.rcParams.update({'ytick.minor.width': 1.})
    plt.rcParams.update({'ytick.minor.size':6})
    plt.rcParams.update({'ytick.direction':'out'})

    return


class ClusterCats:

    def __init__(self, path, cluster_name, redshift, basename=None, bands=None):
        '''
        Make the redshift histograms. For every bandpass in bands, read in 
        joined_ and _annular (shear) catalogs and put them into a Pandas table
        (Yes, I can do this in Pandas itself, I don't care). Then, plot the redshift
        distributions, probably as KDEs since that's cleaner.  

        Inputs:
              path: location of files
              joinedgals_pds: list of pandas objects for joined gals
              sheargals_pds: list of pandas objects for shear gals
              bands: filters for which histograms are desired
        '''
        
        self.path = path
        self.cluster_name = cluster_name
        self.basename = basename
        self.redshift = redshift
        self.bands = bands
        self.min_snr = 5.0
        self.joinedgals_pd = None
        self.selectgals_pd = None

        if bands == None:
            self.bands = ['blue', 'lum', 'shape']
            print(f'\nusing default bands {self.bands}\n')

        if basename == None:
            self.basename = 'forecast'
            print(f'\nusing basename {self.basename}\n')

            
    def _fits_to_panda(self, name):
        '''
        Utility function to open up a fits file and return a pandas object
        There is this weird thing with Filter being saved as a binary, idk what's 
        going on with that, so have a slightly hacky version to 
        '''

        dat = Table.read(name, format='fits')
        
        return dat.to_pandas()

    
    def _prep_pd_files(self, min_snr):
        '''
        Make Pandas table for all bandpasses desired
        '''
        path = self.path
        cluster_name = self.cluster_name
        min_z = float(self.redshift)
        min_snr = self.min_snr
        
        if min_snr == None:
            min_snr = 5.0
            self.snr_cutouff = min_snr
            
        select_gals_tables = []; joined_gals_tables = []
        
        for band in self.bands:
            
            forecast_name = f'{self.basename}_{band}'
            
            joined_gals_name = os.path.join(path, '{cluster_name}_{forecast_name}_gals_joined_master_cat.fits')
            select_gals_name = os.path.join(path, '{cluster_name}_{forecast_name}_annular_gals_master_cat.fits')

            joined_gals_band = Table.read(joined_gals_name.format(cluster_name=cluster_name, forecast_name=forecast_name))
            select_gals_band = Table.read(select_gals_name.format(cluster_name=cluster_name, forecast_name=forecast_name))
        

            wg = (joined_gals_band['snr_win'] > min_snr) & (joined_gals_band['redshift'] > min_z)

            if band == 'u':
                bandname = 'uv'
            else:
                bandname = band

            joined_gals_band.add_column(bandname, name='Filter')
            select_gals_band.add_column(bandname, name='Filter')

            select_gals_tables.append(select_gals_band)
            joined_gals_tables.append(joined_gals_band[wg])
            
        
        all_joined_gals = vstack([tab for tab in joined_gals_tables])
        all_select_gals = vstack([tab for tab in select_gals_tables])

        # Save to file, y not
        outname1 = os.path.join(path, f'{cluster_name}_all_bands_gals_joined_master_cat.fits')
        outname2 = os.path.join(path, f'{cluster_name}_all_bands_annular_gals_master_cat.fits')

        all_joined_gals.write(outname1, format='fits', overwrite=True)
        all_select_gals.write(outname2, format='fits', overwrite=True)

        selectgals_pd = all_select_gals.to_pandas()
        joinedgals_pd = all_joined_gals.to_pandas()

        self.joinedgals_pd = joinedgals_pd
        self.selectgals_pd = selectgals_pd

        return 0


    def prep_pd_files(self,overwrite=False, min_snr=None):
        '''
        Either make new {cluster_name}_all_bands_gals_joined_master_cat.fits and 
        {cluster_name}_all_bands_annular_gals_master_cat.fits or read in existing ones
        '''

        
        path = self.path
        cluster_name = self.cluster_name

        joined_gals_master = os.path.join(path, f'{cluster_name}_all_bands_gals_joined_master_cat.fits')
        select_gals_master = os.path.join(path, f'{cluster_name}_all_bands_annular_gals_master_cat.fits')

        exist_joined = os.path.exists(joined_gals_master)
        exist_select = os.path.exists(select_gals_master)

        '''

        if (exist_joined == False) or (exist_select == False) or (overwrite == True):
            
            self._prep_pd_files(min_snr=min_snr)
            
        else:

            ipdb.set_trace()
            
            selectgals_pd = self._fits_to_panda(joined_gals_master)
            joinedgals_pd = self._fits_to_panda(select_gals_master)

            self.joinedgals_pd = joinedgals_pd
            self.selectgals_pd = selectgals_pd
            
        '''
        
        self._prep_pd_files(min_snr=min_snr)
        
        return

    
    def save_median_redshifts(self):

        joined = self.joinedgals_pd
        select = self.selectgals_pd
        cluster_name = self.cluster_name

        f = open('median_redshifts.csv', 'a', encoding="utf-8")
        
        f.write('# cluster_name, catalog, band, median_z, mean_z, std_z, n_obj\n')

        for band in np.unique(select['Filter']):
            wg = select['Filter']==band
            med_z = np.median(select['redshift'][wg])
            mean_z = np.mean(select['redshift'][wg])
            std_z = np.std(select['redshift'][wg])
            len_z = len(select['redshift'][wg])
            state = f'{cluster_name}, annular_gals, {band}, {med_z:.3f}, {mean_z:.3f}, {std_z:.3f}, {len_z}\n'
            f.write(state)

        for band in np.unique(joined['Filter']):
            wg = joined['Filter']==band
            med_z = np.median(joined['redshift'][wg])
            mean_z = np.mean(joined['redshift'][wg])
            std_z = np.std(joined['redshift'][wg])
            len_z = len(joined['redshift'][wg])
            state = f'{cluster_name}, joined_gals, {band}, {med_z:.3f}, {med_z:.3f}, {std_z:.3f}, {len_z}\n'
            f.write(state)

        f.close()

        return
       
    def _set_palette(self):

        palette = {}

        for band in self.bands:
            if band == 'u':
                palette['uv'] = 'magenta'
            elif (band == 'blue') or (band == 'b'):
                palette['blue'] = 'C0'
            elif band == 'g':
                palette['g'] ='C2'
            elif band == 'lum':
                palette['lum'] = 'C1'
            elif band == 'shape':
                palette['shape'] = 'C3'
            else:
                raise KeyError(f'No palette color defined for band {band}')
        
        return palette

    
    def _make_kdeplot(self, palette):
        '''
        Make a kernel density estimator-type density plot. 
        Left panel: all joined galaxies. Right panel: lensing-selection galaxies
        '''
        
        z = self.redshift
        min_snr = self.min_snr
        joined = self.joinedgals_pd
        select = self.selectgals_pd

        fig, ax = plt.subplots(1,2, tight_layout=True, figsize=(12,5), sharey=True)
        
        sns.kdeplot(joined, x="redshift", hue="Filter",
                        ax=ax[0], multiple="layer", fill=False,
                        palette=palette, clip=[0,4],
                        lw=2, common_norm=False, bw_adjust=1.5,
                        )

        ax[0].set_xlabel('Redshift', fontsize=16)
        ax[0].set_title(f'All galaxies with S/N $>${min_snr} \& z $>${z}', fontsize=16)

        for band in self.bands:                    
            this_band = joined[joined['Filter'] == band]
            ax[0].axvline(np.median(this_band['redshift']),
                              color=palette[band], lw=2, ls='--')

        sns.kdeplot(joined, x="redshift", hue="Filter",
                        ax=ax[1], multiple="layer", fill=False,
                        palette=palette, clip=[0,4],
                        lw=2, common_norm=False, bw_adjust=1.5,
                        )

        ax[1].set_xlabel('Redshift', fontsize=16)
        ax[1].set_title('Lensing sample', fontsize=16)

        for band in self.bands:
            this_band = select[select['Filter'] == band]
            ax[1].axvline(np.median(this_band['redshift']),
                              color=palette[band], lw=2, ls='--')

        fig.savefig(os.path.join(self.path, f'{self.cluster_name}_z_dists.pdf'))
    
        return
    

    def _make_histplot(self, palette):
        '''
        Make a histogram.
        Left panel: all joined galaxies. Right panel: lensing-selection galaxies
        '''
        
        z = self.redshift
        min_snr = self.min_snr
        joined = self.joinedgals_pd
        select = self.selectgals_pd
        
        fig,ax = plt.subplots(1,2, tight_layout=True, figsize=(12,5), sharey=True)


        sns.histplot(joined, x='redshift', hue="Filter", element="step", \
                         bins=30, stat="probability", common_norm=False, \
                         log_scale=[False, False], binrange=[0,4],
                         fill=True, ax=ax[0], palette=palette, multiple="layer"
                         )
        
        ax[0].set_xlabel('Redshift', fontsize=16)
        ax[0].set_title(f'All galaxies with S/N $>${min_snr} \& z $>${z}', fontsize=16)

        for band in self.bands:                    
            this_band = joined[joined['Filter'] == band]
            ax[0].axvline(np.median(this_band['redshift']),
                              color=palette[band], lw=2, ls='--')

        sns.histplot(select, x='redshift', hue="Filter", element="step", \
                         bins=30, stat="probability", common_norm=False, \
                         log_scale=[False, False], binrange=[0,4],
                         fill=True, ax=ax[1], palette=palette, multiple="layer"
                         )
 
        ax[1].set_xlabel('Redshift', fontsize=16)
        ax[1].set_title('Lensing sample', fontsize=16)

        for band in self.bands:
            this_band = select[select['Filter'] == band]
            ax[1].axvline(np.median(this_band['redshift']),
                              color=palette[band], lw=2, ls='--')

        fig.savefig(os.path.join(self.path, f'{self.cluster_name}_z_hists.pdf'))
 
        return

    
    def make_displot(self, distype="kde"):
        '''
        First, define a palette for plots, then make plots
        '''
        
        palette = self._set_palette()
        
        if distype == "kde":
            self._make_kdeplot(palette=palette)
                     
        elif distype == "hist":
            self._make_histplot(palette=palette)
            
        else:
            raise AssertionError('type must be either "kde" or "hist"')
        
        return

    def _classic_lum_func(self):
        '''
        ##############################################################################
        ########### Calculate depth of an observation -- DO NOT CHANGE ###############
        ##############################################################################
        '''

        '''
        # gals is the analysis object catalog
        n,bins=np.histogram(gals['MAG_AUTO'],bins=100)
        midpoints = np.array([(bins[i+1]+bins[i])*0.5 for i in range(len(bins)-1)])

        wg=(midpoints<26.4) & (midpoints>19)

        fit=np.polyfit(midpoints[wg],np.log(n[wg]),1)
        num = fit[0]*midpoints+fit[1]

        plt.hist(gals['MAG_AUTO'],histtype='step',bins=100,label='mock deep 3hr b',log=True)
        plt.plot(midpoints,np.exp(num),'--k')

        # OK, now to estimate 50% completeness
        fraction=np.log(n)/num
        enum=enumerate(fraction)
        l = list(enum)

        # Here you have to pick your point in the resulting l array.
        # In one instance, I used ind=80 for ~100% completeness,
        # used ind=93 for 90% completeness, and np.mean(93,94) for 50% completeness

        complete=midpoints[86]
        complete90=midpoints[90]
        complete50=np.mean([midpoints[97],midpoints[98]])
        '''
        
        return


    def _make_lum_func_plot(self, catalog, outname):
        '''
        Make a fun plot
        '''
        pass
        
        return 
    
    
    def make_luminosity_function(self, catalog=None,
                                     outname='luminosity_function.pdf'):
        '''
        Make luminosity functions.
        '''

        pass

        if catalog is not None:

            self._make_lum_func_plot(catalog=catalog,
                                        outname=f'{self.cluster_name}_lumfunc.pdf')
            
        else:
            
            self._make_lum_func_plot(catalog=self.joinedgals_pd,
                                         outname=f'{self.cluster_name}_joined_lumfunc.pdf')
            self._make_lum_func_plot(catalog=self.selectgals_pd,
                                         outname=f'{self.cluster_name}_select_lumfunc.pdf')
        
            
        return

    
    def run(self, overwrite=False, distype="kde"):

        # Prep the catalogs
        self.prep_pd_files(overwrite=overwrite)

        # Save outputs
        self.save_median_redshifts()

        # Make distribution plots
        self.make_displot(distype=distype)

        # Make luminosity functions
        self.make_luminosity_function()

        return 0

    
def make_latex_table(path):
    '''
    Take that little CSV table and turn it into a latex tab :) 
    '''
    
    median_z = Table.read(os.path.join(path, 'median_redshifts.csv'), format='ascii.csv')

    # Reformat it
    median_z.remove_column('n_obj')

    # There must be a way to do a regexp find/replace with strings in Python
    median_z.write(os.path.join(path, 'median_redshifts_latex.tab'),
                       format='latex')

    # This is cheap but it will do the trick
    head = '# cluster name & catalog & band & median $z$ & mean $z$ & $sigma_z$//'
    cmd = f"sed -i.bak -e 's/_/\\_/g' median_redshifts.csv"
    os.system(cmd)

    
def main(args):
    
    path = args.path
    masses = args.masses
    redshifts = args.redshifts
    bands = args.bands
    distype = str(args.distype)
    overwrite = args.overwrite

    if bands == None:
        bands = ['blue', 'lum', 'shape']
    if redshifts == None:
        redshifts = ['0.059', '0.3', '0.45']
    if masses == None:  
        masses = ['m4.1e14']

    #if overwrite == True:
        
    try:
        os.remove(os.path.join(path, 'median_redshifts.csv'))
        
    except FileNotFoundError:
        pass 
    
    for mass in masses:
        
        for z in redshifts:
            
            cluster_name = f'cl_{mass}_z{z}'
            
            print(f'\nWorking on {cluster_name}\n')
 
            cluster = ClusterCats(path=path, cluster_name=cluster_name,
                                      redshift=z, bands=bands)

            cluster.run(overwrite=overwrite, distype=distype)
            
    # For convenience!
    
    make_latex_table(path)
    
    return 0


if __name__ == '__main__':

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('make_analysis_plots.py has completed succesfully')
    else:
        print(f'make_analysis_plots.py has failed w/ rc={rc}')

