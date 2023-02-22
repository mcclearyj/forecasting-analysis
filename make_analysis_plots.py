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
    parser.add_argument('--lumfunc_only', action='store_true', default=False,
                            help='Make a luminosity function only [default: False]')    
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
        

            wg = (joined_gals_band['snr_win'] > min_snr) #& (joined_gals_band['redshift'] > min_z)

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
            
            selectgals_pd = self._fits_to_panda(joined_gals_master)
            joinedgals_pd = self._fits_to_panda(select_gals_master)

            self.joinedgals_pd = joinedgals_pd
            self.selectgals_pd = selectgals_pd
            
        '''
        
        self._prep_pd_files(min_snr=min_snr)
        
        return

    
    def save_mean_redshifts(self, zcut, catalog=None, catname=None):
        path = self.path
        f = open(os.path.join(path, 'mean_redshifts.csv'), 'a', encoding="utf-8")
        f.write('# cluster_name, catalog, band, median_z, mean_z, std_z, n_obj\n')
        
        if catalog is not None:
            
            if zcut == True:
                catalog = catalog[(catalog['redshift']>float(self.redshift))]  # prob. terrible coding practice
                
            self._save_mean_redshift(catalog, catname, f)
                
        else:
            joined = self.joinedgals_pd
            select = self.selectgals_pd

            if zcut == True:
                joined = joined[(joined['redshift']>float(self.redshift))]  # prob. terrible coding practice

            self._save_mean_redshift(joined, 'joined_gals', f)
            self._save_mean_redshift(select, 'annular_gals', f)        

        f.close()

        return

    def _save_mean_redshift(self, catalog, catname, f):

        cluster_name = self.cluster_name

        for band in np.unique(catalog['Filter']):
            wg = catalog['Filter']==band
            med_z = np.median(catalog['redshift'][wg])
            mean_z = np.mean(catalog['redshift'][wg])
            std_z = np.std(catalog['redshift'][wg])
            len_z = len(catalog['redshift'][wg])
            state = f'{cluster_name}, {catname}, {band}, {med_z:.2f}, {mean_z:.2f}, {std_z:.2f}, {len_z}\n'
            f.write(state)
   
        
    def _set_palette(self):

        palette = {}

        for band in self.bands:
            if (band == 'u'):
                palette['uv'] = 'magenta'
            elif (band in ['b', 'blue']):
                palette['blue'] = 'C0'
            elif (band == 'g'):
                palette['g'] ='C2'
            elif (band == 'lum'):
                palette['lum'] = 'C1'
            elif (band == 'shape'):
                palette['shape'] = 'C3'
            else:
                raise KeyError(f'No palette color defined for band {band}')
        
        return palette

    
    def _make_kdeplot(self, palette, key, zcut, xlabel=None, outname='distplot.png'):
        '''
        Make a kernel density estimator-type density plot. 
        Left panel: all joined galaxies. Right panel: lensing-selection galaxies
        '''

        min_snr = self.min_snr
        joined = self.joinedgals_pd
        select = self.selectgals_pd
        z = self.redshift
        
        if xlabel == None:
            xlabel=val
            
        if zcut == True:
            joined = joined[(joined['redshift']>float(z))]  # prob. terrible coding practice
            join_title = f'All galaxies with S/N $>$ {min_snr} \& z $>$ {z}'
        else:
            join_title = f'All galaxies with S/N $>${min_snr}'
      

        fig, ax = plt.subplots(1,2, tight_layout=True, figsize=(12,5), sharey=True)
        
        sns.kdeplot(joined, x=key, hue="Filter",
                        ax=ax[0], multiple="layer", fill=False,
                        palette=palette, clip=[0,4],
                        lw=2, common_norm=False, bw_adjust=1.5
                        )

        ax[0].set_xlabel(xlabel, fontsize=16)
        ax[0].set_title(join_title, fontsize=16)

        for band in self.bands:                    
            this_band = joined[joined['Filter'] == band]
            ax[0].axvline(np.median(this_band[key]),
                              color=palette[band], lw=2, ls='--')

        sns.kdeplot(select, x=key, hue="Filter",
                        ax=ax[1], multiple="layer", fill=False,
                        palette=palette, clip=[0,4],
                        lw=2, common_norm=False, bw_adjust=1.5
                        )

        ax[1].set_xlabel(xlabel, fontsize=16)
        ax[1].set_title('Lensing sample', fontsize=16)

        for band in self.bands:
            this_band = select[select['Filter'] == band]
            ax[1].axvline(np.median(this_band[key]),
                              color=palette[band], lw=2, ls='--')

        fig.savefig(os.path.join(self.path, outname))
    
        return
    

    def _make_histplot(self, palette, key, zcut, xlabel=None, outname='histplot.png'):
        '''
        Make a histogram.
        Left panel: all joined galaxies. Right panel: lensing-selection galaxies
        '''
        
        min_snr = self.min_snr
        joined = self.joinedgals_pd
        select = self.selectgals_pd
        z = self.redshift
        
        if xlabel == None:
            xlabel=val
            
        if zcut == True:
            joined = joined[(joined['redshift']>float(z))]  # prob. terrible coding practice
            join_title = f'All galaxies with S/N $>$ {min_snr} \& z $>$ {z}'
        else:
            join_title = f'All galaxies with S/N $>${min_snr}'
           
        fig,ax = plt.subplots(1,2, tight_layout=True, figsize=(12,5), sharey=True)


        sns.histplot(joined, x=key, hue="Filter", element="step", \
                         bins=30, stat="probability", common_norm=False, \
                         log_scale=[False, False], binrange=[0,4],
                         fill=True, ax=ax[0], palette=palette, multiple="layer"
                         )
        
        ax[0].set_xlabel(key, fontsize=16)
        ax[0].set_title(join_title, fontsize=16)

        for band in self.bands:                    
            this_band = joined[joined['Filter'] == band]
            ax[0].axvline(np.median(this_band[key]),
                              color=palette[band], lw=2, ls='--')

        sns.histplot(select, x=key, hue="Filter", element="step", \
                         bins=30, stat="probability", common_norm=False, \
                         log_scale=[False, False], binrange=[0,4],
                         fill=True, ax=ax[1], palette=palette, multiple="layer"
                         )
 
        ax[1].set_xlabel(key, fontsize=16)
        ax[1].set_title('Lensing sample', fontsize=16)

        for band in self.bands:
            this_band = select[select['Filter'] == band]
            ax[1].axvline(np.median(this_band[key]),
                              color=palette[band], lw=2, ls='--')

        fig.savefig(os.path.join(self.path, outname))
 
        return

    
    def make_redshift_plot(self, zcut, distype="kde", key="redshift"):
        '''
        First, define a palette for plots, then make plots
        '''
        
        palette = self._set_palette()
       
        xlabel = 'Redshift'
        
        if distype == "kde":
            outname = f'{self.cluster_name}_z_dists.pdf'
            self._make_kdeplot(palette=palette, key=key,
                                   zcut=zcut, xlabel=xlabel, outname=outname)
                     
        elif distype == "hist":
            outname = f'{self.cluster_name}_z_hists.pdf'
            self._make_histplot(palette=palette, key=key,
                                    zcut=zcut, xlabel=xlabel, outname=outname)
            
        else:
            raise AssertionError('type must be either "kde" or "hist"')
        
        return
    

    def make_single_redshift_plot(self, key="redshift"):
        
        min_snr = self.min_snr
        joined = self.joinedgals_pd

        palette = self._set_palette()
       
        join_title = f'All galaxies with S/N $>${min_snr}'
        xlabel='Redshift'

        fig, ax = plt.subplots(1,1, tight_layout=True, figsize=(7, 5))
        
        sns.kdeplot(joined, x=key, hue="Filter",
                        ax=ax, multiple="layer", fill=False,
                        palette=palette, clip=[0,4],
                        lw=2, common_norm=False, bw_adjust=1.5
                        )

        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_title(join_title, fontsize=16)

        for band in self.bands:                    
            this_band = joined[joined['Filter'] == band]
            ax.axvline(np.median(this_band[key]),
                              color=palette[band], lw=2, ls='--')

            
        outname = f'{self.cluster_name}_allgals_zdist.pdf'
        fig.savefig(os.path.join(self.path, outname))
    
        return

    def _make_lum_func(self, catalog, palette, zcut, outname='lumfunc.pdf'):
        '''
        '''

        cluster_name = self.cluster_name
        path = self.path
        
        if zcut == True:
            print(f'Filtering on redshifts {self.redshift}')
            len1 = len(catalog)
            
            catalog = catalog[(catalog['redshift']>float(self.redshift))]  # prob. terrible coding practice
            len2 = len(catalog)

            print(f'Filtered out {len1} - {len2} objects')
            
        fig, ax = plt.subplots(1, 1, figsize=[7.5, 5.], tight_layout=True)
        bins = np.linspace(18.5, 29, 110)
        
        sn = catalog['flux_auto']/catalog['fluxerr_auto']
        sn_bins = np.linspace(0, 500, 1000)
        sn_ind = np.digitize(sn, sn_bins)
        
        f = open(os.path.join(path, 'photometric_depths.csv'), 'a', encoding="utf-8")
        f.write('# cluster_name, min_z, band, median_abmag, sn10_abmag, argmax_abmag\n')

        for bandpass in np.unique(catalog.Filter):

            color=palette[bandpass]
            bandcat = catalog[catalog.Filter == bandpass]
            
            sn = bandcat['flux_auto']/bandcat['fluxerr_auto'] 
            wg_sn10 = (sn > 9.8) & (sn < 10.2)
            
            ab_mag = bandcat.ab_mag[~np.isnan(bandcat.ab_mag)]

            sn10_abmag = np.median(ab_mag[wg_sn10])
            
            n_b, bins_b, _= ax.hist(ab_mag, bins=bins, log=False, density=True,
                                        color=color, histtype='step', label=bandpass, lw=2)

            depth_b = bins_b[n_b.argmax()]
            print(f'max depth b is {depth_b}; SN = {np.median(sn[wg_sn10])} median depth = {sn10_abmag}')
            
            plt.axvline(sn10_abmag, ls='--', color=color)

            # Do the output writing
            
            min_z = np.min(bandcat.redshift)
            med_ab = np.median(ab_mag)
            mean_ab = np.mean(ab_mag)
            
            state = f'{cluster_name}, {min_z:.2f}, {bandpass}, {med_ab:.2f}, {sn10_abmag:.2f}, {depth_b:.2f}\n'
            f.write(state)

        f.close()
                
        ax.set_xlabel(f'ABmag (S/N $>$ {self.min_snr})')
        ax.set_ylabel('Probability Density')
        #ax.set_xlim([18, 28])
        ax.legend(fontsize=14, loc='upper left')
        
            
        fig.savefig(os.path.join(self.path, outname))
 
        return

   
    def _classic_lum_func(self, gals):
        '''
        ##############################################################################
        ########### Calculate depth of an observation -- DO NOT CHANGE ###############
        ##############################################################################
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
        
        
        return

    
    def make_luminosity_function(self, zcut, catalog=None):
        '''
        Make luminosity functions.
        '''

        palette = self._set_palette()

        if catalog is not None:

            self._make_lum_func(catalog, palette=palette, zcut=zcut,
                                    outname=f'{self.cluster_name}_lumfunc.pdf')
        else:
            
            self._make_lum_func(catalog=self.joinedgals_pd, palette=palette, zcut=zcut,
                                    outname=f'{self.cluster_name}_joined_lumfunc.pdf')
            self._make_lum_func(catalog=self.selectgals_pd, palette=palette, zcut=zcut,
                                    outname=f'{self.cluster_name}_select_lumfunc.pdf')
        
        return

    
    def run(self, zcut, overwrite=False, distype="kde"):

        # Prep the catalogs
        self.prep_pd_files(overwrite=overwrite)
        
        # Save outputs
        self.save_mean_redshifts(zcut=zcut)

        # Make distribution plots
        self.make_redshift_plot(distype=distype, zcut=zcut)

        # Make luminosity functions -- only histos
        self.make_luminosity_function(zcut=zcut)
        
        return 0

    def run_all_gal_zs(self, overwrite=False, distype="kde"):

        # Prep the catalogs
        if self.joinedgals_pd is None:
            self.prep_pd_files(overwrite=overwrite)
        
        # Save outputs
        self.save_mean_redshifts(catalog=self.joinedgals_pd, catname='All', zcut=False)

        # Make distribution plots
        self.make_single_redshift_plot()

        # Make luminosity functions -- only histos
        self.make_luminosity_function(catalog=self.joinedgals_pd, zcut=False)
        
        return 0

    
def make_latex_table(path, zcut=False):
    '''
    Take that little CSV table and turn it into a latex tab :) 
    '''
    
    depths = Table.read(os.path.join(path, 'photometric_depths.csv'), format='ascii.csv')
    median_z = Table.read(os.path.join(path, 'mean_redshifts.csv'), format='ascii.csv')    

    # Reformat a little
    median_z.remove_columns(['n_obj', 'std_z'])
    wb = (median_z['band']=='band')
    median_z.remove_rows(wb)

    wb = (depths['band']=='band')
    depths.remove_rows(wb)

    # There must be a way to do a regexp find/replace with strings in Python
    depths.write(os.path.join(path, 'photometric_depths_latex.tab'),
                     format='latex', overwrite=True)
    median_z.write(os.path.join(path, 'mean_redshifts_latex.tab'),
                       format='latex', overwrite=True)

    tex_tab = Table.copy(median_z)
    tex_tab.add_column(depths['argmax_abmag'], name=r'$\argmax({\rm PDF})$ depth')
    tex_tab.add_column(depths['sn10_abmag'], name=r'$SN \sim 10$ depth')

    cluster_names = np.zeros(len(median_z),dtype='<U30')
    galaxy_sample_names = np.zeros(len(median_z),dtype='<U30')
    
    for unique in np.unique(median_z['# cluster_name']):
        here = median_z['# cluster_name'] == unique
        cluster_names[here] = r'\texttt{%s}' % unique
        if (zcut == True):
            redshift = unique.split('z')[1]
            joined_label = r'$z > {redshift}$'.format(redshift=redshift)
        else:
            joined_label = r'All $z$'
        galaxy_sample_names[here] = joined_label

    sel = tex_tab['catalog']=='annular_gals'; galaxy_sample_names[sel]='Lensing'
    zal = tex_tab['catalog']=='All'; galaxy_sample_names[zal]=r'All $z$'
    
    namecol = Table.Column(cluster_names, name='Cluster name')
    catcol = Table.Column(galaxy_sample_names, name='Galaxy sample')

    tex_tab.replace_column('# cluster_name', namecol)
    tex_tab['# cluster_name'].name = 'Cluster name'
    tex_tab.replace_column('catalog', catcol)
    tex_tab['catalog'].name = 'Galaxy sample'

    tex_tab.write(os.path.join(path, 'combined_results.csv'), format='ascii.csv', overwrite=True)
    tex_tab.write(os.path.join(path, 'combined_results_latex.tab'), format='latex', overwrite=True)

    cmd = r"sed -i.bak -e 's/_/\\_/g' %s/combined_results_latex.tab" % path
    print(cmd)
    os.system(cmd)

    
def main(args):
    
    path = args.path
    masses = args.masses
    redshifts = args.redshifts
    bands = args.bands
    lumfunc_only = args.lumfunc_only
    distype = str(args.distype)
    overwrite = args.overwrite

    if bands == None:
        bands = ['blue', 'lum', 'shape']
    if redshifts == None:
        redshifts = ['0.059', '0.3', '0.45']
        #redshifts = ['0.45']
    if masses == None:  
        masses = ['m4.1e14']

    #if overwrite == True:


    zcut = True
    
    try:
        os.remove(os.path.join(path, 'mean_redshifts.csv'))
        os.remove(os.path.join(path, 'photometric_depths.csv'))
        os.remove(os.path.join(path, 'combined_results_latex.tab'))
        os.remove(os.path.join(path, 'combined_results.csv'))

    except FileNotFoundError:
        pass 
    
    for mass in masses:
        
        for z in redshifts:
            
            cluster_name = f'cl_{mass}_z{z}'
            
            print(f'\nWorking on {cluster_name}\n')
 
            cluster = ClusterCats(path=path, cluster_name=cluster_name,
                                      redshift=z, bands=bands)

            cluster.run(overwrite=overwrite, zcut=zcut, distype=distype)
            

    # Then also do it once with no z cut at to show fg distribution
    
    single = cluster.run_all_gal_zs()
    
    # For convenience!
    
    make_latex_table(path, zcut=zcut)
    
    return 0


if __name__ == '__main__':

    args = parse_args()
    
    rc = main(args)
        
    if rc == 0:
        print('make_analysis_plots.py has completed succesfully')
    else:
        print(f'make_analysis_plots.py has failed w/ rc={rc}')

