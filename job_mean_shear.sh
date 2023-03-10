#!/bin/sh
#SBATCH -t 00:59:59
#SBATCH --partition=express
#SBATCH --mem-per-cpu=5GB
#SBATCH --nodes 1
#SBATCH -n 10
#SBATCH -J meanshear
#SBATCH -v 
# Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID
#SBATCH -o meanshear-%J.out
#SBATCH -e meanshear-%J.err

source /work/mccleary_group/miniconda3/etc/profile.d/conda.sh
conda activate sbmcal_139

export BASEDIR=/work/mccleary_group/superbit/mock-data-forecasting
export CODEDIR=/work/mccleary_group/forecasting-analysis
export SBMCAL=/work/mccleary_group/superbit/superbit-metacal

### Uncomment and indent if a loop over clusters is desired too!
#declare -a CLUSTERS=("cl_m4.1e14_z0.059" "cl_m4.1e14_z0.3" "cl_m4.1e14_z0.45" "cl_m7.8e14_z0.3")
#for $CCLUSTER in "${arr[@]}"
export CLUSTER=cl_m4.1e14_z0.45

declare -a BANDS=("blue" "lum" "shape")

for BAND in "${BANDS[@]}"
do
    export RUN_NAME=forecast_$BAND
    export DATADIR=$BASEDIR/$RUN_NAME/$CLUSTER
    export OUTNAME="$DATADIR/$RUN_NAME"_mean_shearprofile_table.fits

    echo "$RUN_NAME"
    echo "$DATADIR"
    echo "$OUTNAME"

    python $CODEDIR/get_mean_shearprofile.py -shear_cats="$DATADIR"/r*/"$RUN_NAME"_shear_profile_cat.fits  -annular_cats="$DATADIR"/r*/"$RUN_NAME"_annular.fits -stackcat_name="$OUTNAME"

    python $SBMCAL/superbit_lensing/analysis/run_analysis.py $BASEDIR/$RUN_NAME -shear_cut=0.06 --vb --overwrite

done

mv meanshear-*.err meanshear-*.out $BASEDIR
