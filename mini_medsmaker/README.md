Assuming directory is cleaned with cleanup.sh, here are the steps:

1. Next run make_minimocks_job_scripts.py (use --help if needed) and submit jobs

Example command:
python make_minimocks_job_scripts.py run_name -base_dir /scratch/vassilakis.g/mock-data-forecasting/

2. Run make_plots.py with config (or copy of config) plotting_george.yaml

