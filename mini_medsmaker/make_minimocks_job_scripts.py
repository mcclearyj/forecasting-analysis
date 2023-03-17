import glob
from argparse import ArgumentParser
import os
import ipdb, pdb
import re
import io

def parse_args():

    parser = ArgumentParser()
    
    parser.add_argument('run_name',type=str,
                        help='Root name of simulations to be run [assumes @sweverett-style organization]')
    parser.add_argument('-js_dir',type=str, default=None,
                        help='Output directory for job scripts [default: ./job_scripts/]')
    parser.add_argument('-base_dir',type=str, default=None,
                        help='Directory containing sim catalogs & images [default: ./]')
    parser.add_argument('-n_reals',type=int, default=None,
                        help='Number of realizations to group together in a batch script [default: 5]')

    return parser.parse_args()



def print_bashfile(run_name, cluster_name, mock_dir_names, start_r, end_r):

    output = io.StringIO()

    print('#!/bin/sh', file=output)
    print('#SBATCH -t 59:59', file=output)
    print('#SBATCH -N 1', file=output)
    print('#SBATCH -n 16', file=output)
    print('#SBATCH --mem-per-cpu=6g', file=output)
    print('#SBATCH --partition=express', file=output)
    print(f'#SBATCH -J {cluster_name}_minimocks_r{start_r}r{end_r}', file=output)
    print('#SBATCH -v', file=output)
    print('#SBATCH --mail-type=ALL', file=output)
    print('#SBATCH --mail-user=jmac.ftw@gmail.com', file=output)
    print(f'#SBATCH -o minimocks_{run_name}_{cluster_name}_r{start_r}r{end_r}.out', file=output)
    print(f'#SBATCH -e minimocks_{run_name}_{cluster_name}_r{start_r}r{end_r}.err', file=output)

    print('', file=output)
    print('', file=output)

    print('source /work/mccleary_group/miniconda3/bin/activate', file=output)
    print('', file=output)
    print('conda activate sbmcal_139', file=output)

    print('echo $PATH', file=output)
    print('echo $PYTHONPATH', file=output)

    print('dirname="slurm_outfiles"', file=output)
    print('if [ ! -d "$dirname" ]', file=output)
    print('then', file=output)
    print('     echo " Directory $dirname does not exist. Creating now"', file=output)
    print('     mkdir -p -- "$dirname"', file=output)
    print('     echo " $dirname created"', file=output)
    print(' else', file=output)
    print('     echo " Directory $dirname exists"', file=output)
    print(' fi', file=output)
    print('', file=output)
    print(' echo "Proceeding with code..."', file=output)

    print('', file=output)
    print('', file=output)


    for mock_dir_name in mock_dir_names:
        job_string = \
        f'python /work/mccleary_group/forecasting-analysis/mini_medsmaker/prep_dirs.py \
    --mock_dir {mock_dir_name}\n' \
        f'python /work/mccleary_group/forecasting-analysis/mini_medsmaker/mini_medsmaker.py \
    --mock_dir {mock_dir_name} \
    -fname_base {run_name} -run_name {run_name}'
        print(job_string, file=output, flush=True)




    print('', file=output)
    print(f'mv minimocks_{run_name}_{cluster_name}_r{start_r}r{end_r}.out minimocks_{run_name}_{cluster_name}_r{start_r}r{end_r}.err "$dirname"',file=output)
    print('', file=output)

    return output.getvalue()



def main(args):

    run_name = args.run_name
    js_dir = args.js_dir
    base_dir = args.base_dir
    n_reals = args.n_reals

    if js_dir == None:
        cwd = os.getcwd()
        js_dir = os.path.join(cwd, 'job_scripts')

    if base_dir == None:
       base_dir = os.getcwd()

    if n_reals == None:
        n_reals = 5

    run_dir = os.path.join(base_dir, run_name)
    run_js_dir = os.path.join(js_dir, run_name)

    # Make sure that the run_dir exists
    if not os.path.isdir(run_dir):
        raise OSError(f'\nrun_dir {run_dir} not found! \nCheck supplied base directory {base_dir} and run name {run_name}\n')

    # Create js_dir if it doesn't exist
    if not os.path.isdir(run_js_dir):
        cmd = 'mkdir -p {run_js_dir}'
        print(f'Saving batch job scripts in {js_dir}/{run_name}')
        os.system(cmd.format(run_js_dir=run_js_dir))
        
        
    # For convenience, write sbatch commands to a text file
    sb = open(f'submit_{run_name}_jobs.txt','w')

    # get cluster names and loop over them
    clusters_glob = glob.glob(f'{run_dir}/*')
    clusters_glob.sort()

    for cluster in clusters_glob:
        testing = []
        cluster_name = cluster.replace(f'{run_dir}/','')
        sb.write(f'\n\n##\n## {cluster_name}\n##\n\n')

        mock_dirs = glob.glob(f'{run_dir}/{cluster_name}/r*')
        mock_dirs.sort(key=lambda x: int(re.search('r[0-9]{1,2}', x.split(f'{cluster_name}')[1]).group()[1:]))

        for i, mock_dir in enumerate(mock_dirs):

            realization = re.search('r[0-9]{1,2}', mock_dir).group()

            if (i % n_reals) == 0:

                dir_list = mock_dirs[i:i+n_reals]
                testing.extend(dir_list)
                start_r = realization[1:]
                end_r   = str(int(realization[1:])+n_reals-1)

                if (int(end_r) >= len(mock_dirs)):
                    end_r = str(len(mock_dirs)-1)

                bash_name = os.path.join(js_dir,
                                        run_name,
                                        f'job_minimocks_{cluster_name}_r{start_r}_r{end_r}.sh'
                                        )

                sb.write(f'sbatch {bash_name}\n')

                with open(bash_name, 'w') as f:

                    script = print_bashfile(run_name, 
                                            cluster_name,
                                            dir_list,
                                            start_r,
                                            end_r
                                            )
                    f.write(script)

            else:
                continue

            if ((i+1) % n_reals == 0):
                sb.write('\n')

    sb.close()

    return 0



if __name__ == '__main__':

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('make_minimocks_job_scripts.py has completed succesfully')
    else:
        print(f'make_minimocks_job_scripts.py has failed w/ rc={rc}')
