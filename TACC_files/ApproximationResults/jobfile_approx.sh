#!/bin/bash
#SBATCH -J powstudcatmod_approx #Job name
#SBATCH -o powstudcatmod_approx.o%j  # Name of outputfile
#SBATCH -e powstudcatmod_approx.e%j  # Name of stderr error file
#SBATCH -p skx		  # Queue (partition) name
#SBATCH -N 4           # Total number of nodes
#SBATCH --ntasks-per-node 48              # Total # of mpi tasks
#SBATCH -t 05:00:00     # Run time (hh:mm:ss) 
#SBATCH --mail-user=bethanyhamilton@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A TG-MTH250014 #my project

module load python/3.9.18 
pip install paramiko

module spider apptainer

module load tacc-apptainer

module load pylauncher/5.1.1 

apptainer pull docker://bethanyhamilton/powstudcatmod_approx:v0

python3 /work2/08147/bethanyh/stampede3/work_files/pmcmlauncher_approx.py
