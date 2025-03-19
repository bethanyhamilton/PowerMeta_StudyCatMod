#!/bin/bash
#SBATCH -J powstudcatmod_test #Job name
#SBATCH -o powstudcatmod_test.o%j  # Name of outputfile
#SBATCH -e powstudcatmod_test.e%j  # Name of stderr error file
#SBATCH -p skx		  # Queue (partition) name
#SBATCH -N 1           # Total number of nodes
#SBATCH -ntasks-per-node 48              # Total # of mpi tasks
#SBATCH -t 01:00:00     # Run time (hh:mm:ss) 
#SBATCH --mail-user=bethanyhamilton@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A powstudcatmod_test #my project


module load python3 # or newer
#pip install paramiko

module spider apptainer

module load tacc-apptainer

apptainer pull docker://bethanyhamilton/powstudcatmod-container_test:v0

python3 pmcmlauncher.py

