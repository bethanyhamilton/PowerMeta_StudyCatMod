#!/bin/bash
#SBATCH -J powermeta_catmod_test #Job name
#SBATCH -o powermeta_catmod_test.o%j  # Name of outputfile
#SBATCH -e powermeta_catmod_test.e%j  # Name of stderr error file
#SBATCH -p skx		  # Queue (partition) name
#SBATCH -N 1           # Total number of nodes
#SBATCH -n 48              # Total # of mpi tasks
#SBATCH -t 01:00:00     # Run time (hh:mm:ss) 
#SBATCH --mail-user=bethanyhamilton@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

### NEED to make docker container!! 

module load python3 # or newer
#pip install paramiko

# load R module
module load Rstats/4.0.3

python pmcmlauncher.py