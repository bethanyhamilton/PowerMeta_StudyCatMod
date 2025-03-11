#!/bin/bash
#SBATCH -J powermeta_catmod1 #Job name
#SBATCH -o powermeta_catmod1.o%j  # Name of outputfile
#SBATCH -e powermeta_catmod1.e%j  # Name of stderr error file
#SBATCH -p skx		  # Queue (partition) name
#SBATCH -N 1           # Total number of nodes
#SBATCH -n 48              # Total # of mpi tasks
#SBATCH -t 01:00:00     # Run time (hh:mm:ss) 
#SBATCH --mail-user=bethanyhamilton@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

module load python3 # or newer
#pip install paramiko

python pmcmlauncher.py