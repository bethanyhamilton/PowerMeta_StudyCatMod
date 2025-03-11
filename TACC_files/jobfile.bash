#!/bin/bash
#SBATCH -J powermeta_catmod1 #Job name
#SBATCH -o powermeta_catmod1.o%j  # Name of outputfile
#SBATCH -e powermeta_catmod1.e%j  # Name of stderr error file
#SBATCH -p skx		  # Queue (partition) name


# For batch 1 we have 3840 commands of which takes an average of 951.168 hours 
# for 1 core
# with 48 cores that is 19.816 hours. 
# Across 4 nodes it becomes 5 hours. Round to 6 for longer conditions

#SBATCH -N 4           # Total number of nodes
#SBATCH -n 48              # Total # of mpi tasks
#SBATCH -t 06:00:00     # Run time (hh:mm:ss) 
#SBATCH --mail-user=bethanyhamilton@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

module load python3 # or newer
#pip install paramiko

# load R module
module load Rstats/4.4.2

python pmcmlauncher.py