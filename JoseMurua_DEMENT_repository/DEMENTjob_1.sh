#!/bin/bash

#SBATCH --job-name=DEMENT
#SBATCH -A allisons_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --error=slurm-%J.err
#SBATCH --array=1-486

module purge
hostname
module load R/4.2.2

# Get array ID
i=${SLURM_ARRAY_TASK_ID}

## Text files to read
Time_stamps="Timestamps.txt"
Directories="Directories.txt"

# Read line #i from the timestamps file
t_stamp=$(sed "${i}q;d" ${Time_stamps})

## Read line #i from the Directories file
Directory=$(sed "${i}q;d" ${Directories})

# Pass line #i of directories and set directory
cd ${Directory}

# Pass line #i of files to a R script
Rscript --no-restore DEMENTBatch.R ${t_stamp} 1