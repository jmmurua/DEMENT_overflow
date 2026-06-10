#!/bin/bash

#SBATCH --job-name=DEMENT
#SBATCH -A allisons_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --error=slurm-%J.err

module purge
hostname
module load R/4.2.2


# Pass line #i of files to a R script
Rscript --no-restore JoseMurua_DEMENT_rep_Figures_HRo.R