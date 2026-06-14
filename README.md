# DEMENT_overflow

## Data and code file overview
This repository contains 24 parameter files, 16 code scripts, 4 txt files containing simulation inputs, and this README document, with the following files

## Files

### Main scripts
1. JoseMurua_DEMENT_rep_Simulation_Setup.R
2. JoseMurua_DEMENT_rep_Processing_HR.R
3. JoseMurua_DEMENT_rep_Processing_HRo.R
4. JoseMurua_DEMENT_rep_Processing_LR.R
5. JoseMurua_DEMENT_rep_Processing_LRo.R
6. JoseMurua_DEMENT_rep_Figures.R
7. DEMENTjob_1.sh
8. DEMENTjob_2.sh
9. DEMENTjob_3.sh
10. DEMENTjob_4.sh
11. Processing_job_HR.sh
12. Processing_job_HRo.sh
13. Processing_job_LR.sh
14. Processing_job_LRo.sh


### Files for simulations

1. Common files
   - DEMENT.0.7.6.R: This file contains the code for the DEMENT model
   - DEMENTBatch.R. This file executes the DEMENT code provided with a timestamp and parameter file ID (single integer) as inputs
   - climate2.txt: This file contains abiotic forcing variables (temperature and moisture)
   - Ea.txt: This file contains the activation energy for the degradation of different polymers
   - inputs.txt: This file specifies substrate inputs per time-step in DEMENT
   - substrates.txt: This file specifies the initial amount of litter substrate in simulations

2. Model parameter files

   Each simulation scenario directory contains four parameter files
   - params1.txt: 100 taxa, CO2 overflow
   - params2.txt: 1 taxon, CO2 overflow
   - params3.txt: 1 taxon, DOC overflow
   - params4.txt: 100 taxa, DOC overflow
  
   These parameter files vary between scenarios in the following aspects:
   - Flexible biomass stoichiometry
   - Enzyme allocation
   - Uptake allocation

## Workflow

There are two possible workflows:
1. Run the simulations from scratch (requires high performance computing)
2. Reproduce the figures from txt simulation outputs

### Run simulations from scratch

1. Run the script "JoseMurua_DEMENT_rep_Simulation_Setup.R". This does the following:
   - Creates additional directories inside each scenario directory. The new directories correspond to different initial litter chemistries (from C:N = 10 to C:N = 90)
   - The corresponding parameter files are put inside of each of these new directories together with files in the "Common_files" directory
   - The intial litter chemistries are calculated and txt files are created for each of these. The files are then placed in the correct directory
   - A list of batch commands is created to run each of the simulations. This list is exported as a txt file and is later used to run the simulations in paralell
   - A list of directories is created and exported as a txt file to run each simulation. The directory specified in the script **should be replaced with your own directory in the HPC environment**
   - A list of timestamps is generated and exported as a txt file. This will be used to run and identify each of the simulations

2. Put the entire directory in a HPC environment
3. Submit "DEMENTjob" files 1 through 4 as jobs to run the simulations in paralell. Model outputs will be generated as an R object containing a nested list of outputs
4. Submit files "Processing_job" HR, HRo, LR, and LRo as jobs to process the R objects containing the simulation outputs. This jobs will use the "Processing" R scripts to generate a series of txt files containing simulation outputs in a more usable format. The subscripts HR and LR stand for "high redundancy" (i.e., simulations with 1 taxon) and "low redundancy" (i.e., simulations with 100 taxa), respectively. If the subscript includes a lowercase "o", this indicates that overflow occurs as dissolved organic matter (i.e., DOC).

### Reproduce figures from simulation outputs

With the simulation outputs now processed and exported as txt files, these can now be used locally to reproduce the figures.

1. Run the script "JoseMurua_DEMENT_rep_Figures.R". This script will read the corresponding txt files and produce all figures.

## Software versions

- R 4.2.1
- tidyverse 1.3.2
- ggplot2 3.5.1
- tibble 3.2.1
- readr 2.1.3
- dplyr 1.1.4
- stringr 1.5.0
- patchwork 1.1.2
- cowplot 1.1.1
- grid 4.2.1
