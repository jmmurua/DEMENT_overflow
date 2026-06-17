# DEMENT_overflow

## Data and code file overview
This repository contains 24 parameter files, 16 code scripts, 4 txt files containing simulation inputs, and this README document, with the following files

## Files

### Data files and variables

Files are named based on the main variable they contain from simulation outputs, followed by a subscript
The subscript represents the following:
- HR: "High redundancy". Contains data from simulations with 1 taxon, where carbon overflow occurs as CO2
- HRo: "High redundancy - overflow". Contains data from simulations with 1 taxon, where carbon overflow occurs as dissolved organic carbon
- LR: "Low redundancy". Contains data from simulations with 100 taxa, where carbon overflow occurs as CO2
- LRo: "Low redundancy - overflow". Contains data from simulations with 100 taxa, where carbon overflow occurs as dissolved organic carbon

All datasets share the following variables:
- Simulation: Simulation identifier that consists of the mechanism and the litter C:N separated by an underscore
- Mechanism: Specifies the mechanism that solves elemental imbalance. Can take the following values:
	- Overflow
	- Stoichiometry
	- Enzyme
	- Enzyme-Uptake
	- All-NoUptake
	- All
- Litter C:N: Specifies the initial litter carbon to nitrogen ratio in that simulation. Corresponds to a number preceded by "CN"
- Days: Specifies the time-step of the simulation as number of days since time 0.

**Note: In all cases, to get the correct units, variables must be divided by the grid size of simulations. In our case, the grid was 100 x 100 = 10,000**

1. Biomass
	- C.biomass: Total microbial carbon biomass (mg cm-2)
	- N.biomass: Total microbial nitrogen biomass (mg cm-2)
	- P.biomass: Total microbial phosphorus biomass (mg cm-2)
	- Biomass: Total microbial biomass (mg cm-2)
	
2. EnzProdSerL
	- Enzyme: Enzyme identifier
	- Abundance: Amount of enzyme (mg cm-2)
	- C.biomass: Total microbial carbon biomass (mg cm-2)
	- N.biomass: Total microbial nitrogen biomass (mg cm-2)
	- P.biomass: Total microbial phosphorus biomass (mg cm-2)
	- Biomass: Total microbial biomass (mg cm-2)
	- Enz.type: Element(s) targeted by the enzyme in question

3. Growth:
	- Element: Carbon, nitrogen, or phosphorus
	- Growth: Amount of microbial biomass growth in units of a specific element (mg cm-2)
	
4. LitterSubsL
	- Litter.C: Carbon present as polymer in litter (mg cm-2)
	- Litter.N: Nitrogen present as polymer in litter (mg cm-2)
	- Substrate: Litter substrate identifier
	- Abundance: Abundance of a specific substrate in carbon units (mg cm-2)
	- Subs.type: Main element contained in a specific substrate
	
5. N_MonomerUptake
	- NH4: Uptake of ammonium in nitrogen units (mg cm-2)
	- PO4: Uptake of phosphate in nitrogen units (mg cm-2)
	- DeadMic: Uptake of microbial necromass monomers in nitrogen units (mg cm-2)
	- DeadEnz: Uptake of inactive enzyme monomers in nitrogen units (mg cm-2)
	- Cellulose: Uptake of cellulose monomers in nitrogen units (mg cm-2)
	- Hemicellulose: Uptake of hemicellulose monomers in nitrogen units (mg cm-2)
	- Starch: Uptake of starch monomers in nitrogen units (mg cm-2)
	- Chitin: Uptake of chitin monomers in nitrogen units (mg cm-2)
	- Lignin: Uptake of lignin monomers in nitrogen units (mg cm-2)
	- Protein1: Uptake of protein1 monomers in nitrogen units (mg cm-2)
	- Protein2: Uptake of protein2 monomers in nitrogen units (mg cm-2)
	- Protein3: Uptake of protein3 monomers in nitrogen units (mg cm-2)
	- OrgP1: Uptake of organic phosphorus1 monomers in nitrogen units (mg cm-2)
	- OrgP2: Uptake of organic phosphorus2 monomers in nitrogen units (mg cm-2)
	- Tot.recycl: Total nitrogen uptake from from necromass plus inactive enzymes (mg cm-2)
	- Tot.litt: Total nitrogen uptake from litter substrates (mg cm-2)
	
6. Necromass
	- DeadMic: Total microbial necromass (mg cm-2)
	- DeadEnz: Total inactive enzymes (mg cm-2)
	
7. RespComp
	- Source: Specific component of respiration
	- Resp: Amount of CO2 produced from a specific component (mg cm-2)
	
8. UptSer
	- Taxon_Uptake_C: Total carbon uptake (mg cm-2)
	- Taxon_Uptake_N: Total nitrogen uptake (mg cm-2)
	- Taxon_Uptake_P: Total phosphorus uptake (mg cm-2)


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
