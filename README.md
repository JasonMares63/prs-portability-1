# PRS Portability

This project attempts to probe the portability of polygenic risk scores (PRS) across and within populations, using continuous genetic distance measures.

## Data

This project uses data from four sources.
* UK Biobank
* 1000 Genomes Project
* Biobank Japan
* Neale Lab UK Biobank GWAS
    * Included for comparisons of GWAS results

## Project overview

`scripts/` contains all the scripts used for computation and organization.
Everything done on a remote server used a script located in this directory. 
The scripts are numbered in the order they must be run to recreate the project.
Some numbers include multiple scripts, and these are distinguished from one another by letters.
Typically, the bash script should be run, and other scripts are short helpers written in R or Python, that are called from within the script only.
Some bash scripts were replaced by Python scripts to quickly queue multiple computations.
These Python scripts typically have a bash script within them, represented as a string.
These should be run using the provided conda environment.

`nb/` contains Jupyter notebooks used for higher-level analyses.
These are also numbered, and they depend on the outputs of the `scripts/` directory files.

