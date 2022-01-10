#!/usr/bin/env zsh

### general job settings
#SBATCH --job-name=For-Philip
#SBATCH --output=output-%j.txt

### resource requirements
#SBATCH --time=24:00:00
### max recommended memory per node: 187200M
#SBATCH --mem-per-cpu=3900M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
### CLAIX2018 skylake nodes
#SBATCH --partition=c18m

### default R version: 3.6.0
### other available R version: 4.0.2
module load MATH r-project/4.0.2
Rscript Experiment_C3_hidden.R

