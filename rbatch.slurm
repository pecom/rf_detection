#!/bin/bash
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH -t 1:00:00
#SBATCH -o /global/homes/p/pecom/logs/rf_logs/rf.%A.out
#SBATCH -e /global/homes/p/pecom/logs/rf_logs/rf.%A.err

module load R
cd $HOME/rf_detection/scripts
# R CMD BATCH rfor_script.R
# Rscript rfor_panic.r 1
Rscript rfor_tests_regression.r

