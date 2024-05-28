#!/bin/bash
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH -t 8:00:00
#SBATCH -o /global/homes/p/pecom/logs/rf_logs/rf.%A.out
#SBATCH -e /global/homes/p/pecom/logs/rf_logs/rf.%A.err

module load R
cd $HOME/rf_detection/scripts/paper
# R CMD BATCH rfor_script.R

for i in $(seq 1 5)
do
	Rscript rfor_err.r $i
	Rscript rfor_10feat.r $i
	Rscript rfor_9feat.r $i
	Rscript rfor_7feat.r $i
	Rscript rfor_6feat.r $i
	Rscript rfor_2class.r $i
	Rscript rfor_regress.r $i
done
