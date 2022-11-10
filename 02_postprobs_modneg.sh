#!/bin/bash

#SBATCH --job-name=copula_new
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25m
#SBATCH	--array=1-18000
#SBATCH --output=/nas/longleaf/home/ethanalt/Projects/Paper2/Results/0Log/slurmLogFiles%a.out
#SBATCH --error=/nas/longleaf/home/ethanalt/Projects/Paper2/Results/0Error/%a.err

## add R module
module add gcc/6.3.0 
module add r/4.1.0

R CMD BATCH --no-restore /nas/longleaf/home/ethanalt/Projects/Paper2/R/02_compute_postprob_modneg.R /nas/longleaf/home/ethanalt/Projects/Paper2/Results/0Rout/R_$SLURM_ARRAY_TASK_ID.Rout