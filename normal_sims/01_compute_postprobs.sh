#!/bin/bash

#SBATCH --job-name=normal_sims
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1g
#SBATCH	--array=1-960
#SBATCH --output=/nas/longleaf/home/ethanalt/Projects/Paper2/normal_sims/Results/Log/slurmLogFiles%a.out
#SBATCH --error=/nas/longleaf/home/ethanalt/Projects/Paper2/normal_sims/Results/Error/%a.err

## add R module
module add gcc/6.3.0 
module add r/4.1.0

R CMD BATCH --no-restore /nas/longleaf/home/ethanalt/Projects/Paper2/normal_sims/R/01_posterior.R /nas/longleaf/home/ethanalt/Projects/Paper2/normal_sims/Results/Rout/R_$SLURM_ARRAY_TASK_ID.Rout