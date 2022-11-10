#!/bin/bash

#SBATCH --job-name=pp_compare
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25m
#SBATCH	--array=1-800
#SBATCH --output=/nas/longleaf/home/ethanalt/Projects/Paper2/normal_approx/Results/Log/slurmLogFiles%a.out
#SBATCH --error=/nas/longleaf/home/ethanalt/Projects/Paper2/normal_approx/Results/Error/%a.err

## add R module
module add gcc/6.3.0 
module add r/4.1.0

R CMD BATCH --no-restore /nas/longleaf/home/ethanalt/Projects/Paper2/normal_approx/R/01_postprobs_ind.R /nas/longleaf/home/ethanalt/Projects/Paper2/normal_approx/Results/Rout/R_$SLURM_ARRAY_TASK_ID.Rout