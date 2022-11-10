#!/bin/bash

#SBATCH --job-name=copula_pos
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1g
#SBATCH	--array=2
#SBATCH --output=/nas/longleaf/home/ethanalt/Projects/Paper2/Results/0Log/slurmLogFiles%a.out
#SBATCH --error=/nas/longleaf/home/ethanalt/Projects/Paper2/Results/0Error/%a.err

## add R module
module add gcc/6.3.0 
module add r/4.1.0

R CMD BATCH --no-restore /nas/longleaf/home/ethanalt/Projects/Paper2/R/03_compute_pos.R /nas/longleaf/home/ethanalt/Projects/Paper2/Results/0Rout/R_$SLURM_ARRAY_TASK_ID.Rout