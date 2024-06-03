#!/bin/bash
#SBATCH --account=rrg-jnorthru
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4
#SBATCH --job-name="tfr_brms"
#SBATCH --output=success.out
module load StdEnv/2020  gcc/9.3.0  r/4.1.2
Rscript tfr_model_ttharvest.R

