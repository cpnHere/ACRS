#!/bin/bash
# Chamara Rajapakshe (2020)
# To run DISORT to generate LUTs
#-----------------------------------
#SBATCH --job-name=LUTorig
#SBATCH --output=LUTorig.out
#SBATCH --partition=batch
#SBATCH --qos=long+
#SBATCH --mem=5000
##SBATCH --–ntasks-per-node=36
##SBATCH --dependency=afterok:7429488
##SBATCH --array=1-500%250

python LUT_runs.py
