#!/bin/bash

#SBATCH --job-name=LESlwp
#SBATCH --output=slurm-LESlwp.out
#SBATCH --partition=batch
#SBATCH --qos=long
#SBATCH --mem=32000

python Dharma_LES_read.py 
