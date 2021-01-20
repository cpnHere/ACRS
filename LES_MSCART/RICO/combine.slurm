#!/bin/bash

#SBATCH --job-name=Combine_3D
#SBATCH --output=slurm-Combine.out
#SBATCH --partition=batch
#SBATCH --qos=short+
#SBATCH --mem=32000

dpath='results/b0p860/' #Data path
fnames='RICO_dharma_005044_b0p860_MSCART_SZA120_SAA000_Nvaa72VZA002_NPH1e4' #Filename without the job number
runs=3 # Number of batches 
pc_format='LES_3D'#LES/step/fractal
action='Dec2020_RICO_3D_multiple_VAA'
stdbuf -oL python combine.py ${dpath} ${fnames} ${runs} ${pc_format} ${action} &> combine.log

