#!/bin/bash

#SBATCH --job-name=Combine_3D
#SBATCH --output=slurm-Combine.out
#SBATCH --partition=batch
#SBATCH --qos=short+
#SBATCH --mem=32000

SZA='120'
band='0p860'
dpath='results/b0p860/'
fnames='RICO_dharma_005044_b0p860_MSCART_SZA120_SAA000_Nvaa004_NPH1e6'
runs=500
pc_format='LES_3D'#LES/step/fractal
action='Dec2020_RICO_3D_multiple_VAA'
stdbuf -oL python combine.py ${dpath} ${fnames} ${runs} ${pc_format} ${action} &> combine0p860_NPH1e6_SZA120.log

