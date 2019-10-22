#!/bin/bash

#SBATCH --job-name=Combine_3D
#SBATCH --output=slurm-Combine.out
#SBATCH --partition=batch
#SBATCH --qos=short
#SBATCH --mem=32000

SZA='160'
band='0p860'
dpath='results/b0p860/'
fnames='ATEXp_dharma_013067_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH1e6'
runs=500
pc_format='LES_3D'#LES/step/fractal
action='May2018_ATEXp_3D'
stdbuf -oL python combine.py ${dpath} ${fnames} ${runs} ${pc_format} ${action} &> combine0p860_NPH1e6_SZA160.log

