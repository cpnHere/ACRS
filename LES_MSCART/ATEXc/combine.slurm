#!/bin/bash

#SBATCH --job-name=Combine_3D
#SBATCH --output=slurm-Combine.out
#SBATCH --partition=batch
#SBATCH --qos=short
#SBATCH --mem=32000

SZA='160'
band='2p13'
dpath='results/b2p13/'
fnames='ATEXc_dharma_007877_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH1e6'
runs=500
pc_format='LES_3D'#LES/step/fractal
action='May2018_ATEXc_3D'
stdbuf -oL python combine.py ${dpath} ${fnames} ${runs} ${pc_format} ${action} &> combine2p13_NPH1e6_SZA160.log

