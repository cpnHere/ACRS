#!/bin/bash

#SBATCH --job-name=Combine_1D
#SBATCH --output=slurm-Combine.out
#SBATCH --partition=batch
#SBATCH --qos=short
#SBATCH --mem=4000

fnames='ATEXp_dharma_013067_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH1e5'
pc_format='LES_1D'#LES/step/fractal
action='May2018_ATEXc'
stdbuf -oL python combine_1D.py ${fnames} ${pc_format} ${action} > combine2p13_NPH1e5_SZA160.log

