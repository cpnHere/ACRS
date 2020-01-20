#!/bin/bash

#SBATCH --job-name=Combine_1D
#SBATCH --output=slurm-Combine.out
#SBATCH --partition=batch
#SBATCH --qos=short+
#SBATCH --mem=4000

fnames='ATEXp_dharma_013067_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH1e6'
pc_format='LES_1D'#LES/step/fractal
action='Jan2020_ATEXp'
stdbuf -oL python combine_1D.py ${fnames} ${pc_format} ${action} > combine0p860_NPH1e6_SZA160.log

