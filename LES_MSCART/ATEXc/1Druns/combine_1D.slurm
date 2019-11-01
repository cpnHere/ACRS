#!/bin/bash

#SBATCH --job-name=Combine_1D
#SBATCH --output=slurm-Combine.out
#SBATCH --partition=batch
#SBATCH --qos=short+
#SBATCH --mem=4000

fnames='ATEXc_dharma_007877_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH1e5'
pc_format='LES_1D'#LES/step/fractal
action='Oct2019_ATEXc_1D'
stdbuf -oL python combine_1D.py ${fnames} ${pc_format} ${action} > combine0p860_NPH1e5_SZA160.log

