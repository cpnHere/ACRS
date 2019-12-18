#!/bin/bash

#SBATCH --job-name=fcSing0p860_60
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --mem=2500
fractal_StCu_ACA_b0p860_x40km_H0p5_AOT0p25_H0p5_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5
MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
job=235
field_path=fracb0p860_bins/
field_pat=fractal_StCu_ACA_b0p860_x40km_H0p5_AOT0p25_H0p5
echo ${job} $HOSTNAME
python LES_MSCART_column_run.py 120 0 ${field_path} ${field_pat} ${job} ${MSCART} ${NPH}
