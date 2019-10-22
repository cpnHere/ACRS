#!/bin/bash

#SBATCH --job-name=ATcolS0p860_160
#SBATCH --output=ATcolS0p860_160.out
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --mem=2500

MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
job=9167
field_path=LESb0p860_bins/
field_pat=ATEXc_dharma_007877_b0p860
n_grid=144 # Number of x/y grids
echo ${job} $HOSTNAME
python LES_MSCART_column_run.py 160 0 ${field_path} ${field_pat} ${job} ${MSCART} ${NPH} ${n_grid}
