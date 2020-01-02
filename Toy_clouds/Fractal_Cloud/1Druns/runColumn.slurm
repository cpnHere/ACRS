#!/bin/bash

#SBATCH --job-name=facaC0p865_20
#SBATCH --partition=batch
#SBATCH --qos=short+
#SBATCH --mem=2500
#SBATCH --array=0-4096%512
MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
job=$SLURM_ARRAY_TASK_ID
field_path=fracb0p865_bins/
field_pat=fractal_cld_b865re12ve05_x40km 
echo ${job} $HOSTNAME
python LES_MSCART_column_run.py 160 0 ${field_path} ${field_pat} ${job} ${MSCART} ${NPH}
