#!/bin/bash

#SBATCH --job-name=fractal0p860
#SBATCH --output=slurm-fractal_cld_b0p865re12ve05_x40km.out
#SBATCH --partition=batch
#SBATCH --qos=short+
#SBATCH --mem=2500
#SBATCH --time=1:00:00
#SBATCH --array=1-100%100

#MSCART='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/MSCART/engine/build/intel-release-nocaf/mscart/MSCART'
MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
band='0p865'
job=$SLURM_ARRAY_TASK_ID

echo ${job} $HOSTNAME
fn=$(python LES_MSCART_setup.py 100 0 fractal_cld_b865re12ve05_x40km.nc 2>&1)
echo ${fn}
${MSCART} 10 ${NPH} 0 ${fn}.nml results/b${band}/${fn}_NPH${NPH}_${job}.nc
#fn=$(python LES_MSCART_setup.py 180 2>&1)
#echo ${fn}
#${MSCART} 10 ${NPH} 0 ${fn}.nml results/${fn}_NPH${NPH}_${job}.nc

#mv slurm.out slurm${NPH}.out
