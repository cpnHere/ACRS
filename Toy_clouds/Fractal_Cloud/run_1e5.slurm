#!/bin/bash

#SBATCH --job-name=f2p13_60
#SBATCH --output=slurm-fractal_StCu_b2p13re12ve05_x40km_H0p5_NPH1e5SZA60.out
#SBATCH --partition=batch
#SBATCH --qos=short
#SBATCH --mem=2500
#SBATCH --array=1-100%100

MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
#MSCART='/umbc/xfs1/zzbatmos/users/charaj1/MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
band='2p13'
job=$SLURM_ARRAY_TASK_ID

echo ${job} $HOSTNAME
fn=$(python LES_MSCART_setup.py 120 0 fractal_StCu_b2p13re12ve05_x40km_H0p5.nc 2>&1)
echo ${fn}
${MSCART} 10 ${NPH} 0 ${fn}.nml results/fracb${band}/${fn}_NPH${NPH}_${job}.nc
#fn=$(python LES_MSCART_setup.py 180 2>&1)
#echo ${fn}
#${MSCART} 10 ${NPH} 0 ${fn}.nml results/${fn}_NPH${NPH}_${job}.nc

#mv slurm.out slurm${NPH}.out
