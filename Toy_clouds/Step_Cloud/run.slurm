#!/bin/bash

#SBATCH --job-name=CR470_60
#SBATCH --output=slurm-CR470SZA60.out
#SBATCH --partition=batch
#SBATCH --qos=short
#SBATCH --mem=2500
#SBATCH --array=0-55%56

MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
#MSCART='/umbc/xfs1/zzbatmos/users/charaj1/MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
bnd=0p860
job=$SLURM_ARRAY_TASK_ID
NmaxScat=5000

echo ${job} $HOSTNAME
fn=$(python LES_MSCART_setup.py 120 0 ${job} ${bnd} 2>&1)
echo ${fn}
${MSCART} 10 ${NPH} 0 ${fn}.nml color_ratio_LUTs/${fn}_NPH${NPH}.nc
#fn=$(python LES_MSCART_setup.py 180 2>&1)
#echo ${fn}
#${MSCART} 10 ${NPH} 0 ${fn}.nml results/${fn}_NPH${NPH}_${job}.nc

#mv slurm.out slurm${NPH}.out
