#!/bin/bash

#SBATCH --job-name=ATp0p860SZA160
#SBATCH --output=slurm-ATp0p860SZA160_1e6.out
#SBATCH --partition=batch
#SBATCH --qos=medium
#SBATCH --mem=5000
##SBATCH --dependency=afterok:7429488
#SBATCH --array=1-500%200

MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e6
band='0p860'
job=$SLURM_ARRAY_TASK_ID

echo ${job} $HOSTNAME
fn=$(python LES_MSCART_setup.py 160 0 ATEXp_dharma_013067_b0p860.nc 2>&1)
echo ${fn}_NPH${NPH}_${job}.nc
${MSCART} 10 ${NPH} 0 ${fn}.nml results/b${band}/${fn}_NPH${NPH}_${job}.nc
#Random delay to avoid I/O overflow
T=$(( ( RANDOM % 10 )  + 1 ))  #Generates number between 1-10
T=$(( 12*$T ))                 #Scales random number
echo "Waiting for "$T          #Display message
sleep $T                       #Process waits/sleeps before finishing
