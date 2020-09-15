#!/bin/bash
# Chamara Rajapakshe (2020)
# MSCART runs for 1D cloud fields (step/fractal)
# ------------------------------------------------
#SBATCH --job-name=step0p865_120
#SBATCH --output=slurm-step0p865SZA120.out
#SBATCH --partition=batch
#SBATCH --qos=short+
#SBATCH --mem=2500
#SBATCH --array=0-10%10

MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
#MSCART='/umbc/xfs1/zzbatmos/users/charaj1/MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e4
bnd=0p865
job=$SLURM_ARRAY_TASK_ID

echo ${job} $HOSTNAME
fn=$(python LES_MSCART_setup.py 120 0 step_cld_b0p865re10ve02_x15km.nc ${bnd} 2>&1)
echo ${fn}
${MSCART} 10 ${NPH} 0 ${fn}.nml results/stepEx/${fn}_NPH${NPH}_${job}.nc
#Random delay to avoid I/O overflow
T=$(( ( RANDOM % 10 )  + 1 ))  #Generates number between 1-10
T=$(( 12*$T ))                 #Scales random number
echo "Waiting for "$T          #Display message
sleep $T                       #Process waits/sleeps before finishing
