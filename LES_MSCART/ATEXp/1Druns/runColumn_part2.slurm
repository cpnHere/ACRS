#!/bin/bash

#SBATCH --job-name=ATpcol2p13_160
#SBATCH --output=slurm-ATpcol2p13_160_1e5.out
#SBATCH --partition=batch
#SBATCH --qos=short
#SBATCH --mem=2500
#SBATCH --array=0-10368%250
MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
job=$(($SLURM_ARRAY_TASK_ID+10368))
field_path=LESb2p13_bins/
field_pat=ATEXp_dharma_013067_b2p13
n_grid=144 # Number of x/y grids
echo ${job} $HOSTNAME
python LES_MSCART_column_run.py 160 0 ${field_path} ${field_pat} ${job} ${MSCART} ${NPH} ${n_grid}
#Random delay to avoid I/O overflow
T=$(( ( RANDOM % 10 )  + 1 ))  #Generates number between 1-10
T=$(( 12*$T ))                 #Scales random number
echo "Waiting for "$T          #Display message
sleep $T                       #Process waits/sleeps before finishing
