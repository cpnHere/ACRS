#!/bin/bash

#SBATCH --job-name=ATcol0p860_160
#SBATCH --output=slurm%J-AT0p860_160_1e5.out
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --mem=5000
#SBATCH --array=0-5184%250
##SBATCH --nodelist=cnode[002-029,031-044]
##SBATCH --exclude=cnode109 
##SBATCH --constraint=hpcf2018
MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
job=$(($SLURM_ARRAY_TASK_ID))
field_path=LESb0p860_bins/
field_pat=ATEXc_dharma_007877_b0p860
n_grid=144 # Number of x/y grids
echo ${job} $HOSTNAME
python LES_MSCART_column_run.py 160 0 ${field_path} ${field_pat} ${job} ${MSCART} ${NPH} ${n_grid}
#Random delay to avoid I/O overflow
T=$(( ( RANDOM % 10 )  + 1 ))  #Generates number between 1-10
T=$(( 12*$T ))                 #Scales random number
echo "Waiting for "$T          #Display message
sleep $T                       #Process waits/sleeps before finishing
