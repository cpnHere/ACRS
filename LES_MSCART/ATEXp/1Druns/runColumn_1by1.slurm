#!/bin/bash

#SBATCH --job-name=LEScol2p13_160
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --mem=2500
#SBATCH --array=1-83%83
MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
#----------Direct input of missing job numbers ------------------
#job_pool=(4317 4461)
#job=${job_pool[$SLURM_ARRAY_TASK_ID]}
#---------Reading missing job index from a file ----------------
job_pool=ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5_missingRuns.dat
ji=$(($SLURM_ARRAY_TASK_ID+1)) #because the actual list of numbers starts from the second line
job=$(sed ${ji}q';'d ${job_pool})

field_path=LESb2p13_bins/
field_pat=ATEXp_dharma_013067_b2p13
n_grid=144 # Number of x/y grids
echo ${job} $ji $HOSTNAME
python LES_MSCART_column_run.py 160 0 ${field_path} ${field_pat} ${job} ${MSCART} ${NPH} ${n_grid} 

#Random delay to avoid I/O overflow
T=$(( ( RANDOM % 10 )  + 1 ))  #Generates number between 1-10
T=$(( 12*$T ))                 #Scales random number
echo "Waiting for "$T          #Display message
sleep $T                       #Process waits/sleeps before finishing
