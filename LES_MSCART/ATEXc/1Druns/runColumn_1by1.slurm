#!/bin/bash

#SBATCH --job-name=testLEScol0p860_140
#SBATCH --output=slurm%J-AT0p860_140_1e5.out
#SBATCH --partition=batch
#SBATCH --qos=long+
#SBATCH --mem=2500
##SBATCH --array=1-8%8
MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
#----------Direct input of missing job numbers ------------------
#job_pool=(4317 4461)
#job=${job_pool[$SLURM_ARRAY_TASK_ID]}
#---------Reading missing job index from a file ----------------
#job_pool=ATEXc_dharma_007877_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5_missingRuns.dat
#ji=$(($SLURM_ARRAY_TASK_ID+1)) #because the actual list of numbers starts from the second line
#job=$(sed ${ji}q';'d ${job_pool})
job=17359
field_path=LESb0p860_bins/
field_pat=ATEXc_dharma_007877_b0p860
n_grid=144 # Number of x/y grids
echo 'job: '${job} 'ji: '$ji 'Host: '$HOSTNAME 
python LES_MSCART_column_run.py 140 0 ${field_path} ${field_pat} ${job} ${MSCART} ${NPH} ${n_grid} 

#Random delay to avoid I/O overflow
T=$(( ( RANDOM % 10 )  + 1 ))  #Generates number between 1-10
T=$(( 12*$T ))                 #Scales random number
echo "Waiting for "$T          #Display message
sleep $T                       #Process waits/sleeps before finishing
