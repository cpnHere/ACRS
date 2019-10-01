#!/bin/bash

#SBATCH --job-name=ATc2p13SZA160RD
#SBATCH --output=slurm-ATc2p13SZA1601e6.out
#SBATCH --partition=batch
#SBATCH --qos=medium
#SBATCH --mem=5000
##SBATCH --dependency=afterok:7429488
#SBATCH --array=0-79%80

MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e6
band='2p13'

job_pool=(  5   9  11  12  13  21  22  24  25  27  52  72  73  74  77  78  79  81  82  84  85  87  88  89  91  93  94  96  97 100 109 112 129 130 131 153 154 156 165 166 167 186 187 188 198 214 216 217 218 219 221 222 223 224 227 228 229 230 232 235 236 237 242 243 245 246 248 251 253 258 259 260 261 267 268 276 280 282 288 289)
job=${job_pool[$SLURM_ARRAY_TASK_ID]}

echo ${job} $HOSTNAME
fn=$(python LES_MSCART_setup.py 160 0 ATEXc_dharma_007877_b2p13.nc 2>&1)
echo ${fn}_NPH${NPH}_${job}.nc
${MSCART} 10 ${NPH} 0 ${fn}.nml results/b${band}/${fn}_NPH${NPH}_${job}.nc
#fn=$(python LES_MSCART_setup.py 180 2>&1)
#echo ${fn}
#${MSCART} 10 ${NPH} 0 ${fn}.nml results/${fn}_NPH${NPH}_${job}.nc

#mv slurm.out slurm${NPH}.out
