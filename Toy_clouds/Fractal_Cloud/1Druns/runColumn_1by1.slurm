#!/bin/bash

#SBATCH --job-name=fcsing0p470
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --mem=2500
#SBATCH --array=0-8%9
MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
job_pool=(573  577  579  916 1078 2140 2149 3147 3189)
job=${job_pool[$SLURM_ARRAY_TASK_ID]}
field_path=fracb0p470_bins/
field_pat=fractal_StCu_ACA_b0p470_x40km_H0p5_AOT0p25_H0p5
echo ${job} $HOSTNAME
python LES_MSCART_column_run.py 140 0 ${field_path} ${field_pat} ${job} ${MSCART} ${NPH}
