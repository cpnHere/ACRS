#!/bin/bash

#SBATCH --job-name=Kokhanovsky
#SBATCH --output=slurm-Kokhanovsky.out
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --mem=5000
##SBATCH --–ntasks-per-node=36
##SBATCH --dependency=afterok:7429488
##SBATCH --array=1-500%250

MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e3
band='2p13'
SZA=120
#job=$SLURM_ARRAY_TASK_ID

echo ${job} $HOSTNAME
#fn=$(python LES_MSCART_setup.py 160 0 RICO_dharma_005044_b0p860.nc 2>&1)
#echo ${fn}_NPH${NPH}_${job}.nc
#${MSCART} 10 ${NPH} 0 ${fn}.nml results/b${band}/${fn}_NPH${NPH}_${job}.nc
fn=$(python field_file_MSCART.py 2>&1)
#fn='dummy_b0p860_MSCART_SSA120_SAA000_VAA000plus'
sed 's/<SrZA>/'$SZA'/g' dummy.nml > tmp1
sed 's/<FFILE>/b'${band}'_fields\/'${fn}'.nc/g' tmp1 > b${band}_fields/${fn}SZA${SZA}.nml
${MSCART} 10 ${NPH} 0 b${band}_fields/${fn}SZA${SZA}.nml results/b${band}/${fn}SZA${SZA}_NPH${NPH}.nc

#Random delay to avoid I/O overflow
#T=$(( ( RANDOM % 10 )  + 1 ))  #Generates number between 1-10
#T=$(( 12*$T ))                 #Scales random number
#echo "Waiting for "$T          #Display message
#sleep $T                       #Process waits/sleeps before finishing
