#!/bin/bash
# Chamara Rajapakshe (2020)
# To run MSCART column runs to generate high resolution LUTs
#---------------------------------------------------------------
#SBATCH --job-name=LUT0p860NPH1e5
#SBATCH --output=slurm-LUT0p860NPH1e5.out
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --mem=5000
##SBATCH --–ntasks-per-node=36
##SBATCH --dependency=afterok:7429488
##SBATCH --array=1-5340%250
#SBATCH --array=1-23%23

ji=$(($SLURM_ARRAY_TASK_ID))
#ji=$(($SLURM_ARRAY_TASK_ID+5340))
#ji=$(($SLURM_ARRAY_TASK_ID+10680))
#ji=$(($SLURM_ARRAY_TASK_ID+16020))

MSCART='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/mscart/build/intel-release-nocaf/mscart/MSCART'
NPH=1e5
band='0p860'
SZA=120
#job=$SLURM_ARRAY_TASK_ID

echo ${job} $HOSTNAME
#---------Reading job id list from a file ----------------
#job_pool=LUT${band}_ve0p05rexCOTx.dat
job_pool=b${band}_missing.dat
fn=$(sed ${ji}q';'d ${job_pool})
#fn='dummy_b0p860_MSCART_SSA120_SAA000_VAA000plus'
sed 's/<SrZA>/'$SZA'/g' dummy.nml > tmp1
sed 's/<FFILE>/b'${band}'_fields\/'${fn}'.nc/g' tmp1 > b${band}_fields/${fn}SZA${SZA}_SAA000_VAA000plus.nml
${MSCART} 10 ${NPH} 0 b${band}_fields/${fn}SZA${SZA}_SAA000_VAA000plus.nml results/b${band}/${fn}SZA${SZA}_SAA000_VAA000plus_NPH${NPH}.nc

Random delay to avoid I/O overflow
T=$(( ( RANDOM % 10 )  + 1 ))  #Generates number between 1-10
T=$(( 12*$T ))                 #Scales random number
echo "Waiting for "$T          #Display message
sleep $T                       #Process waits/sleeps before finishing
