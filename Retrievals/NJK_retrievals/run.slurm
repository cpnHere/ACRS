#!/bin/bash

#SBATCH --job-name=NJK
#SBATCH --output=slurm%J-NJK.out
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --mem=14000
sdir='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Retrievals/NJK_retrievals/data/V1/'
for case_name in 'RICO'
do
	for sza in '120' '140' '160'
	do
		for res in '0p5km' '1p0km' '2p0km' '5p0km'
		do
			echo ${case_name} ${sza}	
			python cpnbispectral_retrieval_LES.py ${case_name} ${sza} ${sdir} ${res}
		done
	done
done
# data/V1/ : Version 1
