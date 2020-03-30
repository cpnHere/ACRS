#!/bin/bash

#SBATCH --job-name=DYPlRet3Dp860Rnt
#SBATCH --output=slurm%J-DYPlRet3Dp860Rnt.out
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --mem=14000
#SBATCH --array=0-2%3
job_pool=(120 140 160)
job=${job_pool[$SLURM_ARRAY_TASK_ID]}
python pol_cloud_ret.py 0p860 $job 3D Breon full _pol_ret_V7 DY nt
#Look the python file to choose different initila guesses.
#lscpu
#mpirun -n 16 python comm.py
#pol_ret_curvefit_V1_sml: small LUT (lesser number of re and ve)
#pol_ret_curvefit_V2: Breon et. al. 2015 with Pol_ret_PM_library_0p860_2p13.hdf5 library
#pol_ret_curvefit_V3: Breon et. al. 2015 with Pol_ret_PM_library_0p860_2p13_V2.hdf5 library   
#pol_ret_curvefit_V4: Breon et. al. 2005 with Pol_ret_PM_library_0p860_2p13_V2.hdf5 'BreonMod'
#pol_ret_mean_V1    : Breon et. al. 2005 with Pol_ret_PM_library_0p860_2p13_V2.hdf5 'Breon'
#_pol_ret_full       : Breon et. al. 2005 with Pol_ret_PM_library_0p860_2p13_V2.hdf5 'Breon'
#_pol_ret	    : First official product
#*MPI		    : with MPI
#_pol_ret_MPIV2	    : Contains the two fitted pattern
#_pol_ret_V3        : Pol_ret_PM_library_0p860_2p13_V2.hdf5 LUT library
#_pol_ret_V4	    : Pol_ret_PM_library_0p860_2p13_V2.hdf5 wiht 'Breon'
#_pol_ret_V5	    : Pol_ret_PM_library_0p860_2p13_V3.hdf5 with 'Breon' including primary bow (** a bug. 4*gemet done twice!!)
#_pol_ret_V6		: Pol_ret_PM_library_0p860_2p13_V3.hdf5 with 'Breon' including primary bow (bug fixed)
#_pol_ret_V7		: Including coarser resolutions. nt results must agree with V6 results.
#sys.argv[4]: 
#	Breon   : F=aP(x)+b*x+c
#	BreonMod: F=aP(x)+b*cos(x)**2+c
#sys.argv[5]:
#	full/mean/pixel level retrievals
#sys.argv[6]: A tail to outputs
#sys.argv[7]: DY-DYCOMS2, RC-RICO, AC-ATEXc, AP-ATEXp
#sys.argv[8]: resolution Ex. 'nt' for native. '0p1km' format 
