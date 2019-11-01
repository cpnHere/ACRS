#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Fri Apr 28 17:22:48 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
LES_MSCART column by columns runs moderator
Dependancies:
    LES_MSCART_setup.py
05/01/2017:
    set the job id of the slurm from 0 to 16383(16384-1)
    This script will create *.nml files for each column and run MSCART.
    Setup the dtector specs in LES_MSCART_setup.py
05/08/2017:
    Adde _NPH in lin 41 to handle some readings in the LES_MSCARTlib
05/15/2018:
    To run LES MSCART 1D
"""
import os,fnmatch,sys

s_the=float(sys.argv[1])
vphi=float(sys.argv[2])
field_path=str(sys.argv[3])
field_pat=str(sys.argv[4])
job=int(sys.argv[5])
MSCART=str(sys.argv[6])
NPH=str(sys.argv[7])
n_grid=int(sys.argv[8])
over_write=True

print('Assume n_xgrd=n_ygrd=%d'%n_grid)
i=str(int(job/n_grid))
j=str(job%n_grid)
fn_pattern=field_pat+'y'+i+'x'+j+'.nc'
for f in os.listdir(field_path):
    if fnmatch.fnmatch(f,fn_pattern):
        os.system("echo "+f)
#        cmnd_ln="fn=$(python LES_MSCART_setup.py "+str(vphi)+" "+field_path+f+\
#        " 2>&1)\n echo ${fn}\n "+MSCART+" 10 "+NPH+" 0 "+\
#        "${fn}.nml results/${fn}_"+NPH+".nc"
#        print(cmnd_ln)
        run_setup="fn=$(python LES_MSCART_column_setup.py "+str(s_the)+" "+\
                        str(vphi)+" "+field_path+f+" 2>&1)\n echo ${fn}\n "
        run_MSCART=MSCART+" 10 "+NPH+" 0 ${fn}.nml results/${fn}_NPH"+NPH+".nc"
        check_n_run="test=$(python checking_previous_runs.py results/${fn}_NPH"+NPH+".nc 2>&1)\n"+\
        "echo checking_previous_runs.py gave ${test}\n"+\
        "if [ $test -eq 1 ]; then\n"+\
        "\t echo Output file already exist\n"+\
        "else\n"+\
        "\t  "+run_MSCART+"\nfi"
        if over_write:
            print('Going to overite the existing files....')
            os.system(run_setup+run_MSCART)
        else:
            #os.system('${fn}')
            print('Only execute if a result .nc file does not exist.')
            os.system(run_setup+check_n_run)

#        os.path.exists
#        sys.exit(out_file)
#echo ${fn}
#${MSCART} 10 ${NPH} 0 ${fn}.nml results/${fn}_NPH${NPH}.nc
