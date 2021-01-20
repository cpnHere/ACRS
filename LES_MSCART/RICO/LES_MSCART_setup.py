#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Mar 27 11:58:56 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
LES_MSCART Moderator
2017/04/03: 
    dharma_MSCART_dummy.nml will be used as a dummy.
    Newly created *.nml file will be a system output
2017/04/29:
    Now expect to handle both single column and 3D cases.
    (LES_MSCART_column_run.py uses this file. So be careful when edit)
2018/05/10:
    To run LES cases

"""
import numpy as np
import sys


#writing *.nml file for LES_MSCART (Dharma)
s_the=float(sys.argv[1])
vza=float(sys.argv[2])
field_file=str(sys.argv[3])
#field_file='OP_dharma_008036_full_3_26_bins/OP_dharma_008036_full_3_26_binsy50x50.nc'
dummy='dharma_MSCART_dummy.nml'
src_the=s_the# Source polar angle (in right handed upward z coordinate system)
src_phi=0.0# Source azimuth angle
n_PRad_xgrd = 144# size of the grids as in field_file
n_PRad_ygrd = 144# size of the grids as in field_file
n_PRad_zlv = 1 #Number of height values for detectors
n_PRad_the = 1 #Number of viewing polar angles
n_PRad_phi = 72 # Number of veiwing azimuth angles
min_PRad_the = vza # Minimum viewing polar angle
max_PRad_the = vza+1.0 # Max viewing polar angle
min_PRad_phi = 0.0 # Minimum viewing azimuth angle
max_PRad_phi = 360.0 # Max viewing azimuth angle
vec_PRad_zlv1 = 3.000000 # position of the detector
NmaxScat=8000

fr=open(dummy,'r')
data=fr.readlines()[:]
fr.close()

SZA=int(src_the)
SAA=int(src_phi)
'''
# When only using 2 azimuth angles (0,180) (Princial plane)
VAA=int(min_PRad_phi)
out_file=field_file[:-3]+'_MSCART_SZA%03d'%SZA+'_SAA%03d'%SAA+\
    '_VAA%03dplus'%VAA
'''
#Multiple VAAs (larger than 2)
out_file=field_file[:-3]+'_MSCART_SZA%03d'%SZA+'_SAA%03d'%SAA+\
    '_Nvaa72VZA%03d'%vza

fw=open(out_file+'.nml',"w")
for i in np.arange(0,23,1): fw.write(data[i])
wstr=' NSMAX  =  '+str(NmaxScat)+',\n'
fw.write(wstr)
for i in np.arange(24,36,1): fw.write(data[i])
wstr=' field_file = "'+field_file+'",\n'
fw.write(wstr)
for i in np.arange(37,42,1): fw.write(data[i])
fw.write(' src_the = '+str(src_the)+',\n')
fw.write(' src_phi = '+str(src_phi)+',\n')
fw.write(data[44])
fw.write(' n_PRad_xgrd = '+str(n_PRad_xgrd)+',\n')
fw.write(' n_PRad_ygrd = '+str(n_PRad_ygrd)+',\n')
fw.write(' n_PRad_zlv = '+str(n_PRad_zlv)+',\n')
fw.write(' n_PRad_the = '+str(n_PRad_the)+',\n')
fw.write(' n_PRad_phi = '+str(n_PRad_phi)+',\n')
fw.write(' min_PRad_the = '+str(min_PRad_the)+',\n')
fw.write(' max_PRad_the = '+str(max_PRad_the)+',\n')
fw.write(' min_PRad_phi = '+str(min_PRad_phi)+',\n')
fw.write(' max_PRad_phi = '+str(max_PRad_phi)+',\n')
fw.write(' vec_PRad_zlv(1) = '+str(vec_PRad_zlv1)+',\n')
fw.write(data[55])

fw.close()
#print(out_file+'.nml SAVED !!')
sys.exit(out_file)
