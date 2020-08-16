#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Tue Feb  4 22:56:26 2020
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To convert *.pkl NJK output to hdf5
Only COT, CER and flags
"""
import sys, os, h5py
from cpnCommonlib import load_obj
if __name__ == "__main__":
    file = sys.argv[1]
    paras = vars(load_obj(file,encoding='latin1'))
    save_name = file+'.hdf5'
    if os.path.isfile(save_name):
        print('hdf5 file already exist!!')
    else:
        f = h5py.File(save_name,'w')
          
        f.attrs['LUT_file']=paras['LUT_file']
        f.attrs['Physics_file']=paras['Physics_file']
        f.attrs['RT865_file']=paras['RT865_file']
        f.attrs['RT2p13_file']=paras['RT2p13_file']
        f.attrs['NJK_case_name']=paras['NJK_case_name']
        
        PCentry=f.create_dataset('tau',data=paras['tau'])
        PCentry.dims[0].label='x'
        PCentry.dims[1].label='y'
        PCentry.attrs['units']='None'
        PCentry.attrs["long_name"]='Cloud_optical_thickness'
        
        PCentry=f.create_dataset('re',data=paras['re'])
        PCentry.dims[0].label='x'
        PCentry.dims[1].label='y'
        PCentry.attrs['units']='um'
        PCentry.attrs["long_name"]='Cloud_effective_radius'
        
        PCentry=f.create_dataset('flag',data=paras['flag'])
        PCentry.dims[0].label='x'
        PCentry.dims[1].label='y'
        PCentry.attrs['units']='None'
        PCentry.attrs["long_name"]='Flags'
        f.close()
        print(save_name+' saved!')
