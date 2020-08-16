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
import matplotlib.pyplot as plt
import numpy as np
from cpnCommonlib import load_obj, savefig
from analysis_lib import vw_rets, S
from cpnbispectral_retrieval_LES import save_to_hdf5
if __name__ == "__main__":
    save_path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Retrievals/NJK_retrievals/data/V1/hdf5/'
    for case_name  in ['DYC','RIC','ATc','ATp']:
        for SZA in [120,140,160]:
            rotate1D='done'#'not_done'
            #case_name='DYC'; SZA=120; 
            res=''
            les = vw_rets(case_name,180-SZA,0)
            les.S=S(SZA,les)
            les.S.loadNJK(SZA,les,rotate1D=rotate1D,res=res)
            save_name1D=les.S.NJK1D.NJK_case_name
            save_name3D=les.S.NJK3D.NJK_case_name
            save_to_hdf5(vars(les.S.NJK1D),save_name1D,path=save_path)
            save_to_hdf5(vars(les.S.NJK3D),save_name3D,path=save_path)

                # fig,ax=plt.subplots()
            # ax.imshow(les.S.NJK1D.tau,origin='lower')
            # ax.set_title('SZA%03d_'%(SZA)+case_name)
            # savefig(fig,'SZA%03d_'%(SZA)+case_name+'.png')
