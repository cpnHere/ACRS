#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Aug 13 16:08:39 2020
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Converting to hdf5 files
"""
import h5py,os
from analysis_lib import S, vw_rets
from  pol_cloud_ret import save_to_hdf5


if __name__=="__main__":
    for SZA in [120,140,160]:
        for casen in ['DYC','RIC','ATp','ATc']:
             #SZA=140;casen='DYC'
             les = vw_rets(casen,180-SZA,0)
             les.S=S(SZA,les)
             pol_tail='_Breon_full_pol_ret_V4_MPI'.replace('V4','V7')
             les.S.loadPol(SZA,les,pol_tail,res='')
             path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Retrievals/Pol_retrievals/results/pol_ret_V7/hdf5/'
             fnm = les.S.get_RT_filename(les.name,str(SZA),'0p860','3D')
             les.S.NPH3D = fnm.rsplit('/',1)[1].split('.',1)[0].split('NPH',1)[1]
             fnm = les.S.get_RT_filename(les.name,str(SZA),'0p860','1D')
             les.S.NPH1D = fnm.rsplit('/',1)[1].split('.',1)[0].split('NPH',1)[1]
             filename3D=les.name+'_'+les.fname.split('.',1)[0]+'_b0p860_MSCART_SZA%d_SAA000_VAA000plus_NPH'%(SZA)+les.S.NPH3D+pol_tail
             filename1D=les.name+'_'+les.fname.split('.',1)[0]+'_b0p860_MSCART_1D_bins_SZA%d_SAA000_VAA000plus_NPH'%(SZA)+les.S.NPH3D+pol_tail
             
             save_to_hdf5(les.S.Pol3D,les,path,filename3D,replace=True)
             save_to_hdf5(les.S.Pol1D,les,path,filename1D,replace=True)
