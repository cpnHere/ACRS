#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Fri Jul 10 20:58:06 2020
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************

- To re-write inputFiles.dat when first PMOM>1
"""

import numpy as np
from scipy.special import legendre
import os, time
from cpnMielib import bulk_Mie,PSDs,MieSet
from cpnRetrievalslib import Bispec_LUT
if __name__ == "__main__":
    
    t0 = time.time()
    inP = './inputFilesHigh2p13/'
    ouP = './inputFilesHigh2p13_2/'
    
    LUT=Bispec_LUT('/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_retrievals/LUTs/',\
                   'MODIS_LUT_extended_SZA%03d_RAA000.nc'%(60))
    LUT.readLUT()
    
    new_re=np.concatenate([np.array(LUT.re),np.array((LUT.re[1:]+LUT.re[0:-1])/2)])
    new_re.sort()
    new_re=np.unique(new_re.round(2))
    new_tau=np.concatenate([np.array(LUT.tau),np.array((LUT.tau[1:]+LUT.tau[0:-1])/2)])
    new_tau.sort()
    new_tau=np.unique(new_tau.round(2))

    re_list = new_re#np.unique(LUT.re.round(2))
    cot_list = new_tau #np.unique(LUT.tau.round(2))
    i,j=0,0
    ni,nj=re_list.size,cot_list.size
    start=time.time()
    for re in re_list:
        j=0
        for COT in cot_list:
            #re = 12.0#np.arange(1,30,0.5) #12.0 # in um
            ve = 0.05
            #COT = 5.00
            SZA = 60.0
            jobid='LUT'+'_ve%0.2fre%0.2fCOT%0.2fSZA%0.4f'%(ve,re,COT,SZA)
            content = np.loadtxt(inP+jobid.replace('.','p')+'_inputFile.dat')
            if content[2] > 0.99999:
                content[2]=0.999999999

            #print('Writing input file....')
            np.savetxt(ouP+jobid.replace('.','p')+'_inputFile.dat',content)
            pc=(i*float(nj)+(j+1.0))/nj/ni*100.0
            tm=time.time()-start
            print("\r%0.2f%% %d seconds remaining ..."%(pc,tm/pc*100-tm),end=' ')#     
            j+=1
        i+=1