#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Jul 30 10:42:39 2020
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Field file generation for LES Look-up tables
"""

import numpy as np
import time,os,sys
import cpnLES_MSCARTlib as lib
from cpnMielib import bulk_Mie,PSDs,MieSet
from cpnRetrievalslib import Bispec_LUT
def cal_bulk_Pmat(bM,band=1):
    '''
    Copied from cpnMielib to avoid computing for,
        - all wavelengths,
        - all Pij s (only P11 is needed)
    band: 1 for 0.86 band and 2 for 2.13 band
    '''
    avP11=np.zeros(bM.Mie.ang.size)
    avP12=np.zeros(bM.Mie.ang.size)
    avP33=np.zeros(bM.Mie.ang.size)
    avP34=np.zeros(bM.Mie.ang.size)
    #print('Computing bulk Pmat........')
    fb=bM.psd.r**2*bM.Mie.qe[:,band].T*bM.Mie.alb[:,band].T*bM.psd.n_N
    ft=np.einsum('i,ij->ij',fb,bM.Mie.P11[:,band,:])
    avP11[:]=np.trapz(ft,bM.psd.r,axis=0)/np.trapz(fb,bM.psd.r)
    ft=np.einsum('i,ij->ij',fb,bM.Mie.P12[:,band,:])
    avP12[:]=np.trapz(ft,bM.psd.r,axis=0)/np.trapz(fb,bM.psd.r)
    
    ft=np.einsum('i,ij->ij',fb,bM.Mie.P33[:,band,:])
    avP33[:]=np.trapz(ft,bM.psd.r,axis=0)/np.trapz(fb,bM.psd.r)
    
    ft=np.einsum('i,ij->ij',fb,bM.Mie.P34[:,band,:])
    avP34[:]=np.trapz(ft,bM.psd.r,axis=0)/np.trapz(fb,bM.psd.r)
    bM.bulk_P11=avP11
    bM.bulk_P12=avP12
    bM.bulk_P33=avP33
    bM.bulk_P34=avP34
def cal_bulk_albnQe(bM,band=1):
    '''
    To compute average Qe and albedo
    Copied from cpnMielib
    band: 1 for 0.86 band and 2 for 2.13 band
    '''
    #print('Averaging for bulk SSA and Qe........')
    for i in range(0,bM.Mie.wvl.size):
        #fb=qe*n_V/r   
        fb=bM.Mie.qe[:,band]*4.0/3.0*np.pi*bM.psd.r**2*bM.psd.n_N
        #ft=fb*alb
        ft=fb*bM.Mie.alb[:,band]
        avalb=np.trapz(ft,bM.psd.r)/np.trapz(fb,bM.psd.r)
        fb=np.pi*bM.psd.r**2*bM.psd.n_N
        ft=bM.Mie.qe[:,band]*fb
        avQe=np.trapz(ft,bM.psd.r)/np.trapz(fb,bM.psd.r)
        return avalb, avQe
def eq_tau(tau860,Qe860,Qebnd):
    '''
    To compute equivalent optical depth (at 0.860 um) (refer 08/30/2017 research log)
    '''
    taubnd=Qebnd/Qe860*tau860
    return taubnd
def progress_bar(i,j,ni,nj,start):
    pc=(i*float(nj)+(j+1.0))/nj/ni*100.0
    tm=time.time()/60-start/60
    print("\r%0.2f%% %d minutes remaining ..."%(pc,tm/pc*100-tm),end=' ')# 
    
if __name__ == "__main__":
    band = 2# 0-0.470, 1-0.86, 2-2.13
    bnd = {0:'0p470',1:'0p860',2:'2p13'}
    ouP = './b'+bnd[band]+'_fields/'
            
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
    
    #print('Reading Mie....')
    mie = MieSet('DYCOM2_dharma_008036_mie_470_860_2p13',path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_simulations/')
    mie.readMie()
    psd = PSDs("high_res_LUT",D=np.asarray(mie.d,dtype=float))
    ni=re_list.size
    nj=cot_list.size
    start=time.time()
    i=0
    for re in re_list:
        j=0
        for COT in cot_list:
            ve = 0.05
            jobid='LUT'+bnd[band]+'_ve%0.2fre%0.2fCOT%0.2f'%(ve,re,COT)
            f_fileN = jobid.replace('.','p')+'.nc' 
            if not(os.path.isfile(ouP+f_fileN)):
                #print("inputFIle not found! Generating....!")
                #Computing scattering properties
                psd.set_mod_gamma_norm(re,ve)
                bM=bulk_Mie("highres_bulkMie_for_LES",psd=psd,Mie=mie) 
                cal_bulk_Pmat(bM,band=band)
                avalb,avQe=cal_bulk_albnQe(bM,band=band)
                _,Qe860=cal_bulk_albnQe(bM,band=1)
                
                #writing MSCART field file
                dy=lib.LES_field('dummy_field.nc')#DYCOM field file for MSCART
                dy.readLES_field()
                #phase matrix
                dy.nang=bM.Mie.ang.size
                PMatVal=np.zeros((dy.npf,dy.nang,dy.nstokes),dtype=float) # (npf,nang,nstokes)
                PMatVal[0,:,0]=bM.bulk_P11
                PMatVal[0,:,1]=bM.bulk_P11
                PMatVal[0,:,2]=bM.bulk_P33
                PMatVal[0,:,3]=bM.bulk_P33
                PMatVal[0,:,4]=bM.bulk_P12
                PMatVal[0,:,5]=bM.bulk_P34
                dy.PMatVal=PMatVal
                dy.scatAng=bM.Mie.ang
                dy.extp3d[:] = eq_tau(COT,Qe860,avQe)
                dy.omgp3d[:] = avalb
        
                dy.fname = ouP+f_fileN 
                dy.writeLES_field(ouP+f_fileN,showmess=False)
                #sys.exit(f_fileN.split('.',1)[0])
            progress_bar(i,j,ni,nj,start)
            j+=1
        i+=1
