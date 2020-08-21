#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Fri Jul 10 20:58:06 2020
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************

- To write DISORT inputFiles for LUT generation
- To run DISORT
"""

import numpy as np
from scipy.special import legendre
import os, time
from cpnMielib import bulk_Mie,PSDs,MieSet
from cpnRetrievalslib import Bispec_LUT
from cpnCommonlib import progress_bar
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
    print('Computing bulk Pmat........')
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
def eq_tau(tau860,Qe860,Qebnd):
    '''
    To compute equivalent optical depth (at 0.860 um) (refer 08/30/2017 research log)
    '''
    taubnd=Qebnd/Qe860*tau860
    return taubnd
if __name__ == "__main__":
    
    t0 = time.time()
    inP = './inputFilesHigh0p860/'
    ouP = './results/LUTsHigh0p860/'
    write_input = True # True if inputFiles do not exist. This takes more time
    NMOM = 500 # Highest Legendre coefficient
    band = 1 # 1-0.86, 2-2.13
    
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
    
    if write_input:
        print('Reading Mie....')
        mie = MieSet('DYCOM2_dharma_008036_mie_470_860_2p13',path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_simulations/')
        mie.readMie()
        psd = PSDs("high_res_LUT",D=np.asarray(mie.d,dtype=float))
        start=time.time()
        ni,nj=re_list.size,cot_list.size
        i,j=0,0
        for re in re_list:
            j=0
            for COT in cot_list:
                #re = 12.0#np.arange(1,30,0.5) #12.0 # in um
                ve = 0.05
                #COT = 5.00
                SZA = 60.0
                jobid='LUT'+'_ve%0.2fre%0.2fCOT%0.2fSZA%0.4f'%(ve,re,COT,SZA)
                if not(os.path.isfile(inP+jobid.replace('.','p')+'_inputFile.dat')):
                    print("inputFIle not found! Generating....!")
                    psd.set_mod_gamma_norm(re,ve)
                    bM=bulk_Mie("highres_bulkMie_for_LES",psd=psd,Mie=mie)  
                    cal_bulk_Pmat(bM,band=band)
                    Ang = bM.Mie.ang
                    P11 = bM.bulk_P11
                    Mu = np.cos(np.radians(Ang))
                    bM.cal_bulk_albnQe()
                    SCA = bM.bulk_alb[band]
                    
                    '''
                    # check if the phase function is normalized
                    C =-0.5*np.trapz(P11,Mu) # negative sign is because the integration should be from -1 to 1, but Mu is from 1 to -1
                    print('normalization',C)
                    if C>0.99999:
                        print("Phase function is normalized.")
                    else:
                        print("Warning!!: Phase function is not normalized!!")
                    '''
                    #print('Computing Legendre polynomials....')
                    N=np.arange(0,NMOM)
                    g = []
                    for n in N:
                        f =  legendre(n)
                        L = f(Mu)
                        g.append(-0.5*np.trapz(P11*L,Mu))
                    g = np.array(g)
                    if g[0] > 0.99999: #replacing if g>1 value
                        g[0]=0.999999999
                    print('Writing input file....')
                    MU0 = np.cos(np.deg2rad(SZA))
                    input_content=np.concatenate([[eq_tau(COT,bM.bulk_Qe[1],bM.bulk_Qe[band])],[SCA],[MU0],g])
                    np.savetxt(inP+jobid.replace('.','p')+'_inputFile.dat',input_content)
                    print(inP+jobid.replace('.','p')+'_inputFile.dat saved!')
                    print("jobid: "+jobid.replace('.','p'))
                    #t1=time.time()
                    #print('%0.2f minutes elapsed!'%((t1-t0)/60))
                else:
                    print("InputFile already exists!. Going to run DISORT directly")
                
                print('Executing DISORT.....')
                os.system('./disort_driver.sh '+jobid.replace('.','p')+' '+inP+' '+ouP)
                #t2=time.time()
                #print('%d seconds!'%(t2-t1))
                progress_bar(i,j,ni,nj,start) 
                j+=1
            i+=1
        t3=time.time()
        print("%0.2f hours elpased!"%(t3/60/60-t0/60/60))
            
    else:
        print('Executing DISORT.....')
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
               
               
                os.system('./disort_driver.sh '+jobid.replace('.','p')+' '+inP+' '+ouP)
                #t2=time.time()
                #print('%d seconds!'%(t2-t1))
                pc=(i*float(nj)+(j+1.0))/nj/ni*100.0
                tm=time.time()/60-start/60
                print("\r%0.2f%% %d minutes remaining ..."%(pc,tm/pc*100-tm),end=' ')# 
                j+=1
            i+=1
            
            
            
