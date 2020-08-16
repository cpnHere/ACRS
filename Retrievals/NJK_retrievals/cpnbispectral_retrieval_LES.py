#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Tue Aug 29 11:02:46 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To do retrievals by using Dans matlab code.
For the fractal cloud case.
09/16/2017:
    My python code also can be used now.
10/17/2017:
    New fractal cloud case retrievals for Zhibo's presentation
12/20/2017:
    Dan puluwan eka eka resolutions walata ekawara retrievals karanta. 
    [d865SZAx,d2p13SZAx] dan naththan [frac1D865,frac1D2p13] kiyala file eka purawatama
        wenas karanta oona issella. 
    Retrievals of the fractal cloud.
********************************************
Created on Fri May 18 11:02:46 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To do NJK retrievals of LES cases

01/18/2018:
    Moved to Taki from Maya as a new version.

"""
import numpy as np
import os, h5py
import matplotlib.pyplot as plt
from textwrap import wrap
from cpnLES_MSCARTlib import LES_case
from cpnRetrievalslib import doNJK_LES, NJK_retrievals
import cpnCommonlib as cpn
import sys
def save_to_hdf5(paras,save_name,replace=False,path=''):
    '''
    paras: vars(cpnRetrievalslib.NJK_retrievals)
    save_name: cpnRetrievalslib.NJK_retrievals.NJK_case_name
    '''
    if not(replace) and os.path.isfile(path+save_name):
        print("File already exists!")
    else:
        f = h5py.File(path+save_name+'.hdf5','w')
          
        f.attrs['RT865_file']=paras['RT865_file']
        f.attrs['RT2p13_file']=paras['RT2p13_file'] 
        f.attrs['Physics_file']=paras['Physics_file'] 
        f.attrs['LUT_file']=paras['LUT_file']
        f.attrs['Physics_file']=paras['Physics_file']
        f.attrs['NJK_case_name']=paras['NJK_case_name']
    
        PCentry=f.create_dataset('VNIR_lut',data=paras['VNIR_lut'])
        PCentry.dims[0].label='reff_lut'
        PCentry.dims[1].label='tau_lut'
        PCentry.attrs['units']='None'
        PCentry.attrs["long_name"]='LUT_VNIR_reflectances'
    
        PCentry=f.create_dataset('SWIR_lut',data=paras['SWIR_lut'])
        PCentry.dims[0].label='reff_lut'
        PCentry.dims[1].label='tau_lut'
        PCentry.attrs['units']='None'
        PCentry.attrs["long_name"]='LUT_SWIR_reflectances'
        
        PCentry=f.create_dataset('reff_lut',data=paras['reff_lut'])
        PCentry.attrs['units']='um'
        PCentry.attrs["long_name"]='Cloud_effective_radius'
    
        PCentry=f.create_dataset('tau_lut',data=paras['tau_lut'])
        PCentry.attrs['units']='None'
        PCentry.attrs["long_name"]='Cloud_optical_thickness'
    
        PCentry=f.create_dataset('obs_SWIR_set',data=paras['obs_SWIR_set'])
        PCentry.dims[0].label='xgrid'
        PCentry.dims[1].label='ygrid'
        PCentry.attrs['units']='None'
        PCentry.attrs["long_name"]='Observed_SWIR_reflectances'
    
        PCentry=f.create_dataset('obs_VNIR_set',data=paras['obs_VNIR_set'])
        PCentry.dims[0].label='xgrid'
        PCentry.dims[1].label='ygrid'
        PCentry.attrs['units']='None'
        PCentry.attrs["long_name"]='Observed_VNIR_reflectances'
    
        PCentry=f.create_dataset('re',data=paras['re'])
        PCentry.dims[0].label='xgrid'
        PCentry.dims[1].label='ygrid'
        PCentry.attrs['units']='um'
        PCentry.attrs["long_name"]='retrieved_cloud_effective_radius'
    
        PCentry=f.create_dataset('tau',data=paras['tau'])
        PCentry.dims[0].label='xgrid'
        PCentry.dims[1].label='ygrid'
        PCentry.attrs['units']='None'
        PCentry.attrs["long_name"]='retrieved_cloud_optical_thickness'
    
        PCentry=f.create_dataset('flag',data=paras['flag'])
        PCentry.dims[0].label='xgrid'
        PCentry.dims[1].label='ygrid'
        PCentry.attrs['units']='None'
        PCentry.attrs["long_name"]='Flags'
        PCentry.attrs['flag_def']=np.string_(['1-Successful retrievals',\
                                                  '2-Reflectances are outside of the LUT space',\
                                                  '3-AOT and COT guesses are outside of the LUT space',\
                                                  '4-AOT and COT retrievals are outside of the LUT space',\
                                                  '5-last two guesses have > 1e-7 difference'])
        PCentry.attrs["long_name"]='Flag_definitions'
        f.close()
def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'figures/')
if __name__=='__main__':
    cnm = str(sys.argv[1])
    sza = str(sys.argv[2])
    sdir = str(sys.argv[3])
    res = str(sys.argv[4])
    if res == 'None':
        res = None
    #cpn.setup_figures(plt)
    vza=0
    #DYCOMS-II-----------------------------------------------------------------
    '''
    les_new_path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/DYCOMS2/'
    les_name='DYCOMS2_dharma_008036';NPH3D='2e6'
    DYC0p860_sza=LES_case(les_name+'_b0p860_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH2e6.hdf5',les_new_path,\
                            RT1Dname=les_name+'_b0p860_MSCART_1D_bins_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5')
    DYC2p13_sza =LES_case(les_name+'_b2p13_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH2e6.hdf5',les_new_path,\
                            RT1Dname=les_name+'_b2p13_MSCART_1D_bins_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5')
    '''
    #ATEXp-----------------------------------------------------------------
    '''
    les_new_path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXp/'
    
    les_name='ATEXp_dharma_013067';NPH3D='1e6';NPH1D='1e5'
    DYC0p860_sza=LES_case(les_name+'_b0p860_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e6.hdf5',les_new_path,\
                            RT1Dname=les_name+'_b0p860_MSCART_1D_bins_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5')
    DYC2p13_sza =LES_case(les_name+'_b2p13_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e6.hdf5',les_new_path,\
                            RT1Dname=les_name+'_b2p13_MSCART_1D_bins_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5')
    '''

    #ATEXc-----------------------------------------------------------------
    '''
    les_new_path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXc/'
    
    les_name='ATEXc_dharma_007877';NPH3D='1e6';NPH1D='1e5'
    DYC0p860_sza=LES_case(les_name+'_b0p860_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e6.hdf5',les_new_path,\
                            RT1Dname=les_name+'_b0p860_MSCART_1D_bins_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5')
    DYC2p13_sza =LES_case(les_name+'_b2p13_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e6.hdf5',les_new_path,\
                            RT1Dname=les_name+'_b2p13_MSCART_1D_bins_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5')
    '''
    #RICO-----------------------------------------------------------------
    #'''
    DYC0p860_sza=LES_case(cnm+'_'+sza+'_b0p860',res=res)
    DYC2p13_sza =LES_case(cnm+'_'+sza+'_b2p13',res=res)
    #'''

    ret_save_dir=sdir

    #New runs
    #'''
    NJK_ret3D=doNJK_LES(DYC0p860_sza,DYC2p13_sza,sza,vza,check_DB=False)
    NJK_ret1D=doNJK_LES(DYC0p860_sza,DYC2p13_sza,sza,vza,RTdim='1D',check_DB=False)
    #Use save_to_hdf5() function to directly save into hdf5 files in future.
    cpn.save_obj(NJK_ret3D,sdir+NJK_ret3D.NJK_case_name)
    cpn.save_obj(NJK_ret1D,sdir+NJK_ret1D.NJK_case_name)
    #'''
    #Load previous runs
    '''
    NJK_ret3D=cpn.load_obj(ret_save_dir+\
            'NJK_retrievals_'+les_name+'_b2p13_b0p860_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH'+NPH3D)
    NJK_ret1D=cpn.load_obj(ret_save_dir+\
            'NJK_retrievals_'+les_name+'_b2p13_b0p860_MSCART_1D_bins_SZA'+sza+'_SAA000_VAA000plus_NPH'+NPH3D)
    '''
    '''
    #Plotting figure
    fig2,ax2=plt.subplots(3,2,figsize=(6,8),subplot_kw={'aspect':'equal'})
    fig2_ttl=NJK_ret3D.NJK_case_name
    fig2.suptitle("\n".join(wrap(fig2_ttl,50)))
    ctf1=ax2[0,0].contourf(DYC0p860_sza.xcens,DYC0p860_sza.ycens,NJK_ret3D.re ,np.linspace(0,25,50),extend='both')
    ctf2=ax2[1,0].contourf(DYC0p860_sza.xcens,DYC0p860_sza.ycens,NJK_ret3D.tau,np.linspace(0,40,50),extend='both')
    ctf3=ax2[2,0].imshow(NJK_ret3D.flag,origin='lower',cmap=plt.cm.get_cmap('jet',NJK_ret3D.flag.max()))    
    ctf4=ax2[0,1].contourf(DYC0p860_sza.xcens,DYC0p860_sza.ycens,NJK_ret1D.re.T ,np.linspace(0,25,50),extend='both')
    ctf5=ax2[1,1].contourf(DYC0p860_sza.xcens,DYC0p860_sza.ycens,NJK_ret1D.tau.T,np.linspace(0,40,50),extend='both')
    ctf6=ax2[2,1].imshow(NJK_ret1D.flag.T,origin='lower',cmap=plt.cm.get_cmap('jet',NJK_ret1D.flag.max()))      
    ax2[0,0].set_title(r'CER$(\mu m)$ on 3D RT')
    ax2[0,1].set_title(r'CER$(\mu m)$ on 1D RT')
    ax2[1,0].set_title(r'COT on 3D RT')
    ax2[1,1].set_title(r'COT on 1D RT')
    ax2[2,0].set_title('Flags')
    ax2[2,1].set_title('Flags')
    fig2.colorbar(ctf1,ax=ax2[0,0],ticks=np.arange(0,25.1,5))
    fig2.colorbar(ctf2,ax=ax2[1,0],ticks=np.arange(0,40.1,10))
    fig2.colorbar(ctf3,ax=ax2[2,0],ticks=[1,2,3,4,5])   
    fig2.colorbar(ctf4,ax=ax2[0,1],ticks=np.arange(0,25.1,5))
    fig2.colorbar(ctf5,ax=ax2[1,1],ticks=np.arange(0,40.1,10))
    fig2.colorbar(ctf6,ax=ax2[2,1],ticks=[1,2,3,4,5])    
    ax2[2,0].xaxis.set_ticks_position('bottom')
    ax2[2,1].xaxis.set_ticks_position('bottom')
    for ticks in ax2[2,0].get_xticklabels():
        ticks.set_rotation(45)
    for ticks in ax2[2,1].get_xticklabels():
        ticks.set_rotation(45)
    fig2.tight_layout(rect=[0,0,1,0.95])
    fig2.show()
    '''
    
    
    
