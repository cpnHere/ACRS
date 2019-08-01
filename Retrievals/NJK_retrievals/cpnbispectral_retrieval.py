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
"""
import os,fnmatch,sys
import scipy.io as sio
import cpnLES_MSCARTlib as lib
import cpnRetrievalslib as ret_lib
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import pickle

sza='120'
d865SZAx=lib.POLCARTdset('d865SZAx','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/')
d865SZAx.readPOLCARThdf5('fractal_cld_b865re12ve05_x40km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                         dpath='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/results/fracb865/')
d2p13SZAx=lib.POLCARTdset('d2p13SZAx','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/')
d2p13SZAx.readPOLCARThdf5('fractal_cld_b2p13re12ve05_x40km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                          dpath='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/results/fracb2p13/')

dfrac=lib.LES_field(d865SZAx.fname.split('_MSCART',2)[0]+'.nc')#DYCOM field file for MSCART
dfrac.readLES_field()

frac1D865=lib.POLCARTdset('frac1D865','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb865_bins/')
frac1D865.readPOLCARThdf5('fractal_cld_b865re12ve05_x40km_1D_bins_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                              dpath='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/results/fracb865_bins/')
frac1D2p13=lib.POLCARTdset('frac1D2p13','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb2p13_bins/')
frac1D2p13.readPOLCARThdf5('fractal_cld_b2p13re12ve05_x40km_1D_bins_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                              dpath='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/results/fracb2p13_bins/')


VZA=0.0 #degrees

#obs_SWIR_set=np.array([0.33,0.33,0.33])
#obs_VNIR_set=np.array([0.5,0.5,0.5])



def NJK_matlab():
    sio.savemat('LUTs'+'.mat',{'LUT_VNIR':LUT_VNIR,\
            'LUT_SWIR':LUT_SWIR,'LUTtau':LUT.tau,\
            'LUTre':LUT.re,'obs_SWIR_set':obs_SWIR_set,\
            'obs_VNIR_set':obs_VNIR_set})
    
    os.system('matlab -nodisplay -nodesktop -r "run biret.m"')
    
    data=sio.loadmat('LUTs_out.mat')
    re_ret=data['re_ret']
    tau_ret=data['tau_ret']
    RFM=data['RFM']
    return re_ret,tau_ret,RFM

def save_HighestRes_pkl(pickle_name,D865,D2p13):
    with open(pickle_name, "wb") as f:
        pickle.dump(NJK_ret, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(D865,f)
        pickle.dump(D2p13,f)
        pickle.dump(dfrac,f)
    print(pickle_name+' SAVED!')


def save_agg_pkl(pickle_name_agg):
    with open(pickle_name_agg, "wb") as f:
        pickle.dump(NJK_ret, f, pickle.HIGHEST_PROTOCOL)
    print(pickle_name_agg+' SAVED!')

fdpath='../LUTs/'
fname='MODIS_LUT_extended_SZA%03d_RAA000.nc'%(180-int(sza))
LUT=ret_lib.Bispec_LUT(fdpath,fname)
LUT.readLUT()
    
mu_ix=LUT.find_mu_ix(VZA)
LUT_VNIR=LUT.I[0,:,1,:,mu_ix].T
LUT_SWIR=LUT.I[1,:,1,:,mu_ix].T


#fig_LUT,axLUT=plotLUT(VNIR_lut,SWIR_lut,reff_lut,tau_lut)

RT865_file='/umbc/xfs1/zzbatmos/users/charaj1/MSCART_RT_data/'+d865SZAx.fname.split('.',1)[0]+'.hdf5'
RT2p13_file='/umbc/xfs1/zzbatmos/users/charaj1/MSCART_RT_data/'+d2p13SZAx.fname.split('.',1)[0]+'.hdf5'
Physics_file='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_retrievals/NJK_retrievals/'+d865SZAx.fname.split('_MSCART',2)[0]+'.nc'
LUT_file='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_retrievals/LUTs/'+'MODIS_LUT_extended_SZA%03d_RAA000.nc'%(180-int(sza))

obVZA_ix=d865SZAx.find_obVZA_ix(VZA)

#Aggregating
res=np.array([10,20,50,100,250,500,1000,2000])
fti=1

for fti in np.arange(0,np.size(res)):
    if fti==0:
        #Highest resolution
        print('----------------------------------------')
        print('Highest resolution')
        obs_SWIR_set=d2p13SZAx.MeanPRad[obVZA_ix,:,0]
        obs_VNIR_set=d865SZAx.MeanPRad[obVZA_ix,:,0]
    else:
        print('----------------------------------------')
        print('%d m resolution'%res[fti])
        obs_SWIR_set=d2p13SZAx.avg_resolution(ftr=res[fti]/10)[obVZA_ix,:,0]
        obs_VNIR_set=d865SZAx.avg_resolution(ftr=res[fti]/10)[obVZA_ix,:,0]
    
    NJK_ret=ret_lib.NJK_retrievals(RT865_file,RT2p13_file,Physics_file,LUT_file)
    NJK_ret.VNIR_lut=LUT_VNIR
    NJK_ret.SWIR_lut=LUT_SWIR
    NJK_ret.reff_lut=LUT.re
    NJK_ret.tau_lut=LUT.tau
    NJK_ret.obs_SWIR_set=obs_SWIR_set
    NJK_ret.obs_VNIR_set=obs_VNIR_set
    NJK_ret.NJK_Python()

    if fti==0:
        pickle_name='fractal_cld_re12ve05_x40km_'+d2p13SZAx.fname.split('km_',1)[1].split('SZA',1)[0]+'SZA'+sza+'_VZA%0.3d'%VZA+'.pkl'
        pickle_highest_res_name=pickle_name
        sv=raw_input('Save as '+pickle_name+' ?')
        if sv=='y':
            save_HighestRes_pkl(pickle_name,d865SZAx,d2p13SZAx)
        else:
            print('Not saved!')
    else:
        pickle_name_agg=pickle_highest_res_name.split('.',1)[0]+'_%dm'%res[fti]+'.pkl'
        save_agg_pkl(pickle_name_agg)    

data = ["NJK_ret", "d865SZAx", "d2p13SZAx", "dfrac"]


#with open("pickle_test.pkl", "rb") as f:
#    NJK=pickle.load(f)
#    dd=pickle.load(f)
#re_ret,tau_ret,flag=retrieve_NJK(obs_VNIR_set[600],obs_SWIR_set[600],VNIR_lut,SWIR_lut,reff_lut,tau_lut)
   
#lt.plot(LUT_SWIR[0::5,50:80:5],LUT_VNIR[0::5,50:80:5])

xcens=(dfrac.xgrd[1:]+dfrac.xgrd[0:-1])/2
fig1,ax1=plt.subplots(2,1,sharex=True)
fig1_ttl=d865SZAx.fname.split('cld_',1)[1].split('re',1)[0]+\
    d2p13SZAx.fname.split('cld',1)[1].split('.',1)[0]+'_retrievals'
fig1.suptitle(fig1_ttl)
ax1[0].plot(xcens,np.squeeze(NJK_ret.re),'g')
ax1[1].plot(xcens,np.squeeze(NJK_ret.tau),'b')
ax1[0].set_ylabel(r'$r_e(\mu m)$')
ax1[0].set_title('Effective radius',size=10)
ax1[1].set_ylabel(r'$\tau$')
ax1[1].set_title('Optical depth',size=10)
ax1[1].set_xlabel('x (km)')

#From Dan's Matlab code
re_ret,tau_ret,RFM=NJK_matlab()
xcens=(dfrac.xgrd[1:]+dfrac.xgrd[0:-1])/2
fig1,ax1=plt.subplots(2,1,sharex=True)
fig1_ttl=d865SZAx.fname.split('cld_',1)[1].split('re',1)[0]+\
    d2p13SZAx.fname.split('cld',1)[1].split('.',1)[0]+'_retrievals'
fig1.suptitle(fig1_ttl)
ax1[0].plot(xcens,np.squeeze(re_ret),'g')
ax1[1].plot(xcens,np.squeeze(tau_ret),'b')
ax1[0].set_ylabel(r'$r_e(\mu m)$')
ax1[0].set_title('Effective radius',size=10)
ax1[1].set_ylabel(r'$\tau$')
ax1[1].set_title('Optical depth',size=10)
ax1[1].set_xlabel('x (km)')

def savefig(fig,fig_ttl):
    fig.savefig('figures/'+fig_ttl+'.png',format='png',dpi=200)
    print('figures/'+fig_ttl+'.png SAVED!')

#savefig(fig1,fig1_ttl)
plt.show()
 
