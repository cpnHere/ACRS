#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu May 25 14:45:42 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To handle MSCART field file(for step cloud case)

08/17/2017:
    Step cloud case (for Dan)
    -So far worked with only 0.86 band. Now going to add 2.13,3.75
    -djm_blk_SS_MODIS_extended.h5 contains single scattering properties of 6 
    bands, 57 re s and 11 ve s. Thus, can be used to 0.865 band too.
    -MODIS_band_0.865_polarimetric_microphysical_LUT.h5 contains single scattering
    properties of 0.865 band for 141 re s and 30 ve s.
08/30/2017:
    -Note that, when we change the band, we have to change (extp3d) extinctions
    to be matched with the reference optical depth (say 0.865 optical depth, refer
    research log 08/30/2017)
01/15/2018:
    - going to generate LUTs for ACA retrievals by using color ratio technique
    (Jethva et al. 2013)
    
"""
import numpy as np
import cpnLES_MSCARTlib as lib
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d
import sys
from cpnMielib import MieSet,PSDs,bulk_Mie
import pickle
class pkl_classes(object):
    def __init__(self,):
        self.class_names=[]

#Optical properties from my mie runs
with open('/umbc/xfs1/zzbatmos/users/charaj1/mie_cpn/aer_cldre5_bulk_SS_for_color_ratio_LUT.pkl', "rb") as f:
    pkl=pickle.load(f)
    bulk_cloud=pickle.load(f)    
    bulk_aero=pickle.load(f)
            
dy=lib.LES_field('OP_dharma_008036_full_3_26.nc')#DYCOM field file for MSCART
dy.readLES_field()

#xcens=(dy.xgrd[0:-1]+dy.xgrd[1:])/2
#ycens=(dy.ygrd[0:-1]+dy.ygrd[1:])/2
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
##X, Y, Z = axes3d.get_test_data(0.05)
#X,Y=np.meshgrid(xcens,ycens)
#Z=dy.omgp3d[0,0,:,:]
#cset = ax.scatter(X, Y, Z)
#ax.clabel(cset, fontsize=9, inline=1)
#

#field dimensions (in km)
xdim=np.linspace(0,1,2,dtype=float)
ydim=np.linspace(0,1,2,dtype=float)
zdim=np.linspace(0,1,3,dtype=float)

MSCART_field=lib.LES_field('Step_cld.nc')
MSCART_field.nxp=len(xdim)
MSCART_field.nyp=len(ydim)
MSCART_field.nzp=len(zdim)
MSCART_field.nx=len(xdim)-1
MSCART_field.ny=len(ydim)-1
MSCART_field.nz=len(zdim)-1
MSCART_field.nz3=2
MSCART_field.iz3l=1
MSCART_field.npar1=2
MSCART_field.npar3=2
MSCART_field.nkd=1
MSCART_field.nstokes=6
MSCART_field.nang=1000
MSCART_field.npf=2
MSCART_field.npsfc=dy.npsfc
        
MSCART_field.xgrd=xdim # xgrd(nxp)
MSCART_field.ygrd=ydim #(nyp)
MSCART_field.zgrd=zdim # (nzp)
MSCART_field.extp1d=dy.extp1d[:,0:2] # (npar1,nz)
MSCART_field.omgp1d=dy.omgp1d[:,0:2] # (npar1,nz)
MSCART_field.jpfp1d=dy.jpfp1d[:,0:2] # (npar1,nz)
MSCART_field.rpfp1d=dy.rpfp1d[:,0:2] # (npar1,nz)

MSCART_field.extp3d=np.zeros((MSCART_field.npar3,MSCART_field.nz3,MSCART_field.ny,MSCART_field.nx)) # (npar3,nz3,ny,nx)
MSCART_field.omgp3d=np.zeros_like(MSCART_field.extp3d) # (npar3,nz3,ny,nx)


#f=h5py.File('../microphysics/MODIS_band_0.865_polarimetric_microphysical_LUT.h5','r')
#from cpnMicrophylib import mphicDset
#mp865=mphicDset('../../microphysics/MODIS_band_0.865_polarimetric_microphysical_LUT.h5')
#mp865.readd()
#mp=mphicDset('../../microphysics/djm_blk_SS_MODIS_extended.h5')
#mp.readd(bands='multiple')

def eq_tau(bnd,tau860,bulk):
    #finding equivalent optical depth (refer 08/30/2017 research log)
    #bulk: cpnMielib.bulk_Mie object
    Qe860=bulk.bulk_Qe[bulk.Mie.wvl==0.860]
    Qebnd=bulk.bulk_Qe[bulk.Mie.wvl==bnd]
    taubnd=Qebnd/Qe860*tau860
    return taubnd

AOT=np.array([0.0,0.25,0.5,1.0,2.0,3.0,5.0])#0.860
COT=np.array([0,2,5,10,20,30,40,50])#0.860
bnd=float(sys.argv[1])
AOTi=int(sys.argv[2])
COTi=int(sys.argv[3])

#-------------------------------------------------------------------------------------
#bottom layer
MSCART_field.extp3d[1,0,:]=0.0#lower bin aerosol
MSCART_field.extp3d[0,0,:]=eq_tau(bnd,COT[COTi],bulk_cloud)/0.5#lower bin cloud
MSCART_field.omgp3d[1,0,:]=0.0#lower bin aerosol
MSCART_field.omgp3d[0,0,:]=bulk_cloud.bulk_alb[bulk_cloud.Mie.wvl==bnd]#lower bin cloud
#-------------------------------------------------------------------------------------
#top layer
MSCART_field.extp3d[0,1,:]=0.0#upper bin clouds
MSCART_field.extp3d[1,1,:]=eq_tau(bnd,AOT[AOTi],bulk_aero)/0.5#upper bin aerosol
MSCART_field.omgp3d[0,1,:]=0.0#upper bin clouds
MSCART_field.omgp3d[1,1,:]=bulk_aero.bulk_alb[bulk_aero.Mie.wvl==bnd]#upper bin aero
#-------------------------------------------------------------------------------------
MSCART_field.jpfp3d=np.zeros_like(MSCART_field.extp3d) # (npar3,nz3,ny,nx)
MSCART_field.jpfp3d[0,:,:,:]=1
MSCART_field.jpfp3d[1,:,:,:]=2
MSCART_field.rpfp3d=np.zeros_like(MSCART_field.extp3d) # (npar3,nz3,ny,nx)
MSCART_field.rpfp3d[0,:,:,:]=1
MSCART_field.rpfp3d[1,:,:,:]=2

MSCART_field.absg1d=np.array([0.0,0.0])# (nkd,nz) used an array instead of a tuple
MSCART_field.wkd=dy.wkd # (nkd)
MSCART_field.scatAng=bulk_aero.Mie.ang # (nang)

PMatVal=np.zeros((MSCART_field.npf,MSCART_field.nang,MSCART_field.nstokes),dtype=float) # (npf,nang,nstokes)
#cloud
PMatVal[0,:,0]=np.squeeze(bulk_cloud.bulk_P11[:,bulk_cloud.Mie.wvl==bnd][:])
PMatVal[0,:,1]=np.squeeze(bulk_cloud.bulk_P11[:,bulk_cloud.Mie.wvl==bnd][:])
PMatVal[0,:,2]=np.squeeze(bulk_cloud.bulk_P33[:,bulk_cloud.Mie.wvl==bnd][:])
PMatVal[0,:,3]=np.squeeze(bulk_cloud.bulk_P33[:,bulk_cloud.Mie.wvl==bnd][:])
PMatVal[0,:,4]=np.squeeze(bulk_cloud.bulk_P12[:,bulk_cloud.Mie.wvl==bnd][:])
PMatVal[0,:,5]=np.squeeze(bulk_cloud.bulk_P34[:,bulk_cloud.Mie.wvl==bnd][:])
#aerosols
PMatVal[1,:,0]=np.squeeze(bulk_aero.bulk_P11[:,bulk_aero.Mie.wvl==bnd][:])
PMatVal[1,:,1]=np.squeeze(bulk_aero.bulk_P11[:,bulk_aero.Mie.wvl==bnd][:])
PMatVal[1,:,2]=np.squeeze(bulk_aero.bulk_P33[:,bulk_aero.Mie.wvl==bnd][:])
PMatVal[1,:,3]=np.squeeze(bulk_aero.bulk_P33[:,bulk_aero.Mie.wvl==bnd][:])
PMatVal[1,:,4]=np.squeeze(bulk_aero.bulk_P12[:,bulk_aero.Mie.wvl==bnd][:])
PMatVal[1,:,5]=np.squeeze(bulk_aero.bulk_P34[:,bulk_aero.Mie.wvl==bnd][:])
MSCART_field.PMatVal=PMatVal

MSCART_field.tsfc2d=np.zeros((MSCART_field.ny,MSCART_field.nx),dtype='float')
MSCART_field.tsfc2d[:]=300.0 # (ny,nx) surface temperature (K)

MSCART_field.jsfc2d=np.zeros_like(MSCART_field.tsfc2d)
MSCART_field.jsfc2d[:]=1.0 # (ny,nx) 

MSCART_field.psfc2d=np.zeros((MSCART_field.ny,MSCART_field.nx,MSCART_field.npsfc),dtype='float')
MSCART_field.psfc2d[:]=0.0 # (ny,nx,npsfc)

MSCART_field.writeLES_field(str('ACA_color_ratio_meth_LUT_field_Cre5_AOT%0.2f_COT%0.2f_bnd%0.3f'%(AOT[AOTi],COT[COTi],bnd)).replace(".","p")+'.nc')
#fig,ax=plt.subplots()
#ax.imshow(dy.extp3d[0,-1,:,:])
#plt.show()
