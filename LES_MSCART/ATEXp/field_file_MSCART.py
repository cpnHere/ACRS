#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu May 7 14:45:42 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To generage MSCART field file from direct DHARMA output *.cdf files. Adopted
from fractal cloud field file generation.

-djm_blk_SS_MODIS_extended.h5 contains single scattering properties of 6 
 bands, 57 re s and 11 ve s. Thus, can be used to 0.865 band too.
-MODIS_band_0.865_polarimetric_microphysical_LUT.h5 contains single scattering
 properties of 0.865 band for 141 re s and 30 ve s.
-Note that, when we change the band, we have to change (extp3d) extinctions
 to be matched with the reference optical depth (say 0.865 optical depth, refer
 research log 08/30/2017)
-poddak parissam වෙන්ට meka කරද්දි. සේව් වෙන files monawada kiyala balanta HONDATA.    
"""
import numpy as np
import pickle, os
from cpnLES_MSCARTlib import DHARMA_onmp,LES_field
from cpnCommonlib import pkl_classes

def edges(x):
    cens     =x
    cens_shif=np.roll(x,-1)
    cens_shif[-1]=cens_shif[-1]*2+cens[-1]
    edges=(cens+cens_shif)/2
    edges=np.append([0.0],edges)
    return edges
            
dy=LES_field('OP_dharma_008036_full_3_26.nc',dpath='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/')#DYCOM field file for MSCART
dy.readLES_field()

cloud=DHARMA_onmp('dharma_013067.cdf','/umbc/xfs1/zzbatmos/users/mild1/Research/LES/LES_Satellite_Simulator_DJM/LES_Simulations/Ackerman_ATEX_highres/CCN600/');
cloud.readDHARMA()
extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=cloud.setupMieOut('DYCOM2_dharma_008036_mie_470_860_2p13',lesCname='ATEXp',mie_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/')

bnd=2.13
save_name='ATEXp_dharma_013067_b%0.2f'%(bnd)
pkl_stamp='ATEXp_dharma_013067_b%0.2f'%(bnd)
save_name=save_name.replace('.','p')+'.nc'
pkl_stamp=pkl_stamp.replace('.','p')

#field dimensions (in km)
xdim                =edges(cloud.x)*1e-3#km
ydim                =edges(cloud.y)*1e-3#km
zdim                =cloud.DHARMA.zw*1e-3#km

MF        =LES_field('MSCART_field.nc')
MF.nxp    =len(xdim)
MF.nyp    =len(ydim)
MF.nzp    =len(zdim)
MF.nx     =len(xdim)-1
MF.ny     =len(ydim)-1
MF.nz     =len(zdim)-1
MF.nz3    =MF.nz #Number of horizontally inhomogeneous layers
MF.iz3l   =1               #Lowest horozontally inhomogeneous layer index
MF.npar1  =2
MF.npar3  =25
MF.nkd    =1
MF.nstokes=6
MF.nang   =1801
MF.npf    =25
MF.npsfc  =7 #Maximum Number of surface BRDF parameters
        
MF.xgrd   =xdim # xgrd(nxp)
MF.ygrd   =ydim #(nyp)
MF.zgrd   =zdim # (nzp)
MF.extp1d =np.zeros((MF.npar1,MF.nz),dtype=float) # (npar1,nz)
MF.omgp1d =np.ones((MF.npar1,MF.nz),dtype=float) # (npar1,nz)
MF.jpfp1d =np.zeros((MF.npar1,MF.nz),dtype=float) # (npar1,nz)
MF.rpfp1d =np.zeros((MF.npar1,MF.nz),dtype=float) # (npar1,nz)

MF.extp3d=np.squeeze(extp3d[np.isclose(cloud_mie.wvl,bnd,atol=0.001),:])#.swapaxes(2,3) # (npar3,nz3,nx,ny)
MF.omgp3d=np.squeeze(omgp3d[np.isclose(cloud_mie.wvl,bnd,atol=0.001),:])#.swapaxes(2,3) # (npar3,nz3,nx,ny)
    
MF.jpfp3d=np.ones_like(MF.extp3d) # (npar3,nz3,nx,ny)
for i in np.arange(0,MF.npar3,1):
    MF.jpfp3d[i,:,:,:]=i+1
MF.rpfp3d=np.ones_like(MF.extp3d) # (npar3,nz3,nx,ny)
for i in np.arange(0,MF.npar3,1):
    MF.rpfp3d[i,:,:,:]=i+1

MF.absg1d=np.zeros((1,MF.nz),dtype=float)# (nkd,nz) used an array instead of a tuple
MF.wkd=np.array([1]) # (nkd)
MF.scatAng=cloud_mie.ang # (nang)
    
PMatVal=np.zeros((MF.npf,MF.nang,MF.nstokes),dtype=float) # (npf,nang,nstokes)
#cloud PM
PMatVal[:,:,0]=P11[np.isclose(cloud_mie.wvl,bnd,atol=0.001)].reshape(MF.npf,MF.nang)
PMatVal[:,:,1]=P11[np.isclose(cloud_mie.wvl,bnd,atol=0.001)].reshape(MF.npf,MF.nang)
PMatVal[:,:,2]=P33[np.isclose(cloud_mie.wvl,bnd,atol=0.001)].reshape(MF.npf,MF.nang)
PMatVal[:,:,3]=P33[np.isclose(cloud_mie.wvl,bnd,atol=0.001)].reshape(MF.npf,MF.nang)
PMatVal[:,:,4]=P12[np.isclose(cloud_mie.wvl,bnd,atol=0.001)].reshape(MF.npf,MF.nang)
PMatVal[:,:,5]=P34[np.isclose(cloud_mie.wvl,bnd,atol=0.001)].reshape(MF.npf,MF.nang)
MF.PMatVal=PMatVal
    
MF.tsfc2d=np.zeros((MF.ny,MF.nx),dtype='float')
MF.tsfc2d[:]=300.0 # (nx,ny) surface temperature (K)
    
MF.jsfc2d=np.zeros_like(MF.tsfc2d)
MF.jsfc2d[:]=1.0 # (nx,ny) 
    
MF.psfc2d=np.zeros((MF.ny,MF.nx,MF.npsfc),dtype='float')
MF.psfc2d[:]=0.0 # (nx,ny,npsfc)
    
#    save_name='fractal_cld_b865re%dve05_x40km_H0p5.nc'%re
if os.path.isfile(save_name):
    print('File name already exist!')
rd=raw_input('Save as '+save_name+' ?')
if rd=='y':
    MF.writeLES_field(save_name)
    #saving pkl files
    pkl=pkl_classes()
    pkl.class_names=['cpnLES_MSCARTlib.LES_field','DHARMA_LES_converter.DHARMA_onmp',\
    'DHARMA_LES_read.DHARMA','cpnMielib.MieSet']
    pkl.object_order=['pkl','MSCART_field','cloud']
    pkl.stamp=pkl_stamp
    with open(save_name.split('.',1)[0]+'.pkl', "wb") as f:
        pickle.dump(pkl,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(MF, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(cloud,f,pickle.HIGHEST_PROTOCOL)
    print(save_name.split('.',1)[0]+'.pkl SAVED!')
else:
    print('File not saved')
#fig,ax=plt.subplots()
#ax.imshow(dy.extp3d[0,-1,:,:])
#plt.show()

