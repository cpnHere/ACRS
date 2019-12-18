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
08/30/2017
    -Note that, when we change the band, we have to change (extp3d) extinctions
    to be matched with the reference optical depth (say 0.865 optical depth, refer
    research log 08/30/2017)
08/31/2017:
    - To setup fractal cloud field file for MSCART'
09/29/2017:
    - Fractal cloud generation part added
12/13/2017:
    poddak parissam වෙන්ට meka කරද්දි. සේව් වෙන files monawada kiyala balanta HONDATA.
03/06/2018:
    Aerosol on top of a fractal cloud
    
"""
import numpy as np
import cpnLES_MSCARTlib as MSCARTlib
import cpnMicrophylib as Miclib
import pickle, os
from cpnMielib import MieSet, bulk_Mie, PSDs

class pkl_classes(object):
    def __init__(self,):
        self.class_names=[]
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d
            
dy=MSCARTlib.LES_field('OP_dharma_008036_full_3_26.nc',dpath='../Step_Cloud/')#DYCOM field file for MSCART
dy.readLES_field()

# Reading Fractal cloud case for Dan's thesis ################################
#fcld=Miclib.frac_physics('physical_scene_case_new.h5')   
#fcld.read()
#fcld.plot_cld()
##############################################################################

fcld=Miclib.frac_physics('fractal_StCu_b865re12ve05_x40km_H0p5.nc')

with open(fcld.fname.split('.',1)[0]+'.pkl', "rb") as f:
    pkl=pickle.load(f)
    _=pickle.load(f)
    fcld=pickle.load(f)

# Aluth cloud generate karanta ##############################################

#lwp=90#gm^{-2}
#re=12.0#microns
#fcld=Miclib.frac_physics('fractal_StCu_b865re12ve05_x40km.nc')
#
#fcld.f=0.37#fractal parameter fn=fc^n
#fcld.generate_fractal(re,lwp,xorder=12,xdist=40.0)
#fcld.plot_lwp_pdf()
#fcld.plot_cld()
#############################################################################


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

   #Optical properties from my mie runs
ss_file='/umbc/xfs1/zzbatmos/users/charaj1/mie_cpn/aer_cld_bulk_SS_for_color_ratio_LUT.pkl'
with open(ss_file, "rb") as f:
    pkl2=pickle.load(f)
    bulk_cloud=pickle.load(f)    
    bulk_aero=pickle.load(f)

re=bulk_cloud.psd.re 
ve=bulk_cloud.psd.ve;bnd=0.470
delH=0.5#km bhawuthika ganakama walakule
AOT=0.00
save_name='fractal_StCu_ACA_b%0.3f_x40km_H0p5_AOT%0.2f_H0p5'%(bnd,AOT)
pkl_stamp='fractal_StCu_ACA_b%0.3fre%dve%0.2f_x40km_H0p5_AOT%0.2f_H0p5'%(bnd,re,ve,AOT)
save_name=save_name.replace('.','p')+'.nc'
pkl_stamp=pkl_stamp.replace('.','p')
#field dimensions (in km)
if not(fcld.x[0]==0):
    xdim=np.append([0],fcld.x)
else:
    xdim=fcld.x
ydim                =np.linspace(0,1,2,dtype=float)
zdim                =np.linspace(0,1,3,dtype=float)

MSCART_field        =MSCARTlib.LES_field('MSCART_field.nc')
MSCART_field.nxp    =len(xdim)
MSCART_field.nyp    =len(ydim)
MSCART_field.nzp    =len(zdim)
MSCART_field.nx     =len(xdim)-1
MSCART_field.ny     =len(ydim)-1
MSCART_field.nz     =len(zdim)-1
MSCART_field.nz3    =2
MSCART_field.iz3l   =1
MSCART_field.npar1  =2
MSCART_field.npar3  =2
MSCART_field.nkd    =1
MSCART_field.nstokes=6
MSCART_field.nang   =1000
MSCART_field.npf    =2
MSCART_field.npsfc  =dy.npsfc
        
MSCART_field.xgrd   =xdim # xgrd(nxp)
MSCART_field.ygrd   =ydim #(nyp)
MSCART_field.zgrd   =zdim # (nzp)
MSCART_field.extp1d =dy.extp1d[:,0:2] # (npar1,nz)
MSCART_field.omgp1d =dy.omgp1d[:,0:2] # (npar1,nz)
MSCART_field.jpfp1d =dy.jpfp1d[:,0:2] # (npar1,nz)
MSCART_field.rpfp1d =dy.rpfp1d[:,0:2] # (npar1,nz)

def eq_tauCld(bnd,tau865):
    #finding equivalent optical depth (refer 08/30/2017 research log)
    Qe865=bulk_cloud.bulk_Qe[np.isclose(bulk_cloud.Mie.wvl,0.865,atol=0.01)][0]
    Qebnd=bulk_cloud.bulk_Qe[np.isclose(bulk_cloud.Mie.wvl,bnd,atol=0.01)][0]
    taubnd=Qebnd/Qe865*tau865
    return taubnd
    
def eq_tauAero(bnd,tau865):
    #finding equivalent optical depth (refer 08/30/2017 research log)
    Qe865=bulk_aero.bulk_Qe[np.isclose(bulk_aero.Mie.wvl,0.865,atol=0.01)][0]
    Qebnd=bulk_aero.bulk_Qe[np.isclose(bulk_aero.Mie.wvl,bnd,atol=0.01)][0]
    taubnd=Qebnd/Qe865*tau865
    return taubnd

MSCART_field.extp3d=np.zeros((MSCART_field.npar3,MSCART_field.nz3,MSCART_field.ny,MSCART_field.nx)) # (npar3,nz3,ny,nx)
MSCART_field.omgp3d=np.zeros_like(MSCART_field.extp3d) # (npar3,nz3,ny,nx)
#-------------------------------------------------------------------------------------
#bottom layer
MSCART_field.extp3d[0,0,0,:]=eq_tauCld(bnd,fcld.tau)/0.5#lower bin cloud
MSCART_field.omgp3d[0,0,0,:]=bulk_cloud.bulk_alb[np.isclose(bulk_cloud.Mie.wvl,bnd,atol=0.01)][0]#lower bin cloud albedo
MSCART_field.extp3d[1,0,0,:]= 0.0#lower bin aerosol
MSCART_field.omgp3d[1,0,0,:]=bulk_aero.bulk_alb[np.isclose(bulk_aero.Mie.wvl,bnd,atol=0.01)][0]#lower bin aerosol albedo
#-------------------------------------------------------------------------------------
#top layer
MSCART_field.extp3d[0,1,0,:]=0.0#upper bin clouds
MSCART_field.omgp3d[0,1,0,:]=bulk_cloud.bulk_alb[np.isclose(bulk_cloud.Mie.wvl,bnd,atol=0.01)][0]#upper bin cloud albedo
MSCART_field.extp3d[1,1,0,:]=eq_tauAero(bnd,AOT)/0.5#upper bin aerosol
MSCART_field.omgp3d[1,1,0,:]=bulk_aero.bulk_alb[np.isclose(bulk_aero.Mie.wvl,bnd,atol=0.01)][0] #upper bin aerosol albedo
#-------------------------------------------------------------------------------------
    
MSCART_field.jpfp3d=np.ones_like(MSCART_field.extp3d) # (npar3,nz3,ny,nx)
MSCART_field.jpfp3d[0,:,:,:]=1
MSCART_field.jpfp3d[1,:,:,:]=2
MSCART_field.rpfp3d=np.ones_like(MSCART_field.extp3d) # (npar3,nz3,ny,nx)
MSCART_field.rpfp3d[0,:,:,:]=1
MSCART_field.rpfp3d[1,:,:,:]=2

MSCART_field.absg1d=np.array([0.0,0.0]).reshape(1,2)# (nkd,nz) used an array instead of a tuple
MSCART_field.wkd=dy.wkd # (nkd)
MSCART_field.scatAng=bulk_cloud.Mie.ang # (nang)
    
PMatVal=np.zeros((MSCART_field.npf,MSCART_field.nang,MSCART_field.nstokes),dtype=float) # (npf,nang,nstokes)
#cloud PM
PMatVal[0,:,0]=bulk_cloud.bulk_P11[:,np.isclose(bulk_cloud.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[0,:,1]=bulk_cloud.bulk_P11[:,np.isclose(bulk_cloud.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[0,:,2]=bulk_cloud.bulk_P33[:,np.isclose(bulk_cloud.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[0,:,3]=bulk_cloud.bulk_P33[:,np.isclose(bulk_cloud.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[0,:,4]=bulk_cloud.bulk_P12[:,np.isclose(bulk_cloud.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[0,:,5]=bulk_cloud.bulk_P34[:,np.isclose(bulk_cloud.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
#aerosol PM
PMatVal[1,:,0]=bulk_aero.bulk_P11[:,np.isclose(bulk_aero.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[1,:,1]=bulk_aero.bulk_P11[:,np.isclose(bulk_aero.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[1,:,2]=bulk_aero.bulk_P33[:,np.isclose(bulk_aero.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[1,:,3]=bulk_aero.bulk_P33[:,np.isclose(bulk_aero.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[1,:,4]=bulk_aero.bulk_P12[:,np.isclose(bulk_aero.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
PMatVal[1,:,5]=bulk_aero.bulk_P34[:,np.isclose(bulk_aero.Mie.wvl,bnd,atol=0.01)].reshape(MSCART_field.nang)
MSCART_field.PMatVal=PMatVal
    
MSCART_field.tsfc2d=np.zeros((MSCART_field.ny,MSCART_field.nx),dtype='float')
MSCART_field.tsfc2d[:]=300.0 # (ny,nx) surface temperature (K)
    
MSCART_field.jsfc2d=np.zeros_like(MSCART_field.tsfc2d)
MSCART_field.jsfc2d[:]=1.0 # (ny,nx) 
    
MSCART_field.psfc2d=np.zeros((MSCART_field.ny,MSCART_field.nx,MSCART_field.npsfc),dtype='float')
MSCART_field.psfc2d[:]=0.0 # (ny,nx,npsfc)
    
#    save_name='fractal_cld_b865re%dve05_x40km_H0p5.nc'%re
if os.path.isfile(save_name):
    print('File name already exist!')
rd=raw_input('Save as '+save_name+' ?')
if rd=='y':
    MSCART_field.writeLES_field(save_name)
    #saving pkl files
    pkl=pkl_classes()
    pkl.class_names=['cpnLES_MSCARTlib.LES_field','cpnMicrophylib.frac_physics']
    pkl.object_order=['pkl','MSCART_field','fcld']
    pkl.SS_properties=ss_file
    pkl.stamp=pkl_stamp
    with open(save_name.split('.',1)[0]+'.pkl', "wb") as f:
        pickle.dump(pkl,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(MSCART_field, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(fcld,f,pickle.HIGHEST_PROTOCOL)
    print(save_name.split('.',1)[0]+'.pkl SAVED!')
else:
    print('File not saved')
#fig,ax=plt.subplots()
#ax.imshow(dy.extp3d[0,-1,:,:])
#plt.show()

