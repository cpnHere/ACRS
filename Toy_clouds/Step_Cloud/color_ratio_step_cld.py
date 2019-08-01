#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Jan 18 11:42:30 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Color ratio method step cloud
1D and 3D reflectances
"""

from cpnLES_MSCARTlib import POLCARTdset, LES_field
import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
import cpnCommonlib as cpn
plt.rc('font', family='serif',weight='bold')


def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma

def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'color_ratio_results/figures/')
    

AOT='0p25'
s_the='120'
nmlpath='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Step_Cloud/'
fdpath='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Step_Cloud/color_ratio_results/'
fdpath_lut='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Step_Cloud/color_ratio_LUTs/'

name860='ACA_color_ratio_meth_step_field_AOT'+AOT+'_COT5p00_10p00_bnd0p860_MSCART_SZA'+s_the+'_SAA000_VAA000plus_NPH1e5_NSmx5000.nc'
name470='ACA_color_ratio_meth_step_field_AOT'+AOT+'_COT5p00_10p00_bnd0p470_MSCART_SZA'+s_the+'_SAA000_VAA000plus_NPH1e5_NSmx5000.nc'

#name860='ACA_color_ratio_meth_step_field_AOT'+AOT+'_COT5p00_10p00_bnd0p860_MSCART_SZA'+s_the+'_SAA000_VAA000plus_NPH1e6.nc'
#name470='ACA_color_ratio_meth_step_field_AOT'+AOT+'_COT5p00_10p00_bnd0p470_MSCART_SZA'+s_the+'_SAA000_VAA000plus_NPH1e6.nc'

d3D860=POLCARTdset('d3D860',nmlpath)
d3D470=POLCARTdset('d3D470',nmlpath)
d3D860.readMSCARTplus(name860,fdpath=fdpath,step=True)
d3D470.readMSCARTplus(name470,fdpath=fdpath,step=True)

#d3D860.savePOLCARThdf5(d3D860.fname.split('.',1)[0]+'.hdf5',dpath='color_ratio_results/',\
#                   pc_format='MSCART_3D',action='color_ratio_ACA')

dstep=LES_field(d3D860.fname.split('_MSCART',2)[0]+'.nc',dpath=nmlpath)#DYCOM field file for MSCART
dstep.readLES_field()

#column runs (1D)
def set1Dstep(AOT,bnd,s_the):
    cfree=POLCARTdset('cfree',nmlpath)
    cfree.readMSCARTplus('ACA_color_ratio_meth_LUT_field_AOT'+AOT+'_COT5p00_bnd0p'+bnd+'_MSCART_SZA'+s_the+'_SAA000_VAA000plus_NPH1e5.nc',\
                           fdpath=fdpath_lut,clm=True,step=True)
    
    cstep=POLCARTdset('cstep',nmlpath)
    cstep.readMSCARTplus('ACA_color_ratio_meth_LUT_field_AOT'+AOT+'_COT10p00_bnd0p'+bnd+'_MSCART_SZA'+s_the+'_SAA000_VAA000plus_NPH1e5.nc',\
                           fdpath=fdpath_lut,clm=True,step=True)
    
    step1Dto3D=POLCARTdset('step1Dto3D','generated_using_'+cfree.fdpath+cfree.fname)
    step1Dto3D.fname=cfree.fname.split('.',1)[0]+'_1Dto3D.nc'
    step1Dto3D.fdpath=cfree.fdpath
    step1Dto3D.VZA=cfree.VZA
    step1Dto3D.ScatA=cfree.ScatA
    step1Dto3D.MeanPRad=np.zeros_like(d3D860.MeanPRad)
    for i in range(0,200):
        step1Dto3D.MeanPRad[:,i,:]=cfree.MeanPRad
    for i in range(200,1200):
        step1Dto3D.MeanPRad[:,i,:]=cstep.MeanPRad
    for i in range(1200,1500):
        step1Dto3D.MeanPRad[:,i,:]=cfree.MeanPRad
    #step1Dto3D.cc3D='step_1D'
    #step1Dto3D.savePOLCARThdf5(step1Dto3D.fname.split('.',1)[0]+'.hdf5',dpath='color_ratio_results/',pc_format='MSCART_1D',action='color_ratio_ACA')
    return step1Dto3D
    
d1D860=set1Dstep(AOT,'860',s_the)
d1D470=set1Dstep(AOT,'470',s_the)

def plotStepRef(d3D,d1D,name):
    xcens=(dstep.xgrd[1:]+dstep.xgrd[0:-1])/2
    fig1,ax1=plt.subplots()
    fig1_ttl=name.split('.',1)[0]+'_reflectance'
    ax1.plot(xcens,np.mean(d3D.MeanPRad[d3D.VZA==0,:,0].T,axis=1),'k-')
    ax1.plot(xcens,np.mean(d1D.MeanPRad[d3D.VZA==0,:,0].T,axis=1),'k--')
    ax1.set_xlabel('x (km)')
    ax1.set_ylabel('I')
    fig1.suptitle("\n".join(wrap(fig1_ttl,70)))
    fig1.show()
    return fig1, fig1_ttl

fig1,fig1_ttl=plotStepRef(d3D470,d1D470,name470)
fig2,fig2_ttl=plotStepRef(d3D860,d1D860,name860)

fig3,ax3=plt.subplots()
fig3_ttl=name860.split('AOT',1)[0]+'AOT'+AOT+'_SZA'+s_the
fig3.suptitle(fig3_ttl)
xcens=(dstep.xgrd[1:]+dstep.xgrd[0:-1])/2
CR=np.mean(d3D470.MeanPRad[d3D470.VZA==0,:,0].T,axis=1)/np.mean(d3D860.MeanPRad[d3D860.VZA==0,:,0].T,axis=1)
CR10m=movingaverage(CR,10)
xcens10m=movingaverage(xcens,10)
CR1D=np.mean(d1D470.MeanPRad[d1D470.VZA==0,:,0].T,axis=1)/np.mean(d1D860.MeanPRad[d1D860.VZA==0,:,0].T,axis=1)
ax3.plot(xcens,CR,'k-',label='1m')
ax3.plot(xcens10m,CR10m,'r-',linewidth=1.0,label='10m (moving avg.)')
ax3.plot(xcens,CR1D,'k--',label='1D')
ax3.set_xlabel('x (km)')
ax3.set_ylabel('$R_{0.470}/R_{0.860}$')
ax3.legend(loc='best',prop={'size':10})
fig3.show()

fig4,ax4=plt.subplots(1,2,figsize=(12,5))
fig4_ttl='CR_and_R865_3'
ax4[0].text(0.01, 0.95, '(b)', transform=ax4[0].transAxes,size=15)
ax4[0].plot(xcens,CR,'k-',label='3D')
ax4[0].plot(xcens,CR1D,'k--',label='1D')
ax4[0].plot(xcens10m,CR10m,'r-',linewidth=1.0,label='3D(10m)')
ax4[0].set_ylabel('$R_{0.470}/R_{0.860}$')
ax4[0].legend(loc='lower right',prop={'size':11})
ax4[1].text(0.01, 0.95, '(c)', transform=ax4[1].transAxes,size=15)
ax4[1].plot(xcens,np.mean(d3D860.MeanPRad[d3D860.VZA==0,:,0].T,axis=1),'k-',label='3D')
ax4[1].plot(xcens,np.mean(d1D860.MeanPRad[d1D860.VZA==0,:,0].T,axis=1),'k--',label='1D')
#ax4[1].set_xlabel('x (km)')
ax4[1].set_ylabel(r'$R_{0.865}$')
ax4[1].legend(loc='lower left',prop={'size':11})
fig4.tight_layout()
fig4.show()
