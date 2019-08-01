#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Wed Mar 20 10:32:54 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
References:
    /umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Step_Cloud/step_cld.py
For #NASA_RSTheory_grant #NSPIRES
Did not do actual retrievals yet
"""

import os,fnmatch,sys
import scipy.io as sio
from cpnLES_MSCARTlib import POLCARTdset, LES_field, scat_ang
from cpnRetrievalslib import Bispec_LUT
import cpnCommonlib as cpn
from textwrap import wrap
import numpy as np
import matplotlib.pyplot as plt
import pickle
def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,path='figures/for_NSPIRES_ZZ/')
def setup1D_3Deq(band,sza,MeanPRad3D):
    '''
    Use two columns step cloud simulation results to generate equivalent 3D MeanPRad
    MeanPRad3D POLCARTdset.MeanPRad
    band = '865' or '2p13'
    sza = '120' ...
    return POLCARTdset,cfree(optically thin part),cstep(optically thick part)
    '''
    DIR='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Step_Cloud/'
    cfree=POLCARTdset('cfree',DIR+'1Druns/')
    cfree.readMSCARTplus('step_cld_b'+band+'re10ve02_x1Dex0p1_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.nc',\
                           fdpath=DIR+'1Druns/results/stepb'+band+'/',clm=True,step=True)
    cstep=POLCARTdset('cstep',DIR+'1Druns/')
    cstep.readMSCARTplus('step_cld_b'+band+'re10ve02_x1Dex10_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.nc',\
                           fdpath=DIR+'1Druns/results/stepb'+band+'/',clm=True,step=True)
    
    step1Dto3Deq=POLCARTdset('step1Dto3Deq','generated_using_'+cfree.fdpath+cfree.fname)
    step1Dto3Deq.fname=cfree.fname.split('.',1)[0]+'_1Dto3D.nc'
    step1Dto3Deq.fdpath=cfree.fdpath
    step1Dto3Deq.VZA=cfree.VZA
    step1Dto3Deq.ScatA=cfree.ScatA
    step1Dto3Deq.MeanPRad=np.zeros_like(MeanPRad3D)
    for i in range(0,200):
        step1Dto3Deq.MeanPRad[:,i,:]=cfree.MeanPRad
    for i in range(200,1200):
        step1Dto3Deq.MeanPRad[:,i,:]=cstep.MeanPRad
    for i in range(1200,1500):
        step1Dto3Deq.MeanPRad[:,i,:]=cfree.MeanPRad

    return step1Dto3Deq,cfree,cstep

if __name__=="__main__":
#    cpn.setup_figures(plt)
    sza='120'
    VZA=0.0
    caseDIR='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Step_Cloud/'
    #3D RT simulations
    dstep865SZAx=POLCARTdset('dstep865SZAx',caseDIR)
    dstep865SZAx.readPOLCARThdf5('step_cld_b865re10ve02_x15km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e6.hdf5',dpath=caseDIR+'results/stepb865/')
    dstep2p13SZAx=POLCARTdset('dstep865SZAx',caseDIR)
    dstep2p13SZAx.readPOLCARThdf5('step_cld_b2p13re10ve02_x15km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e6.hdf5',dpath=caseDIR+'results/stepb2p13/')
    
    dstep=LES_field(dstep865SZAx.fname[0].astype(str).split('_MSCART',2)[0]+'.nc',dpath=caseDIR)#DYCOM field file for MSCART
    dstep.readLES_field()
    
    #1D RT simulations
    dstep865_1D,thin1D_865,thick1D_865=setup1D_3Deq('865' ,'120',dstep865SZAx.MeanPRad)
    dstep2p131D,thin1D2p13,thick1D2p13=setup1D_3Deq('2p13','120',dstep865SZAx.MeanPRad)
    
    #LUTs
    fdpath='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_retrievals/LUTs/'
    fname='MODIS_LUT_extended_SZA%03d_RAA000.nc'%(180-int(sza))
    LUT=Bispec_LUT(fdpath,fname)
    LUT.readLUT()
        
    mu_ix=LUT.find_mu_ix(VZA)
    LUT_VNIR=LUT.I[0,:,1,:,mu_ix].T
    LUT_SWIR=LUT.I[1,:,1,:,mu_ix].T
    re_lines=np.array([4,6,8,10,15,20,30],dtype=float)
    
    #NJK LUT scatter plot for step cloud
    plt.rc('font', family='serif',weight='bold',size=12)
    fig1,ax1=LUT.plotLUT(VZA,re_lines)
    fig1_ttl='dstep'+dstep865SZAx.fname[0].astype(str).split('b865',1)[1].split('.',1)[0]
#    fig1.suptitle(fig1_ttl,size=12)
    ax1.set_xlim([-0.04,0.9])
    ilP=247
    ax1.plot(dstep865_1D.MeanPRad[61,ilP,0],dstep2p131D.MeanPRad[61,ilP,0],'ro',label='_nolegend_')
    ax1.plot(dstep865SZAx.MeanPRad[61,ilP,0],dstep2p13SZAx.MeanPRad[61,ilP,0],'bo',label='_nolegend_')
    sdP=1199
    ax1.plot(dstep865_1D.MeanPRad[61,sdP,0],dstep2p131D.MeanPRad[61,sdP,0],'ro',label='1D')
    ax1.plot(dstep865SZAx.MeanPRad[61,sdP,0],dstep2p13SZAx.MeanPRad[61,sdP,0],'bo',label='3D')
    ax1.legend(loc='upper left',prop={'size':11})
    fig1.show()

    #1D and 3D reflectance (Step cloud illustration)
    fig2,ax2=plt.subplots(figsize=(8,5))
    fig2_ttl='step_cloud_refs_0p865'
    xcens=(dstep.xgrd[:-1]+dstep.xgrd[1:])
    ax2.plot(xcens,dstep865SZAx.MeanPRad[61,:,0],'b-',label='3D')
    ax2.plot(xcens,dstep865_1D.MeanPRad[61,:,0],'r-',label='1D')
    ax2.plot([xcens[ilP],xcens[ilP]],[dstep865SZAx.MeanPRad[61,:,0].min(),dstep865SZAx.MeanPRad[61,:,0].max()],'k--')
    ax2.plot([xcens[sdP],xcens[sdP]],[dstep865_1D.MeanPRad[61,:,0].min(),dstep865_1D.MeanPRad[61,:,0].max()],'k--')
    ax2.set_ylabel('Reflectance')
    ax2.set_xlabel('x (km)')
    fig2.show()    

    fig5,ax5=plt.subplots(figsize=(8,5))
    fig5_ttl='step_cloud_refs_2p13'
    xcens=(dstep.xgrd[:-1]+dstep.xgrd[1:])
    ax5.plot(xcens,dstep2p13SZAx.MeanPRad[61,:,0],'b-',label='3D')
    ax5.plot(xcens,dstep2p131D.MeanPRad[61,:,0],'r-',label='1D')
    ax5.plot([xcens[ilP],xcens[ilP]],[dstep2p13SZAx.MeanPRad[61,:,0].min(),dstep2p13SZAx.MeanPRad[61,:,0].max()],'k--')
    ax5.plot([xcens[sdP],xcens[sdP]],[dstep2p131D.MeanPRad[61,:,0].min(),dstep2p131D.MeanPRad[61,:,0].max()],'k--')
    ax5.set_ylabel('Reflectance')
    ax5.set_xlabel('x (km)')
    fig5.show()    
    
    #Angular pattern of the polarized reflectance at the illuminating edge
    fig3,ax3=plt.subplots(figsize=(8,5))
    fig3_ttl='dstep'+dstep865SZAx.fname[0].astype(str).split('b865',1)[1].split('.',1)[0]+'_Rp'
#    fig3.suptitle(fig3_ttl,size=12)
    ilP=210
    ax3.plot(dstep865SZAx.VZA,-dstep865SZAx.MeanPRad[:,ilP,1],'b-',label='3D(x=%0.2fkm)'%(xcens[ilP]))
    ax3.plot(dstep865_1D.VZA,-dstep865_1D.MeanPRad[:,ilP,1],'r-',label='1D(x=%0.2fkm)'%(xcens[ilP]))
    ax3.set_xlim([-60,60])
    ax3.set_xlabel('VZA')
    ax3.set_ylabel('-Polarized Reflectance')
    ax3t=ax3.twiny()
    ax3tTicks=ax3.get_xticks()
    def tick_fun(VZA):
        V=[np.unique(dstep865SZAx.ScatA[np.where(np.isclose(i,dstep865SZAx.VZA))]) for i in VZA ]
        return ["%d"%z for z in V]
    ax3t.set_xticks(ax3tTicks)
    ax3t.set_xbound(ax3.get_xbound())
    ax3t.set_xticklabels(tick_fun(ax3tTicks))
    ax3.legend(loc='best',prop={'size':11})
    fig3.tight_layout(rect=[0,0,1,0.98])
    fig3.show()

    #Angular pattern of the polarized reflectance at the shadowing edge
    fig4,ax4=plt.subplots(figsize=(8,5))
    fig4_ttl='dstep'+dstep865SZAx.fname[0].astype(str).split('b865',1)[1].split('.',1)[0]+'_Rp_Shad'
#    fig4.suptitle(fig4_ttl,size=12)
    ax4.plot(dstep865SZAx.VZA,-dstep865SZAx.MeanPRad[:,sdP,1],'b-',label='3D(x=%0.2fkm)'%(xcens[sdP]))
    ax4.plot(dstep865_1D.VZA,-dstep865_1D.MeanPRad[:,sdP,1],'r-',label='1D(x=%0.2fkm)'%(xcens[sdP]))
    ax4.set_xlim([-60,60])
    ax4.set_xlabel('VZA')
    ax4.set_ylabel('-Polarized Reflectance')
    ax4t=ax4.twiny()
    ax4tTicks=ax4.get_xticks()
    def tick_fun(VZA):
        V=[np.unique(dstep865SZAx.ScatA[np.where(np.isclose(i,dstep865SZAx.VZA))]) for i in VZA ]
        return ["%d"%z for z in V]
    ax4t.set_xticks(ax4tTicks)
    ax4t.set_xbound(ax4.get_xbound())
    ax4t.set_xticklabels(tick_fun(ax4tTicks))
    ax4.legend(loc='best',prop={'size':11})
    fig4.tight_layout(rect=[0,0,1,0.98])
    fig4.show()

    '''
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
    '''