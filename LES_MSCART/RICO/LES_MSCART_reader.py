#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Fri May 11 06:58:38 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
"""
import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap

from cpnLES_MSCARTlib import LES_case, POLCARTdset, LES_field
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cpnCommonlib as cpn

def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'figures/')
    
def iqu_cb(fig,ctf,ax,ticks=None,orientation='horizontal',label='label',pad=0.2):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=pad)
    fig.colorbar(ctf, cax=cax,ticks=ticks,orientation=orientation,label=label)

def iqu_DYCOMS2(case,VZA=0,SZA=140,vci=None,vc=None,RTdim='3D',case_name='DYCOMS2'): 
    '''
    vc: dictionary to define v and colorbars for contoursf
        vc={'vIR':np.linspace(0   ,1,50) ,'vQR':np.linspace(-.05,0,50)    ,'vUR':np.linspace(0.0,0.1,50) ,\
           'vIe':np.linspace(0.0,2.0,20),'vQe':np.linspace(0.0,1.00,20)  ,'vUe':np.linspace(0.0,1.0,20),\
           'cIR':np.arange(0,1.1,0.25)  ,'cQR':np.arange(-.05,0.011,0.02),'cUR':np.arange(0.0,0.11,0.05),\
           'cIe':np.arange(0.0,2.1,0.5) ,'cQe':np.arange(0.0,1.1,0.5)    ,'cUe':np.arange(0.0,1.1,0.5)},
    vci: Give index (ix) as the following table
        Contourf color bars and v (VC)
        +++++++++++++++++++++++++++++++++++++++
        DYCOMS2
        +++++++++++++++++++++++++++++++++++++++
        ---------------------------------------
        band   |  SZA  |  VZA  |  ix  |  SAA  |
        ---------------------------------------
        0p860  |  140  |  000  |   1  |  030  |
        2p13   |  140  |  000  |   2  |  030  |
        0p860  |  120  |  000  |   3  |  030  |
        2p13   |  120  |  000  |   4  |  030  |
        ---------------------------------------
        2p13   |  140  |  000  |   5  |  000  |
        0p860  |  140  |  000  |   6  |  000  |
        +++++++++++++++++++++++++++++++++++++++
        ATEXc
        +++++++++++++++++++++++++++++++++++++++
        ---------------------------------------
        0p860  |  120  |  000  |   7  |  000  |
        2p13   |  120  |  000  |   8  |  000  |
    RTdim: '1D' or '3D' RT transfer string
    '''
    VC={1:{'vIR':np.linspace(0   ,1,50) ,'vQR':np.linspace(-.05,0,50)    ,'vUR':np.linspace(0.0,0.1,50) ,\
           'vIe':np.linspace(0.0,2.0,20),'vQe':np.linspace(0.0,1.00,20)  ,'vUe':np.linspace(0.0,1.0,20),\
           'cIR':np.arange(0,1.1,0.25)  ,'cQR':np.arange(-.05,0.011,0.02),'cUR':np.arange(0.0,0.11,0.05),\
           'cIe':np.arange(0.0,2.1,0.5) ,'cQe':np.arange(0.0,1.1,0.5)    ,'cUe':np.arange(0.0,1.1,0.5)},
        2:{'vIR':np.linspace(0.1 ,0.3,50) ,'vQR':np.linspace(-.02,0,50)    ,'vUR':np.linspace(0.02,0.03,50) ,\
           'vIe':np.linspace(0.75,2.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.1,0.31,0.1)  ,'cQR':np.arange(-.02,0.001,0.01),'cUR':np.arange(0.02,0.031,0.005),\
           'cIe':np.arange(0.75,2.1,0.25) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,2.1,0.5)},\
        3:{'vIR':np.linspace(0.2 ,0.8,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.2,0.81,0.2)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        4:{'vIR':np.linspace(0.1 ,0.3,50) ,'vQR':np.linspace(-.005,0,50)   ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.1,0.31,0.1)  ,'cQR':np.arange(-.005,.001,.002),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,2.1,0.5)},\
        5:{'vIR':np.linspace(0.2 ,0.4,50) ,'vQR':np.linspace(-.05,0,50) ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)  ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.2,0.41,0.1)  ,'cQR':np.arange(-.05,.01,.02),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)   ,'cUe':np.arange(0,2.1,0.5)},\
        6:{'vIR':np.linspace(0.2 ,0.82,50) ,'vQR':np.linspace(-.1,0,50) ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)  ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.2,0.82,0.2)  ,'cQR':np.arange(-.1,.01,.05),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)   ,'cUe':np.arange(0,2.1,0.5)},\
        7:{'vIR':np.linspace(0.2 ,0.8,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.2,0.81,0.2)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.0,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        8:{'vIR':np.linspace(0.0 ,0.4,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,0.41,0.1)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)}    
        }        
    if VZA==0:
        VZAi=61
    if vc==None and vci==None:
        print('Give either vc or vci')
    elif not(vci==None):
        vc=VC[vci]
    else:
        print('vc given explicitly')
    fig1,ax1=plt.subplots(2,3,figsize=(8,6),subplot_kw={'aspect':'equal'})

    if RTdim=='3D':
        ctfI,ctfQ,ctfU,cteI,cteQ,cteU,fig1_ttl=iqu_LEScase_3D(fig1,case,ax1,VZAi,vc,case_name=case_name)
    elif RTdim=='1D':
        ctfI,ctfQ,ctfU,cteI,cteQ,cteU,fig1_ttl=iqu_LEScase_1D(fig1,case,ax1,VZAi,vc,case_name=case_name)
    ax1[0,1].tick_params(labelleft=False)
    ax1[0,2].tick_params(labelleft=False)
    ax1[1,1].tick_params(labelleft=False)
    ax1[1,2].tick_params(labelleft=False)
    
    iqu_cb(fig1,ctfI,ax1[0,0],ticks=vc['cIR'],label='R_I')
    iqu_cb(fig1,ctfQ,ax1[0,1],ticks=vc['cQR'],label='R_Q')
    iqu_cb(fig1,ctfU,ax1[0,2],ticks=vc['cUR'],label='R_U')
    iqu_cb(fig1,cteI,ax1[1,0],ticks=vc['cIe'],label='RMS%')
    iqu_cb(fig1,cteQ,ax1[1,1],ticks=vc['cQe'],label='RMS%')
    iqu_cb(fig1,cteU,ax1[1,2],ticks=vc['cUe'],label='RMS%')
    fig1.suptitle("\n".join(wrap(fig1_ttl,60)),size=10)   
    return fig1,fig1_ttl,vc
    
def iqu_LEScase_3D(fig1,case,ax1,VZAi,vc,case_name="DYCOMS2"):
    fig1_ttl=case.RT.fname.split('.',1)[0]+'_'+case_name+'_IQU_top_RMSE_bot'    
    ctfI=ax1[0,0].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,0],vc['vIR'],cmap=plt.cm.jet,extend='both')
    ctfQ=ax1[0,1].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,1],vc['vQR'],cmap=plt.cm.jet,extend='both')
    ctfU=ax1[0,2].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,2],vc['vUR'],cmap=plt.cm.jet,extend='both')
    cteI=ax1[1,0].contourf(case.xcens,case.ycens,case.RT.RMSEPRad[VZAi,:,:,0]/case.RT.MeanPRad[VZAi,:,:,0]*100,vc['vIe'],cmap=plt.cm.jet,extend='both')
    cteQ=ax1[1,1].contourf(case.xcens,case.ycens,case.RT.RMSEPRad[VZAi,:,:,1]/case.RT.MeanPRad[VZAi,:,:,1]*100,vc['vQe'],cmap=plt.cm.jet,extend='both')
    cteU=ax1[1,2].contourf(case.xcens,case.ycens,case.RT.RMSEPRad[VZAi,:,:,2]/case.RT.MeanPRad[VZAi,:,:,2]*100,vc['vUe'],cmap=plt.cm.jet,extend='both')
    return ctfI,ctfQ,ctfU,cteI,cteQ,cteU,fig1_ttl
def iqu_LEScase_1D(fig1,case,ax1,VZAi,vc,case_name="DYCOMS2"):    
    fig1_ttl=case.RT1D.fname.split('.',1)[0]+'_'+case_name+'_IQU_top_RMSE_bot'
    ctfI=ax1[0,0].contourf(case.xcens,case.ycens,case.RT1D.MeanPRad[VZAi,:,:,0].T,vc['vIR'],cmap=plt.cm.jet,extend='both')
    ctfQ=ax1[0,1].contourf(case.xcens,case.ycens,case.RT1D.MeanPRad[VZAi,:,:,1].T,vc['vQR'],cmap=plt.cm.jet,extend='both')
    ctfU=ax1[0,2].contourf(case.xcens,case.ycens,case.RT1D.MeanPRad[VZAi,:,:,2].T,vc['vUR'],cmap=plt.cm.jet,extend='both')
    cteI=ax1[1,0].contourf(case.xcens,case.ycens,(case.RT1D.RMSEPRad[VZAi,:,:,0]/case.RT1D.MeanPRad[VZAi,:,:,0]*100).T,vc['vIe'],cmap=plt.cm.jet,extend='both')
    cteQ=ax1[1,1].contourf(case.xcens,case.ycens,(case.RT1D.RMSEPRad[VZAi,:,:,1]/case.RT1D.MeanPRad[VZAi,:,:,1]*100).T,vc['vQe'],cmap=plt.cm.jet,extend='both')
    cteU=ax1[1,2].contourf(case.xcens,case.ycens,(case.RT1D.RMSEPRad[VZAi,:,:,2]/case.RT1D.MeanPRad[VZAi,:,:,2]*100).T,vc['vUe'],cmap=plt.cm.jet,extend='both')
    return ctfI,ctfQ,ctfU,cteI,cteQ,cteU,fig1_ttl
def DYCOMS2_3D_1D_bias(case,VZA,SZA,vci):
    '''
    DYCOMS2 3D and 1D biases plot
    Uses iqu_DYCOMS2()
    '''
    if VZA==0:
        VZAi=61
    _,_,vc=iqu_DYCOMS2(case,VZA,SZA,vci=vci)
    RTbias=POLCARTdset('LES','')
    RTbias.MeanPRad=case.RT.MeanPRad-case.RT1D.MeanPRad.swapaxes(1,2)

    fig7,ax7=plt.subplots(3,2,figsize=(5,7),subplot_kw={'aspect':'equal'})
    fig7_ttl=case.RT.fname.split('.',1)[0]+'_DYCOM2_IQU_3D_1D_biases'
    fig7.suptitle("\n".join(wrap(fig7_ttl,50)),size=10)
    
    ctfI=ax7[0,0].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,0]  ,vc['vIR'],cmap=plt.cm.jet,extend='both');ax7[0,0].set_title('I_3D',size=10)
    ctfQ=ax7[0,1].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,1]  ,vc['vQR'],cmap=plt.cm.jet,extend='both');ax7[0,1].set_title('Q_3D',size=10)
    ctI1=ax7[1,0].contourf(case.xcens,case.ycens,(case.RT1D.MeanPRad[VZAi,:,:,0]).T,vc['vIR'],cmap=plt.cm.jet,extend='both');ax7[1,0].set_title('I_1D',size=10)
    ctQ1=ax7[1,1].contourf(case.xcens,case.ycens,(case.RT1D.MeanPRad[VZAi,:,:,1]).T,vc['vQR'],cmap=plt.cm.jet,extend='both');ax7[1,1].set_title('Q_1D',size=10)
    
    ctfb1=ax7[2,0].contourf(case.xcens,case.ycens,RTbias.MeanPRad[VZAi,:,:,0],np.linspace(-.1,.1,50)  ,cmap=plt.cm.RdBu_r,extend='both');ax7[2,0].set_title('I Bias',size=10)
    ctfb2=ax7[2,1].contourf(case.xcens,case.ycens,RTbias.MeanPRad[VZAi,:,:,1],np.linspace(-.01,.01,50),cmap=plt.cm.RdBu_r,extend='both');ax7[2,1].set_title('Q Bias',size=10)
    
    ax7[0,0].set_ylabel('km',size=10)
    ax7[1,0].set_ylabel('km',size=10)
    ax7[2,0].set_ylabel('km',size=10)
    
    ax7[0,1].tick_params(labelleft=False)
    ax7[1,1].tick_params(labelleft=False)
    ax7[2,1].tick_params(labelleft=False)
    
    iqu_cb(fig7,ctfI,ax7[0,0],ticks=vc['cIR'],label='')
    iqu_cb(fig7,ctfQ,ax7[0,1],ticks=vc['cQR'],label='')
    iqu_cb(fig7,ctI1,ax7[1,0],ticks=vc['cIR'],label='')
    iqu_cb(fig7,ctQ1,ax7[1,1],ticks=vc['cQR'],label='')
    iqu_cb(fig7,ctfb1,ax7[2,0],ticks=np.arange(-.1,.11,.1) ,label='')
    iqu_cb(fig7,ctfb2,ax7[2,1],ticks=np.arange(-.01,.011,.01),label='')
    cpn.sub_labels(ax7)
    return fig7,fig7_ttl

def old_vs_new_SZA140(dynyx,dyo):    
    #Comparing old and new runs
    fig1,ax1=plt.subplots(3,3,figsize=(8,8),subplot_kw={'aspect':'equal'})
    fig1_ttl=dynyx.RT.fname.split('.',1)[0]+'_DYCOM2_IQU_old_vs_new_3D'
    fig1.suptitle("\n".join(wrap(fig1_ttl,60)))
    
    ctfI=ax1[0,0].contourf(dynyx.xcens,dynyx.ycens,dynyx.RT.MeanPRad[61,:,:,0],np.linspace(0   ,1,50) ,cmap=plt.cm.jet)
    ctfQ=ax1[0,1].contourf(dynyx.xcens,dynyx.ycens,dynyx.RT.MeanPRad[61,:,:,1],np.linspace(-.05,0,50) ,cmap=plt.cm.jet)
    ctfU=ax1[0,2].contourf(dynyx.xcens,dynyx.ycens,dynyx.RT.MeanPRad[61,:,:,2],np.linspace(0.0,0.1,50),cmap=plt.cm.jet)
    ax1[1,0].contourf(dyo.xcens,dynyx.ycens  ,dyo.RT.MeanPRad[61,:,:,0],np.linspace(0   ,1,50) ,cmap=plt.cm.jet)
    ax1[1,1].contourf(dyo.xcens,dynyx.ycens  ,dyo.RT.MeanPRad[61,:,:,1]  ,np.linspace(-.05,0,50) ,cmap=plt.cm.jet)
    ax1[1,2].contourf(dyo.xcens,dynyx.ycens  ,dyo.RT.MeanPRad[61,:,:,2]  ,np.linspace(0.0,0.1,50),cmap=plt.cm.jet)
    
    b1,b2,b3=dynyx.RT.MeanPRad[61,:,:,0]-dyo.RT.MeanPRad[61,:,:,0],\
             dynyx.RT.MeanPRad[61,:,:,1]-dyo.RT.MeanPRad[61,:,:,1],\
             dynyx.RT.MeanPRad[61,:,:,2]-dyo.RT.MeanPRad[61,:,:,2]
    ctfb1=ax1[2,0].contourf(dyo.xcens,dynyx.ycens,b1,np.linspace(-.003,.003,50)  ,cmap=plt.cm.RdBu_r,extend='both')
    ctfb2=ax1[2,1].contourf(dyo.xcens,dynyx.ycens,b2,np.linspace(-.002,.002,50),cmap=plt.cm.RdBu_r,extend='both')
    ctfb3=ax1[2,2].contourf(dyo.xcens,dynyx.ycens,b3,np.linspace(-.002,.002,50),cmap=plt.cm.RdBu_r,extend='both')
    
    ax1[0,1].set_title('New',size=10)
    ax1[1,1].set_title('Old',size=10)
    ax1[2,1].set_title('Bias',size=10)
    
    ax1[0,0].set_ylabel('km',size=10)
    ax1[1,0].set_ylabel('km',size=10)
    ax1[2,0].set_ylabel('km',size=10)
    
    ax1[0,1].tick_params(labelleft=False)
    ax1[0,2].tick_params(labelleft=False)
    iqu_cb(fig1,ctfI,ax1[0,0],ticks=np.arange(0,1.1,0.25)     ,label='')
    iqu_cb(fig1,ctfQ,ax1[0,1],ticks=np.arange(-.05,0.011,0.02),label='')
    iqu_cb(fig1,ctfU,ax1[0,2],ticks=np.arange(0.0,0.11,0.05)  ,label='')
    iqu_cb(fig1,ctfI,ax1[1,0],ticks=np.arange(0,1.1,0.25)     ,label='')
    iqu_cb(fig1,ctfQ,ax1[1,1],ticks=np.arange(-.05,0.011,0.02),label='')
    iqu_cb(fig1,ctfU,ax1[1,2],ticks=np.arange(0.0,0.11,0.05)  ,label='')
    iqu_cb(fig1,ctfb1,ax1[2,0],np.arange(-.003,.0031,0.003),label='R_I')
    iqu_cb(fig1,ctfb2,ax1[2,1],np.arange(-.002,.0021,0.002),label='R_Q')
    iqu_cb(fig1,ctfb3,ax1[2,2],np.arange(-.002,.0021,0.002),label='R_U')
    
    fig1.show() 
    return fig1,fig1_ttl
def get_binTicks(data,dc):
    bmx=np.round(data.max(),dc)
    bmn=np.round(data.min(),dc)
    return bmn,bmx


if __name__=='__main__':
    cpn.setup_figures(plt)
    #ATEXp full cases
    '''
    base_dir='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/'    

    ATEXp0p860_120=LES_case('ATEXp_dharma_013067_b0p860_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                          RT1Dname='ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp2p13_120 =LES_case('ATEXp_dharma_013067_b2p13_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                            RT1Dname='ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp0p860_140=LES_case('ATEXp_dharma_013067_b0p860_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                          RT1Dname='ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp2p13_140 =LES_case('ATEXp_dharma_013067_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                            RT1Dname='ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp0p860_160=LES_case('ATEXp_dharma_013067_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                          RT1Dname='ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp2p13_160 =LES_case('ATEXp_dharma_013067_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                            RT1Dname='ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5')

    fig1,fig1_ttl,_=iqu_DYCOMS2(ATEXp0p860_120,0,120,vci=7,case_name='ATEXp');fig1.tight_layout(rect=[0,0,1,0.96]);fig1.show()
    fig2,fig2_ttl,_=iqu_DYCOMS2(ATEXp2p13_120 ,0,120,vci=8,case_name='ATEXp');fig2.tight_layout(rect=[0,0,1,0.96]);fig2.show()
    fig3,fig3_ttl,_=iqu_DYCOMS2(ATEXp0p860_120,0,120,vci=7,case_name='ATEXp',RTdim='1D');fig3.tight_layout(rect=[0,0,1,0.96]);fig3.show()
    fig4,fig4_ttl,_=iqu_DYCOMS2(ATEXp2p13_120 ,0,120,vci=8,case_name='ATEXp',RTdim='1D');fig4.tight_layout(rect=[0,0,1,0.96]);fig4.show()    
    fig5,fig5_ttl  =DYCOMS2_3D_1D_bias(ATEXp0p860_120,0,120,7);fig5.show()
    fig6,fig6_ttl  =DYCOMS2_3D_1D_bias(ATEXp2p13_120,0,120,8);fig6.show()
    
    fig7,fig7_ttl,_=iqu_DYCOMS2(ATEXp0p860_140,0,140,vci=7,case_name='ATEXp');fig7.tight_layout(rect=[0,0,1,0.96]);fig7.show()
    fig8,fig8_ttl,_=iqu_DYCOMS2(ATEXp2p13_140 ,0,140,vci=8,case_name='ATEXp');fig8.show()
    fig9,fig9_ttl,_=iqu_DYCOMS2(ATEXp0p860_140,0,140,vci=7,case_name='ATEXp',RTdim='1D');fig9.show()
    fig10,fig10_ttl,_=iqu_DYCOMS2(ATEXp2p13_140 ,0,140,vci=8,case_name='ATEXp',RTdim='1D');fig10.show()
    fig11,fig11_ttl  =DYCOMS2_3D_1D_bias(ATEXp0p860_140,0,140,7);fig11.show()
    fig12,fig12_ttl  =DYCOMS2_3D_1D_bias(ATEXp2p13_140,0,140,8);fig12.show()

    fig13,fig13_ttl,_=iqu_DYCOMS2(ATEXp0p860_160,0,160,vci=7,case_name='ATEXp');fig13.tight_layout(rect=[0,0,1,0.96]);fig13.show()
    fig14,fig14_ttl,_=iqu_DYCOMS2(ATEXp2p13_160 ,0,160,vci=8,case_name='ATEXp');fig14.show()
    fig15,fig15_ttl,_=iqu_DYCOMS2(ATEXp0p860_160,0,160,vci=7,case_name='ATEXp',RTdim='1D');fig15.show()
    fig16,fig16_ttl,_=iqu_DYCOMS2(ATEXp2p13_160 ,0,160,vci=8,case_name='ATEXp',RTdim='1D');fig16.show()
    fig17,fig17_ttl  =DYCOMS2_3D_1D_bias(ATEXp0p860_160,0,160,7);fig17.show()
    fig18,fig18_ttl  =DYCOMS2_3D_1D_bias(ATEXp2p13_160,0,160,8);fig18.show()


    '''
    cpn.setup_figures(plt)
    base_dir='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/'
    band='0p860'
    MSCARThdf='RICO_dharma_005044_b'+band+'_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5'
    RT=POLCARTdset('LES',base_dir)
    RT.readPOLCARThdf5(MSCARThdf,dpath=base_dir+'results/b'+band+'/')
    RT_field=LES_field(RT.fname.split('_MSCART',1)[0]+'.nc',dpath=base_dir)
    RT_field.readLES_field()
    fig1,ax1=plt.subplots(2,3,figsize=(8,6),subplot_kw={'aspect':'equal'})
    #0p860 120
#    vc={'vIR':np.linspace(0.0 ,0.8,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
#     'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
#     'cIR':np.arange(0.2,0.81,0.2)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(0.0,0.021,0.01),\
#     'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)}    
    #0p860 140
    vc={'vIR':np.linspace(0.0 ,0.8,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
     'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
     'cIR':np.arange(0.2,0.81,0.2)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(0.0,0.021,0.01),\
     'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)}    

    #2p13 160
#    vc={'vIR':np.linspace(0.0 ,0.3,50) ,'vQR':np.linspace(-.005,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
#     'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
#     'cIR':np.arange(0.0,0.31,0.1)  ,'cQR':np.arange(-.005,.001,0.005),'cUR':np.arange(-0.01,0.011,0.01),\
#     'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)}    

    #2p13 140/120
#    vc={'vIR':np.linspace(0.0 ,0.3,50) ,'vQR':np.linspace(-.005,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
#     'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
#     'cIR':np.arange(0.0,0.31,0.1)  ,'cQR':np.arange(-.005,.001,0.005),'cUR':np.arange(-0.01,0.011,0.01),\
#     'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)}    


    xcens=(RT_field.xgrd[1:]+RT_field.xgrd[0:-1])/2
    ycens=(RT_field.ygrd[1:]+RT_field.ygrd[0:-1])/2       
    VZAi=61
    fig1_ttl=RT.fname.split('.',1)[0]+'_RICO_IQU_top_RMSE_bot'

    ctfI=ax1[0,0].contourf(xcens,ycens,RT.MeanPRad[VZAi,:,:,0],vc['vIR'],cmap=plt.cm.jet,extend='both')
    ctfQ=ax1[0,1].contourf(xcens,ycens,RT.MeanPRad[VZAi,:,:,1],vc['vQR'],cmap=plt.cm.jet,extend='both')
    ctfU=ax1[0,2].contourf(xcens,ycens,RT.MeanPRad[VZAi,:,:,2],vc['vUR'],cmap=plt.cm.jet,extend='both')
    cteI=ax1[1,0].contourf(xcens,ycens,RT.RMSEPRad[VZAi,:,:,0]/RT.MeanPRad[VZAi,:,:,0]*100,vc['vIe'],cmap=plt.cm.jet,extend='both')
    cteQ=ax1[1,1].contourf(xcens,ycens,RT.RMSEPRad[VZAi,:,:,1]/RT.MeanPRad[VZAi,:,:,1]*100,vc['vQe'],cmap=plt.cm.jet,extend='both')
    cteU=ax1[1,2].contourf(xcens,ycens,RT.RMSEPRad[VZAi,:,:,2]/RT.MeanPRad[VZAi,:,:,2]*100,vc['vUe'],cmap=plt.cm.jet,extend='both')

    ax1[0,1].tick_params(labelleft=False)
    ax1[0,2].tick_params(labelleft=False)
    ax1[1,1].tick_params(labelleft=False)
    ax1[1,2].tick_params(labelleft=False)
    
    iqu_cb(fig1,ctfI,ax1[0,0],ticks=vc['cIR'],label='R_I')
    iqu_cb(fig1,ctfQ,ax1[0,1],ticks=vc['cQR'],label='R_Q')
    iqu_cb(fig1,ctfU,ax1[0,2],ticks=vc['cUR'],label='R_U')
    iqu_cb(fig1,cteI,ax1[1,0],ticks=vc['cIe'],label='RMS%')
    iqu_cb(fig1,cteQ,ax1[1,1],ticks=vc['cQe'],label='RMS%')
    iqu_cb(fig1,cteU,ax1[1,2],ticks=vc['cUe'],label='RMS%')
    fig1.suptitle("\n".join(wrap(fig1_ttl,60)),size=10)   
    fig1.tight_layout(rect=[0,0,1,0.98])
    fig1.show()
    


