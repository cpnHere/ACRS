#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Dec 14 09:28:24 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Hariyata wade welada kiyala balanta
"""

import numpy as np
import matplotlib.pyplot as plt
from cpnLES_MSCARTlib import POLCARTdset
from cpnMicrophylib import find_bows
import pickle,sys

class pkl_classes(object):
    def __init__(self,):
        self.class_names=[]

def contourf_Iqu(fig,ax,iqu,pb_VZA_ang,sb_VZA_ang,obj_VZA,obj_MeanPRad,v=None,ticks=None,cmap=None):
    #pb_VZA_ang: primarybow VZA
    #xb_VZA_ang: supernumerary VZA
    if iqu==0:
        sgn=1;iq_lab='I'
    elif iqu==1:
        sgn=-1;iq_lab='-Q'
    fig_ttl=iq_lab+'_contours'
    if v==None:
        v = np.linspace(-0.05,0.15,num=256)
    else:
        v=v
    ctf=ax.contourf(xcens,obj_VZA,sgn*obj_MeanPRad[:,:,iqu],v,extend='both',cmap=cmap)
    if ticks==None:
        fig.colorbar(ctf,ax=ax,orientation='vertical',ticks=np.arange(-0.05,0.16,0.05))
    else:
        fig.colorbar(ctf,ax=ax,orientation='vertical',ticks=ticks)
    if np.size(pb_VZA_ang)>1:
        for i in range(0,np.size(pb_VZA_ang)):
            ax.plot(xcens,np.zeros_like(xcens)+pb_VZA_ang[i],'k--')
    else:
        ax.plot(xcens,np.zeros_like(xcens)+pb_VZA_ang,'k--')
    if np.size(sb_VZA_ang)>1:
        for i in range(0,np.size(sb_VZA_ang)):
            ax.plot(xcens,np.zeros_like(xcens)+sb_VZA_ang[i],'k--')
    else:
        ax.plot(xcens,np.zeros_like(xcens)+sb_VZA_ang,'k--')
    ax.set_ylabel('VZA',size=10)
    return fig_ttl


band=str(sys.argv[1])
sza=str(sys.argv[2])

fracSZAx=POLCARTdset('fracSZAx','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/')
#fracSZAx.readMSCARTplus('fractal_cld_b'+band+'re12ve05_x40km_H0p5_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.nc',\
#                        fdpath='results/fracb'+band+'/',step=True)
fracSZAx.readMSCARTplus('fractal_cld_b'+band+'re12ve05_x40km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.nc',\
                        fdpath='results/fracb'+band+'/',step=True)

#with open('fractal_cld_b'+band+'re12ve05_x40km_H0p5.pkl', "rb") as f:
with open('fractal_cld_b'+band+'re12ve05_x40km.pkl', "rb") as f:
    pkl=pickle.load(f)
    MSCART_field=pickle.load(f)
    fcld=pickle.load(f)


#------ spatial and temperal distribution

xcens=fcld.x
if band=='865':
    pb_ix=np.where(np.isclose(fracSZAx.ScatA,142,atol=0.1))
    sb_ix=np.where(np.isclose(fracSZAx.ScatA,153,atol=0.1))
elif band=='2p13':
    pb_ix=np.where(np.isclose(fracSZAx.ScatA,141,atol=0.1))
    sb_ix=np.where(np.isclose(fracSZAx.ScatA,160,atol=0.1))
elif band=='3p75':
    pb_ix=np.where(np.isclose(fracSZAx.ScatA,158,atol=0.1))
    sb_ix=np.where(np.isclose(fracSZAx.ScatA,170,atol=0.1))
else:    
    pb,sb=find_bows(MSCART_field.scatAng,MSCART_field.PMatVal[0,:,0])
    pb_ix=np.where(np.isclose(fracSZAx.ScatA,pb,atol=0.1))
    sb_ix=np.where(np.isclose(fracSZAx.ScatA,sb,atol=0.1))

fig2,ax2=plt.subplots(3,1,figsize=(7,10),sharex=True)
fig2_ttl=contourf_Iqu(fig2,ax2[0],1,fracSZAx.VZA[pb_ix],fracSZAx.VZA[sb_ix],fracSZAx.VZA,fracSZAx.MeanPRad)
fig2.show()
