#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Apr 30 13:09:45 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************

LES simulations from DHARMA

"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import cpnCommonlib as cpn
import h5py
import time
from cpnLES_MSCARTlib import DHARMA

  #dharma_005263.cdf  dharma_007877.cdf  dharma_0.cdf  dharma_012540.cdf

def findLWP(ATEX):
    '''
    Find LWP by integrating alog the each bin
    Approx. time 130 sec.
    ATEX: DHARMA object
    return:
        lwp array (ATEX.x.size,ATEX.y.size)
    '''
    r=ATEX.r_drops#um
    z=ATEX.zw[1:]-ATEX.zw[0:-1]#m
    lwp=np.zeros((ATEX.x.size,ATEX.y.size),dtype=float)
    st=time.time()
    for j in np.arange(0,lwp.shape[0]):
        for k in np.arange(0,lwp.shape[1]):
#            b='%d'%(j/lwp.shape[0]*100)
#            print(b,end='\r')
            t1=np.zeros_like(z,dtype=float)
            for i in np.arange(0,z.size):
                nr=ATEX.dN_drops[:,i,j,k]
                t1[i]=np.sum(4.0/3.0*np.pi*r**3*nr) #um^3/cm^3
                t1[i]=t1[i]*1e-12 #unitless
            rhol=1e6 #g/m^3
            lwp[j,k]=np.sum(t1*rhol*z)#g/m^2
    et=time.time()
    print('%d sec elapsed'%(et-st))
    return lwp

def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,path='figures/')
    
def addcase(f,obj):
    PCentry=f.create_dataset('lwp50m',data=obj.lwp50m)
    PCentry.dims[0].label='x'
    PCentry.dims[1].label='y'
    PCentry.attrs['units']='g/m^2'
    PCentry.attrs["long_name"]='Liquid_water_path_50m_resolution'
    PCentry.attrs["DHARAMA_output_name"]=obj.fname

    PCentry=f.create_dataset('x',data=obj.x)
    PCentry.attrs['units']='m'
    PCentry.attrs['long_name']='x_dimension'
    
    PCentry=f.create_dataset('y',data=obj.y)
    PCentry.attrs['units']='m'
    PCentry.attrs['long_name']='y_dimension'

def calNsave_lwp50m(caseDict,des):
    for i in np.arange(0,np.size(caseDict.keys())):
        caseDict[i].lwp50m=findLWP(caseDict[i])
        fname=des+'_'+caseDict[i].fname.split('.',1)[0]+'_lwp50m'
        f=h5py.File(fname+'.h5','w')
        f.attrs['DHARMA_case_description']=des
        addcase(f,caseDict[i])
        f.close()
def read_lwp50(caseDict,des):
    for i in np.arange(0,np.size(caseDict.keys())):
        f=h5py.File(des+'_'+caseDict[i].fname.split('.',1)[0]+'_lwp50m.h5','r')
        caseDict[i].lwp50m=f['lwp50m'][:]

def add_common_cb(fig,ctf,ts=None,label=None):
    '''
    To add a common color bar for a subplot
    fig: Figure object matplotlib.figure.Figure
    ctf:contourf output instance
    '''
    fig.subplots_adjust(right=0.81)
    cbar_ax = fig.add_axes([0.84, 0.15, 0.01, 0.7])
    if not(ts==None) and not(label==None):
        fig.colorbar(ctf, cax=cbar_ax, ticks=ts,label=label)
    else:
        fig.colorbar(ctf,cax=cbar_ax)

def NineDPlots(fig,ax,caseDict,ttl):
    '''
    ax is 3by3
    '''
    cmap=plt.cm.jet
    fig_ttl=ttl
    fig.suptitle(fig_ttl)
    v =np.linspace(0,300,50)
    ts=np.arange  (0,301,100)
    axs=ax.flat
    i=0
    for n, axn in enumerate(axs):
        ctf=axn.contourf(caseDict[i].x/1e3,caseDict[i].y/1e3,caseDict[i].lwp50m,v,cmap=cmap);
        axn.set_title('%d'%caseDict[i].time+':'+caseDict[i].fname,size=10)
        i+=1
    add_common_cb(fig,ctf,ts,r'LWP ($g/m^2$)')
    ax[0,0].set_ylabel('y(km)')
    ax[1,0].set_ylabel('y(km)')
    ax[2,0].set_ylabel('y(km)')  
    ax[2,0].set_xlabel('x(km)')
    ax[2,1].set_xlabel('x(km)')
    ax[2,2].set_xlabel('x(km)') 
    cpn.sub_labels(ax,'k')
    return fig_ttl

def EightDPlots(fig,ax,caseDict,ttl):
    '''
    ax is 2 by 4
    '''
    cmap=plt.cm.jet
    fig_ttl=ttl
    fig.suptitle(fig_ttl)
    v =np.linspace(0,300,50)
    ts=np.arange  (0,301,100)
    axs=ax.flat
    i=0
    for n, axn in enumerate(axs):
        ctf=axn.contourf(caseDict[i].x/1e3,caseDict[i].y/1e3,caseDict[i].lwp50m,v,cmap=cmap);
        axn.set_title('%d'%caseDict[i].time+':'+caseDict[i].fname,size=10)
        i+=1
    add_common_cb(fig,ctf,ts,r'LWP ($g/m^2$)')
    ax[0,0].set_ylabel('y(km)')
    ax[1,0].set_ylabel('y(km)')  
    ax[1,0].set_xlabel('x(km)')
    ax[1,1].set_xlabel('x(km)')
    ax[1,2].set_xlabel('x(km)')
    ax[1,3].set_xlabel('x(km)')
    cpn.sub_labels(ax,'k')
    return fig_ttl

if __name__=='__main__':
    cpn.setup_figures(plt)
    ATEX040_path  ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/ATEX/Ackerman_ATEX_highres/CCN40/'
    ATEX600_path  ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/ATEX/Ackerman_ATEX_highres/CCN600/'
    ATEXl040_path ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/ATEX/Ackerman_ATEX/CCN40/'
    ATEXl075_path ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/ATEX/Ackerman_ATEX/CCN75/'
    ATEXl600_path ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/ATEX/Ackerman_ATEX/CCN600/'
    DYCOMS2_path  ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/DYCOMS/Ackerman_DYCOMS2_25bins/'
    RICOn_path    ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/RICO/Ackerman_RICO_NEW/'
    ATEX040 ={0:DHARMA('dharma_003122.cdf',ATEX040_path),\
              1:DHARMA('dharma_005263.cdf',ATEX040_path),\
              2:DHARMA('dharma_007877.cdf',ATEX040_path),\
              3:DHARMA('dharma_010376.cdf',ATEX040_path),\
              4:DHARMA('dharma_012540.cdf',ATEX040_path)}
    ATEX600 ={0:DHARMA('dharma_003239.cdf',ATEX600_path),\
              1:DHARMA('dharma_005384.cdf',ATEX600_path),\
              2:DHARMA('dharma_007724.cdf',ATEX600_path),\
              3:DHARMA('dharma_010314.cdf',ATEX600_path),\
              4:DHARMA('dharma_013067.cdf',ATEX600_path)}
    ATEXl040={0:DHARMA('dharma_003314.cdf',ATEXl040_path),\
              1:DHARMA('dharma_003734.cdf',ATEXl040_path),\
              2:DHARMA('dharma_004155.cdf',ATEXl040_path),\
              3:DHARMA('dharma_004576.cdf',ATEXl040_path),\
              4:DHARMA('dharma_004997.cdf',ATEXl040_path),\
              5:DHARMA('dharma_005418.cdf',ATEXl040_path),\
              6:DHARMA('dharma_005839.cdf',ATEXl040_path),\
              7:DHARMA('dharma_006260.cdf',ATEXl040_path),\
              8:DHARMA('dharma_006681.cdf',ATEXl040_path)}
    ATEXl075={0:DHARMA('dharma_003313.cdf',ATEXl075_path),\
              1:DHARMA('dharma_003733.cdf',ATEXl075_path),\
              2:DHARMA('dharma_004154.cdf',ATEXl075_path),\
              3:DHARMA('dharma_004575.cdf',ATEXl075_path),\
              4:DHARMA('dharma_004996.cdf',ATEXl075_path),\
              5:DHARMA('dharma_005417.cdf',ATEXl075_path),\
              6:DHARMA('dharma_005838.cdf',ATEXl075_path),\
              7:DHARMA('dharma_006259.cdf',ATEXl075_path),\
              8:DHARMA('dharma_006680.cdf',ATEXl075_path)}
    ATEXl600={0:DHARMA('dharma_003313.cdf',ATEXl600_path),\
              1:DHARMA('dharma_003733.cdf',ATEXl600_path),\
              2:DHARMA('dharma_004154.cdf',ATEXl600_path),\
              3:DHARMA('dharma_004575.cdf',ATEXl600_path),\
              4:DHARMA('dharma_004996.cdf',ATEXl600_path),\
              5:DHARMA('dharma_005417.cdf',ATEXl600_path),\
              6:DHARMA('dharma_005838.cdf',ATEXl600_path),\
              7:DHARMA('dharma_006259.cdf',ATEXl600_path),\
              8:DHARMA('dharma_006680.cdf',ATEXl600_path)}
    DYCOMS2 ={0:DHARMA('dharma_002304.cdf',DYCOMS2_path),\
              1:DHARMA('dharma_003043.cdf',DYCOMS2_path),\
              2:DHARMA('dharma_003784.cdf',DYCOMS2_path),\
              3:DHARMA('dharma_004637.cdf',DYCOMS2_path),\
              4:DHARMA('dharma_005416.cdf',DYCOMS2_path),\
              5:DHARMA('dharma_006263.cdf',DYCOMS2_path),\
              6:DHARMA('dharma_007138.cdf',DYCOMS2_path),\
              7:DHARMA('dharma_008036.cdf',DYCOMS2_path),\
              8:DHARMA('dharma_008884.cdf',DYCOMS2_path)}
    RICOn   ={0:DHARMA('dharma_001440.cdf',RICOn_path),\
              1:DHARMA('dharma_002160.cdf',RICOn_path),\
              2:DHARMA('dharma_002880.cdf',RICOn_path),\
              3:DHARMA('dharma_003602.cdf',RICOn_path),\
              4:DHARMA('dharma_004324.cdf',RICOn_path),\
              5:DHARMA('dharma_005044.cdf',RICOn_path),\
              6:DHARMA('dharma_005772.cdf',RICOn_path),\
              7:DHARMA('dharma_006501.cdf',RICOn_path)}
    new_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/Ann_Dan/'
    new={0:DHARMA('dharma_001800.cdf',new_path),\
         1:DHARMA('dharma_001800.cdf',new_path)}

#    calNsave_lwp50m(ATEX040 ,'ATEX_higres_CCN040')
#    calNsave_lwp50m(ATEX600 ,'ATEX_higres_CCN600')
#    calNsave_lwp50m(ATEXl040,'ATEX_lowres_CCN040')
#    calNsave_lwp50m(ATEXl075,'ATEX_lowres_CCN075')
#    calNsave_lwp50m(ATEXl600,'ATEX_lowres_CCN600')
#    calNsave_lwp50m(DYCOMS2,'DYCOMS2_25bin')
#    calNsave_lwp50m(RICOn,'RICO')
    
    read_lwp50(ATEX040 ,'ATEX_higres_CCN040')
    read_lwp50(ATEX600 ,'ATEX_higres_CCN600')
    read_lwp50(ATEXl040,'ATEX_lowres_CCN040')    
    read_lwp50(ATEXl075,'ATEX_lowres_CCN075')    
    read_lwp50(ATEXl600,'ATEX_lowres_CCN600')
    read_lwp50(DYCOMS2 ,'DYCOMS2_25bin')
    read_lwp50(RICOn   ,'RICO')
#------------------------------------------------------------------------------
#ATEX high resolution cases
#========================
    fig1,ax1=plt.subplots(1,5,figsize=(16,3),subplot_kw={'aspect':'equal'})
    cmap=plt.cm.jet
    fig1_ttl='ATEX_higres_CCN40_LWP'
    fig1.suptitle(fig1_ttl)
    v =np.linspace(0,300,50)
    ts=np.arange  (0,301,100)
    for i in np.arange(0,np.size(ATEX040.keys())):
        ctf=ax1[i].contourf(ATEX040[i].x/1e3,ATEX040[i].y/1e3,ATEX040[i].lwp50m,v,cmap=cmap);
        ax1[i].set_title('%d'%ATEX040[i].time+':'+ATEX040[i].fname,size=10)
        ax1[i].set_xlabel('x(km)')
        
    add_common_cb(fig1,ctf,ts,r'LWP ($g/m^2$)')
    ax1[0].set_ylabel('y(km)')
    cpn.sub_labels(ax1,'k')    
    fig1.show()

    fig2,ax2=plt.subplots(1,5,figsize=(16,3),subplot_kw={'aspect':'equal'})
    cmap=plt.cm.jet
    fig2_ttl='ATEX_higres_CCN600_LWP'
    fig2.suptitle(fig2_ttl)
    v =np.linspace(0,300,50)
    ts=np.arange  (0,301,100)
    for i in np.arange(0,np.size(ATEX600.keys())):
        ctf=ax2[i].contourf(ATEX600[i].x/1e3,ATEX600[i].y/1e3,ATEX600[i].lwp50m,v,cmap=cmap);
        ax2[i].set_title('%d'%ATEX600[i].time+':'+ATEX600[i].fname,size=10)
        ax2[i].set_xlabel('x(km)')
    add_common_cb(fig2,ctf,ts,r'LWP($g/m^2$)')
    ax2[0].set_ylabel('y(km)')
    cpn.sub_labels(ax2,'k')    
    fig2.show()
    
#--------------------------------------------------------------------------------
#ATEX low resolution cases
#==========================
    fig3,ax3=plt.subplots(3,3,figsize=(8,8),subplot_kw={'aspect':'equal'})
    fig3_ttl=NineDPlots(fig3,ax3,ATEXl040,'ATEX_lowres_CCN040_LWP')
    fig3.show()
    
    fig4,ax4=plt.subplots(3,3,figsize=(8,8),subplot_kw={'aspect':'equal'})
    fig4_ttl=NineDPlots(fig4,ax4,ATEXl075,'ATEX_lowres_CCN075_LWP')
    fig4.show()
    
    fig5,ax5=plt.subplots(3,3,figsize=(8,8),subplot_kw={'aspect':'equal'})
    fig5_ttl=NineDPlots(fig5,ax5,ATEXl600,'ATEX_lowres_CCN600_LWP')
    fig5.show()

#-------------------------------------------------------------------------------
#DYCOMS2
#=======
    fig6,ax6=plt.subplots(3,3,figsize=(8,8),subplot_kw={'aspect':'equal'})
    fig6_ttl=NineDPlots(fig6,ax6,DYCOMS2,'DYCOM2_25bins')
    fig6.show()

#------------------------------------------------------------------------------
#RICO new case
#=============
    fig7,ax7=plt.subplots(2,4,figsize=(12,5),subplot_kw={'aspect':'equal'})
    fig7_ttl=EightDPlots(fig7,ax7,RICOn,'RICO')
    fig7.show()
    
    fig7,ax7=plt.subplots(2,4,figsize=(12,5),subplot_kw={'aspect':'equal'})
    fig7_ttl=EightDPlots(fig7,ax7,new,'new')
    fig7.show()


    