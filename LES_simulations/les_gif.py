#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Tue Apr  2 13:00:55 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To create GIF for LES psudo retrievals.
"""
import numpy as np
import matplotlib.pyplot as plt
import time
from textwrap import wrap
from matplotlib.animation import FuncAnimation
  
from cpnLES_MSCARTlib import DHARMA_onmp
import cpnCommonlib as cpn

def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'GIFs/')
cpn.setup_figures(plt)
if __name__=='__main__':    
    replace='0' #'0' to load from existing files. '1' to replace previous results.
    #DYCOMS2------------------------------------------------------------------------------------------
    DYCOMS2_path  ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/DYCOMS/Ackerman_DYCOMS2_25bins/'
    DYCOMS2=DHARMA_onmp('dharma_008036.cdf',DYCOMS2_path );DYCOMS2.readDHARMA()
    mie_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/'
    mie_name='DYCOM2_dharma_008036_mie_470_860_2p13'
    lesCname='DYCOMS2';lesBinEd=None
    a=1;b=0;SZA=40;VZA=0
    band=1
    wvl={1:'0p860',2:'2p13'}
    DYCOMS2.setup_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band,replace=replace)
    #colorbars
    v={'re':np.linspace(0,25,50),'ve':np.linspace(0,0.2,50) ,'revw':np.linspace(0,25,50),'vevw':np.linspace(0,0.2,50)}
    c={'re':np.arange(0,25.1,5) ,'ve':np.arange(0,0.21,0.1) ,'revw':np.arange(0,25.1,5) ,'vevw':np.arange(0,0.21,0.1)}  
    v['tu']=np.linspace(0,30,50)
    c['tu']=np.arange(0,30.1,10)
    version='V4'

    #DYCOMS2 end ---------------------------------------------------------------------------------------
    '''
    #ATEXc ----------------------------------------------------------------------------------------------
    ATEXc_path  ='/umbc/xfs1/zzbatmos/users/mild1/Research/LES/LES_Satellite_Simulator_DJM/LES_Simulations/Ackerman_ATEX_highres/CCN40/'  
    ATEXc=DHARMA_onmp('dharma_007877.cdf',ATEXc_path );ATEXc.readDHARMA()
    mie_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/'
    mie_name='DYCOM2_dharma_008036_mie_470_860_2p13'
    lesCname='ATEXc'
    a=1;b=0;SZA=40;VZA=0
    band=1
    wvl={1:'0p860',2:'2p13'}
    ATEXc.setup_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)
    #ATEXc end---------------------------------------------------------------------------------------------
    
    #ATEXp -------------------------------------------------------------------------------------------------------------------
    ATEXp_path  ='/umbc/xfs1/zzbatmos/users/mild1/Research/LES/LES_Satellite_Simulator_DJM/LES_Simulations/Ackerman_ATEX_highres/CCN600/'  
    ATEXp=DHARMA_onmp('dharma_013067.cdf',ATEXp_path );ATEXp.readDHARMA()
    mie_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/'
    mie_name='DYCOM2_dharma_008036_mie_470_860_2p13'
    lesCname='ATEXp'
    a=1;b=0;SZA=40;VZA=0;band=1;wvl={1:'0p860',2:'2p13'}
    ATEXp.setup_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)
    #ATEXp end ------------------------------------------------------------------------------------------------------------------------

    #RICO ------------------------------------------------------------------------------------------------------------------------------
    RICO_path='/umbc/xfs1/zzbatmos/users/mild1/Research/LES/LES_Satellite_Simulator_DJM/LES_Simulations/Ackerman_RICO_NEW/'
    RICO=DHARMA_onmp('dharma_005044.cdf',RICO_path );RICO.readDHARMA()
    mie_path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_simulations/'
    mie_name='RICO_dharma_005044_mie_470_860_2p13'
    lesCname='RICO'
    a=1;b=0;SZA=40;VZA=0;band=1;wvl={1:'0p860',2:'2p13'}
    RICO.setup_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)

    #RICO end---------------------------------------------------------------------------------------------------------------------------
    '''
#    v={'re':np.linspace(0,100,50) ,'ve':np.linspace(0,0.5,50) ,'revw':np.linspace(0,25,50),'vevw':np.linspace(0,0.2,50)}
#    c={'re':np.arange(0,250.1,125),'ve':np.arange(0,0.51,0.25),'revw':np.arange(0,25.1,5) ,'vevw':np.arange(0,0.21,0.1)}        
    case=DYCOMS2
    case.setup_dTau(mie_name,lesCname,mie_path=mie_path)
    #re_tau,ve_tau and re_vw,ve_vw with altitude gif
    fig6=plt.figure(figsize=(10,6))
    fig6_ttl=case.VW.fname.split('.',1)[0]+'LES_re_ve_summary'  
    fig6.suptitle("\n".join(wrap(fig6_ttl,80)),size=12)
    ax6p1=plt.subplot2grid((2,10),(0,0),rowspan=2,)
    ax6p2=plt.subplot2grid((2,10),(0,1),colspan=3,aspect='equal')
    ax6p3=plt.subplot2grid((2,10),(1,1),colspan=3,aspect='equal')
    ax6p4=plt.subplot2grid((2,10),(0,4),colspan=3,aspect='equal')
    ax6p5=plt.subplot2grid((2,10),(1,4),colspan=3,aspect='equal')
    ax6p6=plt.subplot2grid((2,10),(0,7),colspan=3,aspect='equal')
    ax6p7=plt.subplot2grid((2,10),(1,7),colspan=3,aspect='equal')
    xcens=case.x*1e-3#km
    ycens=case.x*1e-3#km
    ax6p1.plot(case.VW.tau.sum(axis=(1,2))/case.VW.tau.sum(axis=(1,2)).max(),case.DHARMA.z/1e3)
#    ax6p1.set_xscale('log')
    zi=95
    Re_vw_tau=case.VW.re_tau*case.VW.w_tau
    Ve_vw_tau=case.VW.ve_tau*case.VW.w_tau
    Tu_vw_tau=case.VW.tau*case.VW.w_tau
    line,=ax6p1.plot(np.array([0,1]),case.DHARMA.z[zi]/1e3*np.array([1,1]),'k-',linewidth=3)
    ax6p1.set_ylim(0,case.DHARMA.z[-1]/1e3)
    ctf1=[ax6p2.contourf(xcens,ycens,case.VW.re_tau[zi,:],v['re'],cmap=plt.cm.jet)]
    ctf2=[ax6p3.contourf(xcens,ycens,np.trapz(Re_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),v['revw'],cmap=plt.cm.jet)]
    ctf3=[ax6p4.contourf(xcens,ycens,case.VW.ve_tau[zi,:],v['ve'],cmap=plt.cm.jet)]
    ctf4=[ax6p5.contourf(xcens,ycens,np.trapz(Ve_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),v['vevw'],cmap=plt.cm.jet)]
    ctf5=[ax6p6.contourf(xcens,ycens,case.dTau[band,zi,:],v['tu'],cmap=plt.cm.jet)]
    ctf6=[ax6p7.contourf(xcens,ycens,np.trapz(Tu_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),v['tu'],cmap=plt.cm.jet)]
    cpn.add_cb(fig6,ctf1[0],ax6p2,pad=0.3,ticks=c['re']  ,label=r'$r_e(Z)$ $\mu m$')
    cpn.add_cb(fig6,ctf2[0],ax6p3,pad=0.3,ticks=c['revw'],label=r'$r_e^{vw}(Z)$ $\mu m$')
    cpn.add_cb(fig6,ctf3[0],ax6p4,pad=0.3,ticks=c['ve']  ,label=r'$v_e(Z)$')
    cpn.add_cb(fig6,ctf4[0],ax6p5,pad=0.3,ticks=c['vevw'],label=r'$v_e^{vw}(Z)$')
    cpn.add_cb(fig6,ctf5[0],ax6p6,pad=0.3,ticks=c['tu']  ,label=r'$COT(Z)$')
    cpn.add_cb(fig6,ctf6[0],ax6p7,pad=0.3,ticks=c['tu']  ,label=r'$COT^{vw}(Z)$')
    ax6p1.set_ylabel('Altitude (km)')
    ax6p1.set_xticklabels([])
    fig6.tight_layout(rect=[0,0,1,0.96])
#    fig6.show()
    def updatefig6(zi):
        line.set_ydata(case.DHARMA.z[zi]/1e3*np.array([1,1]))
        ctf1[0]=[ax6p2.contourf(xcens,ycens,case.VW.re_tau[zi,:],v['re'],cmap=plt.cm.jet)]
        ctf2[0]=[ax6p3.contourf(xcens,ycens,np.trapz(Re_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),v['revw'],cmap=plt.cm.jet)]
        ctf3[0]=[ax6p4.contourf(xcens,ycens,case.VW.ve_tau[zi,:],v['ve'],cmap=plt.cm.jet)]
        ctf4[0]=[ax6p5.contourf(xcens,ycens,np.trapz(Ve_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),v['vevw'],cmap=plt.cm.jet)]
        ctf5[0]=[ax6p6.contourf(xcens,ycens,case.dTau[band,zi,:],v['tu'],cmap=plt.cm.jet)]
        ctf6[0]=[ax6p7.contourf(xcens,ycens,np.trapz(Tu_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),v['tu'],cmap=plt.cm.jet)]
        plt.draw()
    t0=time.time()
    print('Generating GIF ...')
    anim = FuncAnimation(fig6, updatefig6, frames=np.arange(0, 96), interval=100)
    anim.save('GIFs/'+fig6_ttl+'dpi200'+version+'.gif', dpi=200, writer='imagemagick')
    t1=time.time()
    print('%d mins elapsed!'%(t1/60-t0/60))

