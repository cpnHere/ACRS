#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Fri Apr  5 15:53:02 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To test and analyze results
"""
import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap

from cpnLES_MSCARTlib import LES_case, DHARMA_onmp, POLCARTdset
from cpnRetrievalslib import  Pol_ret_P12Lib,Pmat,do_fitting

import cpnCommonlib as cpn
def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'figures/')
class S: pass
def check_fit(Q_in2,abc,ret_Re,ret_Ve,Rsq,P12Lib,LEScase,method):
    '''
    To check the fit between observed F and library F
    LEScase: cpnLES_MSCARTlib.LEScase
    P12Lib: cpnRetrievalslib.P12Lib
    '''
    muV=np.cos(np.deg2rad(LEScase.RT.VZA))
    muS=np.cos(np.deg2rad(180-LEScase.RT.SZA))    
    obsSca=LEScase.RT.ScatA
    P=Pmat(P12Lib.re,P12Lib.ve,P12Lib.bulk_Mie_ang,P12Lib.avP12['0p860'],obsSca,method=method)
    P.set_reve(ret_Re,ret_Ve)
    fig6,ax6=plt.subplots(figsize=(8,4))
    ax6.plot(obsSca[P.Q_a1:P.Q_a2],P.imitateF(obsSca[P.Q_a1:P.Q_a2],*abc),'g.--',label='fit')
    ax6.plot(obsSca[P.Q_a1:P.Q_a2],(-Q_in2*4*(muV+muS))[P.Q_a1:P.Q_a2],'k.-')
#    ax6.plot(RT.ScatA[P.Q_a1:P.Q_a2],-P.getP(RT.ScatA[P.Q_a1:P.Q_a2]),'r--')
    ax6.text(0.5,0.6,r'$r_e$ %0.2f, $v_e$ %0.2f, $R^2$ %0.2f'%(ret_Re,ret_Ve,Rsq),transform=ax6.transAxes)
    ax6.set_xlabel('Scattering Angle')
    ax6.set_ylabel(r'$R_p^*$')
    ax6.legend(frameon=False,loc='best')
    return fig6,ax6

def check_fit2(y,x,abc,ret_Re,ret_Ve,Rsq,P12Lib,LEScase,method):
    '''
    To check the fit between observed F and library F
    LEScase: cpnLES_MSCARTlib.LEScase
    P12Lib: cpnRetrievalslib.P12Lib
    '''
    P=Pmat(P12Lib.re,P12Lib.ve,P12Lib.bulk_Mie_ang,P12Lib.avP12['0p860'],LEScase.RT.ScatA,method=method)
    P.set_reve(ret_Re,ret_Ve)
    fig6,ax6=plt.subplots(figsize=(8,4))
    ax6.plot(x,P.imitateF(x,*abc),'g.--',label='fit')
    ax6.plot(x,y,'k.-')
#    ax6.plot(RT.ScatA[P.Q_a1:P.Q_a2],-P.getP(RT.ScatA[P.Q_a1:P.Q_a2]),'r--')
    ax6.text(0.5,0.6,r'$r_e$ %0.2f, $v_e$ %0.2f, $R^2$ %0.2f'%(ret_Re,ret_Ve,Rsq),transform=ax6.transAxes)
    ax6.set_xlabel('Scattering Angle')
    ax6.set_ylabel(r'$R_p^*$')
    ax6.legend(frameon=False,loc='best')
    return fig6,ax6

if __name__=='__main__':
    cpn.setup_figures(plt)
    #Reading Polarimetric retrievals
#    filename='DYCOMS2_dharma_008036_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH2e6_Breon_full_pol_ret_MPIV2'
    filename='DYCOMS2_dharma_008036_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH2e6_Breon_mean_pol_ret_V3'
#    filename='mpi_out_final'
    data=cpn.load_obj(filename,encoding='latin1')
    ret_Re=data['ret_Re']
    ret_Ve=data['ret_Ve']
    abc   =data['abc']
    Qls   =data['Qls']
    Rsq   =data['Rsq']
    yAll  =data['yAll']
    x     =data['x']
    
    #case definitions
    sza='140'
    vza=0
    band='0p860'
    method='Breon' #see Pmat help 
    
    #Loading P12 library
    P12Lib=Pol_ret_P12Lib(fname='Pol_ret_PM_library_0p860_2p13_V2.hdf5')
    P12Lib.loadP12Lib()
    
    #Reading RT simulations
    LEScase=LES_case('DYCOMS2_dharma_008036_b'+band+'_MSCART_SZA'+str(sza)+'_SAA000_VAA000plus_NPH2e6.hdf5',\
                     '/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/',\
                     RT1Dname='DYCOMS2_dharma_008036_b'+band+'_MSCART_1D_bins_SZA'+str(sza)+'_SAA000_VAA000plus_NPH1e5.hdf5')
    #Reading NJK retrievals
    ret_save_dir='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Retrievals/NJK_retrievals/data/'
    NJK=cpn.load_obj(ret_save_dir+'NJK_retrievals_DYCOMS2_dharma_008036_b2p13_b0p860_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH2e6',\
                               encoding='latin1')
    
    #Checking the fit full
    '''
    X,Y=120,120
    Q_in2= LEScase.RT.MeanPRad[:,X,Y,1]
    abc=abc[X,Y,:]
    ret_Re=ret_Re[X,Y]
    ret_Ve=ret_Ve[Y,X]
    Rsq=Rsq[X,Y]
    y=yAll[X,Y]
    '''
    #Checking the fit mean
    Q_in2= LEScase.RT.MeanPRad[:,:,:,1].mean(axis=1).mean(axis=1)
    y=yAll
#    '''
    fig1_ttl=filename+'_checking_fit'
    fig1,ax1=check_fit2(y,x,abc,ret_Re,ret_Ve,Rsq,P12Lib,LEScase,method)
    fig1.suptitle("\n".join(wrap(fig1_ttl,70)),size=12)
    fig1.tight_layout(rect=[0,0,1,0.96])    
    fig1.show()
#    '''
    #Comparing with NJK
#    '''
    fig2,ax2=plt.subplots()
    ax2.plot(NJK.re,ret_Re,'b.')
    cpn.Corr_plot_frames(ax2)
    fig2.show()
#    '''
    #Cost functiona analysis to check the convergence
#    '''
    Rsqs=np.loadtxt('Rsqs_'+filename+'.dat',skiprows=2)
    Rsqs_max=np.loadtxt('Rsqs_max_'+filename+'.dat',skiprows=2)
    fig3,ax3=plt.subplots(1,2,figsize=(11,5))
    fig3_ttl=filename+'_Rsq_convergence'
    ax3[0].set_title('Rsqs',size=10)#-----------------
    ctf=ax3[0].scatter(Rsqs[:,1],Rsqs[:,2],c=Rsqs[:,0],cmap=plt.cm.jet,zorder=2)
    cpn.add_common_cb(fig3,ctf,label=r'$R^2$')
    ax3[1].set_title('Rsqs_max',size=10)#----------------
    ax3[1].plot(Rsqs_max[:,1],Rsqs_max[:,2], 'C3', zorder=1, lw=1,color='k')
    ctf=ax3[1].scatter(Rsqs_max[:,1],Rsqs_max[:,2],c=Rsqs_max[:,0],cmap=plt.cm.jet,zorder=2)
    #-----------------------------------------
    for n, ax in enumerate(ax3.flat):
        ax.set_xlim([0,40])
        ax.set_ylim([0,0.3])
    fig3.suptitle(fig3_ttl,size=12)
    ax3[0].set_xlabel('CER',size=12)
    ax3[1].set_xlabel('CER',size=12)
    ax3[0].set_ylabel('CEV',size=12)
    fig3.show()
#    '''
    
   
    fig5,ax5=plt.subplots()
    fig5_ttl=RT.fname[0].astype(str).split('.',1)[0]+'_whole_domain_pol_ret_curve_fitting'
    fig5.suptitle("\n".join(wrap(fig5_ttl,65)))
    ctf=ax5.contourf(LEScase.xcens,LEScase.ycens,ret_Re,np.linspace(10,20,100),extend='both')
    fig5.colorbar(ctf,ax=ax5,ticks=np.arange(10,20.1,2),label=r'$r_e$ ($\mu m$)')
    ax5.set_xlabel('km')
    ax5.set_ylabel('km')
    fig5.show()
