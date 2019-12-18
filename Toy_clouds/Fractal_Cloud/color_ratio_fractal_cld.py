#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Aug 31 12:11:59 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Fractal cloud (for Dan)
09/29/2017:
    To use different fractal cloud generations.
03/08/2018:
    The ACA study
    Color ration method for ACA above a fractal cloud

"""
import matplotlib.pyplot as plt
import numpy as np
import cpnMicrophylib as Miclib
import cpnLES_MSCARTlib    
import cpnCommonlib as cpn
from textwrap import wrap
import pickle, os
from scipy.fftpack import fft


class pkl_classes(object):
    def __init__(self,):
        self.class_names=[]

def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'results_aca/figures/')

def readOnlyCloud(band,sza):
    #Only cloud case
    onlyCld=cpnLES_MSCARTlib.POLCARTdset('onlyCld','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/')
    onlyCld.readPOLCARThdf5('fractal_StCu_b865re12ve05_x40km_H0p5_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/results/fracb865/')
    onlyCld.dset='Cloud_only_3D'
    onlyCld_1D=cpnLES_MSCARTlib.POLCARTdset('onlyCld_1D','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb'+band+'_bins/')
    onlyCld_1D.readPOLCARThdf5('fractal_StCu_b865re12ve05_x40km_H0p5_1D_bins_H0p5_MSCART_SZA120_SAA000_VAA000plus_NPH1e5.hdf5',\
                                  dpath='/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/results/fracb'+band+'_bins/')
    onlyCld_1D.dset='Cloud_only_1D'
    return onlyCld, onlyCld_1D

def biases(resol,AOT0p25_0p860,AOT0p25_0p860_1D):
    X=resol
    R8653D=AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T
    R8651D=AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T
    if R8651D.ndim==2:
        R8651D=np.mean(R8651D,axis=1)
        R8653D=np.mean(R8653D,axis=1)
    movingXm_MeanPRad3D=cpn.movingaverage(R8653D,X/10)
    movingXm_MeanPRad1D=cpn.movingaverage(R8651D,X/10)
    CR1DXm=cpn.movingaverage(CR1D,X/10)
    CR3DXm=cpn.movingaverage(CR,X/10)
    
    
    R3Dmin1D=movingXm_MeanPRad3D-movingXm_MeanPRad1D
    CR3Dmin1D=CR3DXm-CR1DXm
    xwdth=abs(R3Dmin1D).max()
    ywdth=abs(CR3Dmin1D).max()
    
    return R3Dmin1D,CR3Dmin1D,xwdth,ywdth
    
def plotbiasses(ax,resol,AOT0p25_0p860,AOT0p25_0p860_1D):
    b=biases(resol,AOT0p25_0p860,AOT0p25_0p860_1D) 
    R=b[0];CR=b[1];xw=b[2];yw=b[3]
    ax.plot(R,CR,'.',markersize=2)
    ax.plot([-xw,xw],[0,0],'k-',[0,0],[-yw,yw],'k-')
    ax.set_xlim(-xw,xw)
    ax.set_ylim(-yw,yw)
    ax.set_title('%d m'%resol,size=11)

if __name__=='__main__':
    # Reading Fractal cloud case for Dan's thesis ################################
    #cld=Miclib.frac_physics('physical_scene_case_new.h5')   
    #cld.read() 
    #cld.plot_cld()
    ##############################################################################
    sza='120'
    VZA=0
    #band=str(sys.argv[1])
    #sza=str(sys.argv[2])
    #vza=str(sys.argv[3])#vza=d to go with default values
    #------------------------------------------------------------------------------
    #Reading all information about the ACA scene
    with open('fractal_StCu_ACA_b0p860_x40km_H0p5_AOT0p25_H0p5.pkl', "rb") as f:
        pkl=pickle.load(f)
        MSCART_field=pickle.load(f)
        fcld=pickle.load(f)
    
    ss_file=pkl.SS_properties
    with open(ss_file, "rb") as f:
        pkl2=pickle.load(f)
        bulk_cloud=pickle.load(f)    
        bulk_aero=pickle.load(f)
    #------------------------------------------------------------------------------
    AOT0p25_0p860=cpnLES_MSCARTlib.POLCARTdset('AOT0p25_0p860','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/')
    AOT0p25_0p860.readPOLCARThdf5('fractal_StCu_ACA_b0p860_x40km_H0p5_AOT0p25_H0p5_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5','results_aca/b0p860/')
    AOT0p25_0p860.dset='AOT0p25_0p860_3D'
    frac=cpnLES_MSCARTlib.LES_field(AOT0p25_0p860.fname.split('_MSCART',2)[0]+'.nc')#DYCOM field file for MSCART
    frac.readLES_field()
    xcens=(frac.xgrd[1:]+frac.xgrd[0:-1])/2
    
    AOT0p25_0p470=cpnLES_MSCARTlib.POLCARTdset('AOT0p25_0p470','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/')
    AOT0p25_0p470.readPOLCARThdf5('fractal_StCu_ACA_b0p470_x40km_H0p5_AOT0p25_H0p5_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5','results_aca/b0p470/')
    AOT0p25_0p470.dset='AOT0p25_0p470_3D'


    onlyCld,onlyCld_1D=readOnlyCloud('0p865',sza)
    
    AOT0p25_0p860_1D=cpnLES_MSCARTlib.POLCARTdset('AOT0p25_0p860_1D','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb0p860_bins/')
    AOT0p25_0p860_1D.readPOLCARThdf5('fractal_StCu_ACA_b0p860_x40km_H0p5_AOT0p25_H0p5_MSCART_1D_bins_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                                  dpath='1Druns/results/fracb0p860_bins/')
    AOT0p25_0p860_1D.dset='AOT0p25_0p860_1D'
    
    AOT0p25_0p470_1D=cpnLES_MSCARTlib.POLCARTdset('AOT0p25_0p470_1D','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb0p470_bins/')
    AOT0p25_0p470_1D.readPOLCARThdf5('fractal_StCu_ACA_b0p470_x40km_H0p5_AOT0p25_H0p5_MSCART_1D_bins_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                                  dpath='1Druns/results/fracb0p470_bins/')
    AOT0p25_0p470_1D.dset='AOT0p25_0p470_1D'
    
    AOT0p25_0p860.avg50m_MeanPRad=AOT0p25_0p860.avg_resolution(5)
    AOT0p25_0p860.avg100m_MeanPRad=AOT0p25_0p860.avg_resolution(10)
    AOT0p25_0p860.avg200m_MeanPRad=AOT0p25_0p860.avg_resolution(20)
    AOT0p25_0p860_1D.avg50m_MeanPRad=AOT0p25_0p860_1D.avg_resolution(5)
    AOT0p25_0p860_1D.avg100m_MeanPRad=AOT0p25_0p860_1D.avg_resolution(10)
    AOT0p25_0p860_1D.avg200m_MeanPRad=AOT0p25_0p860_1D.avg_resolution(20)
    
    AOT0p25_0p470.avg50m_MeanPRad=AOT0p25_0p470.avg_resolution(5)
    AOT0p25_0p470.avg100m_MeanPRad=AOT0p25_0p470.avg_resolution(10)
    AOT0p25_0p470.avg200m_MeanPRad=AOT0p25_0p470.avg_resolution(20)
    AOT0p25_0p470_1D.avg50m_MeanPRad=AOT0p25_0p470_1D.avg_resolution(5)
    AOT0p25_0p470_1D.avg100m_MeanPRad=AOT0p25_0p470_1D.avg_resolution(10)
    AOT0p25_0p470_1D.avg200m_MeanPRad=AOT0p25_0p470_1D.avg_resolution(20)

    #%%%%%%%%%%%%%%%%%% figure 01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #I and -Q along x axis.
    print(sza)
    #if (band=='0p865' or band=='2p13') and vza=='d':
    #    if sza=='120':
    #        VZA=-20.0
    #    elif sza=='140':
    #        VZA=1.0
    #    elif sza=='160':
    #        VZA=20.0
    #    else:
    #        print('Check and assign best VZA angle to observe the 3D features in figure 1.')
    #elif band=='3p75' and vza=='d':
    #    if sza=='120':
    #        VZA=-35.0
    #    elif sza=='140':
    #        VZA=-15.0
    #    elif sza=='160':
    #        VZA=1.0
    #    else:
    #        print('Check and assign best VZA angle to observe the 3D features in figure 1.')
    #if vza!='d':
    #    VZA=int(vza)
    
    #VZA=20.0 # note that, there are two set of values for 0.0 VZA
    
    fig1,ax1=plt.subplots(2,1,sharex=True)
    fig1_ttl=AOT0p25_0p860.fname.split('.')[0]+'_VZA_'+str(VZA)
    if VZA==0:
        ax1[0].plot(xcens,np.mean(AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,axis=1),'b-',label=AOT0p25_0p860.dset)
    #    ax1[0].plot(xcens,np.mean(onlyCld.MeanPRad[onlyCld.VZA==VZA,:,0].T,axis=1),'r-',label=onlyCld.dset)
    #    ax1[0].plot(xcens,np.mean(onlyCld_1D.MeanPRad[onlyCld_1D.VZA==VZA,:,0].T,axis=1),'r--',label=onlyCld_1D.dset)
        ax1[0].plot(xcens,np.mean(AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,axis=1),'g--',label=AOT0p25_0p860_1D.dset)
    else:
        ax1[0].plot(xcens,AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,'b-',label=AOT0p25_0p860.dset)
    #    ax1[0].plot(xcens,onlyCld.MeanPRad[onlyCld.VZA==VZA,:,0].T,'r-',label=onlyCld.dset)
    #    ax1[0].plot(xcens,onlyCld_1D.MeanPRad[onlyCld_1D.VZA==VZA,:,0].T,'r--',label=onlyCld_1D.dset)
        ax1[0].plot(xcens,AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,'g--',label=AOT0p25_0p860_1D.dset)
        
    ax1[0].set_ylabel('I',size=10)
    
    ax1[0].legend(loc='best',prop={'size':10})
    if VZA==0:
        ax1[1].plot(xcens,-np.mean(AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,1].T,axis=1),'b-')
    #    ax1[1].plot(xcens,-np.mean(onlyCld.MeanPRad[onlyCld.VZA==VZA,:,1].T,axis=1),'r-')
    #    ax1[1].plot(xcens,-np.mean(onlyCld_1D.MeanPRad[onlyCld_1D.VZA==VZA,:,1].T,axis=1),'r--')
        ax1[1].plot(xcens,-np.mean(AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,1].T,axis=1),'g--')
    else:
        ax1[1].plot(xcens,-AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,1].T,'b-')
    #    ax1[1].plot(xcens,-onlyCld.MeanPRad[onlyCld.VZA==VZA,:,1].T,'r-')
    #    ax1[1].plot(xcens,-onlyCld_1D.MeanPRad[onlyCld_1D.VZA==VZA,:,1].T,'r--')
        ax1[1].plot(xcens,-AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,1].T,'g--')
    ax1[1].set_ylabel('-Q',size=10)
    ax1[1].set_xlabel('x km',size=10)
    fig1.suptitle("\n".join(wrap(fig1_ttl,60)))
    fig1.show()
    
    '''
    ################################################################################
    Figure 4 and 5 : Color ratio (R470/R860) and R860 from 1D and 3D simulations
    '''
    fig4,ax4=plt.subplots(1,2,figsize=(12,5))
    fig4_ttl='CR_and_R865_aca_fractal_v4'
    fig4.suptitle(fig4_ttl)
    ax4[0].text(0.01, 0.95, '(a)', transform=ax4[0].transAxes,size=15)
    #ax4[0].plot(xcens10m,CR10m,'r-',linewidth=1.0,label='3D(10m)')
    ax4[0].set_ylabel('$R_{0.470}/R_{0.860}$')
    ax4[0].legend(loc='lower right',prop={'size':11})
    ax4[1].text(0.01, 0.95, '(b)', transform=ax4[1].transAxes,size=15)
    if VZA==0:
        CR=np.mean(AOT0p25_0p470.MeanPRad[AOT0p25_0p470.VZA==VZA,:,0].T,axis=1)/np.mean(AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==0,:,0].T,axis=1)
        CR1D=np.mean(AOT0p25_0p470_1D.MeanPRad[AOT0p25_0p470_1D.VZA==VZA,:,0].T,axis=1)/np.mean(AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==0,:,0].T,axis=1)
    #    ax4[1].plot(xcens,np.mean(AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,axis=1),'k-',label=AOT0p25_0p860.dset)
    #    ax4[1].plot(xcens,np.mean(AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,axis=1),'b-',label=AOT0p25_0p860_1D.dset)
        AOT0p25_0p860.moving1km_MeanPRad=cpn.movingaverage(np.mean(AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,axis=1),100)
        AOT0p25_0p860_1D.moving1km_MeanPRad=cpn.movingaverage(np.mean(AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,axis=1),100)
    else:
        CR=AOT0p25_0p470.MeanPRad[AOT0p25_0p470.VZA==VZA,:,0].T/AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T
        CR1D=AOT0p25_0p470_1D.MeanPRad[AOT0p25_0p470_1D.VZA==VZA,:,0].T/AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T
    #    ax4[1].plot(xcens,AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,'k-',label=AOT0p25_0p860.dset)
    #    ax4[1].plot(xcens,AOT0p25_0p860_1D.MeanePRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,'b-',label=AOT0p25_0p860_1D.dset)
        AOT0p25_0p860.moving1km_MeanPRad=cpn.movingaverage(AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,100)
        AOT0p25_0p860_1D.moving1km_MeanPRad=cpn.movingaverage(AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,100)    
    #ax4[0].plot(xcens,CR,'k-',label='3D')
    xcens1km=cpn.movingaverage(xcens,100)
    CR1D1km=cpn.movingaverage(CR1D,100)
    CR3D1km=cpn.movingaverage(CR,100)
    
    ax4[1].plot(xcens1km,AOT0p25_0p860.moving1km_MeanPRad,'k-',label=AOT0p25_0p860.dset+'_1km')
    ax4[1].plot(xcens1km,AOT0p25_0p860_1D.moving1km_MeanPRad,'b-',label=AOT0p25_0p860_1D.dset+'_1km')
    ax4[0].plot(xcens1km,CR1D1km,'b-',label='1D_1km')
    ax4[0].plot(xcens1km,CR3D1km,'k-',label='3D_1km')
    
    #ax4[1].set_xlabel('x (km)')
    ax4[1].set_ylabel(r'$R_{0.865}$')
    ax4[1].legend(loc='lower left',prop={'size':11})
    ax4[0].legend(loc='best',prop={'size':11})
    ax4[0].grid()
    ax4[1].grid()
    #fig4.tight_layout()
    fig4.show()

    #----------------------------------------------------------------------------------
    #Figure 5

   
    fig5,ax5=plt.subplots(2,2)
    fig5_ttl='ACA_fractal_CR_and_r860_biases_corr'
    plotbiasses(ax5[0,0],10  ,AOT0p25_0p860,AOT0p25_0p860_1D);
    plotbiasses(ax5[0,1],50  ,AOT0p25_0p860,AOT0p25_0p860_1D)
    plotbiasses(ax5[1,0],500 ,AOT0p25_0p860,AOT0p25_0p860_1D);
    plotbiasses(ax5[1,1],1000,AOT0p25_0p860,AOT0p25_0p860_1D)
    ax5[1,0].set_xlabel(r'$R_{865}^{3D-1D}$')
    ax5[1,1].set_xlabel(r'$R_{865}^{3D-1D}$')
    ax5[0,0].set_ylabel(r'$CR^{3D-1D}$')
    ax5[1,0].set_ylabel(r'$CR^{3D-1D}$')
    fig5.suptitle(fig5_ttl,size=11)
    fig5.show()

    #========3d-1d difference =====================================================
    #def avg_resolution(obj,ftr):
    #    # Averaging radiance over given factor.
    #    # ex. If we have radiance field with 1000 x elements, if ftr=5 we'll have
    #    # a new radiance field with 1000/5 x. Every 5 elements will be averaged.
    #    #obj: POLCARTdset object
    #    #ftr: factor to do averaging.
    #    shape=obj.MeanPRad.shape
    #    newMPRad=np.zeros((shape[0],shape[1]/ftr,shape[2]))
    #    for i in range(0,obj.MeanPRad.shape[0]):
    #        MPRad=np.squeeze(obj.MeanPRad[i,:,:])
    #        rmv_ix=MPRad.shape[0]-MPRad.shape[0]/ftr*ftr
    #        B=MPRad[0:-rmv_ix].reshape(-1,ftr,shape[2])
    #        newMPRad[i,:,:]=np.mean(B,1)
    #    return newMPRad
    
    fig3,ax3=plt.subplots()
    
    
    #------------------------------------------------------------------------------
    if VZA==0:
        del_t=np.squeeze(np.mean(AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,axis=1)-np.mean(AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,axis=1))
    else:
        del_t=np.squeeze(AOT0p25_0p860.MeanPRad[AOT0p25_0p860.VZA==VZA,:,0]-AOT0p25_0p860_1D.MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0])
    #weights = np.ones_like(del_t)/len(del_t)       
    #ax3.hist(abs(del_t),50,weights=weights)
    #del_il=del_t[del_t>0]
    #del_sh=del_t[del_t<0]
    weights = np.ones_like(del_t)/len(del_t)
    ax3.hist(del_t,40,weights=weights,histtype='step',label='10m')
    #-----------------50 m---------------------------------
    if VZA==0:
        del_t=np.squeeze(np.mean(AOT0p25_0p860.avg50m_MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,axis=1)-\
                                 np.mean(AOT0p25_0p860_1D.avg50m_MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,axis=1))
    else:
        del_t=np.squeeze(AOT0p25_0p860.avg50m_MeanPRad[AOT0p25_0p860.VZA==VZA,:,0]-\
                         AOT0p25_0p860_1D.avg50m_MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0])
    weights = np.ones_like(del_t)/len(del_t)
    ax3.hist(del_t,40,weights=weights,histtype='step',label='50m')
    #-----------------100 m---------------------------------
    if VZA==0:
        del_t=np.squeeze(np.mean(AOT0p25_0p860.avg100m_MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,axis=1)-\
                                 np.mean(AOT0p25_0p860_1D.avg100m_MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,axis=1))
    else:
        del_t=np.squeeze(AOT0p25_0p860.avg100m_MeanPRad[AOT0p25_0p860.VZA==VZA,:,0]-\
                         AOT0p25_0p860_1D.avg100m_MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0])
    weights = np.ones_like(del_t)/len(del_t)
    ax3.hist(del_t,40,weights=weights,histtype='step',label='100m')
    #-------------------200m-----------------------------------
    if VZA==0:
        del_t=np.squeeze(np.mean(AOT0p25_0p860.avg200m_MeanPRad[AOT0p25_0p860.VZA==VZA,:,0].T,axis=1)-\
                                 np.mean(AOT0p25_0p860_1D.avg200m_MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0].T,axis=1))
    else:
        del_t=np.squeeze(AOT0p25_0p860.avg200m_MeanPRad[AOT0p25_0p860.VZA==VZA,:,0]-\
                         AOT0p25_0p860_1D.avg200m_MeanPRad[AOT0p25_0p860_1D.VZA==VZA,:,0])
    weights = np.ones_like(del_t)/len(del_t)
    ax3.hist(del_t,40,weights=weights,histtype='step',label='200m')
    
    
    ax3.legend(loc='best')
    fig3_ttl=AOT0p25_0p860.fname.split('.')[0]+'_VZA_'+str(VZA)+'_delta_Reflectance_3Dmin1D'
    ax3.set_xlabel(r'$R_{3D}-R_{1D}$')
    fig3.suptitle("\n".join(wrap(fig3_ttl,60)))
    
    
    
    #=======3d-1d difference end ==================================================

'''
#=======3d-1d difference separate figures =====================================

def refhist(ax):
    fig_ttl=fracSZAx.fname.split('.')[0]+'_VZA_'+str(VZA)+'_delta_Reflectance_3Dmin1D_res10m'
    if VZA==0:
        del_t=np.squeeze(np.mean(fracSZAx.MeanPRad[fracSZAx.VZA==VZA,:,0].T,axis=1)-np.mean(frac1D.MeanPRad[frac1D.VZA==VZA,:,0].T,axis=1))
    else:
        del_t=np.squeeze(fracSZAx.MeanPRad[fracSZAx.VZA==VZA,:,0]-frac1D.MeanPRad[frac1D.VZA==VZA,:,0])
    #weights = np.ones_like(del_t)/len(del_t)       
    #ax3.hist(abs(del_t),50,weights=weights)
    #del_il=del_t[del_t>0]
    #del_sh=del_t[del_t<0]
    weights = np.ones_like(del_t)/len(del_t)
    ax.hist(del_t,40,weights=weights,histtype='step',label='10m')
    
    return fig_ttl


def refhist50(ax):
    fig_ttl=fracSZAx.fname.split('.')[0]+'_VZA_'+str(VZA)+'_delta_Reflectance_3Dmin1D_res50m'
    if VZA==0:
        del_t=np.squeeze(np.mean(fracSZAx.avg50m_MeanPRad[fracSZAx.VZA==VZA,:,0].T,axis=1)-\
                                 np.mean(frac1D.avg50m_MeanPRad[frac1D.VZA==VZA,:,0].T,axis=1))
    else:
        del_t=np.squeeze(fracSZAx.avg50m_MeanPRad[fracSZAx.VZA==VZA,:,0]-\
                         frac1D.avg50m_MeanPRad[frac1D.VZA==VZA,:,0])
    weights = np.ones_like(del_t)/len(del_t)
    ax.hist(del_t,40,weights=weights,histtype='step',label='50m')
    return fig_ttl
    
def refhist100(ax):
    fig_ttl=fracSZAx.fname.split('.')[0]+'_VZA_'+str(VZA)+'_delta_Reflectance_3Dmin1D_res100m'
    if VZA==0:
        del_t=np.squeeze(np.mean(fracSZAx.avg100m_MeanPRad[fracSZAx.VZA==VZA,:,0].T,axis=1)-\
                                 np.mean(frac1D.avg100m_MeanPRad[frac1D.VZA==VZA,:,0].T,axis=1))
    else:
        del_t=np.squeeze(fracSZAx.avg100m_MeanPRad[fracSZAx.VZA==VZA,:,0]-\
                         frac1D.avg100m_MeanPRad[frac1D.VZA==VZA,:,0])
    weights = np.ones_like(del_t)/len(del_t)
    ax.hist(del_t,40,weights=weights,histtype='step',label='100m')
    return fig_ttl
    
def refhist200(ax):
    fig_ttl=fracSZAx.fname.split('.')[0]+'_VZA_'+str(VZA)+'_delta_Reflectance_3Dmin1D_res200m'
    if VZA==0:
        del_t=np.squeeze(np.mean(fracSZAx.avg200m_MeanPRad[fracSZAx.VZA==VZA,:,0].T,axis=1)-\
                                 np.mean(frac1D.avg200m_MeanPRad[frac1D.VZA==VZA,:,0].T,axis=1))
    else:
        del_t=np.squeeze(fracSZAx.avg200m_MeanPRad[fracSZAx.VZA==VZA,:,0]-\
                         frac1D.avg200m_MeanPRad[frac1D.VZA==VZA,:,0])
    weights = np.ones_like(del_t)/len(del_t)
    ax.hist(del_t,40,weights=weights,histtype='step',label='200m')
    return fig_ttl

fracSZAx.avg50m_MeanPRad=fracSZAx.avg_resolution(5)
fracSZAx.avg100m_MeanPRad=fracSZAx.avg_resolution(10)
fracSZAx.avg200m_MeanPRad=fracSZAx.avg_resolution(20)
frac1D.avg50m_MeanPRad=frac1D.avg_resolution(5)
frac1D.avg100m_MeanPRad=frac1D.avg_resolution(10)
frac1D.avg200m_MeanPRad=frac1D.avg_resolution(20)

fig3,ax3=plt.subplots()
ax3.set_xlabel(r'$R_{3D}-R_{1D}$')
fig3_ttl=refhist200(ax3)
fig3.suptitle("\n".join(wrap(fig3_ttl,60)))
#=======3d-1d difference separate figures END =====================================
#============== angular dependancy of Avg. Q ==================================
dom_meanMPRad3D=np.mean(fracSZAx.MeanPRad,axis=1)
dom_meanMPRad1D=np.mean(frac1D.MeanPRad,axis=1)
fig4_ttl=fracSZAx.fname.split('.')[0]+'_dom_avg_-Q'
fig4,ax4=plt.subplots()
ax4.plot(fracSZAx.ScatA,-dom_meanMPRad3D[:,1],'b-',label='3D')
ax4.plot(frac1D.ScatA,-dom_meanMPRad1D[:,1],'r-',label='1D')
ax4.legend(loc='best')
fig4.suptitle("\n".join(wrap(fig4_ttl,60)))
#============== angular dependancy of Avg. Q END ==================================
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

if band=='0p865:
    pb_ix=np.where(np.isclose(fracSZAx.ScatA,142,atol=0.1))
    sb_ix=np.where(np.isclose(fracSZAx.ScatA,153,atol=0.1))
elif band=='2p13':
    pb_ix=np.where(np.isclose(fracSZAx.ScatA,141,atol=0.1))
    sb_ix=np.where(np.isclose(fracSZAx.ScatA,160,atol=0.1))
elif band=='3p75':
    pb_ix=np.where(np.isclose(fracSZAx.ScatA,158,atol=0.1))
    sb_ix=np.where(np.isclose(fracSZAx.ScatA,170,atol=0.1))
else:    
    pb,sb=Miclib.find_bows(frac.scatAng,frac.PMatVal[0,:,0])
    pb_ix=np.where(np.isclose(fracSZAx.ScatA,pb,atol=0.1))
    sb_ix=np.where(np.isclose(fracSZAx.ScatA,sb,atol=0.1))

fig2,ax2=plt.subplots(3,1,figsize=(7,10),sharex=True)
fig2_ttl=contourf_Iqu(fig2,ax2[0],1,fracSZAx.VZA[pb_ix],fracSZAx.VZA[sb_ix],fracSZAx.VZA,fracSZAx.MeanPRad)
contourf_Iqu(fig2,ax2[1],1,frac1D.VZA[pb_ix],frac1D.VZA[sb_ix],frac1D.VZA,frac1D.MeanPRad)
fig2_ttl=contourf_Iqu(fig2,ax2[2],1,fracSZAx.VZA[pb_ix],fracSZAx.VZA[sb_ix],fracSZAx.VZA,\
                      fracSZAx.MeanPRad-frac1D.MeanPRad,v=np.linspace(-0.01,0.01,num=50),\
                        ticks=np.arange(-0.01,0.011,0.005),cmap=plt.cm.coolwarm)
ax2[2].set_xlabel('x (km)',size=10,y=-0.5)
fig2_ttl=fracSZAx.fname.split('.',1)[0]+'_'+fig2_ttl
fig2.suptitle("\n".join(wrap(fig2_ttl,60)))
ax2[0].set_title(r'$-Q_{3D}$')
ax2[1].set_title(r'$-Q_{1D}$')
ax2[2].set_title(r'$-(Q_{3D}-Q_{1D})$')

#=============

#savefig(fig1,fig1_ttl)
#savefig(fig2,fig2_ttl)

#plt.show()
'''
