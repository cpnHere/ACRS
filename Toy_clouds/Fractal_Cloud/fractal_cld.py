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

"""
import matplotlib.pyplot as plt
import numpy as np
import cpnMicrophylib as Miclib
import cpnLES_MSCARTlib    
from textwrap import wrap
import sys
import pickle
from scipy.fftpack import fft


class pkl_classes(object):
    def __init__(self,):
        self.class_names=[]

def savefig(fig,fig_ttl):
    fig.savefig('results_aca/figures/'+fig_ttl+'.png',format='png',dpi=200)
    print('results_aca/figures/'+fig_ttl+'.png SAVED!')

# Reading Fractal cloud case for Dan's thesis ################################
#cld=Miclib.frac_physics('physical_scene_case_new.h5')   
#cld.read() 
#cld.plot_cld()
##############################################################################
band='0p470'
sza='120'
vza='0'
#band=str(sys.argv[1])
#sza=str(sys.argv[2])
#vza=str(sys.argv[3])#vza=d to go with default values
#------------------------------------------------------------------------------
#Reading all information about the ACA scene
with open('fractal_StCu_ACA_b'+band+'_x40km_H0p5_AOT0p25_H0p5.pkl', "rb") as f:
    pkl=pickle.load(f)
    MSCART_field=pickle.load(f)
    fcld=pickle.load(f)

ss_file=pkl.SS_properties
with open(ss_file, "rb") as f:
    pkl2=pickle.load(f)
    bulk_cloud=pickle.load(f)    
    bulk_aero=pickle.load(f)
#------------------------------------------------------------------------------
re=12
xw=40
ve='05'
H='H0p5'

AOT0p25=cpnLES_MSCARTlib.POLCARTdset('AOT0p25','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/')
AOT0p25.readPOLCARThdf5('fractal_StCu_ACA_b'+band+'_x40km_H0p5_AOT0p25_H0p5_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5','results_aca/b'+band+'/')
AOT0p25.dset='AOT0p25_3D'
frac=cpnLES_MSCARTlib.LES_field(AOT0p25.fname.split('_MSCART',2)[0]+'.nc')#DYCOM field file for MSCART
frac.readLES_field()

def readOnlyCloud():
    #Only cloud case
    onlyCld=cpnLES_MSCARTlib.POLCARTdset('onlyCld','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/')
    onlyCld.readPOLCARThdf5('fractal_StCu_b865re12ve05_x40km_H0p5_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5','results/fracb865/')
    onlyCld.dset='Cloud_only_3D'
    onlyCld_1D=cpnLES_MSCARTlib.POLCARTdset('onlyCld_1D','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb'+band+'_bins/')
    onlyCld_1D.readPOLCARThdf5('fractal_StCu_b865re12ve05_x40km_H0p5_1D_bins_H0p5_MSCART_SZA120_SAA000_VAA000plus_NPH1e5.hdf5',\
                                  dpath='1Druns/results/fracb'+band+'_bins/')
    onlyCld_1D.dset='Cloud_only_1D'
    return onlyCld, onlyCld_1D

#onlyCld,onlyCld_1D=readOnlyCloud()

AOT0p25_1D=cpnLES_MSCARTlib.POLCARTdset('AOT0p25_1D','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb'+band+'_bins/')
AOT0p25_1D.readPOLCARThdf5('fractal_StCu_ACA_b'+band+'_x40km_H0p5_AOT0p25_H0p5_1D_bins_H0p5_AOT0p25_H0p5_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                              dpath='1Druns/results/fracb'+band+'_bins/')
AOT0p25_1D.dset='AOT0p25_1D'
#%%%%%%%%%%%%%%%%%% figure 01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#I and -Q along x axis.
print(sza)
if (band=='0p865' or band=='2p13') and vza=='d':
    if sza=='120':
        VZA=-20.0
    elif sza=='140':
        VZA=1.0
    elif sza=='160':
        VZA=20.0
    else:
        print('Check and assign best VZA angle to observe the 3D features in figure 1.')
elif band=='3p75' and vza=='d':
    if sza=='120':
        VZA=-35.0
    elif sza=='140':
        VZA=-15.0
    elif sza=='160':
        VZA=1.0
    else:
        print('Check and assign best VZA angle to observe the 3D features in figure 1.')
if vza!='d':
    VZA=int(vza)
           
#VZA=20.0 # note that, there are two set of values for 0.0 VZA
xcens=(frac.xgrd[1:]+frac.xgrd[0:-1])/2
fig1,ax1=plt.subplots(2,1,sharex=True)
fig1_ttl=AOT0p25.fname.split('.')[0]+'_VZA_'+str(VZA)
if VZA==0:
    ax1[0].plot(xcens,np.mean(AOT0p25.MeanPRad[AOT0p25.VZA==VZA,:,0].T,axis=1),'b-',label=AOT0p25.dset)
    ax1[0].plot(xcens,np.mean(onlyCld.MeanPRad[onlyCld.VZA==VZA,:,0].T,axis=1),'r-',label=onlyCld.dset)
    ax1[0].plot(xcens,np.mean(onlyCld_1D.MeanPRad[onlyCld_1D.VZA==VZA,:,0].T,axis=1),'r--',label=onlyCld_1D.dset)
    ax1[0].plot(xcens,np.mean(AOT0p25_1D.MeanPRad[AOT0p25_1D.VZA==VZA,:,0].T,axis=1),'g--',label=AOT0p25_1D.dset)
else:
    ax1[0].plot(xcens,AOT0p25.MeanPRad[AOT0p25.VZA==VZA,:,0].T,'b-',label=AOT0p25.dset)
    ax1[0].plot(xcens,onlyCld.MeanPRad[onlyCld.VZA==VZA,:,0].T,'r-',label=onlyCld.dset)
    ax1[0].plot(xcens,onlyCld_1D.MeanPRad[onlyCld_1D.VZA==VZA,:,0].T,'r--',label=onlyCld_1D.dset)
    ax1[0].plot(xcens,AOT0p25_1D.MeanPRad[AOT0p25_1D.VZA==VZA,:,0].T,'g--',label=AOT0p25_1D.dset)
    

#ax1[0].plot(xcens,I1D,'r-',label='1D')
ax1[0].set_ylabel('I',size=10)

#ax2 = ax1[0].twinx()
#ax2.plot(xcens,np.squeeze(frac.extp3d),'r',label=r'$\tau$')
#ax2.set_ylabel(r'Optical depth $\tau$',size=10)
ax1[0].legend(loc='best',prop={'size':10})
#ax2.legend(loc='best',prop={'size':10})
#ax2.set_ylim(0,15)
if VZA==0:
    ax1[1].plot(xcens,-np.mean(AOT0p25.MeanPRad[AOT0p25.VZA==VZA,:,1].T,axis=1),'b-')
    ax1[1].plot(xcens,-np.mean(onlyCld.MeanPRad[onlyCld.VZA==VZA,:,1].T,axis=1),'r-')
    ax1[1].plot(xcens,-np.mean(onlyCld_1D.MeanPRad[onlyCld_1D.VZA==VZA,:,1].T,axis=1),'r--')
    ax1[1].plot(xcens,-np.mean(AOT0p25_1D.MeanPRad[AOT0p25_1D.VZA==VZA,:,1].T,axis=1),'g--')
else:
    ax1[1].plot(xcens,-AOT0p25.MeanPRad[AOT0p25.VZA==VZA,:,1].T,'b-')
    ax1[1].plot(xcens,-onlyCld.MeanPRad[onlyCld.VZA==VZA,:,1].T,'r-')
    ax1[1].plot(xcens,-onlyCld_1D.MeanPRad[onlyCld_1D.VZA==VZA,:,1].T,'r--')
    ax1[1].plot(xcens,-AOT0p25_1D.MeanPRad[AOT0p25_1D.VZA==VZA,:,1].T,'g--')
ax1[1].set_ylabel('-Q',size=10)
ax1[1].set_xlabel('x km',size=10)
#ax22 = ax1[1].twinx()
#ax22.plot(xcens,np.squeeze(frac.extp3d),'r')
#ax22.set_ylabel(r'Optical depth $\tau$',size=10)
#ax22.set_ylim(0,15)
fig1.suptitle("\n".join(wrap(fig1_ttl,60)))
fig1.show()

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
fracSZAx.avg50m_MeanPRad=fracSZAx.avg_resolution(5)
fracSZAx.avg100m_MeanPRad=fracSZAx.avg_resolution(10)
fracSZAx.avg200m_MeanPRad=fracSZAx.avg_resolution(20)
frac1D.avg50m_MeanPRad=frac1D.avg_resolution(5)
frac1D.avg100m_MeanPRad=frac1D.avg_resolution(10)
frac1D.avg200m_MeanPRad=frac1D.avg_resolution(20)

#------------------------------------------------------------------------------
if VZA==0:
    del_t=np.squeeze(np.mean(fracSZAx.MeanPRad[fracSZAx.VZA==VZA,:,0].T,axis=1)-np.mean(frac1D.MeanPRad[frac1D.VZA==VZA,:,0].T,axis=1))
else:
    del_t=np.squeeze(fracSZAx.MeanPRad[fracSZAx.VZA==VZA,:,0]-frac1D.MeanPRad[frac1D.VZA==VZA,:,0])
#weights = np.ones_like(del_t)/len(del_t)       
#ax3.hist(abs(del_t),50,weights=weights)
#del_il=del_t[del_t>0]
#del_sh=del_t[del_t<0]
weights = np.ones_like(del_t)/len(del_t)
ax3.hist(del_t,40,weights=weights,histtype='step',label='10m')
#-----------------50 m---------------------------------
if VZA==0:
    del_t=np.squeeze(np.mean(fracSZAx.avg50m_MeanPRad[fracSZAx.VZA==VZA,:,0].T,axis=1)-\
                             np.mean(frac1D.avg50m_MeanPRad[frac1D.VZA==VZA,:,0].T,axis=1))
else:
    del_t=np.squeeze(fracSZAx.avg50m_MeanPRad[fracSZAx.VZA==VZA,:,0]-\
                     frac1D.avg50m_MeanPRad[frac1D.VZA==VZA,:,0])
weights = np.ones_like(del_t)/len(del_t)
ax3.hist(del_t,40,weights=weights,histtype='step',label='50m')
#-----------------100 m---------------------------------
if VZA==0:
    del_t=np.squeeze(np.mean(fracSZAx.avg100m_MeanPRad[fracSZAx.VZA==VZA,:,0].T,axis=1)-\
                             np.mean(frac1D.avg100m_MeanPRad[frac1D.VZA==VZA,:,0].T,axis=1))
else:
    del_t=np.squeeze(fracSZAx.avg100m_MeanPRad[fracSZAx.VZA==VZA,:,0]-\
                     frac1D.avg100m_MeanPRad[frac1D.VZA==VZA,:,0])
weights = np.ones_like(del_t)/len(del_t)
ax3.hist(del_t,40,weights=weights,histtype='step',label='100m')
#-------------------200m-----------------------------------
if VZA==0:
    del_t=np.squeeze(np.mean(fracSZAx.avg200m_MeanPRad[fracSZAx.VZA==VZA,:,0].T,axis=1)-\
                             np.mean(frac1D.avg200m_MeanPRad[frac1D.VZA==VZA,:,0].T,axis=1))
else:
    del_t=np.squeeze(fracSZAx.avg200m_MeanPRad[fracSZAx.VZA==VZA,:,0]-\
                     frac1D.avg200m_MeanPRad[frac1D.VZA==VZA,:,0])
weights = np.ones_like(del_t)/len(del_t)
ax3.hist(del_t,40,weights=weights,histtype='step',label='200m')


ax3.legend(loc='best')
fig3_ttl=fracSZAx.fname.split('.')[0]+'_VZA_'+str(VZA)+'_delta_Reflectance_3Dmin1D'
ax3.set_xlabel(r'$R_{3D}-R_{1D}$')
fig3.suptitle("\n".join(wrap(fig3_ttl,60)))



#=======3d-1d difference end ==================================================
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

