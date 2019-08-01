#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Jun  8 10:12:58 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Step cloud
"""
import matplotlib.pyplot as plt
import numpy as np
import cpnLES_MSCARTlib as lib
from scipy.signal import find_peaks_cwt
from textwrap import wrap
import sys
import cpnCommonlib as cpn
#import netCDF4
#import h5py


band='865'
sza='120'
vza=0
#band=str(sys.argv[1])
#sza=str(sys.argv[2])
#vza=str(sys.argv[3])#vza=d to go with default values

#/umbc/xfs1/zzbatmos/users/charaj1
dstepSZAx=lib.POLCARTdset('dstepSZAx','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Step_Cloud/')
dstepSZAx.readMSCARTplus('step_cld_b'+band+'re10ve02_x15km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e6.nc',fdpath='results/stepb'+band+'/',step=True)
#dstepSZAx.savePOLCARThdf5(dstepSZAx.fdpath+dstepSZAx.fname.split('.',1)[0]+'.hdf5',pc_format='step_cloud',action='For_Dan')
dstep=lib.LES_field(dstepSZAx.fname.split('_MSCART',2)[0]+'.nc')#DYCOM field file for MSCART
dstep.readLES_field()

#column runs (1D)
cfree=lib.POLCARTdset('cfree','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Step_Cloud/1Druns/')
cfree.readMSCARTplus('step_cld_b'+band+'re10ve02_x1Dex0p1_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.nc',\
                       fdpath='1Druns/results/stepb'+band+'/',clm=True,step=True)

cstep=lib.POLCARTdset('cstep','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Step_Cloud/1Druns/')
cstep.readMSCARTplus('step_cld_b'+band+'re10ve02_x1Dex10_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.nc',\
                       fdpath='1Druns/results/stepb'+band+'/',clm=True,step=True)

step1Dto3D=lib.POLCARTdset('step1Dto3D','generated_using_'+cfree.fdpath+cfree.fname)
step1Dto3D.fname=cfree.fname.split('.',1)[0]+'_1Dto3D.nc'
step1Dto3D.fdpath=cfree.fdpath
step1Dto3D.VZA=cfree.VZA
step1Dto3D.ScatA=cfree.ScatA
step1Dto3D.MeanPRad=np.zeros_like(dstepSZAx.MeanPRad)
for i in range(0,200):
    step1Dto3D.MeanPRad[:,i,:]=cfree.MeanPRad
for i in range(200,1200):
    step1Dto3D.MeanPRad[:,i,:]=cstep.MeanPRad
for i in range(1200,1500):
    step1Dto3D.MeanPRad[:,i,:]=cfree.MeanPRad

#%%%%%%%%%%%%%%%%%% figure 01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#I and -Q along x axis.
cpn.setup_figures(plt)
print(sza)
if (band=='865' or band=='2p13') and vza=='d':
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
    
def SS3DIQ_theoretical(thick_mul=1):    
    #Calculating single scattering I and Q for 3D case
    #(Check my 09/11/2017 notes)

#    F0=1
    d_the_scat=dstepSZAx.ScatA[60]
    P11_the_scat=dstep.PMatVal[0,np.where(dstep.scatAng==d_the_scat),0][0][0]
    P21_the_scat=dstep.PMatVal[0,np.where(dstep.scatAng==d_the_scat),4][0][0]
    omega=dstep.omgp3d[0,0,0,10]
    B1=dstep.extp3d[0,0,0,10]
#    print('B1:'+str(B1))
    B2=dstep.extp3d[0,0,0,400]
#    print('B2:'+str(B2))
    mu0=np.cos(np.deg2rad(180-dstepSZAx.SZA))
    mu=np.cos(np.deg2rad(VZA))
     
    KI=omega/4*P11_the_scat/mu0*thick_mul
    KQ=omega/4*P21_the_scat/mu0*thick_mul
    H=1;
    C2=B1/mu0+B2/mu
    C3=B1/mu+B2/mu0
    C1=1/mu0+1/mu
    
    l1=2#left edge
    l2=12#right edge
    tanSZA=np.tan(np.deg2rad(180-dstepSZAx.SZA))
    
    xcens=(dstep.xgrd[1:]+dstep.xgrd[0:-1])/2
#    IssX=np.zeros((xcens.shape[0]+1))
#    QssX=np.zeros((xcens.shape[0]+1))
    
    #left thin
    x1=xcens[xcens<l1]
    I1=KI/C1*(1-np.exp(-B1*C1*H))
    IssXl1=I1
    I1=np.zeros_like(x1)+I1       
    Q1=KQ/C1*(1-np.exp(-B1*C1*H))
    QssXl1=Q1
    Q1=np.zeros_like(x1)+Q1
    ssX=np.append(x1,2)
    IssX=np.append(I1,IssXl1)
    QssX=np.append(Q1,QssXl1)
    
    #illuminating edge
    x2=xcens[(xcens>l1)*(xcens<(l1+H*tanSZA))]
    ssX=np.append(ssX,x2)
    zp=(x2-l1)/tanSZA

    I2=KI*B2*(1/B2/C1*(1-np.exp(-B2*C1*zp))+np.exp(-(B2-B1)*zp/mu0)/C2*(np.exp(-C2*zp)-np.exp(-H*C2)))#*thick_mul
    Q2=KQ*B2*(1/B2/C1*(1-np.exp(-B2*C1*zp))+np.exp(-(B2-B1)*zp/mu0)/C2*(np.exp(-C2*zp)-np.exp(-H*C2)))#*thick_mul
    IssX=np.append(IssX,I2)
    QssX=np.append(QssX,Q2)

    #thick
    x3=xcens[(xcens>(l1+H*tanSZA))*(xcens<l2)]
    I3=KI/C1*(1-np.exp(-B2*C1*H))#*thick_mul
    IssXl2=I3
    I3=np.zeros_like(x3)+I3
    Q3=KQ/C1*(1-np.exp(-B2*C1*H))#*thick_mul
    QssXl2=Q3
    Q3=np.zeros_like(x3)+Q3
    ssX=np.append(ssX,x3)
    ssX=np.append(ssX,l2)
    IssX=np.append(IssX,I3)
    IssX=np.append(IssX,IssXl2)
    QssX=np.append(QssX,Q3)
    QssX=np.append(QssX,QssXl2)
    
    #shadowing
    x4=xcens[((xcens)>l2)*(xcens<l2+H*tanSZA)]
    ssX=np.append(ssX,x4)
    zp=(x4-l2)/tanSZA
    I4=KI*B1*(1/B1/C1*(1-np.exp(-B1*C1*zp))+np.exp((B2-B1)*zp/mu0)/C3*(np.exp(-C3*zp)-np.exp(-H*C3)))
    Q4=KQ*B1*(1/B1/C1*(1-np.exp(-B1*C1*zp))+np.exp((B2-B1)*zp/mu0)/C3*(np.exp(-C3*zp)-np.exp(-H*C3)))
    IssX=np.append(IssX,I4)
    QssX=np.append(QssX,Q4)
    
    #right thin (same as left thin)
    x5=xcens[xcens>l2+H*tanSZA]
    I5=KI/C1*(1-np.exp(-B1*C1*H))
    I5=np.zeros_like(x5)+I5       
    Q5=KQ/C1*(1-np.exp(-B1*C1*H))
    Q5=np.zeros_like(x5)+Q5
    ssX=np.append(ssX,x5)
    IssX=np.append(IssX,I5)
    QssX=np.append(QssX,Q5)
    
    return ssX,IssX,QssX
    
#VZA=20.0 # note that, there are two set of values for 0.0 VZA
xcens=(dstep.xgrd[1:]+dstep.xgrd[0:-1])/2
fig1,ax1=plt.subplots(2,1,sharex=True)
fig1_ttl=dstepSZAx.fname.split('.')[0]+'_VZA_'+str(VZA)
if VZA==0:
#    ax1[0].plot(xcens,np.mean(dstepSZAx.MeanPRad[dstepSZAx.VZA==VZA,:,0].T,axis=1),'b',label=r'$3D$')
    ax1[0].plot(xcens,np.mean(step1Dto3D.MeanPRad[step1Dto3D.VZA==VZA,:,0].T,axis=1),'r',label=r'$1D$')
    ax1[1].plot(xcens,-np.mean(dstepSZAx.MeanPRad[dstepSZAx.VZA==VZA,:,1].T,axis=1),'b')
    ax1[1].plot(xcens,-np.mean(step1Dto3D.MeanPRad[step1Dto3D.VZA==VZA,:,1].T,axis=1),'r')
    
    ssX,IssX,QssX=SS3DIQ_theoretical()  
#    ax1[0].plot(ssX,IssX,'k--',label=r'$3D^{theory}_{SS}$')
    if band=='865':
        m=10;c=0.32
    elif band=='2p13':
        m=10;c=0.18
    elif band=='3p75':
        m=6;c=0.075
#    ax1[0].plot(ssX,IssX*m+c,'k-',label=r'$3D^{theory}_{SS}\times %0.2f+%0.3f$'%(m,c))
#    ax1[1].plot(ssX,-QssX,'k--')
    if band=='865':
        m=4;c=0.002
    elif band=='2p13':
        m=3;c=0.003
    elif band=='3p75':
        m=1;c=0
#    ax1[1].plot(ssX,-QssX*m+c,'k-',label=r'$3D^{theory}_{SS}\times %0.2f+%0.3f$'%(m,c))
else:
    ax1[0].plot(xcens,dstepSZAx.MeanPRad[dstepSZAx.VZA==VZA,:,0].T,'b',label='3D')
    ax1[0].plot(xcens,step1Dto3D.MeanPRad[step1Dto3D.VZA==VZA,:,0].T,'r',label='1D')
    ax1[1].plot(xcens,-dstepSZAx.MeanPRad[dstepSZAx.VZA==VZA,:,1].T,'b',label='3D')
    ax1[1].plot(xcens,-step1Dto3D.MeanPRad[step1Dto3D.VZA==VZA,:,1].T,'r',label='1D')
ax1[0].set_ylabel('I',size=10)
for ln in ax1[0].lines:
    ln.set_linewidth(4.5)
ax1[1].legend(loc='best',prop={'size':9})
ax1[0].legend(loc='best',prop={'size':9})
#ax2 = ax1[0].twinx()
#ax2.plot(xcens,np.squeeze(dstep.extp3d),'r--',linewidth=2.0)
#ax2.set_ylabel(r'Optical depth $\tau$',size=10)
#ax2.set_ylim(0,15)

ax1[1].set_ylabel('-Q',size=10)
ax1[1].set_xlabel('x km',size=10)
for ln in ax1[1].lines:
    ln.set_linewidth(1.5)

#ax22 = ax1[1].twinx()
#ax22.plot(xcens,np.squeeze(dstep.extp3d),'r--',linewidth=2.0)
#ax22.set_ylabel(r'Optical depth $\tau$',size=10)
#ax22.set_ylim(0,15)
fig1.suptitle("\n".join(wrap(fig1_ttl,60)))
fig1.show()
#%%%%%%%%%%%%%%%%%% figure 02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# I and -Q near the step edges w.r.t. SZA and apparent scattering angle
def stepEdge(fig,ax,x,delt,iqu,ap_ttl):
    if iqu==1:
        sgn=-1
    else:
        sgn=1
        
    if x==200:
        neg_del='b';pos_del='r'
    elif x==1200:
        neg_del='r';pos_del='b'
    else:
        print('!!Adjust stepEdge() function according to the step location!!')
    fig_ttl=dstepSZAx.fname.split('.')[0]+ap_ttl
    fig.suptitle(fig_ttl)
    fig.subplots_adjust(top=0.85)
    ax.plot(dstepSZAx.VZA,sgn*dstepSZAx.MeanPRad[:,x-2*delt,iqu],':',color=neg_del,linewidth=1.5,label='x=%0.2f'%((x-2*delt)*1e-2))
    ax.plot(dstepSZAx.VZA,sgn*dstepSZAx.MeanPRad[:,x-delt,iqu],'--',color=neg_del,linewidth=1.5,label='x=%0.2f'%((x-delt)*1e-2))
    ax.plot(dstepSZAx.VZA,sgn*dstepSZAx.MeanPRad[:,x,iqu],color='k',linewidth=2.0,label='x=%0.1f'%(x*1e-2))
    ax.plot(dstepSZAx.VZA,sgn*dstepSZAx.MeanPRad[:,x+delt,iqu],'--',color=pos_del,linewidth=1.5,label='x=%0.2f'%((x+delt)*1e-2))
    ax.plot(dstepSZAx.VZA,sgn*dstepSZAx.MeanPRad[:,x+2*delt,iqu],':',color=pos_del,linewidth=1.5,label='x=%0.2f'%((x+2*delt)*1e-2))
    ax.plot(cfree.VZA,sgn*cfree.MeanPRad[:,iqu],'b',label='1D w/o step')
    ax.plot(cstep.VZA,sgn*cstep.MeanPRad[:,iqu],'r',label='1D step')
    ax.legend(loc='best',prop={'size':10})
    return fig_ttl
    
def iq_err(x1,x2,iqu):
    #x1: illuminating step
    #x2: shadowing step
    if iqu==1:
        sgn=-1;app_ttl='-Q'
    else:
        sgn=1;app_ttl='I'
    fig,ax=plt.subplots(2,1)
    fig_ttl=dstepSZAx.fname.split('.')[0]+'_'+app_ttl+'_error'
    fig.suptitle(fig_ttl)
    fig.subplots_adjust(top=0.85)
    ax[0].errorbar(dstepSZAx.VZA,sgn*dstepSZAx.MeanPRad[:,x1,iqu],sgn*dstepSZAx.RMSEPRad[:,x1,iqu],fmt='r')
    ax[0].xaxis.set_ticks_position('top') 
    ax2=ax[0].twinx()
    ax2.plot(dstepSZAx.VZA,np.abs(dstepSZAx.RMSEPRad[:,x1,iqu]/dstepSZAx.MeanPRad[:,x1,iqu])*100,'g')
    ax[0].set_title('Illuminating step',y=1.08)
    ax[0].set_ylabel(app_ttl)
    ax2.set_ylabel('error(%)')
    ax[1].errorbar(dstepSZAx.VZA,sgn*dstepSZAx.MeanPRad[:,x2,iqu],sgn*dstepSZAx.RMSEPRad[:,x2,iqu],fmt='r')
    ax22=ax[1].twinx()
    ax22.plot(dstepSZAx.VZA,np.abs(dstepSZAx.RMSEPRad[:,x2,iqu]/dstepSZAx.MeanPRad[:,x2,iqu])*100,'g')
    ax[1].set_title('Shadowing step')
    ax[1].set_ylabel(app_ttl)
    ax22.set_ylabel('error(%)')
    ax2.legend(['error(%)'],loc='best',prop={'size':10})
    ax[1].set_xlabel('VZA',size=10)
    
    return fig_ttl,fig
        
def tick_fun(VZA,SZA,SAA):
#    V=180-(SZA+np.cos(np.deg2rad(SAA))*VZA)
    V=np.zeros_like(VZA)
    for i in range(np.size(VZA)):
        if VZA[i]>0:
            V[i]=lib.scat_ang(SZA,VZA[i],SAA-0)
        else:
            V[i]=lib.scat_ang(SZA,-VZA[i],SAA-180)
    
    return ["%d"%z for z in V]

#%%%%%%%%%%%%%%%%%%%% Contour plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        cb=fig.colorbar(ctf,ax=ax,orientation='vertical',ticks=np.arange(-0.05,0.16,0.05))
    else:
        cb=fig.colorbar(ctf,ax=ax,orientation='vertical',ticks=ticks)
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
    return fig_ttl,ctf,cb
def put_scat_ang(ax2,dstepSZA):
    ax2t=ax2.twiny()
    ax2ticks=ax2.get_xticks()
    ax2tticks=ax2ticks
    ax2t.set_xticks(ax2tticks)
    ax2t.set_xbound(ax2.get_xbound())
    ax2t.set_xticklabels(tick_fun(ax2tticks,dstepSZA.SZA,dstepSZA.SAA))
    ax2.set_xlabel('VZA',size=10)
    ax2t.set_xlabel('Apparent Scat. Angle',size=10)
    
#fig2,ax2=plt.subplots()
#fig2_ttl=stepEdge(fig2,ax2,200,10,1,'_Q_atIlluStep')
#put_scat_ang(ax2,dstepSZAx)
#ax2.set_ylabel('-Q',size=10)
#
#fig3,ax3=plt.subplots()
#fig3_ttl=stepEdge(fig3,ax3,200,10,0,'_I_atIlluStep')
#put_scat_ang(ax3,dstepSZAx)
#ax3.set_ylabel('I',size=10)
#
#fig4,ax4=plt.subplots()
#fig4_ttl=stepEdge(fig4,ax4,1200,5,0,'_IatShadStep')
#put_scat_ang(ax4,dstepSZAx)
#ax4.set_ylabel('I',size=10)
#
#fig5,ax5=plt.subplots()
#fig5_ttl=stepEdge(fig5,ax5,1200,5,1,'_QatShadStep')
#put_scat_ang(ax5,dstepSZAx)
#ax5.set_ylabel('-Q',size=10)
#
#fig6_ttl,fig6=iq_err(200,1200,0)
#fig7_ttl,fig7=iq_err(200,1200,1)

def find_bows(ang,P11):
    #ang: scattering angle in degrees
    #P11: phase function
    plt.figure()
    axx=plt.gca();figg=plt.gcf()
    axx.plot(P11)
    axx.set_yscale('log')
    figg.show()
    left=input('Left most index:')
    right=input('Right most index:')
    data=P11[left:right]
    xx=ang[left:right]
    plt.figure()
    ax=plt.gca();fig=plt.gcf()
    ax.plot(xx,data)
#    ax.set_yscale('log')
    ixx=find_peaks_cwt(data,np.arange(1,50))
    ax.plot(xx[ixx],data[ixx],'r+')
    fig.show()
    print(xx[ixx])
    pb=input('Select Primarybow angle:')
    sb=input('Select Supernumerarybow angle:')
    pb_ix=np.where(np.isclose(dstepSZAx.ScatA,pb,atol=0.1))
    sb_ix=np.where(np.isclose(dstepSZAx.ScatA,sb,atol=0.1))
    plt.figure()
    return pb_ix,sb_ix

if band=='865':
    pb_ix=np.where(np.isclose(dstepSZAx.ScatA,142,atol=0.1))
    sb_ix=np.where(np.isclose(dstepSZAx.ScatA,153,atol=0.1))
elif band=='2p13':
    pb_ix=np.where(np.isclose(dstepSZAx.ScatA,141,atol=0.1))
    sb_ix=np.where(np.isclose(dstepSZAx.ScatA,160,atol=0.1))
elif band=='3p75':
    pb_ix=np.where(np.isclose(dstepSZAx.ScatA,158,atol=0.1))
    sb_ix=np.where(np.isclose(dstepSZAx.ScatA,170,atol=0.1))
else:    
    pb_ix,sb_ix=find_bows(dstep.scatAng,dstep.PMatVal[0,:,0])

fig8,ax8=plt.subplots(3,1,figsize=(7,10),sharex=True)

def getQcontourf_vticks(band):
    if band=='3p75':
        v=np.linspace(-0.05,0.05,num=512)
        ticks=np.arange(-0.05,0.051,0.025)
        band_given=True
    else:
        v=0;ticks=0
        band_given=False
        
    return v,ticks,band_given

v,ticks,bg=getQcontourf_vticks(band)
if bg:
    _,ctf,cb=contourf_Iqu(fig8,ax8[0],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,dstepSZAx.MeanPRad,\
                      v=v,ticks=ticks)
    contourf_Iqu(fig8,ax8[1],1,step1Dto3D.VZA[pb_ix],step1Dto3D.VZA[sb_ix],step1Dto3D.VZA,step1Dto3D.MeanPRad,\
                 v=v,ticks=ticks)
else:    
    _,ctf,cb=contourf_Iqu(fig8,ax8[0],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,dstepSZAx.MeanPRad)
    contourf_Iqu(fig8,ax8[1],1,step1Dto3D.VZA[pb_ix],step1Dto3D.VZA[sb_ix],step1Dto3D.VZA,step1Dto3D.MeanPRad)
    
if band=='865' and sza=='120':
    fig8_ttl,_,_=contourf_Iqu(fig8,ax8[2],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,\
                      dstepSZAx.MeanPRad-step1Dto3D.MeanPRad,v=np.linspace(-0.06,0.06,num=512),\
                        ticks=np.arange(-0.06,0.061,0.06),cmap=plt.cm.coolwarm)
elif band=='865' and sza=='140':
    fig8_ttl,_,_=contourf_Iqu(fig8,ax8[2],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,\
                      dstepSZAx.MeanPRad-step1Dto3D.MeanPRad,v=np.linspace(-0.05,0.05,num=512),\
                        ticks=np.arange(-0.05,0.051,0.05),cmap=plt.cm.coolwarm)
elif band=='865' and sza=='160':
    fig8_ttl,_,_=contourf_Iqu(fig8,ax8[2],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,\
                      dstepSZAx.MeanPRad-step1Dto3D.MeanPRad,v=np.linspace(-0.01,0.01,num=512),\
                        ticks=np.arange(-0.01,0.011,0.01),cmap=plt.cm.coolwarm)
elif band=='2p13' and sza=='120':
    fig8_ttl,_,_=contourf_Iqu(fig8,ax8[2],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,\
                      dstepSZAx.MeanPRad-step1Dto3D.MeanPRad,v=np.linspace(-0.1,0.1,num=512),\
                        ticks=np.arange(-0.1,0.11,0.1),cmap=plt.cm.coolwarm)
elif band=='2p13' and sza=='140':
    fig8_ttl,_,_=contourf_Iqu(fig8,ax8[2],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,\
                      dstepSZAx.MeanPRad-step1Dto3D.MeanPRad,v=np.linspace(-0.04,0.04,num=512),\
                        ticks=np.arange(-0.04,0.041,0.04),cmap=plt.cm.coolwarm)
elif band=='2p13' and sza=='160':
    fig8_ttl,_,_=contourf_Iqu(fig8,ax8[2],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,\
                      dstepSZAx.MeanPRad-step1Dto3D.MeanPRad,v=np.linspace(-0.02,0.02,num=512),\
                        ticks=np.arange(-0.02,0.021,0.02),cmap=plt.cm.coolwarm)
elif band=='3p75' and sza=='120':
    fig8_ttl,_,_=contourf_Iqu(fig8,ax8[2],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,\
                      dstepSZAx.MeanPRad-step1Dto3D.MeanPRad,v=np.linspace(-0.05,0.05,num=512),\
                        ticks=np.arange(-0.05,0.051,0.05),cmap=plt.cm.coolwarm)
    
else:
    fig8_ttl,_,_=contourf_Iqu(fig8,ax8[2],1,dstepSZAx.VZA[pb_ix],dstepSZAx.VZA[sb_ix],dstepSZAx.VZA,\
                      dstepSZAx.MeanPRad-step1Dto3D.MeanPRad,v=np.linspace(-0.06,0.06,num=512),\
                        ticks=np.arange(-0.06,0.061,0.06),cmap=plt.cm.coolwarm)
    
ax8[2].set_xlabel('x (km)',size=10,y=-0.5)
fig8_ttl=dstepSZAx.fname.split('.',1)[0]+'_'+fig8_ttl
fig8.suptitle("\n".join(wrap(fig8_ttl,60)))
ax8[0].set_title(r'$-Q_{3D}$')
ax8[1].set_title(r'$-Q_{1D}$')
ax8[2].set_title(r'$-(Q_{3D}-Q_{1D})$')


def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'results/stepb'+band+'/figures/')
#savefig(fig1,fig1_ttl)
#savefig(fig2,fig2_ttl)
#savefig(fig3,fig3_ttl)
#savefig(fig4,fig4_ttl)
#savefig(fig5,fig5_ttl)
#savefig(fig6,fig6_ttl)
#savefig(fig7,fig7_ttl)
#fig8.show()
#savefig(fig8,fig8_ttl)

plt.show()
