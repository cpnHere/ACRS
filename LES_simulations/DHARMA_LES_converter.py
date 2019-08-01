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

05/09/2018:
    To convert DHARMA outputs to do furthur analysis. Especially to generate/convert
    required variables for MSCART inputs.
    

"""
from __future__ import print_function
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import cpnCommonlib as cpn
import h5py
import netCDF4, time
from cpnLES_MSCARTlib import LES_field, DHARMA_onmp
import cpnMielib as Mie
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D

  #dharma_005263.cdf  dharma_007877.cdf  dharma_0.cdf  dharma_012540.cdf

        
def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,path='figures/')
    
def lwcHSlices(fig,ax,case,v,ts):
    '''
    LWC 2D contourf at several altitudes based on the size of ax.
    fig: figure instance
    ax: subplots axeses 
    case: DHARMA_onmp object
    v: np.linspace() for color steps
    ts: np.arange() for colorbar ticks
    #----------------------------------------------------------------------------------------------------------------------------------
    Ex:-
        fig2,ax2=plt.subplots(2,4,figsize=(8,4),subplot_kw={'aspect':'equal'})
        fig2_ttl='Tsc_lwc'
        lwcHSlices(fig2,ax2,ToweringSc,np.linspace(0,0.25,50),np.arange(0,0.26,0.10))
        fig2.suptitle(fig2_ttl)
        fig2.show()
    '''
    axs=ax.flat
    z=case.DHARMA.z
    step=z.size/ax.size
    i=0
    for n, axn in enumerate(axs):
        ctf=axn.contourf(case.x/1e3,case.y/1e3,case.lwc[step,:,:],v,cmap=plt.cm.jet);
        axn.set_title('Altitude: %d m'%z[step],size=10)
        i+=1
        step+=z.size/ax.size
    print(step)
    cpn.add_common_cb(fig,ctf,ts,r'LWC ($g/m^3$)')

def com_w_Dan():
    cpn.setup_figures(plt)
    ATEXl040_path ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/ATEX/Ackerman_ATEX/CCN40/'
    ATEXl600_path ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/ATEX/Ackerman_ATEX/CCN600/'
    DYCOMS2_path  ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/DYCOMS/Ackerman_DYCOMS2_25bins/'
    RICOn_path    ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/RICO/Ackerman_RICO_NEW/'
    
    ATEX040=DHARMA_onmp('dharma_004997.cdf',ATEXl040_path);ATEX040.readDHARMA()
    ATEX600=DHARMA_onmp('dharma_004996.cdf',ATEXl600_path);ATEX600.readDHARMA()
    DYCOMS2=DHARMA_onmp('dharma_008036.cdf',DYCOMS2_path );DYCOMS2.readDHARMA()
    RICO   =DHARMA_onmp('dharma_005044.cdf',RICOn_path   );RICO.readDHARMA()
    
    lesf=LES_field('OP_dharma_008036_full_3_26.nc',dpath='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/')#DYCOM field file for MSCART
    lesf.readLES_field()    
    
    DYD=DHARMA_onmp('dharma_008036.cdf',DYCOMS2_path)#To compare with Dan's calculations
    f=h5py.File(DYCOMS2_path+'Ackerman_DYCOMS2_25bins_dharma_008036_Physics_Optics.h5','r')
    DYD.lwp=f['LWP_50m'][:].swapaxes(0,1)
    DYD.lwc=f['LWC/T008036'][:].swapaxes(1,2)*1e6#in g/m^3
    DYD.tau=f['TAU_W086_total_50m'][:].swapaxes(0,1)
    DYD.DTAU=f['TAU_W086'][:].swapaxes(1,2)
    DYD.x=f['X_50m'][:].reshape(128)
    DYD.y=f['Y_50m'][:].reshape(128)
    f.close()
        
# LWP and tau comparison with Dan's 
    fig1,ax1=plt.subplots(1,2,figsize=(9,4),subplot_kw={'aspect':'equal'})
    cmap=plt.cm.jet
    fig1_ttl='DYCOM2_lwp_comp_w_Dans_right'
    fig1.suptitle(fig1_ttl)
    v =np.linspace(0,300,50)
    ts=np.arange  (0,300.1,100)
    ctf=ax1[0].contourf(DYCOMS2.x/1e3,DYCOMS2.y/1e3,DYCOMS2.lwp,v,cmap=cmap)
    ctf=ax1[1].contourf(DYCOMS2.x/1e3,DYCOMS2.y/1e3,DYD.lwp,v,cmap=cmap)
    cpn.add_common_cb(fig1,ctf,ts,r'LWP(g/m^2)')
    ax1[0].set_ylabel('Y(km)')
    ax1[0].set_xlabel('X(km)')
    ax1[1].set_xlabel('X(km)')
    fig1.show()

#LWC comparison with Dan's
    fig2,ax2=plt.subplots(3,3,figsize=(8,8),subplot_kw={'aspect':'equal'})
    fig2_ttl='DYCOMS2_lwc'
    lwcHSlices(fig2,ax2,DYCOMS2,np.linspace(0,1.25,50),np.arange(0,1.26,0.25))
    fig2.suptitle(fig2_ttl)
    fig2.show()
        
    fig3,ax3=plt.subplots(3,3,figsize=(8,8),subplot_kw={'aspect':'equal'})
    fig3_ttl='DYCOMS2_lwc_Dan'
    lwcHSlices(fig3,ax3,DYD,np.linspace(0,1.25,50),np.arange(0,1.26,0.25))
    fig3.suptitle(fig3_ttl)
    fig3.show()

def Ex1_DYCOMS2_case(runMie=False,MieName='DYCOM2_dharma_008036_mie_470_860_2p13'):
    '''
    Example case to convert DYCOMS-II DHARMA outputs to generate MSCART input files.
    First example
    -----------------------------------------------------
    Previously generated Mie runs' file names
    [1] 'DYCOM2_dharma_008036_mie_470_860_2p13'
    -----------------------------------------------------
    When the Mie run results exist, do as follows to generate properties that are required to 
    produce MSCART input field file.
        extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=DYCOMS2.setupMieOut(MieName,lesCname='DYCOMS2')
    '''
    DYCOMS2_path  ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/DYCOMS/Ackerman_DYCOMS2_25bins/'  
    DYCOMS2=DHARMA_onmp('dharma_008036.cdf',DYCOMS2_path );DYCOMS2.readDHARMA()
    if runMie==False:
        print('Already generated a Mie output: '+MieName)
        print('Required properties to generate MSCART fieldfiles can be obtained from setupMieOut()')
        #extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=DYCOMS2.setupMieOut(MieName,lesCname='DYCOMS2')
    else:
        #----------------clouds-----------------------    
        cloud_mie=Mie.MieSet(MieName)
    
        refs={'lam':[0.443,0.470,0.532,0.550,0.660,0.860,0.865,0.960,1.064,1.24,1.63,2.10,2.13,3.75],\
              'nr':[1.344590,1.341951,1.337115,1.335943,1.330267,1.324481,1.324372,1.322450,1.320509,1.317240,1.308836,1.291837,1.290108,1.351865],\
              'nc':[0.8905030E-09,0.7285308E-09,0.1818175E-08,0.2461249E-08,0.1920720E-07,0.3380264E-06,0.3547385E-06,\
            0.3367559E-05,0.1279514E-05,0.1134782E-04,0.8084244E-04,0.4615043E-03,0.3941974E-03,0.3402400E-02]}
        lam=[refs['lam'][1],refs['lam'][5],refs['lam'][12]]
        nr =[refs[ 'nr'][1],refs[ 'nr'][5],refs[ 'nr'][12]]
        nc =[refs[ 'nc'][1],refs[ 'nc'][5],refs[ 'nc'][12]]
    
        r=np.linspace(DYCOMS2.DHARMA.rbound_drops.min(),DYCOMS2.DHARMA.rbound_drops.max(),10000)#radias um
        D=2*r
    #    ---------------------------------------
        na=1801
        ang=np.linspace(0,180,na)
        ang=np.hstack((np.array([na]),ang))
        np.savetxt('ang.dat',ang,delimiter='\n',fmt='%.2f')
        print('Running Mie code')
        cloud_mie.runMie(lam,nr,nc,D)
    return DYCOMS2
    
def Ex2_ToweringSc_case(runMie=False,MieName='ToweringSc_dharma_001800_mie_470_860_2p13_v2'):
    '''
    Example case to convert ToweringSc DHARMA outputs to generate MSCART input files.
    Second example
    -----------------------------------------------------
    Previously generated Mie runs' file names
    [1] 'ToweringSc_dharma_001800_mie_470_860_2p13_v2'
    -----------------------------------------------------
    When the Mie run results exist, do as follows to generate properties that are required to 
    produce MSCART input field file.
        extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=ToweringSc.setupMieOut('ToweringSc_dharma_001800_mie_470_860_2p13',lesBinEd=lesBinEd)
    ***************************************************************************
    08/14/2018:
        For the ToweringSc case, I have Mie runs that are corresponded to RSP bands from Dan.
        Use DHARMA_onmp.readNsetFromRSPMieOut() directly to get the properties of that are needed to setup MSCART field_file.
    '''
    
    #Ex2_TowSc_case
    ToweringSc=DHARMA_onmp('dharma_001800.cdf','/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/Ann_Dan/')
    ToweringSc.readDHARMA()
    cloud_mie=Mie.MieSet(MieName,path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/')
    if runMie:
        #----------------clouds-----------------------    
    
        refs={'lam':[0.443,0.470,0.532,0.550,0.660,0.860,0.865,0.960,1.064,1.24,1.63,2.10,2.13,3.75],\
              'nr':[1.344590,1.341951,1.337115,1.335943,1.330267,1.324481,1.324372,1.322450,1.320509,1.317240,1.308836,1.291837,1.290108,1.351865],\
              'nc':[0.8905030E-09,0.7285308E-09,0.1818175E-08,0.2461249E-08,0.1920720E-07,0.3380264E-06,0.3547385E-06,\
            0.3367559E-05,0.1279514E-05,0.1134782E-04,0.8084244E-04,0.4615043E-03,0.3941974E-03,0.3402400E-02]}
        lam=[refs['lam'][1],refs['lam'][5],refs['lam'][12]]
        nr =[refs[ 'nr'][1],refs[ 'nr'][5],refs[ 'nr'][12]]
        nc =[refs[ 'nc'][1],refs[ 'nc'][5],refs[ 'nc'][12]]
    
        r=np.concatenate((np.linspace(ToweringSc.DHARMA.rbound_drops.min(),50,500),np.linspace(50,ToweringSc.DHARMA.rbound_drops.max(),9501)[1:]),axis=0)#radias um
        D=2*r
    #    ---------------------------------------
        na=1801
        ang=np.linspace(0,180,na)
        ang=np.hstack((np.array([na]),ang))
        np.savetxt('ang.dat',ang,delimiter='\n',fmt='%.2f')
        print('Running Mie code')
        cloud_mie.runMie(lam,nr,nc,D)
    else:
        print('Already generated a Mie output: '+MieName)
        print('Required properties to generate MSCART fieldfiles can be obtained from setupMieOut(). Input the lesBinEd too.')
        #extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=ToweringSc.setupMieOut('ToweringSc_dharma_001800_mie_470_860_2p13',lesBinEd=lesBinEd)
        
    cloud_mie.readMie()
    r=cloud_mie.d/2
    lesBinEd=np.zeros_like(ToweringSc.DHARMA.rbound_drops,dtype=int)
    for i in np.arange(0,ToweringSc.DHARMA.rbound_drops.size):
        lesBinEd[i]=(abs(r-ToweringSc.DHARMA.rbound_drops[i])).argmin()
    return ToweringSc,lesBinEd

def Ex3_ATEXc_case(runMie=False,MieName='DYCOM2_dharma_008036_mie_470_860_2p13'):
    '''
    Example case to convert ATEXc DHARMA outputs to generate MSCART input files.
    Second example
    -----------------------------------------------------
    Previously generated Mie runs' file names
    [1] 'ToweringSc_dharma_001800_mie_470_860_2p13_v2'
    -----------------------------------------------------
    When the Mie run results exist, do as follows to generate properties that are required to 
    produce MSCART input field file.
        extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=ToweringSc.setupMieOut(MieName,lesBinEd=lesBinEd)

    '''
    ATEXc_path  ='/umbc/xfs1/zzbatmos/users/mild1/Research/LES/LES_Satellite_Simulator_DJM/LES_Simulations/Ackerman_ATEX_highres/CCN40/'  
    ATEXc=DHARMA_onmp('dharma_007877.cdf',ATEXc_path );ATEXc.readDHARMA()
    cloud_mie=Mie.MieSet(MieName,path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/')
    if runMie==False:
        print('Already generated a Mie output: '+MieName)
        print('Required properties to generate MSCART fieldfiles can be obtained from setupMieOut()')
        #extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=ATEXc.setupMieOut(MieName,lesCname='ATEXc')
    else:
        #----------------clouds-----------------------        
        refs={'lam':[0.443,0.470,0.532,0.550,0.660,0.860,0.865,0.960,1.064,1.24,1.63,2.10,2.13,3.75],\
              'nr':[1.344590,1.341951,1.337115,1.335943,1.330267,1.324481,1.324372,1.322450,1.320509,1.317240,1.308836,1.291837,1.290108,1.351865],\
              'nc':[0.8905030E-09,0.7285308E-09,0.1818175E-08,0.2461249E-08,0.1920720E-07,0.3380264E-06,0.3547385E-06,\
            0.3367559E-05,0.1279514E-05,0.1134782E-04,0.8084244E-04,0.4615043E-03,0.3941974E-03,0.3402400E-02]}
        lam=[refs['lam'][1],refs['lam'][5],refs['lam'][12]]
        nr =[refs[ 'nr'][1],refs[ 'nr'][5],refs[ 'nr'][12]]
        nc =[refs[ 'nc'][1],refs[ 'nc'][5],refs[ 'nc'][12]]
    
        r=np.linspace(ATEXc.DHARMA.rbound_drops.min(),ATEXc.DHARMA.rbound_drops.max(),10000)#radias um
        D=2*r
    #    ---------------------------------------
        na=1801
        ang=np.linspace(0,180,na)
        ang=np.hstack((np.array([na]),ang))
        np.savetxt('ang.dat',ang,delimiter='\n',fmt='%.2f')
        print('Running Mie code')
        cloud_mie.runMie(lam,nr,nc,D)

    cloud_mie.readMie()
    r=cloud_mie.d/2
    lesBinEd=np.zeros_like(ATEXc.DHARMA.rbound_drops,dtype=int)
    for i in np.arange(0,ATEXc.DHARMA.rbound_drops.size):
        lesBinEd[i]=(abs(r-ATEXc.DHARMA.rbound_drops[i])).argmin()
    return ATEXc,lesBinEd    

if __name__=='__main__':   
    DYCOMS2=Ex1_DYCOMS2_case(MieName='DYCOM2_dharma_008036_mie_470_860_2p13')
    ToweringSc,lesBinEd=Ex2_ToweringSc_case(MieName='ToweringSc_dharma_001800_mie_470_860_2p13')    
    ATEXc,lesBinEd=Ex3_ATEXc_case(MieName='DYCOM2_dharma_008036_mie_470_860_2p13')
#    lesf=LES_field('OP_dharma_008036_full_3_26.nc',dpath='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/')#DYCOM field file for MSCART
#    lesf.readLES_field()
    '''
    Example case to convert RICO DHARMA outputs to generate MSCART input files.
    Second example
    -----------------------------------------------------
    When the Mie run results exist, do as follows to generate properties that are required to 
    produce MSCART input field file.
        extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=ToweringSc.setupMieOut(MieName,lesBinEd=lesBinEd)

    '''
    runMie=False;MieName='RICO_dharma_005044_mie_470_860_2p13'
    RICO_path='/umbc/xfs1/zzbatmos/users/mild1/Research/LES/LES_Satellite_Simulator_DJM/LES_Simulations/Ackerman_RICO_NEW/'
    RICO=DHARMA_onmp('dharma_005044.cdf',RICO_path );RICO.readDHARMA()
    if runMie==False:
        print('Already generated a Mie output: '+MieName)
        print('Required properties to generate MSCART fieldfiles can be obtained from setupMieOut()')
        #extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=RICO.setupMieOut(MieName,lesCname='RICO')
    else:
        #----------------clouds-----------------------        
        refs={'lam':[0.443,0.470,0.532,0.550,0.660,0.860,0.865,0.960,1.064,1.24,1.63,2.10,2.13,3.75],\
              'nr':[1.344590,1.341951,1.337115,1.335943,1.330267,1.324481,1.324372,1.322450,1.320509,1.317240,1.308836,1.291837,1.290108,1.351865],\
              'nc':[0.8905030E-09,0.7285308E-09,0.1818175E-08,0.2461249E-08,0.1920720E-07,0.3380264E-06,0.3547385E-06,\
            0.3367559E-05,0.1279514E-05,0.1134782E-04,0.8084244E-04,0.4615043E-03,0.3941974E-03,0.3402400E-02]}
        lam=[refs['lam'][1],refs['lam'][5],refs['lam'][12]]
        nr =[refs[ 'nr'][1],refs[ 'nr'][5],refs[ 'nr'][12]]
        nc =[refs[ 'nc'][1],refs[ 'nc'][5],refs[ 'nc'][12]]
    
        r=np.linspace(RICO.DHARMA.rbound_drops.min(),RICO.DHARMA.rbound_drops.max(),10000)#radias um
        D=2*r
    #    ---------------------------------------
        na=1801
        ang=np.linspace(0,180,na)
        ang=np.hstack((np.array([na]),ang))
        np.savetxt('ang.dat',ang,delimiter='\n',fmt='%.2f')
        print('Running Mie code')
        cloud_mie=Mie.MieSet(MieName)
        cloud_mie.runMie(lam,nr,nc,D)

    cloud_mie=Mie.MieSet(MieName,path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_simulations/')
    cloud_mie.readMie()
    r=cloud_mie.d/2
    lesBinEd=np.zeros_like(RICO.DHARMA.rbound_drops,dtype=int)
    for i in np.arange(0,RICO.DHARMA.rbound_drops.size):
        lesBinEd[i]=(abs(r-RICO.DHARMA.rbound_drops[i])).argmin()

 

    extp3d,omgp3d,P11,P33,P12,P34,cloud_mie=RICO.setupMieOut('RICO_dharma_005044_mie_470_860_2p13',\
                            lesCname='RICO',mie_path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_simulations/')
  